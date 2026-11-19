"""Forecast cells (FYRF) from Open-Meteo's open data on AWS (s3://openmeteo, CC BY 4.0).

Reads the latest completed run of each model straight from the bucket's data_spatial layout with
range requests (only the variables we use) and writes one deflated FYRF file per unit (ADR-098):
coastal 1-degree cells (the chart's coast cells) and 5-degree ocean tiles for passages:

  metno_nordic_pp  MET Nordic, 1 km, hourly, ~58 h      -> section 1, cells only, 0.01 x 0.02 degree grid
  ecmwf_ifs        ECMWF IFS HRES, O1280 (~9 km), 10 d  -> section 2, 0.1 degree (cells) / 0.25 (tiles)
  ecmwf_wam        ECMWF WAM waves, O1280, 10 d         -> section 3, as IFS

Output: OUTDIR/<run id>/c_<lat>_<lon>.fyrf.z (cells), OUTDIR/<run id>/o_<lat>_<lon>.fyrf.z (tiles, named by
their south-west corner, multiples of 5) and OUTDIR/index.json (written last: run id, model runs, units).

FYRF v2 (little-endian; the whole file is raw DEFLATE, like the other .z kinds):
  "FYRF" | u8 version=2 | i16 lat | i16 lon | u8 size (degrees: 1 cell, 5 tile) | u8 nSections
  per section: u8 model | u32 reference time (unix s) | u16 nSteps | u16 rows | u16 cols
               | f64 lat0, lon0, dlat, dlon | u8 nFields | u8 field[nFields] | u16 hour[nSteps]
  then per section, per field: rows x cols x nSteps bytes, row 0 = south, point-major, each point's
  steps stored as differences from the step before (mod 256; the first from 0).
Grids include both cell edges (rows = 1/dlat + 1), so a point inside the cell interpolates from one file.
Bytes (255 = no value, e.g. waves on land or outside a regional model):
  1 wind speed 0.5 m/s   2 wind direction (from) 360/255 deg   3 gust 0.5 m/s   4 msl hPa - 900
  5 precipitation since the step before (a 3- or 6-hourly step holds the 3 or 6 h total): 0.1 mm to 10 mm,
    then 100 + 15 * sqrt(mm - 10)
  6 CAPE 20 J/kg   7 wave height 0.1 m   8 wave direction (from) 360/255 deg   9 wave period 0.1 s
  10 temperature 2 m 0.5 C from -40 C   11 current speed 0.02 m/s   12 current direction (towards) 360/255 deg
  13 visibility: 100 m to 10 km, then 100 + km - 10 (to 164 km)

Usage:
  om_forecast.py build OUTDIR [--cells FILE | --box LAT0 LAT1 LON0 LON1] [--tiles all | --tiles-box LAT0 LAT1 LON0 LON1]
  om_forecast.py check OUTDIR [N]      # re-reads N random values per section from the bucket
  om_forecast.py latest                # the run id a build would make now (newest complete IFS run)
  om_forecast.py dump FILE             # header and the first point's series
"""
import argparse, asyncio, json, math, os, random, struct, sys, time, urllib.request, zlib
from datetime import datetime, timedelta, timezone
import numpy as np, fsspec, omfiles
import s3fs

HTTP = 'https://openmeteo.s3.amazonaws.com'
FS = fsspec.filesystem('s3', anon=True)          # listings; data reads use an async client per step
MAGIC, VERSION = b'FYRF', 2
TILE, TILE_D = 5, 0.25                          # ocean tiles: 5 degrees at 0.25 degree
NONE = 255
KN = 1.943844

# MARK: models

class Nordic:
    id, name, dlat, dlon, horizon, reach = 1, 'metno_nordic_pp', 0.01, 0.02, 72, 48
    fields = [1, 2, 3, 4, 5, 10]
    nx, ny, R = 1796, 2321, 6371229.0
    _p1, _l0 = math.radians(63.0), math.radians(15.0)
    _n = math.sin(_p1)
    _F = math.cos(_p1) * math.tan(math.pi / 4 + _p1 / 2) ** _n / _n
    _rho0 = R * _F / math.tan(math.pi / 4 + _p1 / 2) ** _n

    @classmethod
    def _xy(cls, lat, lon):
        rho = cls.R * cls._F / np.tan(np.pi / 4 + np.radians(lat) / 2) ** cls._n
        t = cls._n * (np.radians(lon) - cls._l0)
        return rho * np.sin(t), cls._rho0 - rho * np.cos(t)

    @classmethod
    def index(cls, lat, lon):
        """(row, col) into the (ny, nx) arrays, row 0 = south; -1 outside the grid."""
        x0, y0 = cls._xy(52.302723, 1.918457)        # first grid point (south-west corner)
        x, y = cls._xy(lat, lon)
        j, i = np.round((y - y0) / 1000).astype(int), np.round((x - x0) / 1000).astype(int)
        out = (j < 0) | (j >= cls.ny) | (i < 0) | (i >= cls.nx)
        return np.where(out, -1, j), np.where(out, -1, i)

    variables = ['wind_speed_10m', 'wind_direction_10m', 'wind_gusts_10m', 'pressure_msl', 'precipitation',
                 'temperature_2m']

    @staticmethod
    def values(v):
        return {1: v['wind_speed_10m'], 2: v['wind_direction_10m'], 3: v['wind_gusts_10m'],
                4: v['pressure_msl'], 5: v['precipitation'], 10: v['temperature_2m']}

class O1280:
    """ECMWF's octahedral reduced Gaussian grid: 2560 rows, 4i + 16 points on row i from each pole."""
    _nl = np.array([4 * (k + 1) + 16 if k < 1280 else 4 * (2560 - k) + 16 for k in range(2560)])
    _off = np.concatenate([[0], np.cumsum(_nl)[:-1]])
    _lat = np.degrees(np.arcsin(np.sort(np.polynomial.legendre.leggauss(2560)[0])[::-1]))  # north first

    @classmethod
    def index(cls, lat, lon):
        k = np.abs(cls._lat[None, :] - np.asarray(lat).reshape(-1, 1)).argmin(axis=1).reshape(np.shape(lat))
        nl = cls._nl[k]
        c = np.round((np.asarray(lon) % 360) / (360 / nl)).astype(int) % nl
        return cls._off[k] + c

class IFS(O1280):
    id, name, dlat, dlon, horizon, reach = 2, 'ecmwf_ifs', 0.1, 0.1, 240, 240
    fields = [1, 2, 3, 4, 5, 6, 10, 11, 12, 13]

    variables = ['wind_u_component_10m', 'wind_v_component_10m', 'wind_gusts_10m', 'pressure_msl',
                 'precipitation', 'cape', 'temperature_2m', 'ocean_u_current', 'ocean_v_current', 'visibility']

    @staticmethod
    def values(v):
        u, w = v['wind_u_component_10m'], v['wind_v_component_10m']
        cu, cv = v['ocean_u_current'], v['ocean_v_current']
        return {1: np.hypot(u, w), 2: (np.degrees(np.arctan2(-u, -w)) + 360) % 360, 3: v['wind_gusts_10m'],
                4: v['pressure_msl'] / 100, 5: v['precipitation'], 6: v['cape'], 10: v['temperature_2m'],   # Pa -> hPa
                11: np.hypot(cu, cv), 12: (np.degrees(np.arctan2(cu, cv)) + 360) % 360, 13: v['visibility']}

class WAM(O1280):
    id, name, dlat, dlon, horizon, reach = 3, 'ecmwf_wam', 0.1, 0.1, 240, 240
    fields = [7, 8, 9]

    variables = ['wave_height', 'wave_direction', 'wave_period']

    @staticmethod
    def values(v):
        # Under sea ice WAM gives a wave height of 0 with a period of 20-40 s: no waves have no period, so the period
        # is "no value" wherever the height rounds to 0 m (better no data than wrong data). The direction stays: the
        # app draws a height with its direction, and 0 m in a calm fjord must still read as 0 m, not as a hole.
        h = np.asarray(v['wave_height'], np.float64)
        calm = np.isfinite(h) & (h < 0.05)
        return {7: v['wave_height'], 8: v['wave_direction'], 9: np.where(calm, np.nan, v['wave_period'])}

MODELS = [Nordic, IFS, WAM]

# MARK: quantisation

def quantise(field, a):
    a = np.asarray(a, np.float64)
    bad = ~np.isfinite(a)
    a = np.nan_to_num(a)
    if field in (1, 3):
        q = np.clip(np.round(a / 0.5), 0, 254)
    elif field in (2, 8, 12):
        q = np.round(a * 255 / 360) % 255                 # 0..254; 255 is reserved, 360 folds onto 0
    elif field == 4:
        q = np.clip(np.round(a - 900), 0, 254)
    elif field == 5:
        mm = np.maximum(a, 0)
        q = np.where(mm <= 10, np.round(mm * 10), np.round(100 + 15 * np.sqrt(np.maximum(mm - 10, 0)))).clip(0, 254)
    elif field == 6:
        q = np.clip(np.round(a / 20), 0, 254)
    elif field in (7, 9):
        q = np.clip(np.round(a * 10), 0, 254)
    elif field == 10:
        q = np.clip(np.round((a + 40) * 2), 0, 254)
    elif field == 11:
        q = np.clip(np.round(a / 0.02), 0, 254)
    elif field == 13:
        m = np.maximum(a, 0)
        q = np.where(m <= 10000, np.round(m / 100), np.round(100 + (m - 10000) / 1000)).clip(0, 254)
    else:
        raise ValueError(field)
    return np.where(bad, NONE, q).astype(np.uint8)

def dequantise(field, q):
    q = np.asarray(q, np.float64)
    if field in (1, 3):
        v = q * 0.5
    elif field in (2, 8, 12):
        v = q * 360 / 255
    elif field == 4:
        v = q + 900
    elif field == 5:
        v = np.where(q <= 100, q / 10, 10 + ((q - 100) / 15) ** 2)
    elif field == 6:
        v = q * 20
    elif field == 10:
        v = q / 2 - 40
    elif field == 11:
        v = q * 0.02
    elif field == 13:
        v = np.where(q <= 100, q * 100, 10000 + (q - 100) * 1000)
    else:
        v = q / 10
    return np.where(q == NONE, np.nan, v)

# MARK: runs

def latest_run(model):
    """The newest completed run in data_spatial that reaches model.reach hours:
    (prefix, reference datetime, [valid time strings])."""
    now = datetime.now(timezone.utc)
    for back in range(3):
        day = now - timedelta(days=back)
        prefix = f'openmeteo/data_spatial/{model.name}/{day:%Y/%m/%d}/'
        try:
            runs = sorted(p for p in FS.ls(prefix, refresh=True) if p.endswith('Z'))
        except FileNotFoundError:
            continue
        for run in reversed(runs):
            try:
                meta = json.load(urllib.request.urlopen(f'{HTTP}/{run.split("/", 1)[1]}/meta.json', timeout=30))
            except Exception:
                continue
            if meta.get('completed'):
                ref = datetime.strptime(meta['reference_time'], '%Y-%m-%dT%H:%M:%SZ').replace(tzinfo=timezone.utc)
                last = datetime.strptime(meta['valid_times'][-1], '%Y-%m-%dT%H:%MZ').replace(tzinfo=timezone.utc)
                if (last - ref).total_seconds() / 3600 < model.reach:
                    continue        # ECMWF's 06/18 UTC runs stop at 144 h: take the last 00/12 UTC run instead
                return run, ref, meta['valid_times']
    raise SystemExit(f'{model.name}: no completed run in the last 3 days')

def run_at(model, ref_ts):
    """The run whose reference time is ref_ts (Unix seconds): (prefix, reference datetime, [valid time strings]).
    check() compares against the run a file was built from, so a newer run published meanwhile does not matter."""
    ref = datetime.fromtimestamp(ref_ts, timezone.utc)
    prefix = f'openmeteo/data_spatial/{model.name}/{ref:%Y/%m/%d}/'
    for run in sorted(p for p in FS.ls(prefix, refresh=True) if p.endswith('Z')):
        meta = json.load(urllib.request.urlopen(f'{HTTP}/{run.split("/", 1)[1]}/meta.json', timeout=30))
        if meta.get('reference_time') == f'{ref:%Y-%m-%dT%H:%M:%SZ}':
            return run, ref, meta['valid_times']
    raise SystemExit(f'{model.name}: run {ref:%Y-%m-%d %H}Z is gone from the bucket')

def pick_steps(ref, valid, horizon):
    """Every step the run has, up to the horizon: ECMWF is hourly to 90 h, 3-hourly to 144 h, then 6-hourly."""
    out = []
    for t in valid:
        h = round((datetime.strptime(t, '%Y-%m-%dT%H:%MZ').replace(tzinfo=timezone.utc) - ref).total_seconds() / 3600)
        if h > horizon:
            break
        out.append((h, t.replace(':', '').replace('Z', '')))
    return out

async def _read_step(fs, path, names, sel, shape):
    r = await omfiles.OmFileReaderAsync.from_fsspec(fs, path)

    async def one(name):
        try:
            child = r.get_child_by_name(name)
            if asyncio.iscoroutine(child):
                child = await child
        except ValueError:                           # the step lacks it (gusts at +0 h): byte 255
            return name, np.full(shape, np.nan, np.float32)
        return name, np.asarray(await child.read_array(sel)).reshape(shape)
    return dict(await asyncio.gather(*(one(n) for n in names)))

async def _read_steps(run, stamps, names, sel, shape, reduce, parallel, tries):
    fs = s3fs.S3FileSystem(anon=True, asynchronous=True, skip_instance_cache=True)
    client = await fs.set_session()
    gate = asyncio.Semaphore(parallel)

    async def step(i, stamp):
        async with gate:
            for k in range(tries):
                try:
                    return reduce(await _read_step(fs, f'{run}/{stamp}.om', names, sel, shape), i)
                except (Exception, asyncio.CancelledError) as e:  # omfiles surfaces dropped requests as cancellations
                    if k == tries - 1:
                        raise RuntimeError(f'{run}/{stamp}.om: {e!r}') from e
                    await asyncio.sleep(2 ** k)
    try:
        return await asyncio.gather(*(step(i, t) for i, t in enumerate(stamps)))
    finally:
        await client.close()

def fetch_steps(run, stamps, names, sel, shape, reduce=lambda v, i: v, parallel=8, tries=5):
    """reduce(variables of a step over sel, step index) for every step, in order; read concurrently on one loop."""
    return asyncio.run(_read_steps(run, stamps, names, sel, shape, reduce, parallel, tries))

def fetch_steps_local(run, stamps, names, sel, shape, reduce, parallel=8, tries=5):
    """As fetch_steps, but each step file is downloaded whole and read from local disk (ADR-098: range reads cost
    one round trip each to us-west-2, so a run of ~125 steps x 10 variables is latency-bound; whole files are a few
    large transfers). Files are deleted as soon as they are read; at most `parallel` are on disk at once."""
    import shutil, subprocess, tempfile
    from concurrent.futures import ThreadPoolExecutor
    scratch = tempfile.mkdtemp(prefix='om-')

    def one(item):
        i, stamp = item
        path = f'{scratch}/{i}.om'
        url = f"{HTTP}/{run.split('/', 1)[1]}/{stamp}.om"
        for k in range(tries):
            try:
                subprocess.run(['curl', '-fsS', '--retry', '3', '-o', path, url], check=True, timeout=900)
                reader = omfiles.OmFileReader.from_path(path)
                values = {}
                for name in names:
                    try:
                        child = reader.get_child_by_name(name)
                    except ValueError:                   # the step lacks it (gusts at +0 h): byte 255
                        values[name] = np.full(shape, np.nan, np.float32)
                        continue
                    values[name] = np.asarray(child.read_array(sel)).reshape(shape)
                reader.close()
                return reduce(values, i)
            except Exception as e:
                if k == tries - 1:
                    raise RuntimeError(f'{url}: {e!r}') from e
                time.sleep(2 ** k)
            finally:
                if os.path.exists(path):
                    os.remove(path)
    try:
        with ThreadPoolExecutor(parallel) as pool:
            return list(pool.map(one, enumerate(stamps)))
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

def fetch(run, stamp, names, sel, shape):
    return fetch_steps(run, [stamp], names, sel, shape, parallel=1)[0]

# MARK: build

def unit_name(unit):
    lat, lon, size = unit
    return f"{'c' if size == 1 else 'o'}_{lat}_{lon}"

def unit_grid(model, unit):
    """rows, cols and the point coordinates (row 0 = south, both edges included) of a unit for a model."""
    lat, lon, size = unit
    dlat, dlon = (model.dlat, model.dlon) if size == 1 else (TILE_D, TILE_D)
    rows, cols = int(round(size / dlat)) + 1, int(round(size / dlon)) + 1
    la = lat + dlat * np.arange(rows)
    lo = lon + dlon * np.arange(cols)
    return rows, cols, dlat, dlon, np.repeat(la, cols), np.tile(lo, rows)

def build(outdir, units):
    """One FYRF per unit into OUTDIR/<run id>/, then OUTDIR/index.json. Returns the run id."""
    sections = {u: [] for u in units}
    runs = {}
    for model in MODELS:
        run, ref, valid = latest_run(model)
        steps = pick_steps(ref, valid, model.horizon)
        runs[model.name] = ref
        t0 = time.time()
        used = [u for u in units if model is not Nordic or u[2] == 1]
        grids = {u: unit_grid(model, u) for u in used}
        if model is Nordic:
            idx = {u: Nordic.index(g[4], g[5]) for u, g in grids.items()}
            used = [u for u in used if (idx[u][0] >= 0).mean() >= 0.5]       # cells mostly inside the domain
            if not used:
                print(f'{model.name}: no cell inside the domain'); continue
            J = np.concatenate([idx[u][0] for u in used]); I = np.concatenate([idx[u][1] for u in used])
            inside = J >= 0
            j0, j1, i0, i1 = J[inside].min(), J[inside].max() + 1, I[inside].min(), I[inside].max() + 1
            names, sel, shape = Nordic.variables, (slice(j0, j1), slice(i0, i1)), (j1 - j0, i1 - i0)
            jj, ii = np.maximum(J, 0) - j0, np.maximum(I, 0) - i0

            def pick(w, f):
                return np.where(inside, quantise(f, w[f][jj, ii]), NONE)
        else:
            K = np.concatenate([model.index(grids[u][4], grids[u][5]) for u in used])
            a0, a1 = int(K.min()), int(K.max()) + 1
            names, sel, shape = model.variables, (slice(0, 1), slice(a0, a1)), (a1 - a0,)
            kk = K - a0

            def pick(w, f):
                return quantise(f, w[f][kk])
        npts = sum(grids[u][0] * grids[u][1] for u in used)
        out = {f: np.empty((npts, len(steps)), np.uint8) for f in model.fields}   # all points x steps, per field

        done = [0]

        def reduce(v, i):
            w = model.values(v)
            for f in model.fields:
                out[f][:, i] = pick(w, f)
            done[0] += 1
            if done[0] % 12 == 0 or done[0] == len(steps):
                print(f'  {model.name}: {done[0]}/{len(steps)} steps, {time.time() - t0:.0f} s', flush=True)

        fetch_steps_local(run, [t for _, t in steps], names, sel, shape, reduce)
        at = 0
        for u in used:
            rows, cols, dlat, dlon = grids[u][:4]
            n = rows * cols
            fields = [(f, out[f][at:at + n]) for f in model.fields]
            at += n
            if model is WAM and all((q == NONE).all() for _, q in fields):
                continue                                                     # all land
            sections[u].append((model, ref, steps, rows, cols, dlat, dlon, fields))
        print(f'{model.name}: run {ref:%Y-%m-%d %H}Z, {len(steps)} steps (to +{steps[-1][0]} h), {len(used)} units, '
              f'{npts} points, {time.time() - t0:.0f} s', flush=True)
    run_id = runs['ecmwf_ifs'].strftime('%Y%m%d%H')
    folder = f'{outdir}/{run_id}'
    os.makedirs(folder, exist_ok=True)
    sizes, written = [], {'cells': [], 'tiles': []}
    for u, secs in sections.items():
        if not secs or (u[2] == TILE and not any(m is WAM for m, *_ in secs)):
            continue                                                         # nothing, or a land-only tile
        z = deflate(encode(u, secs))
        path = f'{folder}/{unit_name(u)}.fyrf.z'
        with open(path + '.tmp', 'wb') as fh:
            fh.write(z)
        os.replace(path + '.tmp', path)
        sizes.append(len(z))
        written['cells' if u[2] == 1 else 'tiles'].append(f'{u[0]}_{u[1]}')
    index = {'format': 'FYRF', 'version': VERSION, 'run': run_id,
             'generatedAt': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
             'models': {k: v.strftime('%Y-%m-%dT%H:%M:%SZ') for k, v in runs.items()},
             'tileSize': TILE, 'cells': sorted(written['cells']), 'tiles': sorted(written['tiles']),
             'bytes': int(sum(sizes)),
             'attribution': 'Open-Meteo.com (CC BY 4.0): ECMWF IFS and WAM, MET Norway MET Nordic'}
    with open(f'{outdir}/index.json.tmp', 'w') as fh:
        json.dump(index, fh, separators=(',', ':'))
    os.replace(f'{outdir}/index.json.tmp', f'{outdir}/index.json')
    print(f"{len(written['cells'])} cells + {len(written['tiles'])} tiles, {sum(sizes) / 1e6:.1f} MB, "
          f'median {np.median(sizes) / 1e3:.0f} kB, max {max(sizes) / 1e3:.0f} kB, run {run_id}')
    return run_id

HEAD, SECTION = '<4sBhhBB', '<BIHHHddddB'

def deflate(raw):
    """Raw DEFLATE (wbits -15), like every other .z pack (fyr_formats.py, ChartCore Deflate.swift)."""
    c = zlib.compressobj(9, zlib.DEFLATED, -15)
    return c.compress(raw) + c.flush()

def inflate(z):
    return zlib.decompressobj(-15).decompress(z)

def encode(unit, secs):
    lat, lon, size = unit
    out = [struct.pack(HEAD, MAGIC, VERSION, lat, lon, size, len(secs))]
    body = []
    for model, ref, steps, rows, cols, dlat, dlon, fields in secs:
        out.append(struct.pack(SECTION, model.id, int(ref.timestamp()), len(steps), rows, cols,
                               float(lat), float(lon), dlat, dlon, len(fields)))
        out.append(bytes(f for f, _ in fields))
        out.append(struct.pack(f'<{len(steps)}H', *(h for h, _ in steps)))
        for _, q in fields:
            body.append((np.diff(q.astype(np.int16), axis=1, prepend=0) % 256).astype(np.uint8).tobytes())
    return b''.join(out) + b''.join(body)

def decode(data):
    magic, version, lat, lon, size, n = struct.unpack_from(HEAD, data, 0)
    assert magic == MAGIC and version == VERSION, (magic, version)
    pos, heads = struct.calcsize(HEAD), []
    for _ in range(n):
        mid, ref, ns, rows, cols, la0, lo0, dla, dlo, nf = struct.unpack_from(SECTION, data, pos)
        pos += struct.calcsize(SECTION)
        fields = list(data[pos:pos + nf]); pos += nf
        hours = list(struct.unpack_from(f'<{ns}H', data, pos)); pos += 2 * ns
        heads.append(dict(model=mid, ref=ref, rows=rows, cols=cols, lat0=la0, lon0=lo0, dlat=dla, dlon=dlo,
                          fields=fields, hours=hours))
    for h in heads:
        npts, ns = h['rows'] * h['cols'], len(h['hours'])
        h['data'] = {}
        for f in h['fields']:
            d = np.frombuffer(data, np.uint8, npts * ns, pos).reshape(npts, ns); pos += npts * ns
            h['data'][f] = (np.cumsum(d.astype(np.int64), axis=1) % 256).astype(np.uint8)
    assert pos == len(data), (pos, len(data))
    return (lat, lon, size), heads

# MARK: check

def check(outdir, n):
    """Re-read random values from the bucket by point and compare with the decoded cells."""
    run_id = json.load(open(f'{outdir}/index.json'))['run']
    outdir = f'{outdir}/{run_id}'
    files = sorted(f for f in os.listdir(outdir) if f.endswith('.fyrf.z'))
    by_model = {m.id: m for m in MODELS}
    runs = {}
    bad = total = 0
    rng = random.Random(7)
    for name in rng.sample(files, min(n, len(files))):
        _, heads = decode(inflate(open(f'{outdir}/{name}', 'rb').read()))
        for h in heads:
            model = by_model[h['model']]
            if (model.id, h['ref']) not in runs:
                runs[model.id, h['ref']] = run_at(model, h['ref'])
            run, ref, valid = runs[model.id, h['ref']]
            steps = pick_steps(ref, valid, model.horizon)
            p = rng.randrange(h['rows'] * h['cols']); s = rng.randrange(len(steps))
            plat, plon = h['lat0'] + (p // h['cols']) * h['dlat'], h['lon0'] + (p % h['cols']) * h['dlon']
            if model is Nordic:
                j, i = Nordic.index(np.array([plat]), np.array([plon]))
                if j[0] < 0:
                    continue
                sel = (slice(int(j[0]), int(j[0]) + 1), slice(int(i[0]), int(i[0]) + 1))
                vals = Nordic.values(fetch(run, steps[s][1], Nordic.variables, sel, (1,)))
            else:
                k = int(model.index(np.array([plat]), np.array([plon]))[0])
                vals = model.values(fetch(run, steps[s][1], model.variables, (slice(0, 1), slice(k, k + 1)), (1,)))
            for f in h['fields']:
                want = int(quantise(f, np.atleast_1d(vals[f]))[0])
                got = int(h['data'][f][p, s])
                total += 1
                if want != got:
                    bad += 1
                    print(f'MISMATCH {name} {model.name} field {f} at {plat:.3f},{plon:.3f} +{steps[s][0]} h: {got} != {want}')
    # Units: a field pinned at the edge of its byte range is a unit error (IFS pressure came in Pa), not weather.
    for name in files:
        _, heads = decode(inflate(open(f'{outdir}/{name}', 'rb').read()))
        for h in heads:
            for f, q in h['data'].items():
                valid = q[q != NONE]
                if valid.size and (valid == 254).mean() > 0.2:
                    bad += 1
                    print(f'CLIPPED {name} model {h["model"]} field {f}: {100 * (valid == 254).mean():.0f}% at the top of the scale')
                if valid.size and f in (4, 10) and (valid == 0).mean() > 0.2:
                    bad += 1
                    print(f'CLIPPED {name} model {h["model"]} field {f}: {100 * (valid == 0).mean():.0f}% at the bottom of the scale')
    print(f'checked {total} values, {bad} mismatches')
    if bad == 0 and total > 0:
        print('FYR_GATE_FORECAST_CELLS_OK')

def dump(path):
    (lat, lon, size), heads = decode(inflate(open(path, 'rb').read()))
    print(f'unit {lat}_{lon} size {size}')
    for h in heads:
        name = next(m.name for m in MODELS if m.id == h['model'])
        ref = datetime.fromtimestamp(h['ref'], timezone.utc)
        print(f'  {name}: run {ref:%Y-%m-%d %H}Z, {len(h["hours"])} steps +{h["hours"][0]}..+{h["hours"][-1]} h, '
              f'{h["rows"]} x {h["cols"]} at {h["dlat"]} x {h["dlon"]}, fields {h["fields"]}')

def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    b = sub.add_parser('build'); b.add_argument('outdir'); b.add_argument('--cells'); b.add_argument('--box', nargs=4, type=int)
    b.add_argument('--tiles', choices=['all']); b.add_argument('--tiles-box', nargs=4, type=int)
    c = sub.add_parser('check'); c.add_argument('outdir'); c.add_argument('n', nargs='?', type=int, default=40)
    d = sub.add_parser('dump'); d.add_argument('file')
    sub.add_parser('latest')
    a = ap.parse_args()
    if a.cmd == 'build':
        units = []
        if a.cells:
            units += [tuple(int(x) for x in line.split('_')) + (1,) for line in open(a.cells) if line.strip()]
        elif a.box:
            la0, la1, lo0, lo1 = a.box
            units += [(la, lo, 1) for la in range(la0, la1) for lo in range(lo0, lo1)]
        if a.tiles or a.tiles_box:
            la0, la1, lo0, lo1 = a.tiles_box or (-80, 85, -180, 180)     # WindGrid.latitudeLimit is 85
            units += [(la, lo, TILE) for la in range(la0 - la0 % TILE, la1, TILE) for lo in range(lo0 - lo0 % TILE, lo1, TILE)]
        if not units:
            sys.exit('--cells/--box and/or --tiles/--tiles-box')
        build(a.outdir, units)
    elif a.cmd == 'check':
        check(a.outdir, a.n)
    elif a.cmd == 'latest':
        print(latest_run(IFS)[1].strftime('%Y%m%d%H'))
    else:
        dump(a.file)

if __name__ == '__main__':
    main()

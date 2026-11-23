"""Is what the Fyr app reads there and fresh? Run by .github/workflows/watch.yml; a failed run is the alarm.

The forecast is built from the 00 and 12 UTC runs and published some 8 h after each, so it is 20 h old at most
when all is well. 30 h means a build has been missed (the gate stopped it, the source was late, or the schedule
did not fire). Standard library only.

  python watch.py            checks the live addresses, exit 1 with one line per fault
"""
import json
import sys
import urllib.request
from datetime import datetime, timezone

STALE_HOURS = 30
INDEX = 'https://data.seaqio.com/forecast/v1/index.json'
REACHABLE = [
    'https://data.seaqio.com/v1/manifest.json',
    'https://data.seaqio.com/v1/manifest.json.sig',
    'https://seaqio.com/fyr/',
    'https://seaqio.com/fyr/privacy.html',
]


def forecast_fault(index, now):
    """None when the published forecast is fresh, otherwise what is wrong with it."""
    try:
        made = datetime.fromisoformat(index['generatedAt'].replace('Z', '+00:00'))
    except (KeyError, TypeError, ValueError, AttributeError):
        return 'forecast index has no readable generatedAt'
    age = (now - made).total_seconds() / 3600
    if age >= STALE_HOURS:
        return f"forecast is {age:.0f} h old (run {index.get('run')}, published {index['generatedAt']})"
    if not index.get('cells'):
        return f"forecast run {index.get('run')} lists no cells"
    return None


def fetch(url, limit=None):
    request = urllib.request.Request(url, headers={'User-Agent': 'fyr-watch'})
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read(limit) if limit else response.read()


def main():
    faults = []
    try:
        fault = forecast_fault(json.loads(fetch(INDEX)), datetime.now(timezone.utc))
        if fault:
            faults.append(fault)
    except Exception as error:  # noqa: BLE001 - any failure to read it is the fault
        faults.append(f'{INDEX}: {error}')
    for url in REACHABLE:
        try:
            if not fetch(url, 64):
                faults.append(f'{url}: empty')
        except Exception as error:  # noqa: BLE001
            faults.append(f'{url}: {error}')
    for fault in faults:
        print(f'FYR_WATCH_FAULT {fault}')
    if faults:
        return 1
    print('FYR_WATCH_OK')
    return 0


def selftest():
    now = datetime(2026, 11, 23, 0, 0, tzinfo=timezone.utc)
    assert forecast_fault({'run': 'r', 'generatedAt': '2026-11-22T19:46:08Z', 'cells': ['a']}, now) is None
    assert 'h old' in forecast_fault({'run': 'r', 'generatedAt': '2026-11-21T17:59:00Z', 'cells': ['a']}, now)
    assert forecast_fault({'run': 'r', 'generatedAt': '2026-11-21T18:01:00Z', 'cells': ['a']}, now) is None
    assert 'no cells' in forecast_fault({'run': 'r', 'generatedAt': '2026-11-22T19:46:08Z', 'cells': []}, now)
    assert 'generatedAt' in forecast_fault({}, now)
    assert 'generatedAt' in forecast_fault({'generatedAt': 7}, now)
    print('FYR_WATCH_SELFTEST_OK')


if __name__ == '__main__':
    if sys.argv[1:] == ['selftest']:
        selftest()
    else:
        sys.exit(main())

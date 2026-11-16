#!/usr/bin/env bash
# upload.sh OUTDIR: the run folder first (immutable), then index.json (short cache), then every older run goes:
# the host holds only the newest forecast (the app re-reads the index if a run vanishes under it).
# Needs AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY (an R2 token scoped to the bucket), FYR_S3_ENDPOINT, FYR_BUCKET,
# FYR_PREFIX (for example forecast/v1).
set -euo pipefail
OUT="$1"
RUN=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['run'])" "$OUT/index.json")
# R2 answers a burst with 503 "Reduce your concurrent request rate" now and then (4 of ~8,000 files on 2026-11-16):
# the CLI backs off and retries each request, and sync is re-run so a file that still failed goes up on the next pass.
export AWS_RETRY_MODE=adaptive AWS_MAX_ATTEMPTS=10
S3="aws s3 --endpoint-url $FYR_S3_ENDPOINT --region auto"
DEST="s3://$FYR_BUCKET/$FYR_PREFIX"
for pass in 1 2 3; do
  if $S3 sync --only-show-errors --size-only "$OUT/$RUN" "$DEST/$RUN" \
      --content-type application/octet-stream --cache-control "public, max-age=31536000, immutable"; then
    break
  fi
  [ "$pass" = 3 ] && { echo "upload failed after 3 passes; index.json not published" >&2; exit 1; }
  sleep 30
done
$S3 cp --only-show-errors "$OUT/index.json" "$DEST/index.json" \
  --content-type application/json --cache-control "public, max-age=300"
for old in $($S3 ls "$DEST/" | awk '/PRE [0-9]{10}\// {print $2}' | grep -v "^$RUN/$"); do
  $S3 rm --recursive --only-show-errors "$DEST/$old"
done
echo "published run $RUN"

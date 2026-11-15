#!/usr/bin/env bash
# upload.sh OUTDIR: the run folder first (immutable), then index.json (short cache), then prune to the newest two runs.
# Needs AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY (an R2 token scoped to the bucket), FYR_S3_ENDPOINT, FYR_BUCKET,
# FYR_PREFIX (for example forecast/v1).
set -euo pipefail
OUT="$1"
RUN=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['run'])" "$OUT/index.json")
S3="aws s3 --endpoint-url $FYR_S3_ENDPOINT --region auto"
DEST="s3://$FYR_BUCKET/$FYR_PREFIX"
$S3 cp --recursive --only-show-errors "$OUT/$RUN" "$DEST/$RUN" \
  --content-type application/octet-stream --cache-control "public, max-age=31536000, immutable"
$S3 cp --only-show-errors "$OUT/index.json" "$DEST/index.json" \
  --content-type application/json --cache-control "public, max-age=300"
for old in $($S3 ls "$DEST/" | awk '/PRE [0-9]{10}\// {print $2}' | sort | head -n -2); do
  $S3 rm --recursive --only-show-errors "$DEST/$old"
done
echo "published run $RUN"

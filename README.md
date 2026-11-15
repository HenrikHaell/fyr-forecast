# fyr-forecast

Builds the weather forecast for the [Fyr](https://seaqio.com/fyr/) boating app twice a day from Open-Meteo's open
data, and uploads it for the app.

- `om_forecast.py` reads the newest complete ECMWF IFS / WAM and MET Nordic runs from `s3://openmeteo`
  (range reads, only the variables used) and writes one FYRF v2 file per unit: the chart's coastal 1° cells
  (`cells.txt`) and 5° ocean tiles, then `index.json`. The format is described at the top of the script.
- `.github/workflows/forecast.yml` runs every three hours, builds only when a newer IFS run exists than the one
  published, checks random values against the source (`check`), uploads the run folder, then `index.json`, then
  deletes every older run (only the newest forecast is hosted).
- `upload.sh` talks to the storage over the S3 API. Credentials live in repository secrets, never in the code.

Code: GPL-2.0 (see `LICENSE`; `omfiles` is GPL-2.0). Data: see `NOTICE-DATA.md`.

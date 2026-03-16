# Satellite Endpoint

This folder now contains only the active TSV + navigation RINEX endpoint.

## Files

- `prepare_tsv_nav_inputs.m`: discovers `test.tsv` and the navigation RINEX file under `endpoint_satellite/data`
- `fill_tsv_with_rinex_nav.m`: fills missing satellite position, velocity, clock bias, and clock drift fields
- `run_tsv_nav_endpoint.bat`: Windows launcher
- `run_tsv_nav_endpoint_batch.m`: MATLAB batch entrypoint
- `run_tsv_nav_current_test.m`: simple MATLAB runner for local testing

## Notes

- Existing goGPS engine files under `source/` are untouched.
- This endpoint uses `test.tsv + navigation RINEX` only.
- Output is written to `endpoint_satellite/tsv_nav_output/`.

# Satellite Endpoint

This folder contains only the new endpoint files added for the satellite-position export flow.

## Files

- `main_satellite_endpoint.m`: simple local launcher with editable paths
- `process_and_export_rinex_files_endpoint.m`: endpoint pipeline entry
- `export_processed_mat_to_tsv.m`: MAT to TSV exporter with satellite position and velocity fields
- `run_satellite_endpoint.bat`: Windows batch launcher
- `prepare_endpoint_inputs.m`: auto-discovers observation and nav files under `endpoint_satellite/data`
- `run_rinex_only_endpoint.bat`: RINEX-only launcher that ignores `test.tsv`
- `prepare_rinex_only_inputs.m`: RINEX-only discovery helper

## Notes

- Existing goGPS engine files under `source/` are untouched.
- Existing processing entrypoints such as `main.m` and `process_and_save_rinex_files.m` are reused, not replaced.
- Output is written under the `output_dir` you pass in.
- If no arguments are passed to `run_satellite_endpoint.bat`, it searches `endpoint_satellite/data` automatically.
- If you want to ignore `test.tsv` and use only observation/nav RINEX files, run `run_rinex_only_endpoint.bat`.

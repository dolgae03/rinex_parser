# Satellite Endpoint

This folder now contains only the active TSV + navigation RINEX endpoint.

## Folder Layout

- `data/`: input TSV files only
- `nav/`: user-provided local navigation RINEX files
- `nav_cache/`: automatically downloaded mixed navigation files
- `output/`: filled TSV outputs

## Files

- `prepare_tsv_nav_inputs.m`: discovers `data/*.tsv`, reuses a matching local navigation RINEX from `nav/`, or triggers automatic download
- `resolve_or_download_mixed_nav.m`: infers the TSV day range and downloads cached mixed navigation RINEX files when needed
- `fill_tsv_with_rinex_nav.m`: fills missing satellite position, velocity, clock bias, and clock drift fields
- `run_tsv_nav_endpoint.bat`: Windows launcher
- `run_tsv_nav_endpoint_batch.m`: MATLAB batch entrypoint

## Notes

- Existing goGPS engine files under `source/` are untouched.
- This endpoint uses `test.tsv + navigation RINEX` only.
- Put TSV inputs in `endpoint_satellite/data/`.
- Put local navigation RINEX files in `endpoint_satellite/nav/`.
- If no matching navigation file is found there, the wrapper reads the TSV time span and downloads mixed navigation files into `endpoint_satellite/nav_cache/`.
- Downloaded mixed navigation files are normalized to a goGPS-compatible filename ending in `p` so the legacy parser treats them as mixed nav.
- Mixed navigation files are also filtered down to only the constellation systems actually requested by the TSV.
- SBAS and IRNSS rows are preserved in the TSV, but the current wrapper skips them during navigation loading, so those rows remain unfilled.
- Output is written to `endpoint_satellite/output/` with the suffix `_with_sv_pos.tsv`.

# FastAPI MATLAB Bridge

This service wraps the `endpoint_satellite` MATLAB TSV endpoint behind FastAPI.

## What It Does

- Accepts a `.tsv` upload
- Or accepts a smartphone `.txt` GNSS log and converts it to TSV through the local `gnss_txt_parser`
- Connects to a MATLAB engine
- Calls `endpoint_satellite/run_tsv_nav_endpoint_shared.m`
- Returns a ZIP package as the HTTP response

The main `/process` endpoint does not ask for a navigation file.
It relies on the MATLAB wrapper to reuse or automatically download the required mixed navigation RINEX.
The returned package name is based on the uploaded input name and the processing timestamp, for example
`my_log_20260317_153012_package.zip`.

Each ZIP contains:

- `summary.json`: processing metadata, row counts, time span, constellation list, and MATLAB summary
- `*_original.*`: the uploaded source file
- `*_measurements.tsv`: the measurement-only TSV used as MATLAB input
- `*_with_sv_pos.tsv`: the final filled TSV file

## 1. MATLAB Engine

The API tries to connect to an existing shared MATLAB engine first.
If none is available, it attempts to auto-start a MATLAB engine through the Python MATLAB Engine package.

You can still start a named shared engine manually if you prefer:

```matlab
matlab.engine.shareEngine('gnss_endpoint')
```

Then set:

```powershell
$env:MATLAB_SHARED_ENGINE_NAME = "gnss_endpoint"
```

## 2. Install Python Dependencies

```powershell
pip install -r api_server\requirements.txt
```

Make sure `matlab.engine` is installed in the same Python environment.

## 3. Run the API

```powershell
uvicorn api_server.app:app --host 0.0.0.0 --port 8000
```

or on Windows:

```powershell
.\api_server\run_api_server.bat
```

By default, the batch file activates the `gnss-matlab` conda environment and asks the API to prewarm a MATLAB engine on startup.
You can pass a different conda environment name as the first argument:

```powershell
.\api_server\run_api_server.bat my-env-name
```

## 4. Example Request

TSV directly:

```powershell
curl.exe -X POST "http://127.0.0.1:8000/process" `
  -F "input_file=@C:\path\to\test.tsv"
```

TXT through the bundled smartphone parser:

```powershell
curl.exe -X POST "http://127.0.0.1:8000/process" `
  -F "input_file=@C:\path\to\raw.txt"
```

The `/process` response body is the generated ZIP package itself.

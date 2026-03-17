from __future__ import annotations

import json
import os
import threading
import uuid
from datetime import datetime
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, File, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse

try:
    import matlab.engine  # type: ignore
except ImportError:  # pragma: no cover
    matlab = None
else:
    matlab = matlab.engine

try:
    from gnss_txt_parser import convert_txt_to_tsv as convert_android_txt_to_tsv
except ImportError:  # pragma: no cover
    convert_android_txt_to_tsv = None


REPO_ROOT = Path(__file__).resolve().parents[1]
WORK_ROOT = REPO_ROOT / "api_server_work"

app = FastAPI(title="GNSS MATLAB Endpoint API", version="0.1.0")

_engine_lock = threading.Lock()
_shared_engine = None
_shared_engine_name = None


def _ensure_workdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _connect_engine():
    global _shared_engine, _shared_engine_name

    if matlab is None:
        raise HTTPException(
            status_code=500,
            detail="matlab.engine for Python is not installed in this environment.",
        )

    if _shared_engine is not None:
        return _shared_engine, _shared_engine_name

    requested_name = os.getenv("MATLAB_SHARED_ENGINE_NAME", "").strip()
    try:
        if requested_name:
            engine = matlab.connect_matlab(requested_name)
            engine_name = requested_name
        else:
            engine = matlab.connect_matlab()
            available = list(matlab.find_matlab())
            engine_name = available[0] if available else "auto-started"
    except HTTPException:
        raise
    except Exception as exc:  # pragma: no cover
        raise HTTPException(status_code=503, detail=f"Failed to connect to MATLAB engine: {exc}") from exc

    engine.addpath(engine.genpath(str(REPO_ROOT)), nargout=0)
    _shared_engine = engine
    _shared_engine_name = engine_name
    return _shared_engine, _shared_engine_name


@app.on_event("startup")
def startup_connect_engine() -> None:
    if os.getenv("MATLAB_PREWARM", "1").strip() == "0":
        return
    _connect_engine()


def _save_upload(upload: UploadFile, destination: Path) -> Path:
    _ensure_workdir(destination.parent)
    with destination.open("wb") as fh:
        while True:
            chunk = upload.file.read(1024 * 1024)
            if not chunk:
                break
            fh.write(chunk)
    return destination


def _convert_txt_to_tsv(txt_path: Path, tsv_path: Path) -> Path:
    if convert_android_txt_to_tsv is None:
        raise HTTPException(
            status_code=500,
            detail="gnss_txt_parser is not available in this environment.",
        )

    try:
        summary = convert_android_txt_to_tsv(
            input_path=txt_path,
            output_path=tsv_path,
            filter_valid_pr=True,
        )
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail=f"TXT conversion with gnss_txt_parser failed: {exc}",
        ) from exc

    resolved_output = Path(summary.output_path) if getattr(summary, "output_path", None) else tsv_path
    if not resolved_output.exists():
        raise HTTPException(
            status_code=500,
            detail=f"TXT conversion completed but no TSV output was created: {resolved_output}",
        )
    return resolved_output


def _run_matlab(tsv_path: Path, nav_path: Optional[Path], output_path: Path) -> dict:
    engine, engine_name = _connect_engine()

    nav_arg = str(nav_path) if nav_path else ""
    with _engine_lock:
        try:
            summary_json = engine.run_tsv_nav_endpoint_shared(
                str(tsv_path),
                nav_arg,
                str(output_path),
                nargout=1,
            )
        except Exception as exc:  # pragma: no cover
            raise HTTPException(status_code=500, detail=f"MATLAB execution failed: {exc}") from exc

    try:
        summary = json.loads(summary_json)
    except json.JSONDecodeError as exc:
        raise HTTPException(status_code=500, detail=f"MATLAB returned invalid JSON: {exc}") from exc

    summary["engine_name"] = engine_name
    return summary


def _process_impl(
    input_file: UploadFile,
    nav_file: Optional[UploadFile],
) -> dict:
    suffix = Path(input_file.filename or "").suffix.lower()
    if suffix not in {".txt", ".tsv"}:
        raise HTTPException(status_code=400, detail="input_file must be .txt or .tsv")

    job_id = uuid.uuid4().hex
    job_dir = _ensure_workdir(WORK_ROOT / job_id)
    uploads_dir = _ensure_workdir(job_dir / "uploads")
    output_dir = _ensure_workdir(job_dir / "output")

    input_path = _save_upload(input_file, uploads_dir / (input_file.filename or f"input{suffix}"))
    nav_path = None
    if nav_file is not None and nav_file.filename:
        nav_path = _save_upload(nav_file, uploads_dir / nav_file.filename)

    if suffix == ".txt":
        tsv_path = uploads_dir / f"{input_path.stem}_parsed.tsv"
        tsv_path = _convert_txt_to_tsv(input_path, tsv_path)
        input_kind = "txt"
    else:
        tsv_path = input_path
        input_kind = "tsv"

    output_filename = _build_output_filename(input_path.stem)
    output_path = output_dir / output_filename

    summary = _run_matlab(tsv_path, nav_path, output_path)

    return {
        "job_id": job_id,
        "input_kind": input_kind,
        "input_file": str(input_path),
        "tsv_file": str(tsv_path),
        "nav_file": str(nav_path) if nav_path else None,
        "output_tsv": str(output_path),
        "download_url": f"/outputs/{job_id}/{output_filename}",
        "summary": summary,
    }


def _build_output_response(result: dict) -> FileResponse:
    output_path = Path(result["output_tsv"])
    if not output_path.exists():
        raise HTTPException(
            status_code=500,
            detail=f"MATLAB processing finished but the output file was not created: {output_path}",
        )
    headers = {
        "X-GNSS-Job-Id": result["job_id"],
        "X-GNSS-Input-Kind": result["input_kind"],
    }
    return FileResponse(
        path=output_path,
        filename=output_path.name,
        media_type="text/tab-separated-values",
        headers=headers,
    )


def _build_output_filename(input_stem: str) -> str:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_stem = Path(input_stem).name
    return f"{safe_stem}_{timestamp}_with_sv_pos.tsv"


@app.get("/health")
def health() -> dict:
    available = []
    if matlab is not None:
        try:
            available = list(matlab.find_matlab())
        except Exception:
            available = []

    return {
        "status": "ok",
        "connected_engine": _shared_engine_name,
        "shared_engines": available,
        "txt_parser_available": convert_android_txt_to_tsv is not None,
    }


@app.post("/process")
def process_file(
    input_file: UploadFile = File(...),
) -> FileResponse:
    result = _process_impl(input_file, None)
    return _build_output_response(result)


@app.post("/process_compat")
async def process_file_compat(request: Request) -> dict:
    form = await request.form()

    input_file = form.get("input_file")
    nav_file = form.get("nav_file")
    if not isinstance(input_file, UploadFile):
        raise HTTPException(status_code=400, detail="input_file must be uploaded as a .txt or .tsv file")

    if nav_file == "":
        nav_file = None
    if nav_file is not None and not isinstance(nav_file, UploadFile):
        raise HTTPException(status_code=400, detail="nav_file must be omitted or uploaded as a file")
    result = _process_impl(input_file, nav_file)
    return _build_output_response(result)


@app.get("/outputs/{job_id}/{filename}")
def download_output(job_id: str, filename: str) -> FileResponse:
    target = WORK_ROOT / job_id / "output" / filename
    if not target.exists():
        raise HTTPException(status_code=404, detail="Output file not found")
    return FileResponse(path=target, filename=filename, media_type="text/tab-separated-values")

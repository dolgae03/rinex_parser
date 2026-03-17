from __future__ import annotations

import json
import os
import threading
import uuid
import zipfile
from datetime import datetime, timedelta
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
GPS_EPOCH = datetime(1980, 1, 6)

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


def _convert_txt_to_tsv(txt_path: Path, tsv_path: Path) -> tuple[Path, dict]:
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
    parser_summary = {
        "input_path": str(getattr(summary, "input_path", txt_path)),
        "output_path": str(getattr(summary, "output_path", resolved_output)),
        "total_rows": getattr(summary, "total_rows", None),
        "valid_rows": getattr(summary, "valid_rows", None),
    }
    return resolved_output, parser_summary


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

    txt_parser_summary = None
    if suffix == ".txt":
        tsv_path = uploads_dir / f"{input_path.stem}_parsed.tsv"
        tsv_path, txt_parser_summary = _convert_txt_to_tsv(input_path, tsv_path)
        input_kind = "txt"
    else:
        tsv_path = input_path
        input_kind = "tsv"

    output_filename = _build_output_filename(input_path.stem)
    output_path = output_dir / output_filename

    summary = _run_matlab(tsv_path, nav_path, output_path)
    metadata = _build_job_metadata(
        job_id=job_id,
        input_kind=input_kind,
        input_path=input_path,
        intermediate_tsv=tsv_path,
        output_tsv=output_path,
        nav_path=nav_path,
        matlab_summary=summary,
        txt_parser_summary=txt_parser_summary,
    )
    metadata_path = job_dir / "summary.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=True) + "\n", encoding="utf8")

    bundle_path = job_dir / _build_bundle_filename(input_path.stem)
    _write_bundle_zip(bundle_path, input_path, tsv_path, output_path, metadata_path)

    return {
        "job_id": job_id,
        "input_kind": input_kind,
        "input_file": str(input_path),
        "tsv_file": str(tsv_path),
        "nav_file": str(nav_path) if nav_path else None,
        "output_tsv": str(output_path),
        "bundle_zip": str(bundle_path),
        "download_url": f"/outputs/{job_id}/{bundle_path.name}",
        "summary": summary,
    }


def _build_output_response(result: dict) -> FileResponse:
    bundle_path = Path(result["bundle_zip"])
    if not bundle_path.exists():
        raise HTTPException(
            status_code=500,
            detail=f"Processing finished but the output bundle was not created: {bundle_path}",
        )
    headers = {
        "X-GNSS-Job-Id": result["job_id"],
        "X-GNSS-Input-Kind": result["input_kind"],
    }
    return FileResponse(
        path=bundle_path,
        filename=bundle_path.name,
        media_type="application/zip",
        headers=headers,
    )


def _build_output_filename(input_stem: str) -> str:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_stem = Path(input_stem).name
    return f"{safe_stem}_{timestamp}_with_sv_pos.tsv"


def _build_bundle_filename(input_stem: str) -> str:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_stem = Path(input_stem).name
    return f"{safe_stem}_{timestamp}_package.zip"


def _build_job_metadata(
    job_id: str,
    input_kind: str,
    input_path: Path,
    intermediate_tsv: Path,
    output_tsv: Path,
    nav_path: Optional[Path],
    matlab_summary: dict,
    txt_parser_summary: Optional[dict],
) -> dict:
    return {
        "job_id": job_id,
        "processed_at": datetime.now().isoformat(timespec="seconds"),
        "input_kind": input_kind,
        "input_file": str(input_path),
        "intermediate_tsv": str(intermediate_tsv),
        "final_output_tsv": str(output_tsv),
        "nav_file": str(nav_path) if nav_path else None,
        "txt_parser_summary": txt_parser_summary,
        "intermediate_summary": _summarize_tsv(intermediate_tsv),
        "final_summary": _summarize_tsv(output_tsv) if output_tsv.exists() else None,
        "matlab_summary": matlab_summary,
    }


def _summarize_tsv(tsv_path: Path) -> dict:
    import csv

    row_count = 0
    min_t_sec = None
    max_t_sec = None
    constellations = set()

    with tsv_path.open("r", encoding="utf8", newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        for row in reader:
            row_count += 1
            constellation = row.get("constellation", "").strip()
            if constellation != "":
                constellations.add(constellation)

            t_sec_text = row.get("t_sec", "").strip()
            try:
                t_sec = int(float(t_sec_text))
            except ValueError:
                continue

            min_t_sec = t_sec if min_t_sec is None else min(min_t_sec, t_sec)
            max_t_sec = t_sec if max_t_sec is None else max(max_t_sec, t_sec)

    return {
        "path": str(tsv_path),
        "total_rows": row_count,
        "time_range": {
            "t_sec_start": min_t_sec,
            "t_sec_end": max_t_sec,
            "utc_start": _gps_tsec_to_iso(min_t_sec),
            "utc_end": _gps_tsec_to_iso(max_t_sec),
        },
        "constellations": sorted(constellations),
    }


def _gps_tsec_to_iso(t_sec: Optional[int]) -> Optional[str]:
    if t_sec is None:
        return None
    return (GPS_EPOCH + timedelta(seconds=int(t_sec))).isoformat(timespec="seconds")


def _write_bundle_zip(
    bundle_path: Path,
    input_path: Path,
    intermediate_tsv: Path,
    output_tsv: Path,
    metadata_path: Path,
) -> None:
    intermediate_name = _bundle_member_name(intermediate_tsv.stem, intermediate_tsv.suffix, "measurements")
    output_name = _bundle_member_name(output_tsv.stem, output_tsv.suffix, "with_sv_pos")
    original_name = _bundle_member_name(input_path.stem, input_path.suffix, "original")

    with zipfile.ZipFile(bundle_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.write(metadata_path, arcname="summary.json")
        zf.write(input_path, arcname=original_name)
        zf.write(intermediate_tsv, arcname=intermediate_name)
        zf.write(output_tsv, arcname=output_name)


def _bundle_member_name(file_stem: str, suffix: str, tag: str) -> str:
    safe_stem = Path(file_stem).name
    if safe_stem.endswith(f"_{tag}"):
        return f"{safe_stem}{suffix}"
    return f"{safe_stem}_{tag}{suffix}"


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
async def process_file_compat(request: Request) -> FileResponse:
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
    bundle_target = WORK_ROOT / job_id / filename
    if bundle_target.exists():
        return FileResponse(path=bundle_target, filename=filename, media_type="application/zip")

    tsv_target = WORK_ROOT / job_id / "output" / filename
    if tsv_target.exists():
        return FileResponse(path=tsv_target, filename=filename, media_type="text/tab-separated-values")

    raise HTTPException(status_code=404, detail="Output file not found")

"""Optional FastAPI wrapper, mirroring the other servers in this repo.

    uvicorn wls_outlier_lab.app:app --host 0.0.0.0 --port 8020

POST /analyze runs the full detector-comparison + constellation-ablation
pipeline on a measurement file and writes a report directory. The request only
names *where the data lives*; all analysis is the pure core.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

try:
    from fastapi import FastAPI, HTTPException
    from pydantic import BaseModel
    _FASTAPI = True
except Exception:  # pragma: no cover - server deps optional
    _FASTAPI = False

from .core.detectors import DetectorConfig
from .core.experiment import constellation_ablation, run_experiment
from .core.wls import WeightConfig, WlsConfig
from .reporting import write_report
from .sources.measurement_tsv import load_epochs
from .sources.truth_bestpos import BestposTruthSource
from .sources.truth_columns import ColumnTruthSource


def _run(*, measurements: str, output_dir: str, truth_columns: bool,
         truth_file: Optional[str], truth_bestpos: Optional[str],
         stride: int, max_epochs: Optional[int], detectors: Optional[List[str]],
         min_cn0: float, min_elevation: float, residual_k: float,
         match_tolerance: float, ablation: bool, plots: bool) -> dict:
    if not Path(measurements).exists():
        raise FileNotFoundError(measurements)
    epochs = load_epochs(measurements, stride=stride, max_epochs=max_epochs)
    if not epochs:
        raise ValueError("no epochs loaded")

    if truth_bestpos:
        truth = BestposTruthSource(truth_bestpos).load_track()
    elif truth_file:
        truth = ColumnTruthSource(truth_file).load_track()
    elif truth_columns:
        truth = ColumnTruthSource(measurements).load_track()
    else:
        raise ValueError("no truth source specified")

    wls_cfg = WlsConfig(weights=WeightConfig())
    dcfg = DetectorConfig(min_cn0_dbhz=min_cn0, min_elevation_deg=min_elevation,
                          residual_mad_k=residual_k)
    result = run_experiment(epochs, truth, wls_cfg=wls_cfg, dcfg=dcfg,
                            detectors=detectors, match_tolerance_sec=match_tolerance)
    abl = constellation_ablation(epochs, truth, wls_cfg=wls_cfg,
                                 match_tolerance_sec=match_tolerance) if ablation else None
    meta = {"measurements": measurements, "n_epochs": len(epochs), "stride": stride}
    paths = write_report(output_dir, result, ablation_rows=abl, meta=meta, make_plots=plots)
    return {
        "n_epochs": len(epochs),
        "truth_samples": len(truth.samples),
        "best_detector": result.comparison()[0]["detector"] if result.comparison() else None,
        "comparison": result.comparison(),
        "outputs": paths,
    }


if _FASTAPI:
    app = FastAPI(title="WLS Outlier Lab", version="0.1.0")

    class AnalyzeRequest(BaseModel):
        measurements: str
        output_dir: str = "wls_outlier_lab/results/latest"
        truth_columns: bool = False
        truth_file: Optional[str] = None
        truth_bestpos: Optional[str] = None
        stride: int = 1
        max_epochs: Optional[int] = None
        detectors: Optional[List[str]] = None
        min_cn0: float = 30.0
        min_elevation: float = 10.0
        residual_k: float = 5.0
        match_tolerance: float = 0.5
        ablation: bool = True
        plots: bool = True

    @app.get("/health")
    def health() -> dict:
        return {"status": "ok", "fastapi_available": True}

    @app.post("/analyze")
    def analyze(req: AnalyzeRequest) -> dict:
        try:
            return _run(measurements=req.measurements, output_dir=req.output_dir,
                        truth_columns=req.truth_columns, truth_file=req.truth_file,
                        truth_bestpos=req.truth_bestpos, stride=req.stride,
                        max_epochs=req.max_epochs, detectors=req.detectors,
                        min_cn0=req.min_cn0, min_elevation=req.min_elevation,
                        residual_k=req.residual_k, match_tolerance=req.match_tolerance,
                        ablation=req.ablation, plots=req.plots)
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except Exception as exc:  # noqa: BLE001
            raise HTTPException(status_code=500, detail=repr(exc))
else:  # pragma: no cover
    app = None

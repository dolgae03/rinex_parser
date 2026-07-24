"""Experiment orchestration: baseline vs each detector, and per-constellation ablation.

This is where the study question is answered: *which measurements should we flag,
and does removing them improve WLS against the Novatel/RTK truth?*
Still pure logic — epochs and a truth track go in, comparison tables come out.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple

from ..types import Epoch, TruthTrack, constellation_name
from .detectors import DETECTORS, DetectorConfig
from .metrics import ErrorMetrics, compute_metrics
from .wls import WlsConfig, solve_epoch

SatId = Tuple[int, int]  # (constellation, prn)


@dataclass
class DetectorRun:
    name: str
    metrics: ErrorMetrics
    n_signal_rejections: int = 0
    epochs_with_rejection: int = 0
    rejections_by_constellation: Dict[int, int] = field(default_factory=dict)


@dataclass
class ExperimentResult:
    baseline_name: str
    runs: Dict[str, DetectorRun]

    def comparison(self) -> List[Dict[str, object]]:
        """Rows comparing every detector to the baseline, sorted best-first."""
        base = self.runs[self.baseline_name].metrics
        rows: List[Dict[str, object]] = []
        for name, run in self.runs.items():
            m = run.metrics
            d_h = m.horizontal_rmse_m - base.horizontal_rmse_m
            pct = (100.0 * d_h / base.horizontal_rmse_m
                   if base.horizontal_rmse_m and base.horizontal_rmse_m == base.horizontal_rmse_m
                   else float("nan"))
            rows.append({
                "detector": name,
                "horizontal_rmse_m": m.horizontal_rmse_m,
                "horizontal_p95_m": m.horizontal_p95_m,
                "vertical_rmse_m": m.vertical_rmse_m,
                "cep95_m": m.cep95_m,
                "availability_pct": m.availability_pct,
                "mean_sats_used": m.mean_sats_used,
                "delta_horizontal_rmse_m": d_h,
                "improvement_pct": -pct,  # positive = improvement
                "signal_rejections": run.n_signal_rejections,
                "epochs_with_rejection": run.epochs_with_rejection,
            })
        rows.sort(key=lambda r: (r["horizontal_rmse_m"]
                                 if r["horizontal_rmse_m"] == r["horizontal_rmse_m"]
                                 else float("inf")))
        return rows


def run_experiment(
    epochs: Sequence[Epoch],
    truth: TruthTrack,
    wls_cfg: Optional[WlsConfig] = None,
    dcfg: Optional[DetectorConfig] = None,
    detectors: Optional[Sequence[str]] = None,
    match_tolerance_sec: float = 0.5,
    baseline_name: str = "baseline",
) -> ExperimentResult:
    wls_cfg = wls_cfg or WlsConfig()
    dcfg = dcfg or DetectorConfig()
    names = list(detectors) if detectors else list(DETECTORS.keys())
    if baseline_name not in names:
        names.insert(0, baseline_name)

    runs: Dict[str, DetectorRun] = {}
    for name in names:
        detector = DETECTORS[name]
        solutions = []
        n_rej = 0
        epochs_hit = 0
        rej_by_cons: Dict[int, int] = {}
        for ep in epochs:
            reject: Set[SatId] = detector(ep.obs, ep.t_sec, wls_cfg, dcfg)
            if reject:
                epochs_hit += 1
                n_rej += len(reject)
                for k in reject:
                    rej_by_cons[k[0]] = rej_by_cons.get(k[0], 0) + 1
            solutions.append(solve_epoch(ep.obs, ep.t_sec, wls_cfg, reject=reject))
        metrics = compute_metrics(solutions, truth, label=name,
                                  match_tolerance_sec=match_tolerance_sec,
                                  n_epochs_total=len(epochs))
        runs[name] = DetectorRun(
            name=name, metrics=metrics, n_signal_rejections=n_rej,
            epochs_with_rejection=epochs_hit, rejections_by_constellation=rej_by_cons,
        )
    return ExperimentResult(baseline_name=baseline_name, runs=runs)


def constellation_ablation(
    epochs: Sequence[Epoch],
    truth: TruthTrack,
    wls_cfg: Optional[WlsConfig] = None,
    match_tolerance_sec: float = 0.5,
) -> List[Dict[str, object]]:
    """Solve with each constellation removed (and each alone) to see its effect.

    Answers the 'per satellite constellation' angle: is any constellation
    dragging the WLS solution, and how much does each contribute?
    """
    wls_cfg = wls_cfg or WlsConfig()
    present: Set[int] = set()
    for ep in epochs:
        present.update(o.constellation for o in ep.obs)

    rows: List[Dict[str, object]] = []

    def _solve_subset(keep=None, drop=None) -> ErrorMetrics:
        sols = []
        for ep in epochs:
            reject = {
                (o.constellation, o.prn) for o in ep.obs
                if (keep is not None and o.constellation not in keep)
                or (drop is not None and o.constellation in drop)
            }
            sols.append(solve_epoch(ep.obs, ep.t_sec, wls_cfg, reject=reject))
        return compute_metrics(sols, truth, match_tolerance_sec=match_tolerance_sec,
                               n_epochs_total=len(epochs))

    full = _solve_subset()
    rows.append({"config": "all_constellations", "constellation": "ALL",
                 "horizontal_rmse_m": full.horizontal_rmse_m,
                 "availability_pct": full.availability_pct,
                 "mean_sats_used": full.mean_sats_used})
    for c in sorted(present):
        drop_m = _solve_subset(drop={c})
        only_m = _solve_subset(keep={c})
        rows.append({"config": "drop", "constellation": constellation_name(c),
                     "horizontal_rmse_m": drop_m.horizontal_rmse_m,
                     "availability_pct": drop_m.availability_pct,
                     "mean_sats_used": drop_m.mean_sats_used,
                     "delta_vs_all_m": drop_m.horizontal_rmse_m - full.horizontal_rmse_m})
        rows.append({"config": "only", "constellation": constellation_name(c),
                     "horizontal_rmse_m": only_m.horizontal_rmse_m,
                     "availability_pct": only_m.availability_pct,
                     "mean_sats_used": only_m.mean_sats_used})
    return rows

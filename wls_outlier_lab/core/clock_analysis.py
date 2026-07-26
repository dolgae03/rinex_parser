"""Receiver-clock stability analysis.

Question: does the smartphone's clock-drift instability degrade the (horizontal)
navigation solution, and can we flag anomalies from it?

Each single-epoch WLS estimates a receiver clock bias ``clk(t)`` [m]. From that
series we form:
  - drift(t)        = d(clk)/dt                         [m/s]
  - instability(t)  = | clk(t) - (clk(t-1) + drift(t-1)*dt) |   [m]
    i.e. how far the clock departs from a constant-drift (constant-velocity)
    extrapolation — a jump or a rate change shows up as a spike.

We then flag high-instability epochs and test whether instability correlates
with the horizontal error against truth (Pearson + rank/Spearman), plus a
stable-vs-unstable horizontal-error comparison. Pure logic; the solve/truth come
from the same engine everything else uses.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

from ..types import Epoch, TruthTrack
from . import frames
from .detectors import DETECTORS, DetectorConfig
from .wls import WlsConfig, solve_epoch


@dataclass
class ClockPoint:
    t_sec: float
    dt_s: float
    clk_m: float
    drift_mps: float
    instability_m: float
    is_anomaly: bool
    horizontal_error_m: float


@dataclass
class ClockAnalysis:
    points: List[ClockPoint] = field(default_factory=list)
    n_epochs: int = 0
    n_anomalies: int = 0
    instability_scale_m: float = float("nan")
    instability_threshold_m: float = float("nan")
    median_drift_mps: float = float("nan")
    pearson_r: float = float("nan")
    spearman_r: float = float("nan")
    stable_mean_h_m: float = float("nan")
    unstable_mean_h_m: float = float("nan")
    stable_median_h_m: float = float("nan")
    unstable_median_h_m: float = float("nan")
    bins: List[Dict[str, float]] = field(default_factory=list)
    detector: str = ""

    def summary(self) -> Dict[str, object]:
        return {k: getattr(self, k) for k in (
            "n_epochs", "n_anomalies", "instability_scale_m", "instability_threshold_m",
            "median_drift_mps", "pearson_r", "spearman_r", "stable_mean_h_m",
            "unstable_mean_h_m", "stable_median_h_m", "unstable_median_h_m",
            "bins", "detector")}


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return float("nan")
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return _pearson(rx, ry)


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return float("nan")
    sx, sy = np.std(x), np.std(y)
    if sx == 0 or sy == 0:
        return float("nan")
    return float(np.mean((x - x.mean()) * (y - y.mean())) / (sx * sy))


def _horizontal_error(sol_ecef, truth: TruthTrack, t_sec: float, tol: float) -> float:
    ts = truth.nearest(t_sec, tol)
    if ts is None:
        return float("nan")
    origin = np.asarray(ts.ecef, float)
    enu = frames.enu_rotation_matrix(*frames.ecef_to_lla(*origin)[:2]) @ (
        np.asarray(sol_ecef, float) - origin)
    return float(math.hypot(enu[0], enu[1]))


def analyze_clock_stability(
    epochs: Sequence[Epoch],
    truth: TruthTrack,
    wls_cfg: Optional[WlsConfig] = None,
    dcfg: Optional[DetectorConfig] = None,
    detector: str = "combined",
    match_tolerance_sec: float = 0.5,
    max_dt_s: float = 5.0,
    instability_k: float = 5.0,
) -> ClockAnalysis:
    wls_cfg = wls_cfg or WlsConfig()
    dcfg = dcfg or DetectorConfig()
    det = DETECTORS[detector]

    solved = []
    for ep in epochs:
        reject = det(ep.obs, ep.t_sec, wls_cfg, dcfg)
        s = solve_epoch(ep.obs, ep.t_sec, wls_cfg, reject=reject)
        if s.converged and math.isfinite(s.clock_bias_m) and all(math.isfinite(v) for v in s.ecef):
            solved.append(s)
    solved.sort(key=lambda s: s.t_sec)
    res = ClockAnalysis(detector=detector, n_epochs=len(solved))
    if len(solved) < 4:
        return res

    t = np.array([s.t_sec for s in solved])
    clk = np.array([s.clock_bias_m for s in solved])
    n = len(t)
    dt = np.full(n, float("nan"))
    drift = np.full(n, float("nan"))
    inst = np.full(n, float("nan"))
    for i in range(1, n):
        d = t[i] - t[i - 1]
        if 0 < d <= max_dt_s:
            dt[i] = d
            drift[i] = (clk[i] - clk[i - 1]) / d
    for i in range(2, n):
        if math.isfinite(drift[i]) and math.isfinite(drift[i - 1]) and math.isfinite(dt[i]):
            pred = clk[i - 1] + drift[i - 1] * dt[i]
            inst[i] = abs(clk[i] - pred)

    finite_inst = inst[np.isfinite(inst)]
    if finite_inst.size:
        med = float(np.median(finite_inst))
        scale = max(1.4826 * float(np.median(np.abs(finite_inst - med))), 1e-6)
    else:
        scale = float("nan")
    thr = instability_k * scale if math.isfinite(scale) else float("nan")
    res.instability_scale_m = scale
    res.instability_threshold_m = thr
    res.median_drift_mps = float(np.nanmedian(drift)) if np.isfinite(drift).any() else float("nan")

    herr = np.array([_horizontal_error(s.ecef, truth, s.t_sec, match_tolerance_sec) for s in solved])

    pts: List[ClockPoint] = []
    for i in range(n):
        anom = bool(math.isfinite(inst[i]) and math.isfinite(thr) and inst[i] > thr)
        pts.append(ClockPoint(
            t_sec=float(t[i]), dt_s=float(dt[i]) if math.isfinite(dt[i]) else float("nan"),
            clk_m=float(clk[i]), drift_mps=float(drift[i]) if math.isfinite(drift[i]) else float("nan"),
            instability_m=float(inst[i]) if math.isfinite(inst[i]) else float("nan"),
            is_anomaly=anom, horizontal_error_m=float(herr[i]) if math.isfinite(herr[i]) else float("nan"),
        ))
    res.points = pts
    res.n_anomalies = sum(p.is_anomaly for p in pts)

    # correlation of instability vs horizontal error (where both defined)
    mask = np.isfinite(inst) & np.isfinite(herr)
    if mask.sum() >= 3:
        res.pearson_r = _pearson(inst[mask], herr[mask])
        res.spearman_r = _spearman(inst[mask], herr[mask])

    # stable vs unstable horizontal error
    anom_mask = np.array([p.is_anomaly for p in pts]) & np.isfinite(herr)
    stable_mask = (~np.array([p.is_anomaly for p in pts])) & np.isfinite(herr) & np.isfinite(inst)
    if anom_mask.any():
        res.unstable_mean_h_m = float(np.mean(herr[anom_mask]))
        res.unstable_median_h_m = float(np.median(herr[anom_mask]))
    if stable_mask.any():
        res.stable_mean_h_m = float(np.mean(herr[stable_mask]))
        res.stable_median_h_m = float(np.median(herr[stable_mask]))

    # instability quartile bins vs horizontal error
    if mask.sum() >= 8:
        iv, hv = inst[mask], herr[mask]
        qs = np.quantile(iv, [0, 0.25, 0.5, 0.75, 1.0])
        for q in range(4):
            lo, hi = qs[q], qs[q + 1]
            sel = (iv >= lo) & (iv <= hi) if q == 3 else (iv >= lo) & (iv < hi)
            if sel.any():
                res.bins.append({
                    "quartile": q + 1,
                    "instability_lo_m": round(float(lo), 3),
                    "instability_hi_m": round(float(hi), 3),
                    "n": int(sel.sum()),
                    "mean_h_m": round(float(np.mean(hv[sel])), 3),
                    "median_h_m": round(float(np.median(hv[sel])), 3),
                    "p95_h_m": round(float(np.percentile(hv[sel], 95)), 3),
                })
    return res

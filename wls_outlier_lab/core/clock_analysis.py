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
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..types import Epoch, SatObs, TruthTrack
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
    # confounders + the independent Doppler drift (filled when available)
    n_sats: int = 0
    hdop: float = float("nan")
    clk_sigma_m: float = float("nan")       # formal sigma of the estimated clock
    drift_dop_mps: float = float("nan")     # Doppler-derived clock drift
    drift_dop_change_mps: float = float("nan")  # |d(drift_dop)| between epochs
    is_dop_anomaly: bool = False            # anomaly per the low-noise Doppler drift


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
    # --- is the "instability" real clock motion or just estimation noise? ---
    median_clk_sigma_m: float = float("nan")   # formal sigma of the WLS clock estimate
    noise_expected_inst_m: float = float("nan")  # sqrt(6)*sigma_clk: instability a
    #   pure-noise clock series would already show (2nd difference of white noise)
    noise_fraction: float = float("nan")       # noise_expected / observed scale
    # --- independent Doppler cross-check ---
    doppler_sign: float = float("nan")
    n_doppler_epochs: int = 0
    doppler_resid_rms_mps: float = float("nan")  # fit quality; large => sv_vel/Doppler
    #   inconsistency in the source data, so the Doppler cross-check is not trustworthy
    doppler_rejected_per_epoch: float = float("nan")
    drift_dop_median_mps: float = float("nan")
    drift_dop_std_mps: float = float("nan")
    drift_agreement_rms_mps: float = float("nan")   # RMS(pos-drift - doppler-drift)
    dop_drift_change_scale_mps: float = float("nan")  # robust scale of |d(drift_dop)|
    pos_drift_change_scale_mps: float = float("nan")  # same from position-domain clock
    dop_pearson_r: float = float("nan")   # doppler drift change vs horizontal error
    dop_spearman_r: float = float("nan")
    # Anomalies detected on the Doppler drift — the trustworthy clock-anomaly
    # detector. The position-domain flag is kept for comparison but is dominated
    # by estimation noise, so the two need not agree.
    n_dop_anomalies: int = 0
    dop_instability_threshold_mps: float = float("nan")
    anomaly_agreement: int = 0            # epochs flagged by BOTH detectors
    # --- confounder control / lag ---
    partial_r: float = float("nan")       # instability vs herr, controlling HDOP & n_sats
    lag1_pearson_r: float = float("nan")  # instability(i) vs horizontal error(i+1)
    worst_geometry_pearson_r: float = float("nan")  # within the worst-HDOP quartile

    def summary(self) -> Dict[str, object]:
        return {k: getattr(self, k) for k in (
            "n_epochs", "n_anomalies", "instability_scale_m", "instability_threshold_m",
            "median_drift_mps", "pearson_r", "spearman_r", "stable_mean_h_m",
            "unstable_mean_h_m", "stable_median_h_m", "unstable_median_h_m",
            "median_clk_sigma_m", "noise_expected_inst_m", "noise_fraction",
            "doppler_sign", "n_doppler_epochs", "doppler_resid_rms_mps",
            "doppler_rejected_per_epoch", "drift_dop_median_mps",
            "drift_dop_std_mps", "drift_agreement_rms_mps",
            "dop_drift_change_scale_mps", "pos_drift_change_scale_mps",
            "dop_pearson_r", "dop_spearman_r", "n_dop_anomalies",
            "dop_instability_threshold_mps", "anomaly_agreement",
            "partial_r", "lag1_pearson_r", "worst_geometry_pearson_r",
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


def _partial_r(y: np.ndarray, x: np.ndarray, controls: Sequence[np.ndarray]) -> float:
    """Correlation of x with y after regressing both on the control variables."""
    cols = [c for c in controls if np.isfinite(c).all() and np.std(c) > 0]
    if not cols or len(y) < 5:
        return float("nan")
    A = np.column_stack([np.ones(len(y))] + list(cols))
    try:
        resid = lambda v: v - A @ np.linalg.lstsq(A, v, rcond=None)[0]
    except np.linalg.LinAlgError:
        return float("nan")
    return _pearson(resid(x), resid(y))


def doppler_clock_drift(
    epoch: Epoch,
    rx_ecef: Sequence[float],
    reject: Optional[set] = None,
    sign: float = -1.0,
    min_sats: int = 6,
    mad_k: float = 3.0,
) -> Dict[str, float]:
    """Receiver velocity + **clock drift** from Doppler — independent of the clock
    series obtained by differencing per-epoch position-domain clock estimates.

    Range rate from Doppler: ``rdot = sign * (c / f) * doppler_hz``. Model:
    ``rdot_i = u_i . (v_sv_i - v_rx) + drift`` where ``u_i`` is the line of sight
    and ``drift`` [m/s] is the receiver clock drift, common to all satellites.
    Satellite clock drift (~mm/s) is neglected. Linear least squares in
    ``[v_rx, drift]``, one signal per satellite, with iterated MAD rejection.

    The rejection has to be aggressive: in real processed logs a slice of the
    ``sv_vel`` column can be wrong by hundreds of m/s (seen here for ~18% of the
    BeiDou rows, while GPS/Galileo/QZSS match ``d(sv_pos)/dt`` to 3 mm/s), and
    such rows would otherwise drag the fit. ``n_rejected`` and ``resid_rms_mps``
    are returned so the caller can tell a clean epoch from a contaminated one.
    """
    reject = reject or set()
    rx = np.asarray(rx_ecef, float)
    # one signal per satellite (strongest C/N0): extra bands carry the same
    # range rate, so they only inflate the apparent redundancy
    best: Dict[Tuple[int, int], SatObs] = {}
    for o in epoch.obs:
        if (o.constellation, o.prn) in reject or o.sv_vel is None:
            continue
        if not (math.isfinite(o.doppler_hz) and o.frequency_hz > 0):
            continue
        if not all(math.isfinite(v) for v in o.sv_pos) or not all(math.isfinite(v) for v in o.sv_vel):
            continue
        k = (o.constellation, o.prn)
        cur = best.get(k)
        cn0 = o.cn0_dbhz if math.isfinite(o.cn0_dbhz) else -1.0
        if cur is None or cn0 > (cur.cn0_dbhz if math.isfinite(cur.cn0_dbhz) else -1.0):
            best[k] = o
    rows, rdot = [], []
    for o in best.values():
        d = np.asarray(o.sv_pos, float) - rx
        rng = float(np.linalg.norm(d))
        if rng < 1.0:
            continue
        u = d / rng
        rows.append(np.array([-u[0], -u[1], -u[2], 1.0]))
        rdot.append(sign * frames.C_LIGHT / o.frequency_hz * o.doppler_hz
                    - float(np.dot(u, np.asarray(o.sv_vel, float))))
    out = {"drift_mps": float("nan"), "n_sats": 0, "n_rejected": 0,
           "resid_rms_mps": float("nan"), "speed_mps": float("nan")}
    n0 = len(rows)
    if n0 < min_sats:
        return out
    A = np.array(rows)
    b = np.array(rdot)
    x = v = None
    for _ in range(5):
        try:
            x, *_ = np.linalg.lstsq(A, b, rcond=None)
        except np.linalg.LinAlgError:
            return out
        v = b - A @ x
        s = max(1.4826 * float(np.median(np.abs(v - np.median(v)))), 0.02)
        keep = np.abs(v - np.median(v)) <= mad_k * s
        if keep.all() or keep.sum() < min_sats:
            break
        A, b = A[keep], b[keep]
    out["drift_mps"] = float(x[3])
    out["n_sats"] = int(len(b))
    out["n_rejected"] = int(n0 - len(b))
    out["resid_rms_mps"] = float(np.sqrt(np.mean(v * v)))
    out["speed_mps"] = float(np.linalg.norm(x[:3]))
    return out


def _pick_doppler_sign(pairs: Sequence, n_probe: int = 30) -> float:
    """Doppler sign convention differs by producer; choose it from the data."""
    best, best_rms = -1.0, float("inf")
    for sign in (-1.0, +1.0):
        rms = []
        for ep, s, rej in list(pairs)[:n_probe]:
            r = doppler_clock_drift(ep, s.ecef, reject=rej, sign=sign)
            if math.isfinite(r["resid_rms_mps"]):
                rms.append(r["resid_rms_mps"])
        if rms and float(np.median(rms)) < best_rms:
            best_rms, best = float(np.median(rms)), sign
    return best


def clock_coasting_experiment(
    epochs: Sequence[Epoch],
    truth: TruthTrack,
    wls_cfg: Optional[WlsConfig] = None,
    dcfg: Optional[DetectorConfig] = None,
    detector: str = "combined",
    match_tolerance_sec: float = 0.5,
    sigmas_m: Sequence[float] = (float("inf"), 100.0, 30.0, 10.0, 3.0, 1.0, 0.3),
    max_dt_s: float = 5.0,
    sat_budget: Optional[int] = None,
) -> List[Dict[str, float]]:
    """**Decisive test.** A freely estimated clock absorbs its own instability, so
    zero correlation with horizontal error is the *expected* result and proves
    little. Here we instead make the solution *depend* on clock stability: the
    clock is coasted from the previous epoch (``clk_prev + drift_prev*dt``) and
    imposed as a pseudo-observation with standard deviation ``sigma``. As sigma
    tightens the fix must trust the clock model; a receiver whose clock is
    genuinely unstable degrades, a stable one does not.

    ``sat_budget`` keeps only the N strongest satellites per epoch, which is where
    a clock constraint actually carries weight (with 40 satellites the clock is
    heavily over-determined and the constraint changes nothing).
    """
    wls_cfg = wls_cfg or WlsConfig()
    dcfg = dcfg or DetectorConfig()
    det = DETECTORS[detector]
    eps = sorted(epochs, key=lambda e: e.t_sec)

    # Detector rejections (the expensive part) are computed once and reused.
    rejects = []
    for ep in eps:
        rej = set(det(ep.obs, ep.t_sec, wls_cfg, dcfg))
        if sat_budget:
            cand = sorted(
                {(o.constellation, o.prn) for o in ep.obs} - rej,
                key=lambda sid: -max((o.cn0_dbhz if math.isfinite(o.cn0_dbhz) else -1.0)
                                     for o in ep.obs if (o.constellation, o.prn) == sid))
            rej |= set(cand[sat_budget:])
        rejects.append(rej)

    rows: List[Dict[str, float]] = []
    for sigma in sigmas_m:
        prev_t = prev_clk = prev_drift = None
        herrs, sats, n_coasted = [], [], 0
        for ep, rej in zip(eps, rejects):
            prior = None
            if math.isfinite(sigma) and prev_clk is not None and prev_drift is not None:
                dt = ep.t_sec - prev_t
                if 0 < dt <= max_dt_s:
                    prior = (prev_clk + prev_drift * dt, sigma)
            s = solve_epoch(ep.obs, ep.t_sec, wls_cfg, reject=rej, clock_prior=prior)
            if not (math.isfinite(s.clock_bias_m) and all(math.isfinite(v) for v in s.ecef)):
                continue
            if prior is not None:
                n_coasted += 1
            h = _horizontal_error(s.ecef, truth, s.t_sec, match_tolerance_sec)
            if math.isfinite(h):
                herrs.append(h)
            sats.append(s.n_used)
            if prev_t is not None and 0 < ep.t_sec - prev_t <= max_dt_s:
                prev_drift = (s.clock_bias_m - prev_clk) / (ep.t_sec - prev_t)
            prev_t, prev_clk = ep.t_sec, s.clock_bias_m
        h = np.array(herrs)
        rows.append({
            "clock_sigma_m": float(sigma),
            "mode": "free" if not math.isfinite(sigma) else f"coast_{sigma:g}m",
            "n_epochs": int(len(h)),
            "n_coasted": n_coasted,
            "sat_budget": sat_budget or 0,
            "mean_sats": round(float(np.mean(sats)), 2) if sats else float("nan"),
            "h_rmse_m": round(float(np.sqrt(np.mean(h * h))), 3) if h.size else float("nan"),
            "h_median_m": round(float(np.median(h)), 3) if h.size else float("nan"),
            "h_p95_m": round(float(np.percentile(h, 95)), 3) if h.size else float("nan"),
            "h_max_m": round(float(np.max(h)), 3) if h.size else float("nan"),
        })
    return rows


def analyze_clock_stability(
    epochs: Sequence[Epoch],
    truth: TruthTrack,
    wls_cfg: Optional[WlsConfig] = None,
    dcfg: Optional[DetectorConfig] = None,
    detector: str = "combined",
    match_tolerance_sec: float = 0.5,
    max_dt_s: float = 5.0,
    instability_k: float = 5.0,
    with_doppler: bool = True,
) -> ClockAnalysis:
    wls_cfg = wls_cfg or WlsConfig()
    dcfg = dcfg or DetectorConfig()
    det = DETECTORS[detector]

    pairs = []
    for ep in epochs:
        reject = det(ep.obs, ep.t_sec, wls_cfg, dcfg)
        s = solve_epoch(ep.obs, ep.t_sec, wls_cfg, reject=reject)
        if s.converged and math.isfinite(s.clock_bias_m) and all(math.isfinite(v) for v in s.ecef):
            pairs.append((ep, s, reject))
    pairs.sort(key=lambda p: p[1].t_sec)
    solved = [p[1] for p in pairs]
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

    # --- Is the apparent instability just clock-estimation noise? ---------------
    # The clock is re-estimated every epoch, so its estimate carries noise of
    # sigma_clk. The instability metric is a second difference, which for white
    # noise alone already has RMS sqrt(6)*sigma_clk. If that matches the observed
    # scale, the "instability" says nothing about the physical oscillator.
    clk_sig = np.array([s.clock_sigma_m for s in solved])
    if np.isfinite(clk_sig).any():
        res.median_clk_sigma_m = float(np.nanmedian(clk_sig))
        res.noise_expected_inst_m = math.sqrt(6.0) * res.median_clk_sigma_m
        obs_scale = 1.2533 * scale if math.isfinite(scale) else float("nan")  # MAD->RMS for |N|
        if math.isfinite(obs_scale) and obs_scale > 0:
            res.noise_fraction = res.noise_expected_inst_m / obs_scale

    # --- Independent Doppler-derived clock drift --------------------------------
    drift_dop = np.full(n, float("nan"))
    if with_doppler:
        sign = _pick_doppler_sign(pairs)
        res.doppler_sign = sign
        rr, nrej = [], []
        for i, (ep, s, rej) in enumerate(pairs):
            r = doppler_clock_drift(ep, s.ecef, reject=rej, sign=sign)
            drift_dop[i] = r["drift_mps"]
            if math.isfinite(r["resid_rms_mps"]):
                rr.append(r["resid_rms_mps"])
                nrej.append(r["n_rejected"])
        if rr:
            res.doppler_resid_rms_mps = float(np.median(rr))
            res.doppler_rejected_per_epoch = float(np.mean(nrej))
        ok = np.isfinite(drift_dop)
        res.n_doppler_epochs = int(ok.sum())
        if ok.sum() >= 4:
            res.drift_dop_median_mps = float(np.median(drift_dop[ok]))
            res.drift_dop_std_mps = float(np.std(drift_dop[ok]))
            both = ok & np.isfinite(drift)
            if both.sum() >= 4:
                d = drift[both] - drift_dop[both]
                res.drift_agreement_rms_mps = float(np.sqrt(np.mean((d - np.median(d)) ** 2)))

    # Drift *change* between epochs — the quantity "drift instability" names —
    # measured two ways: from Doppler (low noise) and from the position-domain clock.
    dop_chg = np.full(n, float("nan"))
    for i in range(1, n):
        if math.isfinite(drift_dop[i]) and math.isfinite(drift_dop[i - 1]) and math.isfinite(dt[i]):
            dop_chg[i] = abs(drift_dop[i] - drift_dop[i - 1])
    pos_chg = np.full(n, float("nan"))
    for i in range(2, n):
        if math.isfinite(inst[i]) and math.isfinite(dt[i]) and dt[i] > 0:
            pos_chg[i] = inst[i] / dt[i]
    for name, arr in (("dop_drift_change_scale_mps", dop_chg),
                      ("pos_drift_change_scale_mps", pos_chg)):
        v = arr[np.isfinite(arr)]
        if v.size:
            setattr(res, name, float(np.median(v)))

    # Clock anomalies on the Doppler drift. This is the detector to trust: its fit
    # residual is ~0.01 m/s, whereas the position-domain clock carries metre-level
    # estimation noise that spikes whenever the receiver is moving.
    dv = dop_chg[np.isfinite(dop_chg)]
    dop_thr = float("nan")
    if dv.size:
        dmed = float(np.median(dv))
        dscale = max(1.4826 * float(np.median(np.abs(dv - dmed))), 1e-6)
        dop_thr = dmed + instability_k * dscale
    res.dop_instability_threshold_mps = dop_thr

    pts: List[ClockPoint] = []
    for i in range(n):
        anom = bool(math.isfinite(inst[i]) and math.isfinite(thr) and inst[i] > thr)
        dop_anom = bool(math.isfinite(dop_chg[i]) and math.isfinite(dop_thr)
                        and dop_chg[i] > dop_thr)
        s = solved[i]
        pts.append(ClockPoint(
            t_sec=float(t[i]), dt_s=float(dt[i]) if math.isfinite(dt[i]) else float("nan"),
            clk_m=float(clk[i]), drift_mps=float(drift[i]) if math.isfinite(drift[i]) else float("nan"),
            instability_m=float(inst[i]) if math.isfinite(inst[i]) else float("nan"),
            is_anomaly=anom, horizontal_error_m=float(herr[i]) if math.isfinite(herr[i]) else float("nan"),
            n_sats=int(s.n_used), hdop=float(s.hdop), clk_sigma_m=float(s.clock_sigma_m),
            drift_dop_mps=float(drift_dop[i]), drift_dop_change_mps=float(dop_chg[i]),
            is_dop_anomaly=dop_anom,
        ))
    res.points = pts
    res.n_anomalies = sum(p.is_anomaly for p in pts)
    res.n_dop_anomalies = sum(p.is_dop_anomaly for p in pts)
    res.anomaly_agreement = sum(p.is_anomaly and p.is_dop_anomaly for p in pts)

    # correlation of instability vs horizontal error (where both defined)
    mask = np.isfinite(inst) & np.isfinite(herr)
    if mask.sum() >= 3:
        res.pearson_r = _pearson(inst[mask], herr[mask])
        res.spearman_r = _spearman(inst[mask], herr[mask])

    # Same test with the low-noise Doppler measure of drift instability.
    dmask = np.isfinite(dop_chg) & np.isfinite(herr)
    if dmask.sum() >= 3:
        res.dop_pearson_r = _pearson(dop_chg[dmask], herr[dmask])
        res.dop_spearman_r = _spearman(dop_chg[dmask], herr[dmask])

    # Confounder control: satellite count and geometry drive the horizontal error
    # on their own, so partial out HDOP and n_sats before believing (or dismissing)
    # any instability effect.
    hdop = np.array([p.hdop for p in pts])
    nsat = np.array([float(p.n_sats) for p in pts])
    cmask = mask & np.isfinite(hdop) & np.isfinite(nsat)
    if cmask.sum() >= 8:
        res.partial_r = _partial_r(herr[cmask], inst[cmask], [hdop[cmask], nsat[cmask]])
        q75 = float(np.quantile(hdop[cmask], 0.75))
        worst = cmask & (hdop >= q75)
        if worst.sum() >= 5:
            res.worst_geometry_pearson_r = _pearson(inst[worst], herr[worst])

    # Lag-1: does instability now predict error at the next epoch? (What a filter
    # that propagates the clock would suffer from.)
    if n >= 5:
        a, b = inst[:-1], herr[1:]
        lm = np.isfinite(a) & np.isfinite(b)
        if lm.sum() >= 3:
            res.lag1_pearson_r = _pearson(a[lm], b[lm])

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

"""Per-epoch factor matrix for the factor-vs-horizontal-error correlation study.

Division of labour: this module only *extracts* one row of candidate factors
per solved epoch (it owns the WLS engine, the detector and the truth match).
The statistical analysis itself — speed-bin segmentation, Pearson/Spearman/
partial correlations, block-bootstrap confidence intervals, factor ranking and
every figure — is implemented in MATLAB (``factor_correlation_analysis.m``),
which reads the ``factor_epochs.csv`` this module emits.

Factor families (per epoch, aggregated over the satellites the fix used):

* geometry      — n_used/n_available/n_rejected, H/V/P/G/TDOP, elevation stats,
                  per-constellation counts
* signal        — C/N0 statistics of the used sats and of the whole sky
* consistency   — post-fit residual stats, sigma0, formal clock sigma
* clock         — estimated bias b, position-domain drift b_dot, 1-step
                  instability |clk - (clk_prev + drift_prev*dt)|, plus the
                  independent Doppler drift/its epoch-to-epoch change and the
                  Doppler LS residual (= PR-rate consistency)
* code quality  — Doppler-Code Difference rate (receiver clock drift cancels,
                  so it isolates code multipath/noise rate) and Code-Carrier
                  Divergence rate (2x iono rate + multipath), the two
                  measurement-domain parameters named in the 3rd-year proposal
* dynamics      — truth-derived horizontal speed/acceleration, Doppler speed

The target column is ``h_err_m`` (horizontal error vs truth); vertical error is
exported for completeness but is not an evaluation criterion in this lab.
"""

from __future__ import annotations

import bisect
import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..types import Epoch, SatObs, TruthTrack
from . import frames
from .clock_analysis import _pick_doppler_sign, doppler_clock_drift
from .detectors import DETECTORS, DetectorConfig
from .wls import WlsConfig, solve_epoch

SigKey = Tuple[int, int, int, str]  # (constellation, prn, freq_hz, code_type)


def _sig_key(o: SatObs) -> SigKey:
    return (int(o.constellation), int(o.prn), int(round(o.frequency_hz)), o.code_type)


def _range_rate_mps(o: SatObs, sign: float) -> float:
    if not (math.isfinite(o.doppler_hz) and o.frequency_hz > 0):
        return float("nan")
    return sign * frames.C_LIGHT / o.frequency_hz * o.doppler_hz


def _carrier_range_m(o: SatObs) -> float:
    if not (math.isfinite(o.phase_cycle) and o.frequency_hz > 0) or o.phase_cycle == 0.0:
        return float("nan")
    return frames.C_LIGHT / o.frequency_hz * o.phase_cycle


def _truth_kinematics(truth: TruthTrack) -> Tuple[List[float], List[float], List[float]]:
    """Horizontal speed [m/s] and its rate [m/s^2] at every truth sample
    (central differences in a local ENU frame; endpoints get one-sided)."""
    ts = [s.t_sec for s in truth.samples]
    n = len(ts)
    speed = [float("nan")] * n
    for i in range(n):
        j0, j1 = max(0, i - 1), min(n - 1, i + 1)
        dt = ts[j1] - ts[j0]
        if j1 == j0 or dt <= 0 or dt > 10.0:
            continue
        origin = np.asarray(truth.samples[i].ecef, float)
        rot = frames.enu_rotation_matrix(*frames.ecef_to_lla(*origin)[:2])
        d = rot @ (np.asarray(truth.samples[j1].ecef, float)
                   - np.asarray(truth.samples[j0].ecef, float))
        speed[i] = float(math.hypot(d[0], d[1]) / dt)
    accel = [float("nan")] * n
    for i in range(n):
        j0, j1 = max(0, i - 1), min(n - 1, i + 1)
        dt = ts[j1] - ts[j0]
        if j1 == j0 or dt <= 0 or dt > 10.0:
            continue
        if math.isfinite(speed[j1]) and math.isfinite(speed[j0]):
            accel[i] = (speed[j1] - speed[j0]) / dt
    return ts, speed, accel


def _nearest_idx(ts: Sequence[float], t: float) -> Optional[int]:
    if not ts:
        return None
    i = bisect.bisect_left(ts, t)
    cands = [j for j in (i - 1, i) if 0 <= j < len(ts)]
    return min(cands, key=lambda j: abs(ts[j] - t)) if cands else None


def _agg(vals: List[float]) -> Tuple[float, float, float]:
    """(median|v|, max|v|, rms) over the finite entries; NaN when empty."""
    a = np.asarray([v for v in vals if math.isfinite(v)], float)
    if a.size == 0:
        return float("nan"), float("nan"), float("nan")
    aa = np.abs(a)
    return float(np.median(aa)), float(np.max(aa)), float(np.sqrt(np.mean(a * a)))


def extract_factor_epochs(
    epochs: Sequence[Epoch],
    truth: TruthTrack,
    wls_cfg: Optional[WlsConfig] = None,
    dcfg: Optional[DetectorConfig] = None,
    detector: str = "combined",
    match_tolerance_sec: float = 0.5,
    max_dt_s: float = 5.0,
) -> List[Dict[str, float]]:
    """One dict of candidate factors per solved-and-truth-matched epoch."""
    wls_cfg = wls_cfg or WlsConfig()
    dcfg = dcfg or DetectorConfig()
    det = DETECTORS[detector]

    solved = []
    for ep in epochs:
        reject = det(ep.obs, ep.t_sec, wls_cfg, dcfg)
        s = solve_epoch(ep.obs, ep.t_sec, wls_cfg, reject=reject)
        if s.converged and all(math.isfinite(v) for v in s.ecef):
            solved.append((ep, s, reject))
    solved.sort(key=lambda p: p[1].t_sec)
    if len(solved) < 4:
        return []

    dop_sign = _pick_doppler_sign(solved)
    truth_ts, truth_speed, truth_accel = _truth_kinematics(truth)

    rows: List[Dict[str, float]] = []
    prev: Optional[Tuple[Epoch, object, set]] = None
    prev_clk: Optional[float] = None
    prev_drift: Optional[float] = None
    prev_t: Optional[float] = None
    prev_drift_dop: Optional[float] = None

    for ep, sol, reject in solved:
        # --- truth match -> target ---------------------------------------
        ts = truth.nearest(ep.t_sec, match_tolerance_sec)
        if ts is None:
            prev = (ep, sol, reject)
            continue
        origin = np.asarray(ts.ecef, float)
        rot = frames.enu_rotation_matrix(*frames.ecef_to_lla(*origin)[:2])
        enu = rot @ (np.asarray(sol.ecef, float) - origin)
        row: Dict[str, float] = {
            "t_sec": ep.t_sec,
            "h_err_m": float(math.hypot(enu[0], enu[1])),
            "v_err_m": float(enu[2]),
        }

        # --- dynamics (truth) ---------------------------------------------
        ti = _nearest_idx(truth_ts, ep.t_sec)
        row["speed_mps"] = truth_speed[ti] if ti is not None else float("nan")
        row["accel_mps2"] = truth_accel[ti] if ti is not None else float("nan")

        # --- geometry -------------------------------------------------------
        row.update({
            "n_used": sol.n_used, "n_available": sol.n_available,
            "n_rejected": len(reject),
            "hdop": sol.hdop, "vdop": sol.vdop, "pdop": sol.pdop,
            "gdop": sol.gdop, "tdop": sol.tdop,
        })
        elevs = [sol.elevation_deg[k] for k in sol.used_keys if k in sol.elevation_deg]
        elevs = [e for e in elevs if math.isfinite(e)]
        row["elev_min_deg"] = float(min(elevs)) if elevs else float("nan")
        row["elev_mean_deg"] = float(np.mean(elevs)) if elevs else float("nan")
        row["n_low_elev"] = int(sum(1 for e in elevs if e < 20.0))
        used_sat = {(k[0], k[1]) for k in sol.used_keys}
        for cid, cname in ((0, "gps"), (1, "gal"), (2, "bds"), (4, "qzs")):
            row[f"n_{cname}"] = sum(1 for c, _ in used_sat if c == cid)

        # --- signal (C/N0) ----------------------------------------------------
        by_sat_cn0: Dict[Tuple[int, int], float] = {}
        for o in ep.obs:
            if math.isfinite(o.cn0_dbhz):
                k = (o.constellation, o.prn)
                by_sat_cn0[k] = max(by_sat_cn0.get(k, -1.0), o.cn0_dbhz)
        cn0_all = np.asarray(list(by_sat_cn0.values()), float)
        cn0_used = np.asarray([v for k, v in by_sat_cn0.items() if k in used_sat], float)
        row["cn0_mean_used"] = float(np.mean(cn0_used)) if cn0_used.size else float("nan")
        row["cn0_min_used"] = float(np.min(cn0_used)) if cn0_used.size else float("nan")
        row["cn0_std_used"] = float(np.std(cn0_used)) if cn0_used.size else float("nan")
        row["cn0_mean_all"] = float(np.mean(cn0_all)) if cn0_all.size else float("nan")
        row["cn0_frac_below30"] = (float(np.mean(cn0_all < 30.0))
                                   if cn0_all.size else float("nan"))

        # --- consistency (post-fit) ----------------------------------------
        res = np.asarray([sol.residuals_m[k] for k in sol.used_keys
                          if k in sol.residuals_m and math.isfinite(sol.residuals_m[k])], float)
        row["sigma0_hat"] = sol.sigma0_hat
        row["resid_rms_m"] = float(np.sqrt(np.mean(res * res))) if res.size else float("nan")
        row["resid_max_m"] = float(np.max(np.abs(res))) if res.size else float("nan")
        row["clk_sigma_m"] = sol.clock_sigma_m

        # --- clock -----------------------------------------------------------
        row["clk_m"] = sol.clock_bias_m
        dt = (sol.t_sec - prev_t) if prev_t is not None else float("nan")
        drift = float("nan")
        inst = float("nan")
        if math.isfinite(dt) and 0 < dt <= max_dt_s and prev_clk is not None:
            drift = (sol.clock_bias_m - prev_clk) / dt
            if prev_drift is not None and math.isfinite(prev_drift):
                inst = abs(sol.clock_bias_m - (prev_clk + prev_drift * dt))
        row["clk_drift_mps"] = drift
        row["clk_inst_m"] = inst

        dop = doppler_clock_drift(ep, sol.ecef, reject=reject, sign=dop_sign)
        row["drift_dop_mps"] = dop["drift_mps"]
        row["dop_resid_rms_mps"] = dop["resid_rms_mps"]
        row["dop_n_rejected"] = dop["n_rejected"]
        row["speed_dop_mps"] = dop["speed_mps"]
        dchg = float("nan")
        if (math.isfinite(dop["drift_mps"]) and prev_drift_dop is not None
                and math.isfinite(prev_drift_dop) and math.isfinite(dt) and 0 < dt <= max_dt_s):
            dchg = abs(dop["drift_mps"] - prev_drift_dop) / dt
        row["drift_dop_change_mps"] = dchg

        # --- code quality: DCD & CCD rates vs previous epoch ----------------
        dcd_vals: List[float] = []
        ccd_vals: List[float] = []
        n_phase = 0
        n_loi = 0
        if prev is not None and math.isfinite(dt) and 0 < dt <= max_dt_s:
            prev_ep = prev[0]
            prev_by_key = {_sig_key(o): o for o in prev_ep.obs}
            for o in ep.obs:
                if (o.constellation, o.prn) not in used_sat:
                    continue
                po = prev_by_key.get(_sig_key(o))
                if po is None:
                    continue
                if math.isfinite(o.pseudorange_m) and math.isfinite(po.pseudorange_m):
                    dpr = (o.corrected_pseudorange_m - po.corrected_pseudorange_m) / dt
                    rr0, rr1 = _range_rate_mps(po, dop_sign), _range_rate_mps(o, dop_sign)
                    if math.isfinite(rr0) and math.isfinite(rr1):
                        # trapezoidal Doppler so range acceleration cancels;
                        # rx/sv clock drift is common to both terms and drops out
                        dcd_vals.append(dpr - 0.5 * (rr0 + rr1))
                cr0, cr1 = _carrier_range_m(po), _carrier_range_m(o)
                if math.isfinite(cr0) and math.isfinite(cr1):
                    n_phase += 1
                    if o.loi or po.loi:
                        n_loi += 1
                    else:
                        ccd_vals.append(((o.pseudorange_m - cr1)
                                         - (po.pseudorange_m - cr0)) / dt)
        row["dcd_med_mps"], row["dcd_max_mps"], row["dcd_rms_mps"] = _agg(dcd_vals)
        row["ccd_med_mps"], row["ccd_max_mps"], row["ccd_rms_mps"] = _agg(ccd_vals)
        row["frac_loi"] = (n_loi / n_phase) if n_phase else float("nan")
        row["dt_s"] = dt

        rows.append(row)
        prev = (ep, sol, reject)
        prev_clk, prev_t = sol.clock_bias_m, sol.t_sec
        prev_drift = drift if math.isfinite(drift) else prev_drift
        prev_drift_dop = dop["drift_mps"] if math.isfinite(dop["drift_mps"]) else prev_drift_dop

    return rows


FACTOR_COLUMNS: List[str] = [
    "t_sec", "dt_s", "h_err_m", "v_err_m",
    "speed_mps", "accel_mps2", "speed_dop_mps",
    "n_used", "n_available", "n_rejected",
    "hdop", "vdop", "pdop", "gdop", "tdop",
    "elev_min_deg", "elev_mean_deg", "n_low_elev",
    "n_gps", "n_gal", "n_bds", "n_qzs",
    "cn0_mean_used", "cn0_min_used", "cn0_std_used", "cn0_mean_all", "cn0_frac_below30",
    "sigma0_hat", "resid_rms_m", "resid_max_m", "clk_sigma_m",
    "clk_m", "clk_drift_mps", "clk_inst_m",
    "drift_dop_mps", "drift_dop_change_mps", "dop_resid_rms_mps", "dop_n_rejected",
    "dcd_med_mps", "dcd_max_mps", "dcd_rms_mps",
    "ccd_med_mps", "ccd_max_mps", "ccd_rms_mps", "frac_loi",
]

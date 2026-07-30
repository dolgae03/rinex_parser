"""Measurement-quality calibration from truth.

With a known truth position we can measure *the measurements themselves*, not
just the position solution. For every observation we form the prefit residual at
truth and remove the per-epoch, per-constellation median (which absorbs the
receiver clock + inter-system bias). What remains is the per-signal measurement
error — noise + multipath + NLOS + unmodelled elevation-dependent atmosphere.

Aggregating those residuals answers, from data rather than assumption:
 - the empirical per-constellation sigma (and hence a data-driven weight scale),
 - how sigma varies with C/N0 and elevation (validates the weight model),
 - the gross-outlier rate per constellation.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..types import Epoch, TruthTrack, constellation_name
from . import frames
from .wls import WlsConfig, _pick_signals


@dataclass
class MeasResidual:
    t_sec: float
    constellation: int
    prn: int
    freq_hz: float
    cn0_dbhz: float
    elevation_deg: float
    azimuth_deg: float
    residual_m: float       # detrended (per-epoch, per-constellation median removed)
    raw_offset_m: float     # L = corrected_pr - range(truth), before detrend


def measurement_residuals(
    epochs: Sequence[Epoch],
    truth: TruthTrack,
    wls_cfg: Optional[WlsConfig] = None,
    match_tolerance_sec: float = 0.5,
) -> List[MeasResidual]:
    """Prefit residuals at truth, detrended per (epoch, constellation)."""
    wls_cfg = wls_cfg or WlsConfig()
    out: List[MeasResidual] = []
    for ep in epochs:
        ts = truth.nearest(ep.t_sec, match_tolerance_sec)
        if ts is None:
            continue
        rx = np.asarray(ts.ecef, float)
        lat, lon, _ = frames.ecef_to_lla(*rx)
        R = frames.enu_rotation_matrix(lat, lon)
        recs: List[Tuple] = []
        by_cons: Dict[int, List[float]] = defaultdict(list)
        for o in _pick_signals(ep.obs, set()):
            sv = np.asarray(o.sv_pos, float)
            if wls_cfg.apply_sagnac:
                sv = frames.sagnac_corrected_sv(sv, float(np.linalg.norm(sv - rx)) / frames.C_LIGHT)
            diff = sv - rx
            rng = float(np.linalg.norm(diff))
            enu = R @ diff
            el = math.degrees(math.atan2(enu[2], math.hypot(enu[0], enu[1])))
            az = math.degrees(math.atan2(enu[0], enu[1])) % 360.0
            iono = o.iono_delay_m if (wls_cfg.apply_iono and math.isfinite(o.iono_delay_m)) else 0.0
            L = o.corrected_pseudorange_m - iono - rng
            recs.append((o, el, az, L))
            by_cons[o.constellation].append(L)
        med = {c: float(np.median(v)) for c, v in by_cons.items()}
        for o, el, az, L in recs:
            out.append(MeasResidual(
                t_sec=ep.t_sec, constellation=o.constellation, prn=o.prn,
                freq_hz=o.frequency_hz, cn0_dbhz=o.cn0_dbhz, elevation_deg=el,
                azimuth_deg=az, residual_m=L - med[o.constellation], raw_offset_m=L,
            ))
    return out


def outlier_catalog(residuals: Sequence[MeasResidual], k: float = 5.0,
                    mad_floor_m: float = 6.0) -> Dict:
    """Flag truth-referenced blunders and return a catalog + summary.

    ``residual_m`` is already de-medianed per (epoch, constellation), so a single
    global robust scale (MAD over all residuals) gives each observation a z-score;
    |z| > k marks a blunder. Returns per-record rows (for the sky plot) and
    aggregate counts by constellation and by satellite.
    """
    if not residuals:
        return {"rows": [], "n": 0, "n_outliers": 0, "scale_m": float("nan"),
                "by_constellation": {}, "worst_satellites": []}
    vals = np.array([r.residual_m for r in residuals], float)
    med = float(np.median(vals))
    scale = max(1.4826 * float(np.median(np.abs(vals - med))), mad_floor_m)
    rows, by_c, by_sat = [], defaultdict(lambda: [0, 0]), defaultdict(lambda: [0, 0])
    for r in residuals:
        z = abs(r.residual_m - med) / scale
        is_out = z > k
        rows.append({
            "t_sec": r.t_sec, "constellation": constellation_name(r.constellation),
            "prn": r.prn, "freq_mhz": round(r.freq_hz / 1e6, 2),
            "elevation_deg": round(r.elevation_deg, 2), "azimuth_deg": round(r.azimuth_deg, 2),
            "cn0_dbhz": round(r.cn0_dbhz, 1) if r.cn0_dbhz == r.cn0_dbhz else "",
            "residual_m": round(r.residual_m, 3), "z": round(z, 2),
            "is_outlier": int(is_out),
        })
        cname = constellation_name(r.constellation)
        by_c[cname][0] += 1
        by_c[cname][1] += int(is_out)
        key = f"{cname}-{r.prn}"
        by_sat[key][0] += 1
        by_sat[key][1] += int(is_out)
    by_constellation = {c: {"n": n, "n_outliers": no, "rate": round(no / n, 4)}
                        for c, (n, no) in sorted(by_c.items())}
    worst = sorted(({"satellite": s, "n": n, "n_outliers": no, "rate": round(no / n, 4)}
                    for s, (n, no) in by_sat.items() if no > 0),
                   key=lambda d: (-d["n_outliers"], -d["rate"]))[:15]
    n_out = sum(r["is_outlier"] for r in rows)
    return {"rows": rows, "n": len(rows), "n_outliers": n_out,
            "scale_m": round(scale, 3), "threshold_k": k,
            "by_constellation": by_constellation, "worst_satellites": worst}


def _robust_stats(values: np.ndarray, k: float = 5.0) -> Dict[str, float]:
    """Median/MAD plus an outlier-excluded std and outlier rate."""
    a = np.asarray(values, float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return {"n": 0}
    m = float(np.median(a))
    mad = 1.4826 * float(np.median(np.abs(a - m)))
    if mad > 0:
        keep = np.abs(a - m) <= k * mad
    else:
        keep = np.ones(a.shape, bool)
    clean = a[keep]
    return {
        "n": int(a.size),
        "median_m": m,
        "mad_sigma_m": mad,
        "clean_std_m": float(np.std(clean)) if clean.size else float("nan"),
        "rms_m": float(np.sqrt(np.mean(a * a))),
        "outlier_rate": float(1.0 - keep.mean()),
    }


@dataclass
class CalibrationResult:
    by_constellation: Dict[str, Dict[str, float]] = field(default_factory=dict)
    sigma_scale_by_constellation: Dict[str, float] = field(default_factory=dict)
    by_cn0_bin: List[Dict[str, float]] = field(default_factory=list)
    by_elevation_bin: List[Dict[str, float]] = field(default_factory=list)
    reference_constellation: str = ""
    n_residuals: int = 0

    def to_dict(self) -> Dict:
        return {
            "reference_constellation": self.reference_constellation,
            "n_residuals": self.n_residuals,
            "by_constellation": self.by_constellation,
            "sigma_scale_by_constellation": self.sigma_scale_by_constellation,
            "by_cn0_bin": self.by_cn0_bin,
            "by_elevation_bin": self.by_elevation_bin,
        }


def calibrate(
    residuals: Sequence[MeasResidual],
    reference: str = "GPS",
    cn0_edges: Sequence[float] = (0, 25, 30, 35, 40, 45, 99),
    elev_edges: Sequence[float] = (0, 15, 30, 45, 60, 90),
    scale_metric: str = "clean_std_m",
) -> CalibrationResult:
    """Empirical per-constellation / per-C/N0 / per-elevation noise from residuals.

    ``sigma_scale_by_constellation`` = (constellation sigma) / (reference sigma),
    i.e. the data-driven version of ``WeightConfig.sigma_scale_by_constellation``.
    """
    res = CalibrationResult(reference_constellation=reference, n_residuals=len(residuals))

    by_c: Dict[int, List[float]] = defaultdict(list)
    for r in residuals:
        by_c[r.constellation].append(r.residual_m)
    for c, vals in sorted(by_c.items()):
        res.by_constellation[constellation_name(c)] = _robust_stats(np.array(vals))

    ref_sigma = res.by_constellation.get(reference, {}).get(scale_metric, float("nan"))
    for name, stats in res.by_constellation.items():
        s = stats.get(scale_metric, float("nan"))
        res.sigma_scale_by_constellation[name] = (
            round(s / ref_sigma, 2) if (ref_sigma and math.isfinite(ref_sigma) and s == s) else float("nan")
        )

    def _bin(values_key, edges):
        rows = []
        arr = np.array([(getattr(r, values_key), r.residual_m) for r in residuals], dtype=float)
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = arr[(arr[:, 0] >= lo) & (arr[:, 0] < hi)]
            stats = _robust_stats(sel[:, 1]) if sel.size else {"n": 0}
            rows.append({"bin_lo": lo, "bin_hi": hi, **stats})
        return rows

    res.by_cn0_bin = _bin("cn0_dbhz", cn0_edges)
    res.by_elevation_bin = _bin("elevation_deg", elev_edges)
    return res

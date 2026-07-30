"""Outlier detectors. Each returns the set of satellites to reject for one epoch.

The lab's chosen strategy is *multi-detector comparison*: run several detectors,
remove what each flags, re-solve WLS, and compare the improvement against truth.
Detectors never mutate observations and never touch truth; they only *select*.

Rejection unit
--------------
Detectors reject whole satellites ``(constellation, prn)``, not individual
signals, because a bad broadcast ephemeris corrupts every frequency of a
satellite (and the WLS uses one signal per satellite anyway).

Robustness
----------
A single gross blunder (e.g. a broken ephemeris off by ~10^6 m) masks itself in
an ordinary least-squares fix: it drags the solution so far that no single
residual stands out. So the residual/elevation detectors judge observations from
a *robust* (Huber-IRLS) position that down-weights blunders, then flag satellites
whose residual deviates from the robust centre by more than ``k`` robust scales
(a data-driven MAD, so unmodelled atmosphere does not trigger false alarms).

A detector is a callable ``detect(obs, t_sec, wls_cfg, dcfg) -> set[SatId]``.
Register new ones in ``DETECTORS`` and they appear in the experiment automatically.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, replace
from statistics import median
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from ..types import SatObs
from . import frames
from .wls import WlsConfig, _pick_signals, solve_epoch

SatId = Tuple[int, int]  # (constellation, prn) — the rejection unit


@dataclass
class DetectorConfig:
    min_cn0_dbhz: float = 30.0
    min_elevation_deg: float = 10.0
    residual_mad_k: float = 5.0     # flag |residual - median| > k * robust_scale
    mad_floor_m: float = 6.0        # floor on the robust scale [m]
    min_obs_for_residual: int = 6


def _robust_solve(obs, t_sec, wls_cfg: WlsConfig, reject: Optional[Set[SatId]] = None):
    # Detection uses de-median residuals where common atmosphere cancels, so the
    # per-satellite atmospheric refinement is skipped here (keeps the many
    # detector solves single-pass and fast). Atmosphere is applied in the final
    # metric solves and in calibration.
    cfg = replace(wls_cfg, robust_huber=True, apply_tropo=False, apply_iono=False)
    return solve_epoch(obs, t_sec, cfg, reject=reject)


def _usable(sol) -> bool:
    # Detection only needs a finite position; the Huber IRLS 'converged' flag can
    # stay False near the solution even when the estimate is good.
    return all(math.isfinite(v) for v in sol.ecef)


def _residuals_by_sat(sol) -> Dict[SatId, float]:
    return {(k[0], k[1]): v for k, v in sol.residuals_m.items()}


def _flag_by_mad(res: Dict[SatId, float], dcfg: DetectorConfig) -> Set[SatId]:
    if len(res) < dcfg.min_obs_for_residual:
        return set()
    vals = list(res.values())
    med = median(vals)
    scale = max(median(abs(v - med) for v in vals) * 1.4826, dcfg.mad_floor_m)
    thr = dcfg.residual_mad_k * scale
    return {sat for sat, v in res.items() if abs(v - med) > thr}


# ---------------------------------------------------------------------------
# Detectors
# ---------------------------------------------------------------------------

def detect_none(obs, t_sec, wls_cfg, dcfg) -> Set[SatId]:
    return set()


def detect_cn0(obs: Sequence[SatObs], t_sec, wls_cfg, dcfg) -> Set[SatId]:
    # gate on the signal actually used per satellite (strongest C/N0)
    return {(o.constellation, o.prn) for o in _pick_signals(obs, set())
            if math.isfinite(o.cn0_dbhz) and o.cn0_dbhz < dcfg.min_cn0_dbhz}


def detect_elevation(obs: Sequence[SatObs], t_sec, wls_cfg, dcfg) -> Set[SatId]:
    sol = _robust_solve(obs, t_sec, wls_cfg)
    if not _usable(sol):
        return set()
    pos = np.asarray(sol.ecef, float)
    rej: Set[SatId] = set()
    for o in _pick_signals(obs, set()):
        elev, _ = frames.elevation_azimuth_deg(pos, np.asarray(o.sv_pos, float))
        if elev < dcfg.min_elevation_deg:
            rej.add((o.constellation, o.prn))
    return rej


def detect_residual(obs, t_sec, wls_cfg, dcfg) -> Set[SatId]:
    """Robust-fit residual snooping across all constellations (global scale)."""
    sol = _robust_solve(obs, t_sec, wls_cfg)
    if not _usable(sol):
        return set()
    return _flag_by_mad(_residuals_by_sat(sol), dcfg)


def detect_residual_per_constellation(obs, t_sec, wls_cfg, dcfg) -> Set[SatId]:
    """Robust residual snooping run independently within each constellation (>=5 sats)."""
    picked = _pick_signals(obs, set())
    by_cons: Dict[int, int] = defaultdict(int)
    for o in picked:
        by_cons[o.constellation] += 1
    rej: Set[SatId] = set()
    for cons, n in by_cons.items():
        if n < 5:
            continue
        drop = {(o.constellation, o.prn) for o in picked if o.constellation != cons}
        sol = _robust_solve(obs, t_sec, wls_cfg, reject=drop)
        if _usable(sol):
            rej |= _flag_by_mad(_residuals_by_sat(sol), dcfg)
    return rej


def detect_combined(obs, t_sec, wls_cfg, dcfg) -> Set[SatId]:
    """C/N0 + elevation gating, unioned with robust residual snooping."""
    return (detect_cn0(obs, t_sec, wls_cfg, dcfg)
            | detect_elevation(obs, t_sec, wls_cfg, dcfg)
            | detect_residual(obs, t_sec, wls_cfg, dcfg))


Detector = Callable[..., Set[SatId]]

DETECTORS: Dict[str, Detector] = {
    "baseline": detect_none,
    "cn0": detect_cn0,
    "elevation": detect_elevation,
    "residual": detect_residual,
    "residual_per_constellation": detect_residual_per_constellation,
    "combined": detect_combined,
}

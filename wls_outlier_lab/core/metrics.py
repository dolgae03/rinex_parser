"""Positioning-error metrics against a truth track. Pure logic.

Errors are expressed in the local ENU frame at each epoch's truth point, which
keeps horizontal/vertical separation meaningful over a moving trajectory.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..types import TruthTrack, WlsSolution
from . import frames


@dataclass
class EpochError:
    t_sec: float
    east_m: float
    north_m: float
    up_m: float
    n_used: int
    hdop: float
    pdop: float

    @property
    def horizontal_m(self) -> float:
        return math.hypot(self.east_m, self.north_m)

    @property
    def d3_m(self) -> float:
        return math.sqrt(self.east_m ** 2 + self.north_m ** 2 + self.up_m ** 2)


@dataclass
class ErrorMetrics:
    label: str = ""
    n_epochs_total: int = 0
    n_epochs_solved: int = 0
    n_epochs_matched: int = 0
    availability_pct: float = float("nan")
    mean_sats_used: float = float("nan")
    horizontal_rmse_m: float = float("nan")
    horizontal_mean_m: float = float("nan")
    horizontal_p50_m: float = float("nan")
    horizontal_p95_m: float = float("nan")
    horizontal_max_m: float = float("nan")
    cep50_m: float = float("nan")
    cep95_m: float = float("nan")
    drms_2d_m: float = float("nan")
    vertical_rmse_m: float = float("nan")
    vertical_mean_m: float = float("nan")
    d3_rmse_m: float = float("nan")
    mean_hdop: float = float("nan")
    per_epoch: List[EpochError] = field(default_factory=list)

    def to_summary(self) -> Dict[str, float]:
        d = {k: getattr(self, k) for k in (
            "label", "n_epochs_total", "n_epochs_solved", "n_epochs_matched",
            "availability_pct", "mean_sats_used", "horizontal_rmse_m",
            "horizontal_mean_m", "horizontal_p50_m", "horizontal_p95_m",
            "horizontal_max_m", "cep50_m", "cep95_m", "drms_2d_m",
            "vertical_rmse_m", "vertical_mean_m", "d3_rmse_m", "mean_hdop",
        )}
        return d


def _rms(values: np.ndarray) -> float:
    return float(math.sqrt(np.mean(values ** 2))) if len(values) else float("nan")


def compute_metrics(
    solutions: Sequence[WlsSolution],
    truth: TruthTrack,
    label: str = "",
    match_tolerance_sec: float = 0.5,
    n_epochs_total: Optional[int] = None,
) -> ErrorMetrics:
    """Match each solved epoch to truth and aggregate error statistics."""
    per_epoch: List[EpochError] = []
    n_solved = 0
    sats_used: List[int] = []
    for sol in solutions:
        if not sol.converged or not all(math.isfinite(v) for v in sol.ecef):
            continue
        n_solved += 1
        sats_used.append(sol.n_used)
        ts = truth.nearest(sol.t_sec, match_tolerance_sec)
        if ts is None:
            continue
        origin = np.asarray(ts.ecef, float)
        enu = frames.enu_rotation_matrix(*frames.ecef_to_lla(*origin)[:2]) @ (
            np.asarray(sol.ecef, float) - origin
        )
        per_epoch.append(
            EpochError(sol.t_sec, float(enu[0]), float(enu[1]), float(enu[2]),
                       sol.n_used, sol.hdop, sol.pdop)
        )

    m = ErrorMetrics(label=label)
    m.per_epoch = per_epoch
    m.n_epochs_total = n_epochs_total if n_epochs_total is not None else len(solutions)
    m.n_epochs_solved = n_solved
    m.n_epochs_matched = len(per_epoch)
    if m.n_epochs_total:
        m.availability_pct = 100.0 * n_solved / m.n_epochs_total
    if sats_used:
        m.mean_sats_used = float(np.mean(sats_used))
    if not per_epoch:
        return m

    horiz = np.array([e.horizontal_m for e in per_epoch])
    east = np.array([e.east_m for e in per_epoch])
    north = np.array([e.north_m for e in per_epoch])
    up = np.array([e.up_m for e in per_epoch])
    d3 = np.array([e.d3_m for e in per_epoch])
    hdops = np.array([e.hdop for e in per_epoch if math.isfinite(e.hdop)])

    m.horizontal_rmse_m = _rms(horiz)
    m.horizontal_mean_m = float(np.mean(horiz))
    m.horizontal_p50_m = float(np.percentile(horiz, 50))
    m.horizontal_p95_m = float(np.percentile(horiz, 95))
    m.horizontal_max_m = float(np.max(horiz))
    m.cep50_m = float(np.percentile(horiz, 50))
    m.cep95_m = float(np.percentile(horiz, 95))
    # 2DRMS = 2 * sqrt(var_E + var_N) about the error mean
    m.drms_2d_m = 2.0 * float(math.sqrt(np.var(east) + np.var(north)))
    m.vertical_rmse_m = _rms(up)
    m.vertical_mean_m = float(np.mean(up))
    m.d3_rmse_m = _rms(d3)
    if len(hdops):
        m.mean_hdop = float(np.mean(hdops))
    return m

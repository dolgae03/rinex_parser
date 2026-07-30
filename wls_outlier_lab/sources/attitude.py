"""Attitude source — v1 hook only.

Phone attitude was chosen as a *future* quality covariate: the interface is
defined now so the WLS/outlier study proceeds unblocked, and a real parser
(phone ``OrientationDeg``/IMU, or Novatel ``INSPVA`` vehicle attitude) can drop
in later without touching the core. ``AttitudeSeries`` is the shape the future
covariate analysis will consume: time-indexed roll/pitch/yaw in degrees.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from .base import AttitudeSource


@dataclass(frozen=True)
class AttitudeSample:
    t_sec: float
    roll_deg: float
    pitch_deg: float
    yaw_deg: float


@dataclass
class AttitudeSeries:
    name: str
    samples: List[AttitudeSample] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.samples = sorted(self.samples, key=lambda s: s.t_sec)
        self._times = [s.t_sec for s in self.samples]

    def nearest(self, t_sec: float, tolerance_sec: float = 0.5) -> Optional[AttitudeSample]:
        if not self.samples:
            return None
        idx = bisect.bisect_left(self._times, t_sec)
        cands = []
        if idx < len(self.samples):
            cands.append(self.samples[idx])
        if idx > 0:
            cands.append(self.samples[idx - 1])
        best = min(cands, key=lambda s: abs(s.t_sec - t_sec))
        return best if abs(best.t_sec - t_sec) <= tolerance_sec else None


class NullAttitudeSource(AttitudeSource):
    """Default: no attitude available (v1)."""

    def load_series(self) -> Optional[AttitudeSeries]:
        return None

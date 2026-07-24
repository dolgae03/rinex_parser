"""Core data model for the WLS outlier lab.

These dataclasses are deliberately I/O-free. Data-access adapters in
``wls_outlier_lab.sources`` produce them; the business logic in
``wls_outlier_lab.core`` consumes them. Nothing here reads a file or knows
where the data came from.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

Vec3 = Tuple[float, float, float]

# Constellation integer encoding used by the canonical measurement TSV/CSV
# (matches gnss_txt_parser.data_class.Constellation).
CONSTELLATION_NAMES: Dict[int, str] = {
    0: "GPS",
    1: "GALILEO",
    2: "BEIDOU",
    3: "GLONASS",
    4: "QZSS",
    5: "SBAS",
    6: "IRNSS",
    7: "UNKNOWN",
}


def constellation_name(code: int) -> str:
    return CONSTELLATION_NAMES.get(int(code), f"CONS_{int(code)}")


@dataclass(frozen=True)
class SatObs:
    """One satellite-signal measurement at one epoch.

    ``pseudorange_m`` and ``sv_clock_bias_m`` follow the convention verified
    against the MATLAB endpoint output: the clock-corrected observation is
    ``pseudorange_m + sv_clock_bias_m + pr_correction_m`` (all metres), which
    equals ``geometric_range + c * rx_clock_bias (+ inter-system bias)``.
    """

    constellation: int
    prn: int
    frequency_hz: float
    code_type: str
    pseudorange_m: float
    sv_pos: Vec3
    sv_clock_bias_m: float
    sv_vel: Optional[Vec3] = None
    pr_correction_m: float = 0.0
    doppler_hz: float = float("nan")
    cn0_dbhz: float = float("nan")
    loi: bool = False
    iono_delay_m: float = 0.0

    @property
    def key(self) -> Tuple[int, int, int]:
        """Stable identity of this signal: (constellation, prn, freq_hz)."""
        return (int(self.constellation), int(self.prn), int(round(self.frequency_hz)))

    @property
    def corrected_pseudorange_m(self) -> float:
        return self.pseudorange_m + self.sv_clock_bias_m + self.pr_correction_m


@dataclass
class Epoch:
    """All observations sharing one measurement time."""

    t_sec: float
    obs: List[SatObs]
    gps_week: int = -1
    tow_sec: float = float("nan")

    def constellations(self) -> List[int]:
        return sorted({o.constellation for o in self.obs})


@dataclass(frozen=True)
class TruthSample:
    t_sec: float
    ecef: Vec3


@dataclass
class TruthTrack:
    """Time-indexed ground-truth trajectory in ECEF metres."""

    name: str
    samples: List[TruthSample]
    source_type: str = "unknown"

    def __post_init__(self) -> None:
        self.samples = sorted(self.samples, key=lambda s: s.t_sec)
        self._times = [s.t_sec for s in self.samples]

    def nearest(self, t_sec: float, tolerance_sec: float = 0.5) -> Optional[TruthSample]:
        """Nearest truth sample within ``tolerance_sec`` (None if none/out of range)."""
        if not self.samples:
            return None
        idx = bisect.bisect_left(self._times, t_sec)
        candidates = []
        if idx < len(self.samples):
            candidates.append(self.samples[idx])
        if idx > 0:
            candidates.append(self.samples[idx - 1])
        best = min(candidates, key=lambda s: abs(s.t_sec - t_sec))
        if abs(best.t_sec - t_sec) > tolerance_sec:
            return None
        return best


@dataclass
class WlsSolution:
    """Single-epoch weighted-least-squares fix and its diagnostics."""

    t_sec: float
    ecef: Vec3
    clock_bias_m: float
    isb_m: Dict[int, float] = field(default_factory=dict)
    n_used: int = 0
    n_available: int = 0
    residuals_m: Dict[Tuple[int, int, int], float] = field(default_factory=dict)
    std_residuals: Dict[Tuple[int, int, int], float] = field(default_factory=dict)
    used_keys: List[Tuple[int, int, int]] = field(default_factory=list)
    elevation_deg: Dict[Tuple[int, int, int], float] = field(default_factory=dict)
    gdop: float = float("nan")
    pdop: float = float("nan")
    hdop: float = float("nan")
    vdop: float = float("nan")
    tdop: float = float("nan")
    sigma0_hat: float = float("nan")  # a-posteriori standard deviation of unit weight
    converged: bool = False
    iterations: int = 0
    reason: str = ""  # populated when the epoch could not be solved

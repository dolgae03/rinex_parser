from __future__ import annotations

import csv
import io
import math
from dataclasses import asdict, dataclass, field, is_dataclass
from datetime import datetime, timedelta, timezone
from enum import Enum
from typing import Any, Dict, List, Union


class Constellation(Enum):
    GPS = 0
    GALILEO = 1
    BEIDOU = 2
    GLONASS = 3
    QZSS = 4
    SBAS = 5
    IRNSS = 6
    UNKNOWN = 7


_CONS_ALIAS = {
    "G": Constellation.GPS,
    "E": Constellation.GALILEO,
    "C": Constellation.BEIDOU,
    "R": Constellation.GLONASS,
    "J": Constellation.QZSS,
    "S": Constellation.SBAS,
    "I": Constellation.IRNSS,
    "GPS": Constellation.GPS,
    "GALILEO": Constellation.GALILEO,
    "BEIDOU": Constellation.BEIDOU,
    "BDS": Constellation.BEIDOU,
    "GLONASS": Constellation.GLONASS,
    "QZSS": Constellation.QZSS,
    "SBAS": Constellation.SBAS,
    "IRNSS": Constellation.IRNSS,
    "NAVIC": Constellation.IRNSS,
    1: Constellation.GPS,
    2: Constellation.SBAS,
    3: Constellation.GLONASS,
    4: Constellation.QZSS,
    5: Constellation.BEIDOU,
    6: Constellation.GALILEO,
    7: Constellation.IRNSS,
    0: Constellation.GPS,
}


def to_constellation(x: Union[str, int, Constellation, None]) -> int:
    if isinstance(x, Constellation):
        return x.value
    if isinstance(x, str):
        return _CONS_ALIAS.get(x.strip().upper(), Constellation.UNKNOWN).value
    if isinstance(x, int):
        return _CONS_ALIAS.get(int(x), Constellation.UNKNOWN).value
    return Constellation.UNKNOWN.value


def _leaf_only_unique_keys(
    data: Dict[str, Any],
    out: Dict[str, Any],
    origins: Dict[str, str],
    path: str = "",
) -> None:
    for key, value in data.items():
        here = f"{path}.{key}" if path else key
        if isinstance(value, dict):
            _leaf_only_unique_keys(value, out, origins, here)
            continue

        if key in out:
            raise AssertionError(
                f"Duplicate leaf key '{key}' at '{here}' (already from '{origins[key]}')."
            )
        out[key] = value
        origins[key] = here


def _dataclass_to_leaf_row(obj: Any) -> Dict[str, Any]:
    if is_dataclass(obj):
        data = asdict(obj)
    elif isinstance(obj, dict):
        data = obj
    else:
        data = dict(getattr(obj, "__dict__", {}))

    out: Dict[str, Any] = {}
    origins: Dict[str, str] = {}
    _leaf_only_unique_keys(data, out, origins)
    return out


def _pad3(seq, fill=math.nan):
    seq = list(seq) if seq is not None else []
    return (seq + [fill, fill, fill])[:3]


GPS_EPOCH = datetime(1980, 1, 6, tzinfo=timezone.utc)
SECONDS_IN_WEEK = 7 * 24 * 3600


@dataclass
class TimeTag:
    t_sec: int = math.nan
    gps_week: int = -1
    tow_sec: int = math.nan

    @staticmethod
    def from_sec(t: float) -> "TimeTag":
        if not isinstance(t, (int, float)) or math.isnan(t):
            return TimeTag(t_sec=math.nan, gps_week=-1, tow_sec=math.nan)

        dt = GPS_EPOCH + timedelta(seconds=float(t))
        total_sec = (dt - GPS_EPOCH).total_seconds()
        week = int(total_sec // SECONDS_IN_WEEK)
        tow = float(total_sec - week * SECONDS_IN_WEEK)
        return TimeTag(t_sec=int(t), gps_week=week, tow_sec=int(tow))


@dataclass
class SatelliteInfo:
    sv_pos: List[float] = field(default_factory=lambda: [math.nan, math.nan, math.nan])
    sv_vel: List[float] = field(default_factory=lambda: [math.nan, math.nan, math.nan])
    sv_clock_bias: float = math.nan
    sv_clock_drift: float = math.nan
    pr_correction: float = math.nan
    iono_coff: List[float] = field(default_factory=lambda: [0.0] * 8)
    constellation: Constellation = Constellation.UNKNOWN
    prn: int = -1
    frequency_hz: float = math.nan
    code_type: str = ""

    def validate(self) -> None:
        if len(self.sv_pos) != 3 or len(self.sv_vel) != 3:
            raise ValueError("sv_pos/sv_vel must be length-3 lists.")
        if self.prn < -1:
            raise ValueError("prn must be >= -1.")


@dataclass
class MeasurementValue:
    pseudorange_m: float = math.nan
    phase_cycle: float = math.nan
    doppler_hz: float = math.nan
    snr_dbhz: float = math.nan
    loi: bool = False

    def validate_basic(self) -> None:
        if not math.isnan(self.snr_dbhz) and not (0.0 <= self.snr_dbhz <= 70.0):
            raise ValueError(f"SNR out of plausible range: {self.snr_dbhz}")


@dataclass
class Measurement:
    time: TimeTag = field(default_factory=TimeTag)
    sat: SatelliteInfo = field(default_factory=SatelliteInfo)
    value: MeasurementValue = field(default_factory=MeasurementValue)
    ground_truth: List[float] = field(default_factory=lambda: [math.nan, math.nan, math.nan])

    def validate(self) -> None:
        self.sat.validate()
        self.value.validate_basic()
        if self.time.gps_week < -1:
            raise ValueError("gps_week must be >= -1.")

    @staticmethod
    def headers() -> List[str]:
        return [
            "t_sec",
            "gps_week",
            "tow_sec",
            "constellation",
            "prn",
            "sv_pos_x",
            "sv_pos_y",
            "sv_pos_z",
            "sv_vel_x",
            "sv_vel_y",
            "sv_vel_z",
            "iono_a0",
            "iono_a1",
            "iono_a2",
            "iono_a3",
            "iono_b0",
            "iono_b1",
            "iono_b2",
            "iono_b3",
            "sv_clock_bias",
            "sv_clock_drift",
            "pr_correction",
            "frequency_hz",
            "code_type",
            "pseudorange_m",
            "phase_cycle",
            "doppler_hz",
            "snr_dbhz",
            "loi",
            "gt_pos_x",
            "gt_pos_y",
            "gt_pos_z",
        ]

    def _to_row_dict(self) -> Dict[str, Any]:
        leaf = _dataclass_to_leaf_row(self)

        sv_pos_x, sv_pos_y, sv_pos_z = _pad3(leaf.pop("sv_pos", []))
        sv_vel_x, sv_vel_y, sv_vel_z = _pad3(leaf.pop("sv_vel", []))
        gt_x_m, gt_y_m, gt_z_m = _pad3(leaf.pop("ground_truth", []))
        iono = list(leaf.pop("iono_coff", []))
        iono = (iono + [0.0] * 8)[:8]

        if "loi" in leaf and isinstance(leaf["loi"], bool):
            leaf["loi"] = int(leaf["loi"])

        leaf["constellation"] = to_constellation(leaf.get("constellation", None))
        leaf.update(
            {
                "sv_pos_x": sv_pos_x,
                "sv_pos_y": sv_pos_y,
                "sv_pos_z": sv_pos_z,
                "sv_vel_x": sv_vel_x,
                "sv_vel_y": sv_vel_y,
                "sv_vel_z": sv_vel_z,
                "iono_a0": iono[0],
                "iono_a1": iono[1],
                "iono_a2": iono[2],
                "iono_a3": iono[3],
                "iono_b0": iono[4],
                "iono_b1": iono[5],
                "iono_b2": iono[6],
                "iono_b3": iono[7],
                "gt_pos_x": gt_x_m,
                "gt_pos_y": gt_y_m,
                "gt_pos_z": gt_z_m,
            }
        )
        return leaf

    def to_csv(self, sep: str = "\t") -> str:
        row = self._to_row_dict()
        headers = self.headers()
        missing = [header for header in headers if header not in row]
        if missing:
            raise KeyError(f"Missing columns for CSV: {missing}")

        fields = [row[header] for header in headers]
        buf = io.StringIO()
        writer = csv.writer(buf, delimiter=sep, lineterminator="")
        writer.writerow(fields)
        return buf.getvalue()


def measurements_to_csv(meas_list: list[Measurement], sep: str = "\t") -> str:
    return "\n".join(measurement.to_csv(sep) for measurement in meas_list)

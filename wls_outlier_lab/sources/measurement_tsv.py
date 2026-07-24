"""Measurement source: canonical measurement TSV / processed CSV.

Reads the schema produced by ``gnss_txt_parser`` + the MATLAB SV-position
endpoint (32-column ``*_with_sv_pos*.tsv``) and the 30-column
``data/*_processed*.csv`` variant. Columns are matched by header name, so
column order and the presence of extra columns (``dop_correction``,
``gt_pos_*``) do not matter.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from ..types import Epoch, SatObs
from .base import MeasurementSource


def _f(value: Optional[str]) -> float:
    if value is None:
        return float("nan")
    value = value.strip()
    if value == "" or value.lower() == "nan":
        return float("nan")
    try:
        return float(value)
    except ValueError:
        return float("nan")


def _sniff_delimiter(header_line: str) -> str:
    return "\t" if "\t" in header_line else ","


def _row_to_obs(row: Dict[str, str]) -> Optional[SatObs]:
    try:
        constellation = int(_f(row.get("constellation")))
        prn = int(_f(row.get("prn")))
    except ValueError:
        return None
    if not math.isfinite(_f(row.get("constellation"))):
        return None
    pr_corr = _f(row.get("pr_correction"))
    sv_vel = None
    if "sv_vel_x" in row:
        sv_vel = (_f(row.get("sv_vel_x")), _f(row.get("sv_vel_y")), _f(row.get("sv_vel_z")))
    loi_raw = _f(row.get("loi"))
    return SatObs(
        constellation=constellation,
        prn=prn,
        frequency_hz=_f(row.get("frequency_hz")),
        code_type=(row.get("code_type") or "").strip(),
        pseudorange_m=_f(row.get("pseudorange_m")),
        sv_pos=(_f(row.get("sv_pos_x")), _f(row.get("sv_pos_y")), _f(row.get("sv_pos_z"))),
        sv_clock_bias_m=_f(row.get("sv_clock_bias")),
        sv_vel=sv_vel,
        pr_correction_m=pr_corr if math.isfinite(pr_corr) else 0.0,
        doppler_hz=_f(row.get("doppler_hz")),
        cn0_dbhz=_f(row.get("snr_dbhz")),
        loi=bool(int(loi_raw)) if math.isfinite(loi_raw) else False,
    )


def load_epochs(
    path: str | Path,
    *,
    stride: int = 1,
    max_epochs: Optional[int] = None,
    t_start: Optional[float] = None,
    t_end: Optional[float] = None,
    constellations: Optional[Sequence[int]] = None,
) -> List[Epoch]:
    """Stream the file and group consecutive rows sharing ``t_sec`` into epochs.

    ``stride`` keeps every k-th epoch; ``max_epochs`` caps the count; the time
    window and constellation filter subset the data. The file is assumed sorted
    by ``t_sec`` (as emitted by the parser/endpoint).
    """
    path = Path(path)
    cons_filter = set(constellations) if constellations is not None else None
    epochs: List[Epoch] = []
    kept_epoch_index = -1

    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        first = fh.readline()
        delimiter = _sniff_delimiter(first)
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter=delimiter)

        cur_t: Optional[float] = None
        cur_week = -1
        cur_tow = float("nan")
        cur_obs: List[SatObs] = []

        def flush() -> bool:
            nonlocal kept_epoch_index
            if cur_t is None or not cur_obs:
                return True
            if t_start is not None and cur_t < t_start:
                return True
            if t_end is not None and cur_t > t_end:
                return True
            kept_epoch_index += 1
            if kept_epoch_index % stride != 0:
                return True
            epochs.append(Epoch(t_sec=cur_t, obs=list(cur_obs),
                                gps_week=cur_week, tow_sec=cur_tow))
            if max_epochs is not None and len(epochs) >= max_epochs:
                return False
            return True

        for row in reader:
            t = _f(row.get("t_sec"))
            if not math.isfinite(t):
                continue
            if cur_t is None:
                cur_t = t
            if t != cur_t:
                if not flush():
                    return epochs
                cur_t = t
                cur_obs = []
            cur_week = int(_f(row.get("gps_week"))) if math.isfinite(_f(row.get("gps_week"))) else -1
            cur_tow = _f(row.get("tow_sec"))
            obs = _row_to_obs(row)
            if obs is None:
                continue
            if cons_filter is not None and obs.constellation not in cons_filter:
                continue
            cur_obs.append(obs)
        flush()
    return epochs


class TsvMeasurementSource(MeasurementSource):
    """MeasurementSource wrapper around :func:`load_epochs`."""

    def __init__(self, path: str | Path, **options) -> None:
        self.path = Path(path)
        self.options = options

    def load_epochs(self) -> List[Epoch]:
        return load_epochs(self.path, **self.options)

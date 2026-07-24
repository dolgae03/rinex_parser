"""Truth source: Novatel BESTPOS ASCII (#BESTPOSA ... or <BESTPOS abbrev).

The Novatel field logs supply an RTK (NARROW_INT) BESTPOS stream that is the
ground truth for the smartphone-under-test analysis. This is the one piece of
Novatel handling that lived only in MATLAB (`parse_gps_best_pos.m`); here it is
in Python, producing a ``TruthTrack`` and, optionally, the RTK CSV schema the
existing ``api_server`` truth pipeline already understands.

Record layout (ASCII), header and body split on ';':
    #BESTPOSA,port,seq,idle,timestatus,WEEK,TOW,...;
        sol_status,pos_type,LAT,LON,HGT_HAE,undulation,datum,lat_std,lon_std,...
GPS seconds since 1980-01-06 = WEEK*604800 + TOW (GPS time base, matches the
measurement TSV ``t_sec``).
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from ..core.frames import lla_to_ecef
from ..types import TruthSample, TruthTrack
from .base import TruthSource

SECONDS_PER_WEEK = 604800


@dataclass
class BestposRecord:
    t_sec: float
    week: int
    tow: float
    sol_status: str
    pos_type: str
    lat_deg: float
    lon_deg: float
    height_m: float
    lat_std: float
    lon_std: float
    height_std: float


def _clean(token: str) -> str:
    return token.split("*")[0].strip()  # strip trailing *checksum


def parse_bestpos(
    path: str | Path,
    require_sol_status: str = "SOL_COMPUTED",
    require_pos_types: Optional[Sequence[str]] = None,
) -> List[BestposRecord]:
    """Parse #BESTPOSA lines. ``require_pos_types`` e.g. ('NARROW_INT',) for RTK-fixed only."""
    path = Path(path)
    pos_type_filter = {p.upper() for p in require_pos_types} if require_pos_types else None
    records: List[BestposRecord] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if "BESTPOS" not in line or ";" not in line:
                continue
            header, _, body = line.partition(";")
            htok = header.split(",")
            btok = body.split(",")
            if len(htok) < 7 or len(btok) < 8:
                continue
            try:
                week = int(_clean(htok[5]))
                tow = float(_clean(htok[6]))
            except ValueError:
                continue
            sol_status = _clean(btok[0]).upper()
            pos_type = _clean(btok[1]).upper()
            if require_sol_status and sol_status != require_sol_status.upper():
                continue
            if pos_type_filter is not None and pos_type not in pos_type_filter:
                continue
            try:
                lat = float(_clean(btok[2]))
                lon = float(_clean(btok[3]))
                hgt = float(_clean(btok[4]))
                lat_std = float(_clean(btok[7])) if len(btok) > 7 else float("nan")
                lon_std = float(_clean(btok[8])) if len(btok) > 8 else float("nan")
                hgt_std = float(_clean(btok[9])) if len(btok) > 9 else float("nan")
            except ValueError:
                continue
            records.append(BestposRecord(
                t_sec=week * SECONDS_PER_WEEK + tow, week=week, tow=tow,
                sol_status=sol_status, pos_type=pos_type,
                lat_deg=lat, lon_deg=lon, height_m=hgt,
                lat_std=lat_std, lon_std=lon_std, height_std=hgt_std,
            ))
    return records


def load_bestpos_track(
    path: str | Path,
    name: str = "novatel_bestpos",
    require_pos_types: Optional[Sequence[str]] = ("NARROW_INT", "INS_RTKFIXED"),
) -> TruthTrack:
    """BESTPOS -> ECEF TruthTrack. Defaults to RTK-fixed (NARROW_INT) samples.

    Duplicate timestamps (dual-antenna / multiple ports) are averaged.
    """
    records = parse_bestpos(path, require_pos_types=require_pos_types)
    if not records:
        # fall back to any solution type if no RTK-fixed rows exist
        records = parse_bestpos(path, require_pos_types=None)
    by_t: Dict[float, List[BestposRecord]] = {}
    for r in records:
        by_t.setdefault(r.t_sec, []).append(r)
    samples: List[TruthSample] = []
    for t, recs in by_t.items():
        r = recs[0]
        ecef = lla_to_ecef(r.lat_deg, r.lon_deg, r.height_m)
        samples.append(TruthSample(t, (float(ecef[0]), float(ecef[1]), float(ecef[2]))))
    pos_types = sorted({r.pos_type for r in records})
    return TruthTrack(name=name, samples=samples,
                      source_type=f"novatel_bestpos[{','.join(pos_types)}]")


def export_rtk_csv(path: str | Path, out_csv: str | Path,
                   require_pos_types: Optional[Sequence[str]] = ("NARROW_INT", "INS_RTKFIXED")) -> Path:
    """Write the RTK CSV schema (gps_sec, latitude_deg, longitude_deg, height_m)
    consumed by ``api_server.truth_alignment`` — a bridge to the existing tooling."""
    records = parse_bestpos(path, require_pos_types=require_pos_types)
    out_csv = Path(out_csv)
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["gps_sec", "latitude_deg", "longitude_deg", "height_m"])
        for r in records:
            w.writerow([f"{r.t_sec:.3f}", f"{r.lat_deg:.11f}",
                        f"{r.lon_deg:.11f}", f"{r.height_m:.4f}"])
    return out_csv


class BestposTruthSource(TruthSource):
    def __init__(self, path: str | Path, name: str = "novatel_bestpos",
                 require_pos_types: Optional[Sequence[str]] = ("NARROW_INT", "INS_RTKFIXED")) -> None:
        self.path = Path(path)
        self.name = name
        self.require_pos_types = require_pos_types

    def load_track(self) -> TruthTrack:
        return load_bestpos_track(self.path, self.name, self.require_pos_types)

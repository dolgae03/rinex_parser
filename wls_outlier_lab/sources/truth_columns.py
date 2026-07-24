"""Truth source: ground-truth ECEF already embedded as gt_pos_* columns.

The API server's truth pipeline writes ``gt_pos_x/y/z`` into the measurement
TSV (``*_with_gt_*.tsv``). This adapter lifts those into a ``TruthTrack`` so the
lab can be validated end-to-end without re-parsing the original truth file.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Dict, List

from ..types import TruthSample, TruthTrack
from .base import TruthSource
from .measurement_tsv import _f, _sniff_delimiter


def load_truth_from_columns(path: str | Path, name: str = "embedded_gt") -> TruthTrack:
    path = Path(path)
    by_t: Dict[float, List[tuple]] = {}
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        first = fh.readline()
        delimiter = _sniff_delimiter(first)
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames or "gt_pos_x" not in reader.fieldnames:
            raise ValueError(f"{path} has no gt_pos_x column")
        for row in reader:
            t = _f(row.get("t_sec"))
            x, y, z = _f(row.get("gt_pos_x")), _f(row.get("gt_pos_y")), _f(row.get("gt_pos_z"))
            if not all(math.isfinite(v) for v in (t, x, y, z)):
                continue
            by_t.setdefault(t, []).append((x, y, z))
    samples = [
        TruthSample(t, (
            sum(p[0] for p in pts) / len(pts),
            sum(p[1] for p in pts) / len(pts),
            sum(p[2] for p in pts) / len(pts),
        ))
        for t, pts in by_t.items()
    ]
    return TruthTrack(name=name, samples=samples, source_type="embedded_gt_columns")


class ColumnTruthSource(TruthSource):
    def __init__(self, path: str | Path, name: str = "embedded_gt") -> None:
        self.path = Path(path)
        self.name = name

    def load_track(self) -> TruthTrack:
        return load_truth_from_columns(self.path, self.name)

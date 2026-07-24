"""Synthetic epoch generator for tests — perfect geometry with a known truth."""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import numpy as np

from wls_outlier_lab.core import frames
from wls_outlier_lab.types import Epoch, SatObs, TruthSample, TruthTrack

RX_LLA = (37.2557, 127.0553, 40.0)


def rx_ecef() -> np.ndarray:
    return frames.lla_to_ecef(*RX_LLA)


def _sat_at(rx: np.ndarray, az_deg: float, el_deg: float, radius: float = 2.2e7) -> np.ndarray:
    R = frames.enu_rotation_matrix(*frames.ecef_to_lla(*rx)[:2])
    az, el = math.radians(az_deg), math.radians(el_deg)
    los_enu = np.array([math.cos(el) * math.sin(az), math.cos(el) * math.cos(az), math.sin(el)])
    return rx + radius * (R.T @ los_enu)


def make_epoch(
    t_sec: float = 1_459_488_809.0,
    clk_m: float = 1234.5,
    n_gps: int = 8,
    isb: Optional[Dict[int, float]] = None,
    n_per_extra: int = 6,
    blunders: Optional[Dict[Tuple[int, int], float]] = None,
    noise_m: float = 0.0,
    seed: int = 0,
) -> Tuple[Epoch, np.ndarray]:
    """Build one epoch of clean pseudoranges (+ optional per-sat blunders).

    Returns (epoch, true_rx_ecef). Pseudorange = geometric range + clk + isb,
    so a solver with sagnac/iono disabled must recover rx and clk exactly.
    """
    rng = np.random.default_rng(seed)
    rx = rx_ecef()
    isb = isb or {}
    blunders = blunders or {}
    obs: List[SatObs] = []

    def add(cons: int, count: int) -> None:
        for i in range(count):
            az = (cons * 37 + i * (360.0 / max(count, 1))) % 360.0
            el = 12.0 + (i * 71.0) % 73.0
            sat = _sat_at(rx, az, el)
            geo = float(np.linalg.norm(sat - rx))
            bias = isb.get(cons, 0.0)
            err = blunders.get((cons, i + 1), 0.0)
            n = rng.normal(0, noise_m) if noise_m > 0 else 0.0
            pr = geo + clk_m + bias + err + n
            obs.append(SatObs(
                constellation=cons, prn=i + 1, frequency_hz=1575420000.0,
                code_type="C", pseudorange_m=pr, sv_pos=(sat[0], sat[1], sat[2]),
                sv_clock_bias_m=0.0, cn0_dbhz=45.0,
            ))

    add(0, n_gps)
    for cons in isb:
        if cons != 0:
            add(cons, n_per_extra)
    return Epoch(t_sec=t_sec, obs=obs), rx


def make_track(rx: np.ndarray, t_sec: float) -> TruthTrack:
    return TruthTrack(name="synthetic", samples=[TruthSample(t_sec, (rx[0], rx[1], rx[2]))])

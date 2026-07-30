"""Receiver clock-drift stability analysis (synthetic)."""

from __future__ import annotations

from wls_outlier_lab.core.clock_analysis import analyze_clock_stability
from wls_outlier_lab.core.wls import WlsConfig
from wls_outlier_lab.types import TruthSample, TruthTrack
from wls_outlier_lab.tests.synthetic import make_epoch

EXACT = WlsConfig(apply_sagnac=False, apply_iono=False, apply_tropo=False)
T0 = 1_459_488_000.0
DRIFT = 50.0        # m per 1 s epoch
JUMP_EPOCH = 15
JUMP_M = 300.0


def _dataset(n=30):
    epochs, samples = [], []
    for i in range(n):
        clk = 1000.0 + DRIFT * i + (JUMP_M if i == JUMP_EPOCH else 0.0)
        ep, rx = make_epoch(t_sec=T0 + i, clk_m=clk, n_gps=9, noise_m=1.0, seed=i)
        epochs.append(ep)
        samples.append(TruthSample(ep.t_sec, (rx[0], rx[1], rx[2])))
    return epochs, TruthTrack("synthetic", samples)


def test_clock_drift_and_jump_detection():
    epochs, truth = _dataset()
    ca = analyze_clock_stability(epochs, truth, wls_cfg=EXACT, match_tolerance_sec=1.0)
    assert ca.n_epochs >= 25
    # constant-velocity drift recovered (~50 m/s at 1 Hz)
    assert abs(ca.median_drift_mps - DRIFT) < 5.0
    # the injected clock jump is flagged as an anomaly
    assert ca.n_anomalies >= 1
    worst = max(ca.points, key=lambda p: (p.instability_m if p.instability_m == p.instability_m else -1))
    assert abs(worst.t_sec - (T0 + JUMP_EPOCH)) < 1.5


def test_clock_does_not_wreck_horizontal():
    # despite a drifting/jumping clock, per-epoch WLS keeps horizontal error small
    epochs, truth = _dataset()
    ca = analyze_clock_stability(epochs, truth, wls_cfg=EXACT, match_tolerance_sec=1.0)
    assert ca.stable_mean_h_m < 5.0

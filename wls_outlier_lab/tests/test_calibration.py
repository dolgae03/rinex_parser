"""Measurement-noise calibration from truth (synthetic)."""

from __future__ import annotations

from wls_outlier_lab.core.calibration import calibrate, measurement_residuals
from wls_outlier_lab.core.wls import WlsConfig
from wls_outlier_lab.types import TruthSample, TruthTrack
from wls_outlier_lab.tests.synthetic import make_epoch

EXACT = WlsConfig(apply_sagnac=False, apply_iono=False)


def _dataset(n_epochs=30, noise_m=4.0, blunders=None):
    epochs, samples = [], []
    for i in range(n_epochs):
        ep, rx = make_epoch(t_sec=1_459_488_000.0 + i, n_gps=9, isb={2: 25.0},
                            noise_m=noise_m, blunders=blunders, seed=i)
        epochs.append(ep)
        samples.append(TruthSample(ep.t_sec, (rx[0], rx[1], rx[2])))
    return epochs, TruthTrack("synthetic", samples)


def test_calibration_recovers_noise_level():
    epochs, truth = _dataset(noise_m=4.0)
    res = measurement_residuals(epochs, truth, EXACT, match_tolerance_sec=1.0)
    assert len(res) > 0
    cal = calibrate(res)
    gps = cal.by_constellation["GPS"]
    # detrended residual std should be near the injected 4 m (within a factor ~2)
    assert 1.5 < gps["clean_std_m"] < 8.0
    assert cal.sigma_scale_by_constellation["GPS"] == 1.0


def test_calibration_flags_blunder_constellation():
    # BeiDou gets a persistent gross blunder on one satellite
    epochs, truth = _dataset(noise_m=3.0, blunders={(2, 2): 4000.0})
    cal = calibrate(measurement_residuals(epochs, truth, EXACT, match_tolerance_sec=1.0))
    assert cal.by_constellation["BEIDOU"]["outlier_rate"] > cal.by_constellation["GPS"]["outlier_rate"]

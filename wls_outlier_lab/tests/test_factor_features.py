"""Factor-matrix exporter: shape, and that the factors mean what they claim.

The MATLAB side (factor_correlation_analysis.m) trusts these columns blindly,
so the invariants worth pinning are semantic: a clean static scenario must give
near-zero DCD (Doppler-Code Difference) and near-zero truth speed, and every
column in FACTOR_COLUMNS must actually be emitted.
"""

from __future__ import annotations

import math

from wls_outlier_lab.core.factor_features import FACTOR_COLUMNS, extract_factor_epochs
from wls_outlier_lab.core.wls import WlsConfig
from wls_outlier_lab.types import TruthSample, TruthTrack

from .synthetic import make_epoch


def _static_scene(n_epochs: int = 8, drift_mps: float = -250.0):
    epochs = []
    rx = None
    t0 = 1_459_488_809.0
    clk0 = 1234.5
    for k in range(n_epochs):
        ep, rx = make_epoch(t_sec=t0 + k, clk_m=clk0 + drift_mps * k,
                            n_gps=9, drift_mps=drift_mps, noise_m=0.05, seed=k)
        epochs.append(ep)
    truth = TruthTrack(name="synthetic", samples=[
        TruthSample(t0 + k, (rx[0], rx[1], rx[2])) for k in range(n_epochs)])
    return epochs, truth


def test_exports_every_declared_column():
    epochs, truth = _static_scene()
    cfg = WlsConfig(apply_sagnac=False, apply_tropo=False, apply_iono=False)
    rows = extract_factor_epochs(epochs, truth, wls_cfg=cfg)
    assert len(rows) == len(epochs)
    missing = [c for c in FACTOR_COLUMNS if c not in rows[-1]]
    assert not missing, f"exporter dropped declared columns: {missing}"


def test_static_clean_scene_yields_consistent_factors():
    epochs, truth = _static_scene()
    cfg = WlsConfig(apply_sagnac=False, apply_tropo=False, apply_iono=False)
    rows = extract_factor_epochs(epochs, truth, wls_cfg=cfg)
    last = rows[-1]
    # truth is static -> speed ~ 0; solution is clean -> tiny horizontal error
    assert abs(last["speed_mps"]) < 0.01
    assert last["h_err_m"] < 1.0
    # code and Doppler agree by construction -> DCD stays at the noise level:
    # d(PR)/dt noise = 0.05 m * sqrt(2) over 1 s, far below 1 m/s
    assert math.isfinite(last["dcd_med_mps"]) and abs(last["dcd_med_mps"]) < 1.0
    # position-domain clock drift recovers the injected -250 m/s
    assert abs(last["clk_drift_mps"] - (-250.0)) < 1.0
    assert abs(last["drift_dop_mps"] - (-250.0)) < 1.0
    # no phase in the synthetic scene -> CCD must be NaN, not a fake zero
    assert not math.isfinite(last["ccd_med_mps"])


def test_blunder_epoch_raises_residual_factors():
    epochs, truth = _static_scene()
    cfg = WlsConfig(apply_sagnac=False, apply_tropo=False, apply_iono=False)
    clean = extract_factor_epochs(epochs, truth, wls_cfg=cfg)

    dirty_epochs = []
    t0 = epochs[0].t_sec
    for k in range(len(epochs)):
        ep, _ = make_epoch(t_sec=t0 + k, clk_m=1234.5 - 250.0 * k, n_gps=9,
                           drift_mps=-250.0, noise_m=0.05, seed=k,
                           blunders={(0, 3): 25.0})
        dirty_epochs.append(ep)
    dirty = extract_factor_epochs(dirty_epochs, truth, wls_cfg=cfg)

    # a 25 m blunder on one satellite must show up in the residual factors
    # (the combined detector may or may not evict it; either way sigma0 of the
    # clean scene stays below the dirty scene unless the sat was rejected)
    clean_s = [r["resid_max_m"] for r in clean if math.isfinite(r["resid_max_m"])]
    dirty_s = [max(r["resid_max_m"], 25.0 * (r["n_rejected"] > 0))
               for r in dirty if math.isfinite(r["resid_max_m"])]
    assert max(dirty_s) > max(clean_s)

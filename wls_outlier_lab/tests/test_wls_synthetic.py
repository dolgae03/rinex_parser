"""WLS correctness on synthetic geometry + outlier-detection behaviour."""

from __future__ import annotations

import math

import numpy as np

from wls_outlier_lab.core.detectors import DetectorConfig, detect_residual
from wls_outlier_lab.core.experiment import run_experiment
from wls_outlier_lab.core.wls import WlsConfig, solve_epoch
from wls_outlier_lab.tests.synthetic import make_epoch, make_track

# perfect-data config: no sagnac / iono / tropo so the synthetic ranges are exact
EXACT = WlsConfig(apply_sagnac=False, apply_iono=False, apply_tropo=False)


def test_recovers_position_single_constellation():
    ep, rx = make_epoch(clk_m=1234.5, n_gps=9)
    sol = solve_epoch(ep.obs, ep.t_sec, EXACT)
    assert sol.converged
    err = np.linalg.norm(np.array(sol.ecef) - rx)
    assert err < 1e-3, err
    assert abs(sol.clock_bias_m - 1234.5) < 1e-3


def test_recovers_inter_system_bias():
    ep, rx = make_epoch(clk_m=500.0, n_gps=8, isb={2: 40.0, 1: -15.0}, n_per_extra=6)
    sol = solve_epoch(ep.obs, ep.t_sec, EXACT)
    assert sol.converged
    assert np.linalg.norm(np.array(sol.ecef) - rx) < 1e-3
    assert abs(sol.isb_m[2] - 40.0) < 1e-3
    assert abs(sol.isb_m[1] - (-15.0)) < 1e-3


def test_residual_detects_gross_blunder():
    # one satellite off by +600 m
    ep, rx = make_epoch(n_gps=9, blunders={(0, 4): 600.0})
    rej = detect_residual(ep.obs, ep.t_sec, EXACT, DetectorConfig())
    assert (0, 4) in rej
    clean = solve_epoch(ep.obs, ep.t_sec, EXACT, reject=rej)
    assert np.linalg.norm(np.array(clean.ecef) - rx) < 1e-2


def test_baseline_stays_finite_with_blunder():
    # a huge blunder must not make the solver diverge (line search keeps it finite)
    ep, rx = make_epoch(n_gps=9, blunders={(0, 3): 1.0e6})
    sol = solve_epoch(ep.obs, ep.t_sec, EXACT)
    assert all(math.isfinite(v) for v in sol.ecef)


def test_experiment_residual_beats_baseline():
    ep, rx = make_epoch(n_gps=9, isb={2: 30.0}, blunders={(0, 5): 800.0, (2, 2): 5000.0})
    track = make_track(rx, ep.t_sec)
    result = run_experiment([ep], track, wls_cfg=EXACT,
                            detectors=["baseline", "residual"], match_tolerance_sec=1.0)
    base = result.runs["baseline"].metrics.horizontal_rmse_m
    resid = result.runs["residual"].metrics.horizontal_rmse_m
    assert resid < base
    assert resid < 1.0

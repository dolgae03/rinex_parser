"""Clock pseudo-observation, Doppler drift, and the coasting experiment."""

from __future__ import annotations

import math

import numpy as np

from wls_outlier_lab.core.clock_analysis import (
    clock_coasting_experiment,
    doppler_clock_drift,
)
from wls_outlier_lab.core.wls import WlsConfig, solve_epoch
from wls_outlier_lab.types import TruthSample, TruthTrack
from wls_outlier_lab.tests.synthetic import make_epoch

EXACT = WlsConfig(apply_sagnac=False, apply_iono=False, apply_tropo=False)
T0 = 1_459_488_000.0
DRIFT = 40.0        # clock drift [m/s] == [m per 1 s epoch]


def _track(n, jump_at=None, jump_m=0.0, noise_m=0.5):
    epochs, samples = [], []
    for i in range(n):
        clk = 500.0 + DRIFT * i + (jump_m if (jump_at is not None and i >= jump_at) else 0.0)
        ep, rx = make_epoch(t_sec=T0 + i, clk_m=clk, n_gps=9, noise_m=noise_m, seed=i)
        epochs.append(ep)
        samples.append(TruthSample(ep.t_sec, (rx[0], rx[1], rx[2])))
    return epochs, TruthTrack("synthetic", samples)


def test_clock_prior_consistent_is_harmless_wrong_is_not():
    ep, rx = make_epoch(clk_m=1234.5, n_gps=9)
    free = solve_epoch(ep.obs, ep.t_sec, EXACT)
    assert abs(free.clock_bias_m - 1234.5) < 1e-3

    good = solve_epoch(ep.obs, ep.t_sec, EXACT, clock_prior=(1234.5, 1.0))
    assert np.linalg.norm(np.array(good.ecef) - rx) < 0.1

    # a tight but wrong clock prior must corrupt the fix (that is the mechanism
    # by which a bad clock model damages a coasted solution)
    bad = solve_epoch(ep.obs, ep.t_sec, EXACT, clock_prior=(1234.5 + 500.0, 0.01))
    assert np.linalg.norm(np.array(bad.ecef) - rx) > 10.0


def test_clock_prior_buys_one_satellite():
    ep, rx = make_epoch(clk_m=1234.5, n_gps=3)      # 3 sats: underdetermined for 4 unknowns
    assert not math.isfinite(solve_epoch(ep.obs, ep.t_sec, EXACT).ecef[0])
    ok = solve_epoch(ep.obs, ep.t_sec, EXACT, clock_prior=(1234.5, 0.1))
    assert np.linalg.norm(np.array(ok.ecef) - rx) < 1.0


def test_doppler_recovers_clock_drift():
    ep, rx = make_epoch(clk_m=100.0, n_gps=10, drift_mps=-250.0)
    r = doppler_clock_drift(ep, rx, sign=-1.0)
    assert r["n_sats"] >= 8
    assert abs(r["drift_mps"] - (-250.0)) < 1e-3
    assert r["speed_mps"] < 1e-3          # static receiver


def test_coasting_stable_clock_is_fine_but_a_jump_hurts():
    stable, truth = _track(25)
    rows = clock_coasting_experiment(stable, truth, wls_cfg=EXACT, detector="baseline",
                                     match_tolerance_sec=1.0, sigmas_m=(float("inf"), 1.0))
    free = next(r for r in rows if r["mode"] == "free")
    tight = next(r for r in rows if r["mode"] != "free")
    # a perfectly linear clock can be coasted at no cost
    assert tight["h_rmse_m"] < free["h_rmse_m"] * 3.0 + 1.0

    jumpy, jtruth = _track(25, jump_at=12, jump_m=300.0)
    jrows = clock_coasting_experiment(jumpy, jtruth, wls_cfg=EXACT, detector="baseline",
                                      match_tolerance_sec=1.0, sigmas_m=(float("inf"), 1.0))
    jfree = next(r for r in jrows if r["mode"] == "free")
    jtight = next(r for r in jrows if r["mode"] != "free")
    # with an unmodelled 300 m clock jump, trusting the clock model must cost accuracy
    assert jtight["h_rmse_m"] > jfree["h_rmse_m"] * 2.0

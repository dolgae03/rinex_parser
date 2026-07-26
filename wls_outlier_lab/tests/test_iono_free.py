"""Dual-frequency ionosphere-free combination."""

from __future__ import annotations

from wls_outlier_lab.core.iono_free import form_iono_free_epoch
from wls_outlier_lab.types import Epoch, SatObs

F1 = 1575.42e6   # L1
F2 = 1176.45e6   # L5


def test_iono_free_removes_ionosphere():
    geom_plus_clk = 22_000_000.0 + 1234.0
    iono1 = 15.0
    iono2 = iono1 * (F1 / F2) ** 2  # ionosphere scales as 1/f^2
    o1 = SatObs(0, 4, F1, "C", pseudorange_m=geom_plus_clk + iono1,
                sv_pos=(1.0, 2.0, 3.0), sv_clock_bias_m=0.0, cn0_dbhz=45.0)
    o2 = SatObs(0, 4, F2, "C", pseudorange_m=geom_plus_clk + iono2,
                sv_pos=(1.0, 2.0, 3.0), sv_clock_bias_m=0.0, cn0_dbhz=44.0)
    ie = form_iono_free_epoch(Epoch(t_sec=100.0, obs=[o1, o2]))
    assert len(ie.obs) == 1
    assert ie.obs[0].code_type == "IF"
    assert abs(ie.obs[0].corrected_pseudorange_m - geom_plus_clk) < 1e-3


def test_iono_free_drops_single_band_satellite():
    o1 = SatObs(0, 4, F1, "C", pseudorange_m=2.2e7, sv_pos=(1.0, 2.0, 3.0),
                sv_clock_bias_m=0.0, cn0_dbhz=45.0)
    ie = form_iono_free_epoch(Epoch(t_sec=100.0, obs=[o1]))
    assert len(ie.obs) == 0  # cannot form IF from one band

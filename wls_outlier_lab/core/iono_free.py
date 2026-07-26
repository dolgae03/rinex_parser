"""Dual-frequency ionosphere-free pseudorange combination.

When a satellite is tracked on an upper band (L1/E1/B1 ~1.5 GHz) and a lower band
(L5/E5a/B2a ~1.2 GHz), the ionosphere-free combination

    PR_IF = (f1^2 * PR1 - f2^2 * PR2) / (f1^2 - f2^2)

removes ~all first-order ionospheric delay — no broadcast model, no reliance on
the (here broken) Klobuchar coefficients. Cost: ~3x measurement noise and losing
single-frequency-only satellites. This transforms an epoch's observations into
one IF observation per dual-band satellite; the WLS then runs unchanged with the
ionosphere already gone.
"""

from __future__ import annotations

from collections import defaultdict
from typing import List, Sequence

from ..types import Epoch, SatObs

BAND_SPLIT_HZ = 1_400_000_000.0  # upper band > this, lower band <= this


def form_iono_free_epoch(epoch: Epoch) -> Epoch:
    by_sat = defaultdict(list)
    for o in epoch.obs:
        by_sat[(o.constellation, o.prn)].append(o)

    new_obs: List[SatObs] = []
    for (cons, prn), sigs in by_sat.items():
        upper = [o for o in sigs if o.frequency_hz > BAND_SPLIT_HZ]
        lower = [o for o in sigs if 0 < o.frequency_hz <= BAND_SPLIT_HZ]
        if not upper or not lower:
            continue  # need both bands to form IF

        def _cn0(o):
            return o.cn0_dbhz if o.cn0_dbhz == o.cn0_dbhz else -1.0

        o1 = max(upper, key=_cn0)
        o2 = max(lower, key=_cn0)
        f1, f2 = o1.frequency_hz, o2.frequency_hz
        if f1 == f2:
            continue
        denom = f1 * f1 - f2 * f2
        # IF of the corrected pseudoranges (clock/pr_correction are common to both
        # bands and, since the coefficients sum to 1, pass through unchanged).
        cpr_if = (f1 * f1 * o1.corrected_pseudorange_m - f2 * f2 * o2.corrected_pseudorange_m) / denom
        new_obs.append(SatObs(
            constellation=cons, prn=prn, frequency_hz=f1, code_type="IF",
            pseudorange_m=cpr_if, sv_pos=o1.sv_pos, sv_clock_bias_m=0.0,
            sv_vel=o1.sv_vel, pr_correction_m=0.0,
            doppler_hz=o1.doppler_hz, cn0_dbhz=min(_cn0(o1), _cn0(o2)),
            loi=bool(o1.loi or o2.loi), iono_delay_m=0.0,
        ))
    return Epoch(t_sec=epoch.t_sec, obs=new_obs, gps_week=epoch.gps_week, tow_sec=epoch.tow_sec)


def form_iono_free(epochs: Sequence[Epoch]) -> List[Epoch]:
    return [form_iono_free_epoch(ep) for ep in epochs]

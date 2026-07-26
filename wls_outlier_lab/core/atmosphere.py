"""Atmospheric delay models: Klobuchar ionosphere + Saastamoinen troposphere.

Both return a one-way slant delay in metres to be *subtracted* from the
pseudorange. Pure functions; no I/O.

Klobuchar uses the broadcast alpha/beta coefficients (present in the measurement
TSV as iono_a0..a3 / iono_b0..b3) and follows IS-GPS-200. The L1 delay is scaled
by (f_L1 / f)^2 for other frequency bands.

Saastamoinen uses a standard-atmosphere approximation of pressure/temperature/
humidity at the receiver height with a simple 1/cos(zenith) mapping — the
well-worn RTKLIB `tropmodel` form, good to a few cm at moderate elevation.
"""

from __future__ import annotations

import math
from typing import Sequence

C_LIGHT = 299792458.0
F_L1_HZ = 1575.42e6


def klobuchar_iono_delay_m(lat_deg: float, lon_deg: float, elev_deg: float, azim_deg: float,
                           gps_tow_sec: float, alpha: Sequence[float], beta: Sequence[float],
                           freq_hz: float = F_L1_HZ) -> float:
    """Broadcast (Klobuchar) ionospheric slant delay [m]. 0 if coeffs absent."""
    if elev_deg <= 0.0 or freq_hz <= 0.0:
        return 0.0
    if not any(alpha) and not any(beta):
        return 0.0

    phi_u = lat_deg / 180.0          # semicircles
    lam_u = lon_deg / 180.0
    el_sc = elev_deg / 180.0         # elevation in semicircles
    az = math.radians(azim_deg)

    psi = 0.0137 / (el_sc + 0.11) - 0.022                  # earth-centred angle (sc)
    phi_i = phi_u + psi * math.cos(az)
    phi_i = max(min(phi_i, 0.416), -0.416)
    lam_i = lam_u + psi * math.sin(az) / math.cos(phi_i * math.pi)
    phi_m = phi_i + 0.064 * math.cos((lam_i - 1.617) * math.pi)

    t = 43200.0 * lam_i + gps_tow_sec
    t = t % 86400.0
    if t < 0:
        t += 86400.0

    amp = sum(alpha[n] * phi_m ** n for n in range(4))
    if amp < 0:
        amp = 0.0
    per = sum(beta[n] * phi_m ** n for n in range(4))
    if per < 72000.0:
        per = 72000.0

    x = 2.0 * math.pi * (t - 50400.0) / per
    slant_factor = 1.0 + 16.0 * (0.53 - el_sc) ** 3
    if abs(x) < 1.57:
        delay_s = slant_factor * (5.0e-9 + amp * (1.0 - x * x / 2.0 + x ** 4 / 24.0))
    else:
        delay_s = slant_factor * 5.0e-9

    delay_m = delay_s * C_LIGHT                    # L1 delay
    return delay_m * (F_L1_HZ / freq_hz) ** 2      # scale by 1/f^2 to this band


def saastamoinen_tropo_delay_m(lat_deg: float, height_m: float, elev_deg: float,
                               humidity: float = 0.7) -> float:
    """Saastamoinen tropospheric slant delay [m] (standard atmosphere)."""
    if elev_deg <= 0.0:
        return 0.0
    hgt = min(max(height_m, 0.0), 5000.0)
    el = max(elev_deg, 3.0)  # floor the mapping so low sats don't explode
    pres = 1013.25 * (1.0 - 2.2557e-5 * hgt) ** 5.2568
    temp = 15.0 - 6.5e-3 * hgt + 273.16
    e = 6.108 * humidity * math.exp((17.15 * temp - 4684.0) / (temp - 38.45))
    z = math.radians(90.0 - el)
    lat = math.radians(lat_deg)
    trph = 0.0022768 * pres / (1.0 - 0.00266 * math.cos(2.0 * lat) - 2.8e-7 * hgt) / math.cos(z)
    trpw = 0.002277 * (1255.0 / temp + 0.05) * e / math.cos(z)
    return trph + trpw

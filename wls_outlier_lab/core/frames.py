"""Coordinate-frame math (WGS84). Pure numpy, no I/O.

Kept self-contained so the lab does not depend on the rest of the repo, but the
conventions match ``smartphone_ekf_api.geo`` and ``api_server.truth_alignment``.
"""

from __future__ import annotations

import math
from typing import Tuple

import numpy as np

# WGS84
WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)
WGS84_B = WGS84_A * (1.0 - WGS84_F)

C_LIGHT = 299792458.0
OMEGA_EARTH = 7.2921151467e-5  # rad/s, WGS84 earth rotation rate


def lla_to_ecef(lat_deg: float, lon_deg: float, alt_m: float) -> np.ndarray:
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    sin_lat, cos_lat = math.sin(lat), math.cos(lat)
    n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)
    x = (n + alt_m) * cos_lat * math.cos(lon)
    y = (n + alt_m) * cos_lat * math.sin(lon)
    z = (n * (1.0 - WGS84_E2) + alt_m) * sin_lat
    return np.array([x, y, z], dtype=float)


def ecef_to_lla(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """ECEF -> (lat_deg, lon_deg, alt_m) via Bowring fixed-point iteration."""
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    if p < 1e-9:  # at a pole
        lat = math.copysign(math.pi / 2.0, z)
        alt = abs(z) - WGS84_B
        return math.degrees(lat), math.degrees(lon), alt
    lat = math.atan2(z, p * (1.0 - WGS84_E2))
    for _ in range(8):
        sin_lat = math.sin(lat)
        n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)
        alt = p / math.cos(lat) - n
        lat = math.atan2(z, p * (1.0 - WGS84_E2 * n / (n + alt)))
    sin_lat = math.sin(lat)
    n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)
    alt = p / math.cos(lat) - n
    return math.degrees(lat), math.degrees(lon), alt


def enu_rotation_matrix(lat_deg: float, lon_deg: float) -> np.ndarray:
    """3x3 matrix R such that enu = R @ (ecef_point - ecef_origin)."""
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    sl, cl = math.sin(lat), math.cos(lat)
    so, co = math.sin(lon), math.cos(lon)
    return np.array(
        [
            [-so, co, 0.0],
            [-sl * co, -sl * so, cl],
            [cl * co, cl * so, sl],
        ],
        dtype=float,
    )


def ecef_to_enu(point_ecef: np.ndarray, origin_ecef: np.ndarray) -> np.ndarray:
    lat, lon, _ = ecef_to_lla(*origin_ecef)
    return enu_rotation_matrix(lat, lon) @ (np.asarray(point_ecef, float) - np.asarray(origin_ecef, float))


def elevation_azimuth_deg(rx_ecef: np.ndarray, sv_ecef: np.ndarray) -> Tuple[float, float]:
    """Elevation and azimuth (deg) of a satellite as seen from the receiver."""
    lat, lon, _ = ecef_to_lla(*rx_ecef)
    los = np.asarray(sv_ecef, float) - np.asarray(rx_ecef, float)
    enu = enu_rotation_matrix(lat, lon) @ los
    e, n, u = enu
    horiz = math.hypot(e, n)
    elev = math.degrees(math.atan2(u, horiz))
    azim = math.degrees(math.atan2(e, n)) % 360.0
    return elev, azim


def sagnac_corrected_sv(sv_ecef: np.ndarray, travel_time_s: float) -> np.ndarray:
    """Rotate the satellite ECEF position into the reception-time frame.

    Compensates earth rotation during signal travel time (Sagnac effect).
    """
    theta = OMEGA_EARTH * travel_time_s
    ct, st = math.cos(theta), math.sin(theta)
    x, y, z = sv_ecef
    return np.array([ct * x + st * y, -st * x + ct * y, z], dtype=float)

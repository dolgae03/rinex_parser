"""Coordinate-frame round trips."""

from __future__ import annotations

import numpy as np

from wls_outlier_lab.core import frames


def test_ecef_lla_roundtrip():
    for lla in [(37.2557, 127.0553, 40.0), (0.0, 0.0, 0.0), (-33.9, 151.2, 58.0),
                (64.1, -21.9, 12.0)]:
        ecef = frames.lla_to_ecef(*lla)
        back = frames.ecef_to_lla(*ecef)
        assert abs(back[0] - lla[0]) < 1e-7
        assert abs(back[1] - lla[1]) < 1e-7
        assert abs(back[2] - lla[2]) < 1e-3


def test_zenith_satellite_elevation():
    rx = frames.lla_to_ecef(37.2557, 127.0553, 40.0)
    up = rx / np.linalg.norm(rx)
    sat = rx + up * 2.0e7  # straight up
    elev, _ = frames.elevation_azimuth_deg(rx, sat)
    assert elev > 89.5


def test_horizon_satellite_elevation():
    rx = frames.lla_to_ecef(37.2557, 127.0553, 40.0)
    R = frames.enu_rotation_matrix(*frames.ecef_to_lla(*rx)[:2])
    east = R.T @ np.array([1.0, 0.0, 0.0])
    sat = rx + east * 2.0e7
    elev, azim = frames.elevation_azimuth_deg(rx, sat)
    assert abs(elev) < 1.0
    assert abs(azim - 90.0) < 1.0

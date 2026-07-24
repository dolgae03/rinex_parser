"""Data-access adapters: measurement TSV grouping and BESTPOS parsing."""

from __future__ import annotations

from pathlib import Path

from wls_outlier_lab.sources.measurement_tsv import load_epochs
from wls_outlier_lab.sources.truth_bestpos import load_bestpos_track, parse_bestpos

TSV = (
    "t_sec\tgps_week\ttow_sec\tconstellation\tprn\tsv_pos_x\tsv_pos_y\tsv_pos_z\t"
    "sv_clock_bias\tpr_correction\tfrequency_hz\tcode_type\tpseudorange_m\t"
    "doppler_hz\tsnr_dbhz\tloi\n"
    "100\t2413\t100\t0\t4\t-2.6e7\t4.0e6\t3.8e6\t15000\t\t1575420000\tC\t23459773\t100\t45\t0\n"
    "100\t2413\t100\t0\t14\t-1.2e7\t4.0e6\t3.8e6\t222000\t\t1575420000\tC\t21406411\t90\t46\t0\n"
    "101\t2413\t101\t2\t22\t-2.1e7\t4.0e6\t3.8e6\t134000\t\t1561098000\tC\t23569049\t80\t44\t1\n"
)

BESTPOS = (
    "#BESTPOSA,FILE,0,46.5,FINESTEERING,2425,435981.000,02000800,cdba,18018;"
    "SOL_COMPUTED,NARROW_INT,37.25569551166,127.05525940823,34.3433,23.2,WGS84,"
    "0.013,0.0136,0.031,\"3183\",1.0,0.0,45,44,44,43,00,21,3f,37*0eaf3c56\n"
    "#BESTPOSA,FILE,0,17.5,FINESTEERING,2425,435982.000,02000800,cdba,18018;"
    "SOL_COMPUTED,SINGLE,37.2,127.0,34.0,23.2,WGS84,1.0,1.0,2.0,\"\",1.0,0.0,20,20*abcd1234\n"
)


def test_load_epochs_groups_by_time(tmp_path: Path):
    f = tmp_path / "m.tsv"
    f.write_text(TSV, encoding="utf-8")
    epochs = load_epochs(f)
    assert len(epochs) == 2
    assert len(epochs[0].obs) == 2
    o = epochs[0].obs[0]
    assert o.constellation == 0 and o.prn == 4
    assert abs(o.pseudorange_m - 23459773) < 1
    assert abs(o.sv_clock_bias_m - 15000) < 1
    assert o.pr_correction_m == 0.0  # blank -> 0
    assert epochs[1].obs[0].loi is True


def test_load_epochs_stride_and_constellation_filter(tmp_path: Path):
    f = tmp_path / "m.tsv"
    f.write_text(TSV, encoding="utf-8")
    only_bds = load_epochs(f, constellations=[2])
    # epoch 100 has no BeiDou -> dropped; epoch 101 kept with the BeiDou sat
    assert len(only_bds) == 1
    assert only_bds[0].obs[0].constellation == 2


def test_parse_bestpos_filters_pos_type(tmp_path: Path):
    f = tmp_path / "b.ascii"
    f.write_text(BESTPOS, encoding="utf-8")
    narrow = parse_bestpos(f, require_pos_types=("NARROW_INT",))
    assert len(narrow) == 1
    r = narrow[0]
    assert r.week == 2425 and abs(r.tow - 435981.0) < 1e-3
    assert abs(r.t_sec - (2425 * 604800 + 435981.0)) < 1e-3
    assert abs(r.lat_deg - 37.25569551166) < 1e-9

    track = load_bestpos_track(f, require_pos_types=("NARROW_INT", "SINGLE"))
    assert len(track.samples) == 2  # both rows kept, distinct timestamps

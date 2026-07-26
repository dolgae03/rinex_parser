"""CLI entry point.

Example (validation dataset, truth embedded as gt columns):

    python -m wls_outlier_lab.cli \
        --measurements .../with_sv_pos_with_gt_rtk.tsv --truth-columns \
        --output-dir results/2026_04_06 --stride 5

Example (Novatel BESTPOS as truth for a phone measurement TSV):

    python -m wls_outlier_lab.cli \
        --measurements phone_with_sv_pos.tsv \
        --truth-bestpos .../NMND..._BESTPOS.ASCII \
        --output-dir results/session
"""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import replace
from pathlib import Path

from .core.calibration import calibrate, measurement_residuals, outlier_catalog
from .core.clock_analysis import analyze_clock_stability
from .core.detectors import DetectorConfig
from .core.experiment import compare_ionosphere, constellation_ablation, run_experiment
from .core.iono_free import form_iono_free
from .core.wls import WeightConfig, WlsConfig
from .reporting import write_iono_comparison, write_report
from .sources.measurement_tsv import load_epochs
from .sources.truth_bestpos import BestposTruthSource
from .sources.truth_columns import ColumnTruthSource


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Define GNSS outliers and measure the WLS improvement vs truth.")
    p.add_argument("--measurements", required=True, help="measurement TSV/CSV with SV positions")
    p.add_argument("--output-dir", required=True)

    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--truth-columns", action="store_true",
                   help="use gt_pos_* columns embedded in the measurements file")
    g.add_argument("--truth-file", help="separate file carrying gt_pos_* columns")
    g.add_argument("--truth-bestpos", help="Novatel BESTPOS ASCII file")

    p.add_argument("--bestpos-pos-types", default="NARROW_INT,INS_RTKFIXED",
                   help="comma list of BESTPOS pos types to accept (RTK-fixed default)")
    p.add_argument("--stride", type=int, default=1, help="keep every k-th epoch")
    p.add_argument("--max-epochs", type=int, default=None)
    p.add_argument("--match-tolerance", type=float, default=0.5, help="truth match window [s]")
    p.add_argument("--detectors", default=None,
                   help="comma list; default = all registered detectors")

    p.add_argument("--min-cn0", type=float, default=30.0)
    p.add_argument("--min-elevation", type=float, default=10.0)
    p.add_argument("--residual-k", type=float, default=5.0,
                   help="robust residual snooping threshold (MAD multiplier)")
    p.add_argument("--mad-floor", type=float, default=6.0, help="floor on robust residual scale [m]")

    p.add_argument("--base-sigma", type=float, default=3.0)
    p.add_argument("--no-elevation-weighting", action="store_true")
    p.add_argument("--no-cn0-weighting", action="store_true")
    p.add_argument("--no-sagnac", action="store_true")
    p.add_argument("--no-tropo", action="store_true", help="disable Saastamoinen troposphere")
    p.add_argument("--iono", action="store_true",
                   help="enable broadcast Klobuchar ionosphere (crude single-freq; see README)")
    p.add_argument("--iono-free", action="store_true",
                   help="run on the dual-frequency (L1/L5) ionosphere-free combination")
    p.add_argument("--compare-iono", action="store_true",
                   help="compare the nav solution under none/tropo/klobuchar/iono-free")
    p.add_argument("--clock-analysis", action="store_true",
                   help="analyze receiver clock-drift instability vs horizontal error")

    p.add_argument("--no-ablation", action="store_true")
    p.add_argument("--no-calibrate", action="store_true",
                   help="skip the truth-referenced measurement-noise calibration")
    # Plotting is done in MATLAB (plot_wls_results.m) from the CSV outputs.
    # This flag additionally emits the optional matplotlib PNGs.
    p.add_argument("--matplotlib-plots", action="store_true",
                   help="also write matplotlib PNGs (default: CSV only; plot in MATLAB)")
    return p


def main(argv=None) -> int:
    args = build_arg_parser().parse_args(argv)
    t0 = time.time()

    print(f"[load] measurements: {args.measurements}")
    epochs = load_epochs(args.measurements, stride=args.stride, max_epochs=args.max_epochs)
    print(f"[load] {len(epochs)} epochs")
    if not epochs:
        print("no epochs loaded", file=sys.stderr)
        return 2

    if args.truth_bestpos:
        pos_types = tuple(s.strip() for s in args.bestpos_pos_types.split(",") if s.strip())
        truth = BestposTruthSource(args.truth_bestpos, require_pos_types=pos_types).load_track()
        truth_desc = f"bestpos:{args.truth_bestpos}"
    elif args.truth_file:
        truth = ColumnTruthSource(args.truth_file).load_track()
        truth_desc = f"columns:{args.truth_file}"
    else:
        truth = ColumnTruthSource(args.measurements).load_track()
        truth_desc = "columns:measurements"
    print(f"[load] truth: {truth.name} ({truth.source_type}) {len(truth.samples)} samples")

    wls_cfg = WlsConfig(
        weights=WeightConfig(
            base_sigma_m=args.base_sigma,
            use_elevation=not args.no_elevation_weighting,
            use_cn0=not args.no_cn0_weighting,
        ),
        apply_sagnac=not args.no_sagnac,
        apply_tropo=not args.no_tropo,
        apply_iono=args.iono,
    )
    dcfg = DetectorConfig(
        min_cn0_dbhz=args.min_cn0,
        min_elevation_deg=args.min_elevation,
        residual_mad_k=args.residual_k,
        mad_floor_m=args.mad_floor,
    )
    detectors = [s.strip() for s in args.detectors.split(",")] if args.detectors else None

    iono_cmp = None
    if args.compare_iono:
        print("[run] ionosphere comparison (none / tropo / klobuchar / iono_free) ...")
        iono_cmp = compare_ionosphere(epochs, truth, wls_cfg=wls_cfg, dcfg=dcfg,
                                      match_tolerance_sec=args.match_tolerance)

    if args.iono_free:
        epochs = form_iono_free(epochs)
        wls_cfg = replace(wls_cfg, apply_iono=False)
        n_obs = sum(len(e.obs) for e in epochs)
        print(f"[iono] iono-free (L1/L5): {n_obs} IF obs across {len(epochs)} epochs")

    print("[run] detector comparison ...")
    result = run_experiment(epochs, truth, wls_cfg=wls_cfg, dcfg=dcfg,
                            detectors=detectors, match_tolerance_sec=args.match_tolerance)

    ablation = None
    if not args.no_ablation:
        print("[run] constellation ablation ...")
        ablation = constellation_ablation(epochs, truth, wls_cfg=wls_cfg,
                                          match_tolerance_sec=args.match_tolerance)

    calibration = None
    catalog = None
    if not args.no_calibrate:
        print("[run] measurement-noise calibration + blunder catalog (vs truth) ...")
        residuals = measurement_residuals(epochs, truth, wls_cfg=wls_cfg,
                                          match_tolerance_sec=args.match_tolerance)
        calibration = calibrate(residuals)
        catalog = outlier_catalog(residuals, k=args.residual_k, mad_floor_m=args.mad_floor)

    clock = None
    if args.clock_analysis:
        print("[run] clock-drift stability analysis ...")
        clock = analyze_clock_stability(epochs, truth, wls_cfg=wls_cfg, dcfg=dcfg,
                                        match_tolerance_sec=args.match_tolerance)

    meta = {
        "measurements": str(args.measurements),
        "truth": truth_desc,
        "n_epochs": len(epochs),
        "stride": args.stride,
        "wls": {"base_sigma_m": args.base_sigma,
                "elevation_weighting": not args.no_elevation_weighting,
                "cn0_weighting": not args.no_cn0_weighting,
                "sagnac": not args.no_sagnac,
                "tropo": not args.no_tropo,
                "iono_klobuchar": args.iono,
                "iono_free": args.iono_free},
        "detector_config": {"min_cn0_dbhz": args.min_cn0,
                            "min_elevation_deg": args.min_elevation,
                            "residual_mad_k": args.residual_k,
                            "mad_floor_m": args.mad_floor},
        "runtime_sec": None,
    }
    meta["runtime_sec"] = round(time.time() - t0, 2)
    paths = write_report(args.output_dir, result, ablation_rows=ablation, meta=meta,
                         make_plots=args.matplotlib_plots, calibration=calibration,
                         catalog=catalog, clock=clock)
    if iono_cmp is not None:
        paths.update(write_iono_comparison(args.output_dir, iono_cmp))

    print("\n=== Detector comparison (best first) ===")
    print(f"{'detector':<24}{'hRMSE[m]':>10}{'h95[m]':>9}{'vRMSE[m]':>10}"
          f"{'avail%':>8}{'sats':>6}{'rej':>7}{'impr%':>8}")
    for r in result.comparison():
        print(f"{r['detector']:<24}{_n(r['horizontal_rmse_m']):>10}{_n(r['horizontal_p95_m']):>9}"
              f"{_n(r['vertical_rmse_m']):>10}{_n(r['availability_pct'],1):>8}"
              f"{_n(r['mean_sats_used'],1):>6}{r['signal_rejections']:>7}{_n(r['improvement_pct'],1):>8}")
    if clock is not None and clock.n_epochs:
        ppm = clock.median_drift_mps / 299792458.0 * 1e6 if clock.median_drift_mps == clock.median_drift_mps else float("nan")
        print("\n=== Clock-drift stability vs horizontal error ===")
        print(f"  median drift {_n(clock.median_drift_mps,1)} m/s (~{_n(ppm,3)} ppm); "
              f"{clock.n_anomalies}/{clock.n_epochs} anomalies (instability > {_n(clock.instability_threshold_m,1)} m)")
        print(f"  corr(instability, horizontal err): pearson={_n(clock.pearson_r)}  spearman={_n(clock.spearman_r)}")
        print(f"  horizontal error   stable={_n(clock.stable_mean_h_m)} m   unstable={_n(clock.unstable_mean_h_m)} m")
        sr = clock.spearman_r
        verdict = ("weak/no correlation -> clock-drift instability does NOT drive horizontal error"
                   if not (sr == sr) or abs(sr) < 0.2 else
                   "correlation present -> clock instability may affect the horizontal solution")
        print(f"  verdict: {verdict}")

    if iono_cmp is not None:
        print("\n=== Horizontal nav solution by ionosphere treatment (detector=combined) ===")
        print(f"{'mode':<12}{'hRMSE[m]':>10}{'hCEP95[m]':>11}{'vRMSE[m]':>10}{'vMean[m]':>10}{'avail%':>8}{'sats':>6}")
        for mode, m in iono_cmp.items():
            print(f"{mode:<12}{_n(m.horizontal_rmse_m):>10}{_n(m.cep95_m):>11}{_n(m.vertical_rmse_m):>10}"
                  f"{_n(m.vertical_mean_m):>10}{_n(m.availability_pct,1):>8}{_n(m.mean_sats_used,1):>6}")

    if calibration is not None:
        print("\n=== Measurement noise vs truth (empirical, ref=GPS) ===")
        print(f"{'constellation':<14}{'n':>7}{'clean_std[m]':>13}{'outlier%':>10}{'sigma_scale':>12}")
        for name, s in calibration.by_constellation.items():
            if s.get("n", 0) == 0:
                continue
            print(f"{name:<14}{s['n']:>7}{_n(s.get('clean_std_m'),1):>13}"
                  f"{_n(100*s.get('outlier_rate',0),1):>10}"
                  f"{_n(calibration.sigma_scale_by_constellation.get(name),2):>12}")

    if catalog is not None and catalog.get("n"):
        print(f"\n=== Blunder catalog (truth-referenced, |z|>{catalog['threshold_k']}, "
              f"scale={catalog['scale_m']} m) ===")
        print(f"  {catalog['n_outliers']}/{catalog['n']} obs flagged "
              f"({100*catalog['n_outliers']/catalog['n']:.1f}%)")
        for c, s in catalog["by_constellation"].items():
            print(f"    {c:<9} {s['n_outliers']:>5}/{s['n']:<6} ({100*s['rate']:.1f}%)")
        if catalog["worst_satellites"]:
            print("  worst satellites:", ", ".join(
                f"{w['satellite']}({w['n_outliers']})" for w in catalog["worst_satellites"][:8]))

    print(f"\n[done] outputs in {args.output_dir}  ({meta['runtime_sec']}s)")
    for k, v in paths.items():
        print(f"  {k}: {v}")
    print("\n[plot] in MATLAB:")
    print(f"  addpath('wls_outlier_lab'); plot_wls_results('{args.output_dir}')")
    return 0


def _n(x, nd: int = 3) -> str:
    try:
        if x != x:  # NaN
            return "nan"
        return f"{x:.{nd}f}"
    except (TypeError, ValueError):
        return str(x)


if __name__ == "__main__":
    raise SystemExit(main())

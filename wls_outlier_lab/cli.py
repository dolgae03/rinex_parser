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
from pathlib import Path

from .core.calibration import calibrate, measurement_residuals
from .core.detectors import DetectorConfig
from .core.experiment import constellation_ablation, run_experiment
from .core.wls import WeightConfig, WlsConfig
from .reporting import write_report
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
    )
    dcfg = DetectorConfig(
        min_cn0_dbhz=args.min_cn0,
        min_elevation_deg=args.min_elevation,
        residual_mad_k=args.residual_k,
        mad_floor_m=args.mad_floor,
    )
    detectors = [s.strip() for s in args.detectors.split(",")] if args.detectors else None

    print("[run] detector comparison ...")
    result = run_experiment(epochs, truth, wls_cfg=wls_cfg, dcfg=dcfg,
                            detectors=detectors, match_tolerance_sec=args.match_tolerance)

    ablation = None
    if not args.no_ablation:
        print("[run] constellation ablation ...")
        ablation = constellation_ablation(epochs, truth, wls_cfg=wls_cfg,
                                          match_tolerance_sec=args.match_tolerance)

    calibration = None
    if not args.no_calibrate:
        print("[run] measurement-noise calibration (vs truth) ...")
        residuals = measurement_residuals(epochs, truth, wls_cfg=wls_cfg,
                                          match_tolerance_sec=args.match_tolerance)
        calibration = calibrate(residuals)

    meta = {
        "measurements": str(args.measurements),
        "truth": truth_desc,
        "n_epochs": len(epochs),
        "stride": args.stride,
        "wls": {"base_sigma_m": args.base_sigma,
                "elevation_weighting": not args.no_elevation_weighting,
                "cn0_weighting": not args.no_cn0_weighting,
                "sagnac": not args.no_sagnac},
        "detector_config": {"min_cn0_dbhz": args.min_cn0,
                            "min_elevation_deg": args.min_elevation,
                            "residual_mad_k": args.residual_k,
                            "mad_floor_m": args.mad_floor},
        "runtime_sec": None,
    }
    meta["runtime_sec"] = round(time.time() - t0, 2)
    paths = write_report(args.output_dir, result, ablation_rows=ablation, meta=meta,
                         make_plots=args.matplotlib_plots, calibration=calibration)

    print("\n=== Detector comparison (best first) ===")
    print(f"{'detector':<24}{'hRMSE[m]':>10}{'h95[m]':>9}{'vRMSE[m]':>10}"
          f"{'avail%':>8}{'sats':>6}{'rej':>7}{'impr%':>8}")
    for r in result.comparison():
        print(f"{r['detector']:<24}{_n(r['horizontal_rmse_m']):>10}{_n(r['horizontal_p95_m']):>9}"
              f"{_n(r['vertical_rmse_m']):>10}{_n(r['availability_pct'],1):>8}"
              f"{_n(r['mean_sats_used'],1):>6}{r['signal_rejections']:>7}{_n(r['improvement_pct'],1):>8}")
    if calibration is not None:
        print("\n=== Measurement noise vs truth (empirical, ref=GPS) ===")
        print(f"{'constellation':<14}{'n':>7}{'clean_std[m]':>13}{'outlier%':>10}{'sigma_scale':>12}")
        for name, s in calibration.by_constellation.items():
            if s.get("n", 0) == 0:
                continue
            print(f"{name:<14}{s['n']:>7}{_n(s.get('clean_std_m'),1):>13}"
                  f"{_n(100*s.get('outlier_rate',0),1):>10}"
                  f"{_n(calibration.sigma_scale_by_constellation.get(name),2):>12}")

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

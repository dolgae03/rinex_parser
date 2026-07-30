"""Write experiment outputs: JSON summary, CSV tables, and optional plots.

The only module besides ``sources`` that touches disk. Plots are best-effort:
if matplotlib is missing the tables are still written.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Dict, List, Optional

from .core.experiment import ExperimentResult
from .core.metrics import ErrorMetrics
from .types import constellation_name


def _metrics_dict(m: ErrorMetrics) -> Dict[str, object]:
    return m.to_summary()


def _r(x, nd: int = 2):
    if x is None:
        return ""
    try:
        if x != x:  # NaN
            return ""
        return round(float(x), nd)
    except (TypeError, ValueError):
        return x


def write_report(
    output_dir: str | Path,
    result: ExperimentResult,
    ablation_rows: Optional[List[Dict[str, object]]] = None,
    meta: Optional[Dict[str, object]] = None,
    make_plots: bool = True,
    calibration: Optional[object] = None,
    catalog: Optional[Dict[str, object]] = None,
    clock: Optional[object] = None,
    coasting: Optional[List[Dict[str, object]]] = None,
) -> Dict[str, str]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths: Dict[str, str] = {}

    comparison = result.comparison()
    best = comparison[0]["detector"] if comparison else result.baseline_name

    # --- summary.json -----------------------------------------------------
    summary = {
        "meta": meta or {},
        "baseline_detector": result.baseline_name,
        "best_detector": best,
        "measurement_calibration": calibration.to_dict() if calibration is not None else None,
        "outlier_catalog_summary": (
            {k: v for k, v in catalog.items() if k != "rows"} if catalog is not None else None
        ),
        "clock_stability": clock.summary() if clock is not None else None,
        "clock_coasting": coasting,
        "detector_comparison": comparison,
        "detectors": {
            name: {
                "metrics": _metrics_dict(run.metrics),
                "signal_rejections": run.n_signal_rejections,
                "epochs_with_rejection": run.epochs_with_rejection,
                "rejections_by_constellation": {
                    constellation_name(c): n
                    for c, n in sorted(run.rejections_by_constellation.items())
                },
            }
            for name, run in result.runs.items()
        },
        "constellation_ablation": ablation_rows or [],
    }
    summary_path = out / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    paths["summary"] = str(summary_path)

    # --- detector_comparison.csv -----------------------------------------
    if comparison:
        comp_path = out / "detector_comparison.csv"
        with comp_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(comparison[0].keys()))
            w.writeheader()
            w.writerows(comparison)
        paths["detector_comparison"] = str(comp_path)

    # --- constellation_ablation.csv --------------------------------------
    if ablation_rows:
        keys: List[str] = []
        for r in ablation_rows:
            for k in r:
                if k not in keys:
                    keys.append(k)
        abl_path = out / "constellation_ablation.csv"
        with abl_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(ablation_rows)
        paths["constellation_ablation"] = str(abl_path)

    # --- measurement calibration (empirical noise vs truth) --------------
    if calibration is not None:
        cal = calibration
        cons_path = out / "calibration_by_constellation.csv"
        with cons_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["constellation", "n", "median_m", "mad_sigma_m", "clean_std_m",
                        "rms_m", "outlier_rate", "sigma_scale_vs_ref"])
            for name, s in cal.by_constellation.items():
                w.writerow([name, s.get("n", 0), _r(s.get("median_m")), _r(s.get("mad_sigma_m")),
                            _r(s.get("clean_std_m")), _r(s.get("rms_m")), _r(s.get("outlier_rate"), 4),
                            cal.sigma_scale_by_constellation.get(name)])
        paths["calibration_by_constellation"] = str(cons_path)
        for tag, rows in (("cn0", cal.by_cn0_bin), ("elevation", cal.by_elevation_bin)):
            bp = out / f"calibration_by_{tag}.csv"
            with bp.open("w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh)
                w.writerow(["bin_lo", "bin_hi", "n", "clean_std_m", "mad_sigma_m", "rms_m", "outlier_rate"])
                for r in rows:
                    w.writerow([r.get("bin_lo"), r.get("bin_hi"), r.get("n", 0),
                                _r(r.get("clean_std_m")), _r(r.get("mad_sigma_m")),
                                _r(r.get("rms_m")), _r(r.get("outlier_rate"), 4)])
            paths[f"calibration_by_{tag}"] = str(bp)

    # --- blunder catalog (truth-referenced, all obs + is_outlier flag) ---
    if catalog is not None and catalog.get("rows"):
        rows = catalog["rows"]
        cat_path = out / "blunder_catalog.csv"
        with cat_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        paths["blunder_catalog"] = str(cat_path)

    # --- clock stability (per-epoch clock, drift, instability) -----------
    if clock is not None and clock.points:
        clk_path = out / "clock_analysis.csv"
        with clk_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["t_sec", "dt_s", "clk_m", "drift_mps", "instability_m",
                        "is_anomaly", "horizontal_error_m", "n_sats", "hdop",
                        "clk_sigma_m", "drift_dop_mps", "drift_dop_change_mps",
                        "is_dop_anomaly"])
            for p in clock.points:
                w.writerow([f"{p.t_sec:.3f}", _r(p.dt_s, 3), _r(p.clk_m, 3), _r(p.drift_mps, 4),
                            _r(p.instability_m, 4), int(p.is_anomaly), _r(p.horizontal_error_m, 4),
                            p.n_sats, _r(p.hdop, 3), _r(p.clk_sigma_m, 4),
                            _r(p.drift_dop_mps, 4), _r(p.drift_dop_change_mps, 4),
                            int(p.is_dop_anomaly)])
        paths["clock_analysis"] = str(clk_path)

    # --- clock coasting sweep (does the fix suffer when it trusts the clock?) ---
    if coasting:
        co_path = out / "clock_coasting.csv"
        with co_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(coasting[0].keys()))
            w.writeheader()
            w.writerows(coasting)
        paths["clock_coasting"] = str(co_path)

    # --- per-epoch errors for baseline and best --------------------------
    for tag in {result.baseline_name, best}:
        run = result.runs.get(tag)
        if run is None:
            continue
        pe_path = out / f"per_epoch_{tag}.csv"
        with pe_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["t_sec", "east_m", "north_m", "up_m", "horizontal_m",
                        "d3_m", "n_used", "hdop", "pdop"])
            for e in run.metrics.per_epoch:
                w.writerow([f"{e.t_sec:.3f}", f"{e.east_m:.4f}", f"{e.north_m:.4f}",
                            f"{e.up_m:.4f}", f"{e.horizontal_m:.4f}", f"{e.d3_m:.4f}",
                            e.n_used, f"{e.hdop:.3f}", f"{e.pdop:.3f}"])
        paths[f"per_epoch_{tag}"] = str(pe_path)

    if make_plots:
        try:
            _write_plots(out, result, best, paths)
        except Exception as exc:  # plotting must never break the report
            paths["plot_error"] = repr(exc)

    return paths


def write_factor_epochs(output_dir: str | Path, rows: List[Dict[str, float]]) -> Dict[str, str]:
    """Per-epoch factor matrix consumed by MATLAB (factor_correlation_analysis.m)."""
    from .core.factor_features import FACTOR_COLUMNS

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "factor_epochs.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FACTOR_COLUMNS, extrasaction="ignore",
                           restval="")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if isinstance(v, float) and v != v else v)
                        for k, v in r.items()})
    return {"factor_epochs": str(path)}


def write_iono_comparison(output_dir: str | Path, metrics_by_mode: Dict[str, object]) -> Dict[str, str]:
    """Write the ionosphere-treatment comparison table + per-mode per-epoch errors."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths: Dict[str, str] = {}
    comp_path = out / "iono_comparison.csv"
    with comp_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["mode", "horizontal_rmse_m", "horizontal_cep95_m", "horizontal_p95_m",
                    "vertical_rmse_m", "vertical_mean_m", "availability_pct", "mean_sats_used"])
        for mode, m in metrics_by_mode.items():
            w.writerow([mode, _r(m.horizontal_rmse_m, 3), _r(m.cep95_m, 3), _r(m.horizontal_p95_m, 3),
                        _r(m.vertical_rmse_m, 3), _r(m.vertical_mean_m, 3),
                        _r(m.availability_pct, 2), _r(m.mean_sats_used, 2)])
    paths["iono_comparison"] = str(comp_path)
    for mode, m in metrics_by_mode.items():
        pe = out / f"per_epoch_iono_{mode}.csv"
        with pe.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["t_sec", "east_m", "north_m", "up_m", "horizontal_m"])
            for e in m.per_epoch:
                w.writerow([f"{e.t_sec:.3f}", f"{e.east_m:.4f}", f"{e.north_m:.4f}",
                            f"{e.up_m:.4f}", f"{e.horizontal_m:.4f}"])
        paths[f"per_epoch_iono_{mode}"] = str(pe)
    return paths


def _write_plots(out: Path, result: ExperimentResult, best: str, paths: Dict[str, str]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    comparison = result.comparison()
    base = result.baseline_name

    import numpy as np

    # 1) detector comparison bar chart (log scale: 2 m vs 100s of km both visible)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    names = [r["detector"] for r in comparison]
    vals = [r["horizontal_rmse_m"] for r in comparison]
    bars = ax.bar(names, vals, color=["#c0392b" if n == best else "#7f8c8d" for n in names])
    ax.set_yscale("log")
    ax.set_ylabel("Horizontal RMSE [m] (log)")
    ax.set_title("WLS horizontal RMSE by outlier detector (vs truth)")
    ax.tick_params(axis="x", rotation=30)
    for b, v in zip(bars, vals):
        if v == v:
            ax.annotate(f"{v:.2f}" if v < 100 else f"{v/1000:.0f} km",
                        (b.get_x() + b.get_width() / 2, v), ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    p = out / "detector_comparison.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    paths["plot_comparison"] = str(p)

    # 2) Horizontal + vertical error of the best (cleaned) solution vs truth
    xr = result.runs.get(best)
    if xr and xr.metrics.per_epoch:
        ep = xr.metrics.per_epoch
        t0 = min(e.t_sec for e in ep)
        t = np.array([e.t_sec - t0 for e in ep])
        east = np.array([e.east_m for e in ep])
        north = np.array([e.north_m for e in ep])
        up = np.array([e.up_m for e in ep])
        horiz = np.hypot(east, north)

        fig, axs = plt.subplots(2, 2, figsize=(13, 9))
        m = xr.metrics

        axs[0, 0].scatter(east, north, s=7, alpha=0.5, color="#c0392b")
        axs[0, 0].axhline(0, color="k", lw=0.5); axs[0, 0].axvline(0, color="k", lw=0.5)
        axs[0, 0].set_aspect("equal", adjustable="datalim")
        axs[0, 0].set_xlabel("East error [m]"); axs[0, 0].set_ylabel("North error [m]")
        axs[0, 0].set_title(f"Horizontal scatter  (RMSE {m.horizontal_rmse_m:.2f} m, "
                            f"CEP95 {m.cep95_m:.2f} m)  [{best}]")

        axs[0, 1].plot(t, horiz, lw=0.8, color="#c0392b")
        axs[0, 1].axhline(m.horizontal_rmse_m, color="k", ls="--", lw=0.8,
                          label=f"RMSE {m.horizontal_rmse_m:.2f} m")
        axs[0, 1].set_xlabel("time since start [s]"); axs[0, 1].set_ylabel("horizontal error [m]")
        axs[0, 1].set_title("Horizontal error vs time"); axs[0, 1].legend()

        axs[1, 0].plot(t, up, lw=0.8, color="#2c6fbb")
        axs[1, 0].axhline(0, color="k", lw=0.5)
        axs[1, 0].axhline(m.vertical_mean_m, color="k", ls="--", lw=0.8,
                          label=f"mean {m.vertical_mean_m:.2f} m, RMSE {m.vertical_rmse_m:.2f} m")
        axs[1, 0].set_xlabel("time since start [s]"); axs[1, 0].set_ylabel("vertical (up) error [m]")
        axs[1, 0].set_title("Vertical error vs time"); axs[1, 0].legend()

        for label, data, color in (("horizontal", horiz, "#c0392b"),
                                    ("vertical |U|", np.abs(up), "#2c6fbb")):
            sd = np.sort(data)
            axs[1, 1].plot(sd, np.linspace(0, 100, len(sd)), label=label, color=color)
        axs[1, 1].set_xlabel("error [m]"); axs[1, 1].set_ylabel("percentile [%]")
        axs[1, 1].set_title("Error CDF"); axs[1, 1].grid(alpha=0.3); axs[1, 1].legend()

        fig.suptitle(f"Positioning error vs truth — detector '{best}'", fontsize=13)
        fig.tight_layout()
        p2 = out / "error_horizontal_vertical.png"
        fig.savefig(p2, dpi=130)
        plt.close(fig)
        paths["plot_error"] = str(p2)

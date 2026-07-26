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


def write_report(
    output_dir: str | Path,
    result: ExperimentResult,
    ablation_rows: Optional[List[Dict[str, object]]] = None,
    meta: Optional[Dict[str, object]] = None,
    make_plots: bool = True,
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

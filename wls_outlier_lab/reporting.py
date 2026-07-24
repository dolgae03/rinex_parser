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

    # 1) detector comparison bar chart
    fig, ax = plt.subplots(figsize=(8, 4.5))
    names = [r["detector"] for r in comparison]
    vals = [r["horizontal_rmse_m"] for r in comparison]
    ax.bar(names, vals, color=["#c0392b" if n == best else "#7f8c8d" for n in names])
    ax.set_ylabel("Horizontal RMSE [m]")
    ax.set_title("WLS horizontal RMSE by outlier detector (vs truth)")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    p = out / "detector_comparison.png"
    fig.savefig(p, dpi=130)
    plt.close(fig)
    paths["plot_comparison"] = str(p)

    # 2) EN scatter + horizontal-error time series (baseline vs best)
    br = result.runs.get(base)
    xr = result.runs.get(best)
    if br and xr and br.metrics.per_epoch and xr.metrics.per_epoch:
        fig, (a0, a1) = plt.subplots(1, 2, figsize=(12, 5))
        a0.scatter([e.east_m for e in br.metrics.per_epoch],
                   [e.north_m for e in br.metrics.per_epoch],
                   s=6, alpha=0.4, label=f"{base}", color="#7f8c8d")
        a0.scatter([e.east_m for e in xr.metrics.per_epoch],
                   [e.north_m for e in xr.metrics.per_epoch],
                   s=6, alpha=0.4, label=f"{best}", color="#c0392b")
        a0.axhline(0, color="k", lw=0.5)
        a0.axvline(0, color="k", lw=0.5)
        a0.set_aspect("equal", adjustable="datalim")
        a0.set_xlabel("East error [m]")
        a0.set_ylabel("North error [m]")
        a0.set_title("Horizontal error scatter")
        a0.legend()

        a1.plot([e.t_sec for e in br.metrics.per_epoch],
                [e.horizontal_m for e in br.metrics.per_epoch],
                lw=0.8, label=base, color="#7f8c8d")
        a1.plot([e.t_sec for e in xr.metrics.per_epoch],
                [e.horizontal_m for e in xr.metrics.per_epoch],
                lw=0.8, label=best, color="#c0392b")
        a1.set_xlabel("GPS time [s]")
        a1.set_ylabel("Horizontal error [m]")
        a1.set_title("Horizontal error vs time")
        a1.legend()
        fig.tight_layout()
        p2 = out / "error_baseline_vs_best.png"
        fig.savefig(p2, dpi=130)
        plt.close(fig)
        paths["plot_error"] = str(p2)

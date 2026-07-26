# wls_outlier_lab

Define GNSS **measurement outliers** and quantify the **WLS positioning gain** from
removing them, against a **Novatel BESTPOS / RTK truth**.

> Study question: *by satellite constellation or overall, which measurements should
> we flag as outliers, and how much does removing them improve a weighted
> least-squares fix versus the truth?*

This is a fresh engine, not an extension of `smartphone_ekf_api` (that package
smooths the phone's already-solved `Fix` positions with an EKF; it has no
pseudorange model, no WLS, and no outlier logic). Here we work at the **raw
measurement** level.

## Design: data acquisition is separated from business logic

Many datasets will flow through this over time, so *reading the data* is kept
strictly apart from *the analysis*:

```
sources/   data-access adapters   ->  Epoch / TruthTrack        (all I/O lives here)
core/      WLS · detectors · metrics · experiment               (pure, no I/O)
reporting  summary.json · CSV tables · PNG plots
cli / app  thin entry points (CLI + optional FastAPI)
```

A new dataset = a new `MeasurementSource` / `TruthSource`; the engine never changes.

### Sources (I/O)
- `measurement_tsv` — canonical measurement TSV (`*_with_sv_pos*.tsv`, 32 cols) and
  the `data/*_processed*.csv` variant produced by `gnss_txt_parser` + the MATLAB
  SV-position endpoint. Columns matched by name; SV position, clock and iono are
  already filled, so each row is WLS-ready.
- `truth_columns` — truth lifted from `gt_pos_x/y/z` columns embedded in a
  `*_with_gt_*.tsv` (used for end-to-end validation).
- `truth_bestpos` — **Novatel BESTPOS ASCII → TruthTrack** (RTK-fixed by default).
  This is the only Novatel handling that previously lived only in MATLAB
  (`parse_gps_best_pos.m`). Also exports the RTK CSV schema
  (`gps_sec, latitude_deg, longitude_deg, height_m`) that `api_server` already reads.
- `attitude` — **v1 hook only**. Phone attitude was chosen as a *future* quality
  covariate; the interface (`AttitudeSeries`, roll/pitch/yaw) is defined so a real
  parser (phone `OrientationDeg`/IMU, or Novatel `INSPVA`) drops in without touching
  the core.

### Core (pure logic)
- `wls` — single-epoch weighted least squares on pseudoranges. Per-constellation
  inter-system bias, elevation + C/N0 weighting, Sagnac correction, Gauss-Newton
  with a backtracking line search (so a gross blunder biases but never diverges the
  fix). An optional **Huber-IRLS robust mode** is used by the detectors. Observation
  model, verified numerically against the endpoint output:
  `pseudorange_m + sv_clock_bias_m (+ pr_correction_m) = range + c·clk (+ ISB)`.
- `detectors` — the multi-detector comparison: `baseline` (none), `cn0`,
  `elevation`, `residual` (robust-fit residual snooping, global scale),
  `residual_per_constellation`, `combined`. Rejection is **per satellite**
  `(constellation, prn)` — a bad ephemeris corrupts every frequency.
- `metrics` — ENU error vs truth: horizontal/vertical RMSE, p50/p95, CEP50/95,
  2DRMS, 3D RMSE, availability, mean satellites used, HDOP.
- `experiment` — `run_experiment` (baseline vs each detector, improvement table) and
  `constellation_ablation` (drop / keep-only each constellation).

## Usage

```powershell
pip install -r wls_outlier_lab\requirements.txt   # numpy (+ matplotlib, fastapi optional)

# Validation dataset — truth embedded as gt_pos_* columns
python -m wls_outlier_lab.cli `
  --measurements <..._with_sv_pos_with_gt_rtk.tsv> --truth-columns `
  --output-dir results\2026_04_06 --stride 5

# Phone measurement TSV with a Novatel BESTPOS truth
python -m wls_outlier_lab.cli `
  --measurements phone_with_sv_pos.tsv `
  --truth-bestpos <..._BESTPOS.ASCII> `
  --output-dir results\session
```

Optional API: `uvicorn wls_outlier_lab.app:app --port 8020` → `POST /analyze`.

### Outputs & plotting

The Python pipeline writes **data** only by default: `summary.json`,
`detector_comparison.csv`, `constellation_ablation.csv`, `per_epoch_<detector>.csv`.

**Plots are drawn in MATLAB** (cleaner, and matches the rest of this repo) by
`plot_wls_results.m`, which reads those CSVs:

```matlab
addpath('wls_outlier_lab');
plot_wls_results('results/2026_04_06')   % the --output-dir you used
```

It produces, into the same folder:
- `matlab_error_horizontal_vertical.png` — horizontal scatter, horizontal-vs-time,
  vertical-vs-time, and a horizontal/vertical error CDF for the best detector;
- `matlab_detector_comparison.png` — horizontal RMSE per detector (log scale);
- `matlab_constellation_ablation.png` — per-constellation drop/only RMSE.

(`--matplotlib-plots` still emits equivalent matplotlib PNGs if you have no MATLAB.)

## Validation result (2026-04-06 driving log, RTK truth, 216 epochs)

| detector | horizontal RMSE | vertical RMSE | avail. | note |
|---|---|---|---|---|
| `residual` / `combined` | **2.15 m** | 52 m | 100% | removes the blunder |
| `cn0` | 271.7 km | — | 100% | misses it |
| `elevation` | 271.8 km | — | 100% | misses it |
| `baseline` (naive) | 272.1 km | — | 100% | — |
| `residual_per_constellation` | 272.4 km | — | 100% | masks it |

Finding: a single **broken QZSS broadcast ephemeris** (PRN off by ~10⁶ m) throws
naive WLS ~270 km off. It carries **strong C/N0 (~53 dB-Hz) at high elevation
(~79°)**, so C/N0 and elevation gating do **not** catch it — only residual-based
detection does. Per-constellation scaling *masks* it (QZSS has too few satellites
for an internal robust scale); a **global** robust scale is required.
`constellation_ablation` independently confirms it: dropping QZSS recovers 2.18 m,
while GPS-only / GALILEO-only / BeiDou-only are each 2.5–5.3 m.

## Measurement calibration — data-driven sigma (not guessed)

With truth we measure the noise instead of assuming it (`core/calibration.py`):
each observation's prefit residual at truth, detrended by the per-epoch,
per-constellation median, is the per-signal error. Aggregated it yields the
empirical per-constellation sigma, and how sigma varies with C/N0 and elevation.
On the validation log:

| grouping | empirical clean sigma |
|---|---|
| GPS / GALILEO / BEIDOU | 21 / 18 / 16 m (≈ equal) |
| QZSS | ~1.4 × 10⁶ m (the broken-ephemeris satellite) |
| elevation 0–15° → 30–90° | 21 m → ~5 m |
| C/N0 <25 → >45 dB-Hz | 53 m → 10 m |

Takeaways: the real quality drivers are **elevation and C/N0** (the weight model
is validated), GPS/Galileo/BeiDou are ~equal here, and no GLONASS/SBAS/IRNSS are
even present — so the per-constellation sigma multipliers now **default to a
neutral 1.0** rather than baked-in guesses. Derive real ones per receiver and
apply them: `WeightConfig().with_constellation_scales(cal.sigma_scale_by_constellation)`.
The CLI writes `calibration_by_{constellation,cn0,elevation}.csv` by default
(`--no-calibrate` to skip); `plot_wls_results.m` draws them.

## Applying it to the `samsung_3rd` Novatel session

`samsung_3rd/21-sample_novatel_log/` supplies the **truth** (RTK BESTPOS,
`load_bestpos_track` verified: 11 k fixed samples over 2 h 16 m at Suwon). The
Samsung phone raw log for that same drive is not in the repo yet; once it is parsed
to a measurement TSV, run with `--truth-bestpos <...BESTPOS.ASCII>` and the analysis
runs unchanged.

## v1 limitations / next steps
- **No atmospheric model** (iono/tropo): the clean solution keeps a ~50 m vertical
  bias. Adding Klobuchar iono (coefficients are in the data) + Saastamoinen tropo
  will tighten vertical and sharpen residual detection. Outlier detection already
  works because the robust MAD scale adapts to the unmodelled spread.
- **Attitude** is an interface only (`sources/attitude.py`).
```

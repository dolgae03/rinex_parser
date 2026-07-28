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
- `clock_analysis` — receiver-clock stability: drift and instability from the
  per-epoch clock, an independent Doppler drift estimate, confounder-controlled
  correlation against the horizontal error, and the `clock_coasting_experiment`
  stress test. See the clock section below.

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

# Receiver clock study: instability analysis + the coasting stress test
python -m wls_outlier_lab.cli `
  --measurements <..._with_gt_ppp_aligned.tsv> --truth-columns `
  --output-dir results\clock --detectors combined --no-ablation `
  --clock-analysis --clock-coasting

# ... and the weak-geometry version (one system, minimum satellites)
python -m wls_outlier_lab.cli `
  --measurements <...tsv> --truth-columns --constellations GPS `
  --output-dir results\clock_gpsonly --detectors combined --no-ablation `
  --no-calibrate --clock-coasting --coast-sat-budget 5
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

### Blunder catalog & sky plot

The same truth-referenced residuals are turned into a **blunder catalog**
(`blunder_catalog.csv`: every observation with elevation, azimuth, C/N0,
residual, z-score, and an `is_outlier` flag), plus per-constellation and worst-
satellite tallies in `summary.json`. `plot_wls_results.m` draws a **sky plot**
(`matlab_skyplot.png`) — azimuth/elevation polar, coloured by |residual|, with
flagged outliers marked — so you can see *when and from which direction* bad
measurements arrive. On the validation log 6.8% of obs are flagged; QZSS is 67%
(satellites QZSS-3 and QZSS-4, the broken ephemerides), the rest are low-elevation
BeiDou near the horizon.

## Clock-drift stability vs horizontal error (`core/clock_analysis.py`)

Does the phone's clock-drift instability hurt the horizontal solution? A plain
correlation test is **not** enough to answer this, for two reasons, so the
pipeline attacks it from four sides.

### 1. Correlation, with the confounders removed (`--clock-analysis`)

From the per-epoch WLS clock bias we form drift `d(clk)/dt` and an instability
metric (departure from a constant-drift extrapolation, in metres), flag anomalies,
and correlate instability against the horizontal error vs truth — plainly, and
also partialling out HDOP and satellite count (which drive the error on their
own), at lag 1, and inside the worst-geometry quartile.

### 2. Is the "instability" even real? — noise floor + independent Doppler

**First trap:** the clock is *estimated* each epoch, so its series carries
estimation noise, and the instability metric is a second difference which
amplifies it. The analysis reports the formal clock sigma and the instability a
purely noisy clock would already show (`sqrt(6)·sigma_clk`). Second, it estimates
the clock drift **independently from Doppler** (`doppler_clock_drift`: range rate
`-(c/f)·doppler` against `u·(v_sv − v_rx) + drift`, one signal per satellite,
iterated MAD rejection) and compares the two.

### 3. The decisive test — clock coasting (`--clock-coasting`)

**Second trap, the important one:** a freely estimated clock *absorbs its own
instability by construction*, so `r ≈ 0` is the theoretically expected result and
proves almost nothing. So we make the solution *depend* on clock stability: the
clock is predicted as `clk_prev + drift_prev·dt` and imposed as a
pseudo-observation with standard deviation `sigma` (`solve_epoch(clock_prior=…)`,
which also makes an epoch solvable with one satellite fewer). Sweeping `sigma`
from ∞ (free) down to 0.3 m finds the point where trusting the clock starts to
cost accuracy — that point *is* the receiver's clock predictability.
`--coast-sat-budget N` thins the sky to N satellites, and `--constellations GPS`
restricts to one system, for the weak-geometry case where a clock constraint
actually carries weight.

The test is verified to be *sensitive* before being believed: unit tests confirm
that a tight but wrong clock prior corrupts the fix, and that an injected 300 m
clock jump makes the coasted solution >2× worse.

### Result on the 2026-04-06 log (1079 epochs at 1 Hz, PPP truth)

| measure | value |
|---|---|
| clock drift | −242.6 m/s (≈ −0.81 ppm), smooth ramp, **no ms jumps** |
| Doppler-derived drift | −242.46 m/s (agrees to **0.14 m/s**), fit residual **0.012 m/s** |
| real drift instability (Doppler) | **0.41 m/s per second** |
| same from position-domain clock | 1.63 m/s — i.e. **4× inflated by estimation noise** |
| clock estimate sigma | 3.2 m → a pure-noise clock would show 7.9 m of "instability", **4× more than observed** |
| corr(instability, horizontal error) | Pearson 0.02, Spearman 0.02, partial (HDOP+n_sats) 0.04, lag-1 −0.00, worst-HDOP quartile 0.04, Doppler-based −0.00 / −0.05 |
| horizontal error, stable vs unstable epochs | 2.80 vs 2.68 m |
| 18/1079 instability anomalies | horizontal error 0.5–5.4 m, i.e. unremarkable (overall p95 is 4.5 m) |

Coasting sweep (full sky, horizontal RMSE): free 3.00 m → σ=100/30/10 m all
3.00 m → σ=3 m 3.01 m → σ=1 m 3.23 m → **σ=0.3 m 8.10 m**.

The same sweep **restricted to GPS only with a 5-satellite budget** (weak geometry,
minimum redundancy — so the result cannot be dismissed as an artifact of having 36
satellites) has the identical shape, in median horizontal error: free 8.00 m →
σ=100 m 8.01 → σ=30 m 7.98 → σ=10 m 8.11 → σ=3 m 9.13 → σ=1 m 11.48 →
σ=0.3 m 32.65 m. Coasting is free down to σ ≈ 10 m; the error at 5 satellites is
set by geometry, not by the clock.

### Which clock anomalies are real?

The position-domain instability flags 18/1079 epochs — but they land
**where the car is moving** (67% of the flags, vs 28% of all epochs; median speed
6.4 vs 0.0 m/s) and the independent Doppler drift barely moves at those epochs
(0.79 vs 0.41 m/s median). So most of them are **estimation-noise false alarms
from vehicle dynamics/multipath, not oscillator events**. The analysis therefore
also flags anomalies on the Doppler drift (`is_dop_anomaly`, threshold from its
own MAD) and reports how many the two detectors agree on: **35 Doppler anomalies
(drift jerk > 2.18 m/s per second) vs 18 position-domain, only 7 in common.**
**Use the Doppler flag** — its fit residual is 0.013 m/s versus metre-level noise
in the position-domain clock.

Note the two results are consistent, not contradictory: there *are* ~35 genuine
clock-drift jerks, they are simply too small (≈ 2 m of clock motion) to matter
next to 10–13 m pseudorange noise, and the free clock re-estimate absorbs them.
So the answer to "can we find clock-drift anomalies?" is **yes, from Doppler** —
and separately, "do they wreck the horizontal fix?" is **no**.

### The event-based test — the strongest form of the answer

The Doppler drift curve shows a genuine **clock disturbance episode from ~740 s to
~1010 s**: the drift swings over −225 … −258 m/s and the jerk p95 reaches
**8.05 m/s per second, 21× the quiet median (0.38 m/s)**. Satellite count and HDOP
are identical inside and outside it (35.7 sats, HDOP 0.46 vs 0.48), so it is a
clean natural experiment. Horizontal error:

| | epochs | mean | median | p95 |
|---|---|---|---|---|
| quiet clock | 808 | 2.98 m | 3.04 m | 4.60 m |
| **clock disturbed** | 271 | **2.26 m** | **2.17 m** | **4.12 m** |

The error during the disturbance is not merely equal, it is **lower** (Welch
t = −9.6). The clock is obviously not *improving* the fix — that difference comes
from a confound (the car was driving in the open during that stretch rather than
parked near obstructions). The point is that **a real 21× clock-drift disturbance
produced zero horizontal degradation**, which is a much stronger statement than a
null correlation.

**Verdict: this phone's clock-drift instability does not degrade the horizontal
solution — and now we know why, quantitatively.** Its clock is unpredictable by
only ~0.4 m over one second, roughly **25× less than the pseudorange noise
(10–13 m)**, so it is irrelevant to the fix; and the per-epoch WLS re-estimates
the clock anyway, absorbing it. The coasting sweep locates the limit: the clock
can be coasted at σ ≈ 1 m at under 10% cost (σ = 3 m is free), and only a
sub-metre demand breaks it — and the two independent routes agree, since the
Doppler-measured 0.41 m/s of 1-second unpredictability is exactly the σ at which
the constraint starts to bite.

Scope of the claim, stated in WLS terms: "harmless" means **harmless while the
clock is a free per-epoch parameter**. The coasting sweep is the same engine with
the clock constrained instead, and it shows where that stops being true — tighten
the constraint past σ ≈ 1 m and this clock does hurt. Everything here, including
the σ = 0.3 m row, is single-epoch WLS; nothing in this package propagates state
between epochs.

Outputs: `clock_analysis.csv` (now also n_sats, HDOP, clock sigma, Doppler drift),
`clock_coasting.csv`, `clock_stability` + `clock_coasting` in `summary.json`;
`plot_wls_results.m` draws `matlab_clock_stability.png` (drift from both sources,
instability against its noise floor, scatter, quartiles) and
`matlab_clock_coasting.png`.

### Two data-quality findings this turned up

- **BeiDou `sv_vel` is wrong in ~18% of rows** of this processed TSV — off by
  hundreds of m/s, while GPS/Galileo/QZSS match `d(sv_pos)/dt` to 3 mm/s. It does
  not affect pseudorange positioning (which ignores velocity) but it wrecks any
  Doppler/velocity work, and it is why the Doppler estimator needs aggressive
  rejection. Worth fixing upstream.
- **The older `data/..._with_gt_rtk.tsv` truth for this same drive is offset
  +25.8 m in height** (std 0.11 m) relative to this PPP-aligned truth — an
  orthometric-vs-ellipsoidal datum mismatch (Korea's geoid undulation).
  Horizontal agrees to ~1 m. So earlier vertical RMSE figures (~45–52 m) were
  mostly this offset; with the PPP truth vertical RMSE is **19.8 m**. Horizontal
  conclusions are unaffected.

## Applying it to the `samsung_3rd` Novatel session

`samsung_3rd/21-sample_novatel_log/` supplies the **truth** (RTK BESTPOS,
`load_bestpos_track` verified: 11 k fixed samples over 2 h 16 m at Suwon). The
Samsung phone raw log for that same drive is not in the repo yet; once it is parsed
to a measurement TSV, run with `--truth-bestpos <...BESTPOS.ASCII>` and the analysis
runs unchanged.

## Atmospheric correction (`core/atmosphere.py`)

Delays are applied in a second WLS pass: solve coarsely, compute the per-satellite
iono + tropo slant delay at that fix, subtract, refine. Effect on the best
detector (validation log):

| config | horizontal RMSE | vertical RMSE |
|---|---|---|
| none | 2.10 m | 51.98 m |
| **tropo only (default)** | **2.25 m** | **45.29 m** |
| iono only (Klobuchar) | 4.73 m | 27.28 m |
| tropo + iono | 4.43 m | 20.77 m |

- **Troposphere (Saastamoinen)** is a clean win — vertical down, horizontal flat —
  so it is **on by default**.
- **Ionosphere (broadcast Klobuchar)** more than halves the vertical bias but
  **degrades horizontal (2.1 → 4.7 m)**, so it is **off by default** (`--iono` to
  enable). Two reasons: (1) single-frequency Klobuchar is a coarse model; (2) the
  `iono_b` (beta / period) coefficients in this processed TSV are the wrong scale
  (~1e-7 vs the expected ~1e5), so its diurnal term is degenerate — an upstream
  parsing issue worth fixing.
- **Dual-frequency ionosphere-free** (`core/iono_free.py`, `--iono-free`) forms the
  L1/L5 combination that removes ~all first-order ionosphere without any model.
  `--compare-iono` runs all four treatments with the same detector:

  | mode | horizontal RMSE | vertical RMSE | mean sats |
  |---|---|---|---|
  | none | 2.15 m | 51.9 m | 35.7 |
  | tropo (default) | 2.29 m | 45.2 m | 35.7 |
  | klobuchar | 4.46 m | 20.7 m | 35.7 |
  | iono_free | 9.81 m | 38.4 m | **15.9** |

  Finding: on this phone iono-free is **worse**, not better — L5 is sparse so it
  keeps only dual-band satellites (**satellites halved → weaker geometry**) and the
  combination amplifies noise ~2.6×. Ionosphere removal is real but outweighed here.
  So **troposphere-only is the default**; iono handling stays opt-in and, for this
  receiver/environment, best left off horizontally. `plot_wls_results.m` draws the
  comparison (`matlab_iono_comparison.png`: EN scatter, horizontal CDF, RMSE bars).

## v1 limitations / next steps
- **Ionosphere**: horizontal is best with tropo-only here; a denser-L5 dataset (or
  SBAS/PPP corrections) would let iono-free pay off. Both are wired in.
- **Attitude** is an interface only (`sources/attitude.py`).
- **Clock**: the conclusion (clock instability is harmless) is established for a
  freely estimated per-epoch clock. If anything downstream ever constrains the
  clock across epochs, the coasting sweep already gives the budget: keep the
  constraint looser than the measured ~0.4 m/s per second of drift
  unpredictability (σ ≳ 1 m at 1 Hz).
- **BeiDou `sv_vel`** is corrupt for ~18% of rows upstream; fix it in the parser
  and the Doppler estimator can drop its aggressive rejection.
```

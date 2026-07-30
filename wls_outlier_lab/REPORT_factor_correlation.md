# 삼성 스마트폰 GNSS: 에폭별 환경 요인 vs WLS 수평 항법해 오차 상관분석

**데이터**: `samsung_3rd/gnss_log_2026_04_06_14_12_10_..._with_gt_ppp_aligned.tsv`
(1079 에폭, PPP truth, GPS/GAL/BDS/QZSS, 1 Hz) ·
**엔진**: detector=combined + tropo-only WLS (수평 RMSE 2.996 m) ·
**분석 구현**: MATLAB `factor_correlation_analysis.m` ·
**결과물**: `results/samsung3rd_factors/` (CSV 3종 + `matlab_factor_*.png` 4장 + `factor_summary.txt`)

---

## 1. 질문과 파이프라인

특정 에폭의 환경 요소(DOP, C/N0, PR rate, 잔차, 클락 b·ḃ, 동역학 등)를 최대한 많이
뽑아, **속도 구간별로** WLS 항법해-truth 수평오차와의 상관을 재고, 이 폰에서 실제로
취약하게 작용하는 요인을 가려낸다. 2026-07-31 미팅 PPT의 3차년도 제안(측정치 이상
파라미터: **C/N0, Doppler-Code Difference, Code-Carrier Divergence** + 동적 상황
오차 모델링)과 문헌 슬라이드의 피처 목록(고도각/방위각, residual, innovation,
PR consistency, carrier smoothness 등)을 피처 설계에 반영했다.

```
Python (검증된 WLS 엔진 재사용)                MATLAB (분석 본체 구현)
─────────────────────────────────            ─────────────────────────────────
cli --factor-export                          factor_correlation_analysis.m
  → 에폭별 44개 피처 추출                       → 속도 구간 분할 (4구간)
  → factor_epochs.csv (1079행)                → Spearman/Pearson/부분상관
                                             → moving-block bootstrap 95% CI
                                             → 랭킹, CSV/TXT, PNG 4장
```

## 2. 피처 (에폭당 44개, 위성별 값은 사용-위성 집계)

| 계열 | 피처 | 비고 |
|---|---|---|
| 기하 | `n_used/available/rejected`, `hdop/vdop/pdop/gdop/tdop`, `elev_min/mean`, `n_low_elev`, `n_gps/gal/bds/qzs` | |
| 신호 | `cn0_mean/min/std(used)`, `cn0_mean_all`, `cn0_frac_below30` | |
| 정합성 | `sigma0_hat`, `resid_rms/max`, `clk_sigma_m` | 사후 잔차 |
| 클락 | `clk_m`(b), `clk_drift_mps`(ḃ), `clk_inst_m`, `drift_dop_mps`, `drift_dop_change`, `dop_resid_rms` | 기존 클락 분석 지표 포함 |
| 코드품질 | `dcd_med/max/rms` (**Doppler-Code Difference rate**: dPR/dt − Doppler, 수신기 클락 소거됨), `ccd_med/max/rms` (**Code-Carrier Divergence rate**: d(PR−λφ)/dt), `frac_loi` | PPT 3차년도 파라미터. `phase_cycle` 컬럼을 이번에 로더에 추가 |
| 동역학 | `speed_mps`·`accel_mps2` (truth 중앙차분), `speed_dop_mps` (Doppler LS) | |
| 파생(MATLAB) | `clk_detrend_m`, `clk_drift_dev_mps`, `drift_dop_dev_mps`, `abs_accel_mps2` | b는 단조 드리프트라 원시값은 시간추세 프록시 |

## 3. 통계 방법

- **Spearman을 주지표**로(단조·이상치 강건), Pearson 병기.
- **부분상관 2종**: 기하 통제(HDOP+n_used), 기하+속도 통제 — "동역학 프록시인가"를 판별.
- 오차 시계열은 자기상관이 커서 소박한 p-value는 과신 → **moving-block bootstrap**
  (블록 ≤30 s, B=1000) 95% CI로 유의성 판정(`*` = CI가 0을 제외).
- 속도 구간: truth 수평속도로 **정지(<0.5) / 저속(0.5–3) / 중속(3–8) / 고속(>8 m/s)**.

## 4. 결과

### 구간별 수평오차 — 정지가 가장 나쁘다

| 구간 | n | 중앙값 [m] | p95 [m] | RMSE [m] |
|---|---|---|---|---|
| 전체 | 1079 | 2.875 | 4.517 | 2.996 |
| **정지(<0.5)** | **754** | **2.989** | **4.575** | **3.096** |
| 저속(0.5–3) | 35 | 2.253 | 3.983 | 2.566 |
| 중속(3–8) | 78 | 2.467 | 4.515 | 2.838 |
| 고속(>8) | 212 | 2.498 | 4.202 | 2.745 |

세션 구조: 0–640 s 정지 → 이후 두 번의 주행. **정지 구간이 이동 구간보다 오차
중앙값이 ~0.5 m 나쁘다** — 정지 시 멀티패스가 지속 바이어스로 고이고 이동 시
평균화된다는, 스마트폰에서 알려진 패턴과 일치. `speed_mps` r=−0.16\*,
`abs_accel` r=−0.17\* (동역학이 클수록 오차 **감소**).

### 전체 랭킹 상위 (Spearman, `*`=bootstrap 유의, p.g+v=기하+속도 통제 부분상관)

| 순위 | 피처 | r | p.g+v | 판정 |
|---|---|---|---|---|
| 1 | `cn0_std_used` | **+0.374\*** | +0.345 | **사용 위성 간 C/N0 편차 — 최강 요인. 통제 후에도 유지** |
| 2 | `cn0_mean_all` | +0.321\* | +0.265 | 평균 C/N0가 높을수록 오차↑ (역설, §5) |
| 3 | `clk_m` | +0.317\* | +0.262 | 단조 드리프트 → 시간추세 프록시(인과 아님) |
| 4–7 | `resid_rms/max`, `sigma0_hat`, `cn0_mean_used` | +0.25~0.27\* | +0.18~0.24 | 사후 잔차계열 — 온라인 품질지표로 유효 |
| 8 | `elev_min_deg` | −0.238\* | −0.211 | 최저 고도각 위성이 낮을수록 오차↑ |
| 9–10 | `clk_drift_mps`, `drift_dop_mps` | −0.22\* | −0.17 | 느린 공변(온도 등 추정), 아래 클락 판정 참조 |
| 11–12 | `cn0_min_used`, `cn0_frac_below30` | ±0.20\* | ±0.16 | 약한 위성의 존재 자체가 해로움 |
| 13 | `n_bds` | +0.187\* | +0.116 | **BeiDou 사용 수 ↑ → 오차 ↑** (blunder 카탈로그와 정합) |
| — | `dcd/ccd_*` | −0.06~−0.12 | ≈−0.05 | 에폭 집계로는 예측력 낮음 (§5) |
| — | `clk_inst_m` | +0.017 | +0.041 | **무상관 — 기존 클락 결론 재확인** |
| — | `n_used`, `hdop` 등 기하 | ≈0.00~0.10 | — | 하늘이 풍부해(평균 36위성, HDOP 0.47) 기하는 병목 아님 |

### 속도 구간별 상위 요인

- **정지(n=754)**: `cn0_std_used` +0.28, `elev_min` −0.26, `cn0_min_used` −0.26,
  `resid_rms` +0.25\* — 신호품질 요인이 지배.
- **중속(n=78)**: `n_gal` +0.43\*, `sigma0_hat` +0.40\*, `cn0_mean_all` +0.37\* —
  n이 작고 특정 경로 구간과 얽혀 있어 참고 수준.
- **고속(n=212)**: `clk_detrend_m` +0.49\*, **`cn0_std_used` +0.475\***,
  `cn0_mean_used` +0.42 — 주행 중에도 C/N0 편차가 1위군.
  (`clk_detrend`는 시간구조 잔존 가능성이 높음 — 진짜 클락 이상 지표인
  `clk_inst_m`/`drift_dop_change`가 전 구간 ≈0이므로 인과로 보지 않음.)

### 그림

`matlab_factor_heatmap.png`(피처×구간 r 히트맵), `matlab_factor_rank.png`(부호별
랭킹+CI, 구간별 오차), `matlab_factor_top_scatter.png`(상위 6피처 산점도+십분위
중앙값 추세), `matlab_factor_timeseries.png`(오차·속도·상위 2피처 시계열).

## 5. 해석 — 이 폰의 취약점은 무엇인가

1. **비대칭 신호환경(멀티패스)이 1순위 취약점.** 위성 수·DOP가 아니라 **사용 위성 간
   C/N0 편차**가 전 구간에서 가장 강하고, 기하·속도를 통제해도 살아남는다(+0.345).
   일부 위성만 약해지는 부분 차폐/반사 상황이 가중치로 다 흡수되지 못하고 해에 남는다.
2. **정지 상태가 오히려 위험하다.** 정지 시 멀티패스 기하가 고정되어 바이어스가
   지속되고(중앙값 2.99 m), 주행은 이를 평균화한다(2.50 m). "동적 상황 오차 모델링"
   과제에서 **정지-멀티패스 체류가 별도 모드로 다뤄져야 함**을 시사.
3. **`cn0_mean`의 양(+) 상관 역설**: 반사파가 보강간섭으로 C/N0를 올리면서 코드
   바이어스를 유발하는 정지 지점 특성으로 해석. **"C/N0 높음=측정치 좋음" 가정은 이
   데이터에서 성립하지 않음** — soft weighting 학습 피처로 평균보다 *분산*이 낫다.
4. **BeiDou가 통계적으로 해롭다.** `n_bds` +0.19\* — 기존 truth-referenced blunder
   카탈로그(BDS 최다 flag)·`sv_vel` 손상(18.5% 행)과 3중으로 정합. BDS 가중 하향
   또는 선별 사용 실험 가치 있음.
5. **잔차계열(σ₀, resid RMS/max)은 truth 없이 계산 가능한 온라인 품질지표**로서
   오차와 유의하게 동행 — integrity/가중 조정의 입력 후보.
6. **클락은 취약점이 아니다.** b·ḃ의 느린 공변은 시간추세와 얽혀 있으나, 진짜 이상
   지표(`clk_inst_m`, `drift_dop_change_mps`)는 전 구간 무상관 — 4중 검증(코스팅
   스윕 포함)으로 확정한 기존 결론과 일치.
7. **DCD/CCD는 에폭 집계로는 약하다**(통제 후 ≈−0.05). 원인: 에폭 중앙값/최대로
   뭉개면 위성 단위 정보가 소실. PPT의 취지대로 **위성(측정치) 단위 이상 검출
   피처**로 쓰는 것이 맞고, 그 평가는 후속(blunder 카탈로그와 측정치 단위 결합).

## 6. 질문 대신 자체 판단한 사항

| 결정 | 근거 |
|---|---|
| 분석·플롯은 MATLAB, WLS·피처 추출은 기존 Python 엔진 | "MATLAB으로 구현" 지시 이행 + 2주간 검증된 탐지기/ISB/robust 로직 재작성 위험 회피. 통계(순위상관·부분상관·bootstrap)와 그림은 전부 MATLAB 구현 |
| 속도 경계 0.5/3/8 m/s | 세션 속도 분포(정지 70%, 주행 최고 21 m/s)에서 해석 가능한 4구간; 저속 구간(n=35)은 표본 부족으로 상관 억제(n<60) |
| 오차 목표 = 수평오차만 | 기존 컨벤션(수직 truth 신뢰도 문제 포함) |
| b, ḃ 원시값 포함 + 파생(detrend/dev) 추가 | 요청대로 포함하되, 단조 드리프트가 시간추세 프록시가 되는 문제를 파생 피처와 캐비앗으로 처리 |
| Spearman 주지표 + block bootstrap | 오차 분포 꼬리·자기상관 대응. 통계 함수는 toolbox 없이 로컬 구현 |
| DCD 부호/기준 | Doppler sign은 데이터에서 자동 결정(−1), 사다리꼴 평균으로 가속 성분 상쇄, 위성시계는 corrected PR 사용으로 소거 |

## 7. 한계

- 단일 세션·단일 폰: 정지 구간이 한 지점이라 "정지 취약"과 "그 지점의 멀티패스"가
  완전히는 분리 안 됨(구간 내 상관으로 부분 확인). 세션 추가 시 같은 파이프라인 재사용 가능.
- 중속/저속 구간 표본 부족(78/35 에폭).
- 상관≠인과: 개입 실험(예: BDS 제외 재해석, C/N0-분산 기반 가중)을 후속으로 제안.

**재현**: `python -m wls_outlier_lab.cli --measurements <tsv> --truth-columns
--output-dir results/samsung3rd_factors --factor-export --clock-analysis
--detectors combined --no-ablation --no-calibrate` →
MATLAB `factor_correlation_analysis('results/samsung3rd_factors')`

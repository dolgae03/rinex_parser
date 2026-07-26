"""Single-epoch weighted-least-squares pseudorange positioning.

Pure business logic: given a list of ``SatObs`` it returns a ``WlsSolution``.
No file I/O, no plotting, no outlier policy (detectors live in ``detectors.py``
and only ever *select* which satellites are handed to ``solve_epoch``).

Model
-----
For each observation ``i`` the clock-corrected observation is

    L_i = pseudorange_m + sv_clock_bias_m + pr_correction_m - iono_delay_m

and the predicted observation is

    P_i = |sv_pos_i - rx| + clk + isb[constellation_i]

``clk`` is the receiver clock bias (metres, common to all signals) and ``isb``
is a per-constellation inter-system bias relative to the reference constellation
(GPS when present). Unknowns are solved by Gauss-Newton with a diagonal weight
matrix from elevation and C/N0, and a backtracking line search so a gross
blunder cannot make the iteration diverge (it just biases the estimate — that is
what the detectors are for). The per-iteration design build is vectorised.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

from ..types import SatObs, WlsSolution
from . import frames

Key = Tuple[int, int, int]      # signal identity (constellation, prn, freq_hz)
SatId = Tuple[int, int]         # satellite identity (constellation, prn) — rejection unit


@dataclass
class WeightConfig:
    """Diagonal observation-weight model. Variance = base^2 * f_elev * f_cn0.

    Only *relative* weights affect the position estimate; the absolute scale
    sets the a-priori sigma used for the residual-snooping test statistic.
    """

    base_sigma_m: float = 3.0
    use_elevation: bool = True
    use_cn0: bool = True
    cn0_ref_dbhz: float = 40.0
    min_elevation_deg: float = 5.0  # floor to keep 1/sin(el) finite for low sats
    # Per-constellation base-sigma multipliers (relative to base_sigma_m). Default
    # is NEUTRAL (1.0 for every constellation): we do not bake in unmeasured
    # guesses. Derive real values per receiver from truth with
    # ``core.calibration`` and set them here (see WeightConfig.with_constellation_scales).
    use_constellation: bool = True
    sigma_scale_by_constellation: Dict[int, float] = field(default_factory=dict)
    default_constellation_scale: float = 1.0

    def with_constellation_scales(self, scales_by_name_or_id: Dict) -> "WeightConfig":
        """Return a copy with per-constellation sigma scales (accepts names or ids)."""
        from ..types import CONSTELLATION_NAMES
        name_to_id = {v: k for k, v in CONSTELLATION_NAMES.items()}
        resolved = {}
        for key, val in scales_by_name_or_id.items():
            cid = name_to_id.get(key, key)
            resolved[int(cid)] = float(val)
        import dataclasses
        return dataclasses.replace(self, sigma_scale_by_constellation=resolved)

    def scale_for(self, constellation: int) -> float:
        if not self.use_constellation:
            return 1.0
        return self.sigma_scale_by_constellation.get(int(constellation),
                                                     self.default_constellation_scale)

    def sigma_m(self, elevation_deg: float, cn0_dbhz: float, constellation: int = 0) -> float:
        var = (self.base_sigma_m * self.scale_for(constellation)) ** 2
        if self.use_elevation:
            s = math.sin(math.radians(max(elevation_deg, self.min_elevation_deg)))
            var *= 1.0 / (s * s)
        if self.use_cn0 and cn0_dbhz == cn0_dbhz and cn0_dbhz > 0:
            var *= 10.0 ** (-(cn0_dbhz - self.cn0_ref_dbhz) / 10.0)
        return math.sqrt(var)

    def variance_array(self, elevation_deg: np.ndarray, cn0_dbhz: np.ndarray,
                       cons_scale: Optional[np.ndarray] = None) -> np.ndarray:
        base = self.base_sigma_m ** 2
        var = np.full(elevation_deg.shape, base, dtype=float)
        if cons_scale is not None:
            var = var * (cons_scale * cons_scale)
        if self.use_elevation:
            el = np.radians(np.maximum(elevation_deg, self.min_elevation_deg))
            s = np.sin(el)
            var = var / (s * s)
        if self.use_cn0:
            cn0 = np.where(np.isfinite(cn0_dbhz) & (cn0_dbhz > 0), cn0_dbhz, self.cn0_ref_dbhz)
            var = var * 10.0 ** (-(cn0 - self.cn0_ref_dbhz) / 10.0)
        return var


@dataclass
class WlsConfig:
    weights: WeightConfig = None
    apply_sagnac: bool = True
    apply_iono: bool = True  # subtract SatObs.iono_delay_m when present
    max_iterations: int = 15
    convergence_m: float = 1e-2
    reference_constellation: Optional[int] = 0  # GPS; falls back to most-populous
    # Robust (Huber IRLS) mode — used by detectors to get a blunder-resistant
    # position from which outliers stand out. The naive 'baseline' leaves it off.
    robust_huber: bool = False
    huber_c: float = 1.5
    huber_floor_m: float = 5.0

    def __post_init__(self) -> None:
        if self.weights is None:
            self.weights = WeightConfig()


def _pick_signals(obs: Sequence[SatObs], reject: Set[SatId]) -> List[SatObs]:
    """Validity filter + one signal per satellite (highest C/N0).

    Rows with non-finite satellite position / pseudorange are dropped, as are
    satellites in ``reject``. When a satellite has several frequencies we keep
    the strongest so the design matrix is not inflated by correlated signals.
    Rejection is per-satellite (constellation, prn) because a bad broadcast
    ephemeris corrupts every frequency of that satellite.
    """
    best: Dict[SatId, SatObs] = {}
    for o in obs:
        if (o.constellation, o.prn) in reject:
            continue
        sp = o.sv_pos
        if not all(math.isfinite(v) for v in sp):
            continue
        if not math.isfinite(o.pseudorange_m) or not math.isfinite(o.sv_clock_bias_m):
            continue
        if abs(o.pseudorange_m) > 1e8 or o.pseudorange_m <= 0:
            # near-raw pseudoranges are ~2e7 m; guard against decode garbage
            continue
        sat = (o.constellation, o.prn)
        cur = best.get(sat)
        cn0 = o.cn0_dbhz if math.isfinite(o.cn0_dbhz) else -1.0
        if cur is None or cn0 > (cur.cn0_dbhz if math.isfinite(cur.cn0_dbhz) else -1.0):
            best[sat] = o
    return list(best.values())


def solve_epoch(
    obs: Sequence[SatObs],
    t_sec: float,
    config: Optional[WlsConfig] = None,
    reject: Optional[Iterable[SatId]] = None,
    x0: Optional[np.ndarray] = None,
) -> WlsSolution:
    """Solve one epoch. ``reject`` names satellites (constellation, prn) to drop."""
    cfg = config or WlsConfig()
    reject_set: Set[SatId] = set(reject or ())
    signals = _pick_signals(obs, reject_set)
    n_avail = len({(o.constellation, o.prn) for o in obs})

    sol = WlsSolution(t_sec=t_sec, ecef=(float("nan"),) * 3, clock_bias_m=float("nan"),
                      n_available=n_avail)
    if not signals:
        sol.reason = "no_valid_signals"
        return sol

    cons_present = sorted({o.constellation for o in signals})
    ref = cfg.reference_constellation
    if ref not in cons_present:
        counts: Dict[int, int] = {}
        for o in signals:
            counts[o.constellation] = counts.get(o.constellation, 0) + 1
        ref = max(counts, key=counts.get)
    isb_cons = [c for c in cons_present if c != ref]
    isb_index = {c: 4 + i for i, c in enumerate(isb_cons)}
    n_unknown = 4 + len(isb_cons)
    if len(signals) < n_unknown:
        sol.reason = f"underdetermined ({len(signals)} obs < {n_unknown} unknowns)"
        sol.n_used = len(signals)
        return sol

    # Static per-observation arrays (constant across iterations).
    n = len(signals)
    sv0 = np.array([o.sv_pos for o in signals], dtype=float)          # (n,3)
    L_obs = np.array([o.corrected_pseudorange_m for o in signals], dtype=float)
    if cfg.apply_iono:
        iono = np.array([o.iono_delay_m if math.isfinite(o.iono_delay_m) else 0.0
                         for o in signals], dtype=float)
        L_obs = L_obs - iono
    cn0 = np.array([o.cn0_dbhz for o in signals], dtype=float)
    cons_scale = np.array([cfg.weights.scale_for(o.constellation) for o in signals], dtype=float)
    isb_col = np.array([isb_index.get(o.constellation, -1) for o in signals], dtype=int)
    keys = [o.key for o in signals]

    def _build(state_vec):
        """Vectorised design matrix G, residual r, weight vector w, elevations."""
        rx = state_vec[:3]
        clk = state_vec[3]
        lat, lon, _ = frames.ecef_to_lla(*rx)
        R = frames.enu_rotation_matrix(lat, lon)

        sv = sv0
        if cfg.apply_sagnac:
            travel = np.linalg.norm(sv0 - rx, axis=1) / frames.C_LIGHT
            theta = frames.OMEGA_EARTH * travel
            ct, st = np.cos(theta), np.sin(theta)
            sv = np.column_stack((ct * sv0[:, 0] + st * sv0[:, 1],
                                  -st * sv0[:, 0] + ct * sv0[:, 1],
                                  sv0[:, 2]))
        diff = sv - rx
        rng = np.linalg.norm(diff, axis=1)
        unit = diff / rng[:, None]

        enu = diff @ R.T
        elev = np.degrees(np.arctan2(enu[:, 2], np.hypot(enu[:, 0], enu[:, 1])))

        isb_per_obs = np.where(isb_col >= 0, state_vec[np.clip(isb_col, 0, None)], 0.0)
        predicted = rng + clk + isb_per_obs
        r = L_obs - predicted

        G = np.zeros((n, n_unknown))
        G[:, :3] = -unit
        G[:, 3] = 1.0
        has_isb = isb_col >= 0
        G[np.arange(n)[has_isb], isb_col[has_isb]] = 1.0

        w = 1.0 / cfg.weights.variance_array(elev, cn0, cons_scale)
        if cfg.robust_huber and len(r) > 4:
            # IRLS: down-weight observations far from the robust residual centre
            s = max(1.4826 * float(np.median(np.abs(r - np.median(r)))), cfg.huber_floor_m)
            a = np.abs(r - np.median(r)) / s
            w = w * np.where(a <= cfg.huber_c, 1.0, cfg.huber_c / np.maximum(a, 1e-9))
        return G, r, w, elev

    def _wssr(r, w):
        return float(np.sum(w * r * r))

    if x0 is not None:
        state = np.array([x0[0], x0[1], x0[2], 0.0] + [0.0] * len(isb_cons), dtype=float)
    else:
        mean_sv = sv0.mean(axis=0)
        r0 = np.linalg.norm(mean_sv)
        guess = mean_sv / r0 * frames.WGS84_A if r0 > 0 else np.array([frames.WGS84_A, 0.0, 0.0])
        state = np.array([guess[0], guess[1], guess[2], 0.0] + [0.0] * len(isb_cons), dtype=float)

    # Gauss-Newton with backtracking line search.
    iterations = 0
    converged = False
    for iterations in range(1, cfg.max_iterations + 1):
        G, resid, w, _ = _build(state)
        try:
            N = G.T @ (G * w[:, None])
            dstate = np.linalg.solve(N, G.T @ (w * resid))
        except np.linalg.LinAlgError:
            sol.reason = "singular_normal_matrix"
            return sol
        ssr0 = _wssr(resid, w)
        alpha, accepted = 1.0, False
        for _ in range(8):
            cand = state + alpha * dstate
            _, rc, wc, _ = _build(cand)
            if _wssr(rc, wc) <= ssr0:
                state, accepted = cand, True
                break
            alpha *= 0.5
        if not accepted:
            state = state + dstate
        if np.linalg.norm(alpha * dstate[:3]) < cfg.convergence_m:
            converged = True
            break

    G, resid, w, elev = _build(state)
    W = np.diag(w)
    rx = state[:3]

    v = resid - G @ np.linalg.solve(G.T @ W @ G, G.T @ W @ resid)
    dof = n - n_unknown
    Ninv = np.linalg.inv(G.T @ W @ G)
    Qvv = np.diag(1.0 / w) - G @ Ninv @ G.T
    sigma0_hat = math.sqrt(max(float(v.T @ W @ v) / dof, 0.0)) if dof > 0 else float("nan")

    sol.ecef = (float(rx[0]), float(rx[1]), float(rx[2]))
    sol.clock_bias_m = float(state[3])
    sol.isb_m = {c: float(state[isb_index[c]]) for c in isb_cons}
    sol.n_used = n
    sol.used_keys = keys
    sol.residuals_m = {k: float(vi) for k, vi in zip(keys, v)}
    sol.elevation_deg = {k: float(e) for k, e in zip(keys, elev)}
    sol.iterations = iterations
    sol.converged = converged
    sol.sigma0_hat = sigma0_hat
    # Baarda w-test with a-priori unit variance: w_i = |v_i| / sqrt(Qvv_ii).
    # Using the a-priori sigma (not the contaminated sigma0_hat) stops a single
    # gross blunder from masking itself.
    std_res: Dict[Key, float] = {}
    for i, k in enumerate(keys):
        q = Qvv[i, i]
        std_res[k] = abs(v[i]) / math.sqrt(q) if q > 1e-12 else float("nan")
    sol.std_residuals = std_res

    _fill_dop(sol, sv0, rx)
    return sol


def _fill_dop(sol: WlsSolution, sv0: np.ndarray, rx: np.ndarray) -> None:
    """Dilution of precision from the position+clock geometry (unweighted)."""
    diff = sv0 - rx
    rng = np.linalg.norm(diff, axis=1)
    good = rng > 1.0
    if good.sum() < 4:
        return
    unit = diff[good] / rng[good][:, None]
    A = np.column_stack((unit, np.ones(good.sum())))
    try:
        Q = np.linalg.inv(A.T @ A)
    except np.linalg.LinAlgError:
        return
    lat, lon, _ = frames.ecef_to_lla(*rx)
    R = frames.enu_rotation_matrix(lat, lon)
    Qenu = R @ Q[:3, :3] @ R.T
    sol.gdop = float(math.sqrt(max(np.trace(Q), 0.0)))
    sol.pdop = float(math.sqrt(max(np.trace(Q[:3, :3]), 0.0)))
    sol.hdop = float(math.sqrt(max(Qenu[0, 0] + Qenu[1, 1], 0.0)))
    sol.vdop = float(math.sqrt(max(Qenu[2, 2], 0.0)))
    sol.tdop = float(math.sqrt(max(Q[3, 3], 0.0)))

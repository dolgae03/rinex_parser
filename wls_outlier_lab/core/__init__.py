"""Business logic: WLS, outlier detectors, error metrics, experiments. No I/O."""

from .detectors import DETECTORS, DetectorConfig
from .experiment import (
    DetectorRun,
    ExperimentResult,
    compare_ionosphere,
    constellation_ablation,
    run_experiment,
)
from .iono_free import form_iono_free, form_iono_free_epoch
from .metrics import EpochError, ErrorMetrics, compute_metrics
from .wls import WeightConfig, WlsConfig, solve_epoch

__all__ = [
    "solve_epoch",
    "WlsConfig",
    "WeightConfig",
    "DETECTORS",
    "DetectorConfig",
    "compute_metrics",
    "ErrorMetrics",
    "EpochError",
    "run_experiment",
    "constellation_ablation",
    "compare_ionosphere",
    "form_iono_free",
    "form_iono_free_epoch",
    "ExperimentResult",
    "DetectorRun",
]

"""Business logic: WLS, outlier detectors, error metrics, experiments. No I/O."""

from .detectors import DETECTORS, DetectorConfig
from .experiment import (
    DetectorRun,
    ExperimentResult,
    constellation_ablation,
    run_experiment,
)
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
    "ExperimentResult",
    "DetectorRun",
]

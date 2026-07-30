"""wls_outlier_lab — define GNSS measurement outliers and quantify the WLS gain.

Pipeline for the study: *which measurements (per constellation or overall) should
be flagged as outliers, and how much does removing them improve a weighted
least-squares fix versus a Novatel BESTPOS / RTK truth?*

Layering (data acquisition is fully separated from business logic):

    sources/   data-access adapters  -> Epoch / TruthTrack           (I/O)
    core/      WLS, detectors, metrics, experiment                   (pure)
    reporting  write CSV/JSON/plots
    cli / app  thin entry points (CLI + optional FastAPI)

A new dataset is a new source class; the engine never changes.
"""

from .types import Epoch, SatObs, TruthSample, TruthTrack, WlsSolution, constellation_name

__all__ = [
    "SatObs",
    "Epoch",
    "TruthSample",
    "TruthTrack",
    "WlsSolution",
    "constellation_name",
]

__version__ = "0.1.0"

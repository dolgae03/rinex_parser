"""Abstract data-access interfaces.

Everything that *reads bytes from somewhere* implements one of these. The core
(WLS, detectors, metrics, experiment) only ever sees ``Epoch`` / ``TruthTrack``
objects, so a new dataset = a new source class, with zero changes to the logic.
This is the separation the project asked for: many data programs, one engine.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List, Optional

from ..types import Epoch, TruthTrack


class MeasurementSource(ABC):
    """Yields per-epoch satellite observations for the receiver under test."""

    @abstractmethod
    def load_epochs(self) -> List[Epoch]:
        ...


class TruthSource(ABC):
    """Yields a reference trajectory (e.g. Novatel BESTPOS RTK) in ECEF."""

    @abstractmethod
    def load_track(self) -> TruthTrack:
        ...


class AttitudeSource(ABC):
    """Optional attitude covariate (phone IMU orientation or INS truth).

    v1 ships the interface only; concrete parsers are a later extension so the
    outlier/WLS study can proceed without them.
    """

    @abstractmethod
    def load_series(self) -> "Optional[object]":
        ...

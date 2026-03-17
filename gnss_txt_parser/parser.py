from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from gnss_txt_parser.data_class import Measurement, MeasurementValue, SatelliteInfo, TimeTag


C = 299_792_458.0
GPS_UNIX_OFFSET_MILLIS = 315_964_800_000
WEEK_NANOS = 604_800 * 1_000_000_000

STATE_CODE_LOCK = 1
STATE_TOW_DECODED = 8
STATE_MSEC_AMBIGUOUS = 16
STATE_TOW_KNOWN = 16_384

CONSTELLATION_ANDROID = {
    0: "unknown",
    1: "gps",
    2: "sbas",
    3: "glonass",
    4: "qzss",
    5: "beidou",
    6: "galileo",
    7: "irnss",
}

CODE_TYPE_ANDROID = {
    "gps": {"C": "l1", "Q": "l5", "X": "l5"},
    "glonass": {"C": "l1"},
    "qzss": {"C": "l1", "X": "l5"},
    "galileo": {"C": "e1", "Q": "e5a", "X": "e5a"},
    "beidou": {"I": "b1", "X": "l5"},
}

QZSS_PRN_SVN = {
    193: 1,
    194: 2,
    199: 3,
    195: 4,
    196: 5,
}

RAW_NUMERIC_COLUMNS: dict[str, type] = {
    "utcTimeMillis": int,
    "TimeNanos": int,
    "LeapSecond": int,
    "TimeUncertaintyNanos": float,
    "FullBiasNanos": int,
    "BiasNanos": float,
    "BiasUncertaintyNanos": float,
    "DriftNanosPerSecond": float,
    "DriftUncertaintyNanosPerSecond": float,
    "HardwareClockDiscontinuityCount": int,
    "Svid": int,
    "TimeOffsetNanos": float,
    "State": int,
    "ReceivedSvTimeNanos": int,
    "ReceivedSvTimeUncertaintyNanos": float,
    "Cn0DbHz": float,
    "PseudorangeRateMetersPerSecond": float,
    "PseudorangeRateUncertaintyMetersPerSecond": float,
    "AccumulatedDeltaRangeState": int,
    "AccumulatedDeltaRangeMeters": float,
    "AccumulatedDeltaRangeUncertaintyMeters": float,
    "CarrierFrequencyHz": float,
    "CarrierCycles": float,
    "CarrierPhase": float,
    "CarrierPhaseUncertainty": float,
    "MultipathIndicator": float,
    "SnrInDb": float,
    "ConstellationType": int,
    "AgcDb": float,
    "BasebandCn0DbHz": float,
    "FullInterSignalBiasNanos": float,
    "FullInterSignalBiasUncertaintyNanos": float,
    "SatelliteInterSignalBiasNanos": float,
    "SatelliteInterSignalBiasUncertaintyNanos": float,
    "ChipsetElapsedRealtimeNanos": int,
    "IsFullTracking": int,
}


@dataclass
class ParseSummary:
    input_path: Path
    total_rows: int
    valid_rows: int
    output_path: Path | None = None


def _parse_value(value: str, caster: type) -> int | float | str | None:
    text = value.strip()
    if text == "":
        return None
    try:
        return caster(text)
    except ValueError:
        return None


def _safe_int(value: object, default: int = 0) -> int:
    if value is None:
        return default
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if math.isnan(value):
            return default
        return int(value)
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _safe_float(value: object, default: float = math.nan) -> float:
    if value is None:
        return default
    if isinstance(value, (int, float)):
        return float(value)
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _truncate_frequency_for_output(value: object) -> float:
    freq = _safe_float(value)
    if math.isnan(freq):
        return math.nan
    integer_value = int(freq)
    text = str(integer_value)
    if len(text) <= 6:
        return float(integer_value)
    return float(text[:6] + ("0" * (len(text) - 6)))


def _load_raw_rows(input_path: Path) -> list[dict[str, object]]:
    header: list[str] | None = None
    rows: list[dict[str, object]] = []

    with input_path.open("r", encoding="utf8", newline="") as infile:
        reader = csv.reader(infile)
        for raw_row in reader:
            if not raw_row:
                continue

            first = raw_row[0].strip()
            if first == "# Raw":
                header = ["record_type"] + [item.strip() for item in raw_row[1:]]
                continue
            if first != "Raw":
                continue
            if header is None:
                raise ValueError(f"Found 'Raw' row before '# Raw' header in {input_path}")

            row: dict[str, object] = {}
            padded = raw_row + [""] * (len(header) - len(raw_row))
            for column, value in zip(header, padded):
                if column == "record_type":
                    row[column] = value
                elif column in RAW_NUMERIC_COLUMNS:
                    row[column] = _parse_value(value, RAW_NUMERIC_COLUMNS[column])
                else:
                    row[column] = value.strip()
            rows.append(row)

    if header is None:
        raise ValueError(f"Could not find '# Raw' header in {input_path}")
    if not rows:
        raise ValueError(f"Could not find any 'Raw' rows in {input_path}")
    return rows


def _normalize_qzss_prn(gnss_id: str, svid: object) -> int:
    sv_id = _safe_int(svid, default=-1)
    if gnss_id == "qzss":
        return QZSS_PRN_SVN.get(sv_id, sv_id)
    return sv_id


def _compute_gps_millis(row: dict[str, object]) -> float:
    utc_millis = _safe_float(row.get("utcTimeMillis"))
    leap_seconds = _safe_float(row.get("LeapSecond"), default=18.0)
    if math.isnan(utc_millis):
        return math.nan
    return utc_millis - GPS_UNIX_OFFSET_MILLIS + leap_seconds * 1000.0


def _compute_measurement_rows(raw_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    segment_base_full_bias: int | None = None
    segment_base_bias_nanos = 0.0
    previous_full_bias: int | None = None

    for raw in raw_rows:
        current_full_bias = _safe_int(raw.get("FullBiasNanos"))
        current_bias_nanos = _safe_float(raw.get("BiasNanos"), default=0.0)

        if previous_full_bias is None or abs(current_full_bias - previous_full_bias) > 50_000_000:
            segment_base_full_bias = current_full_bias
            segment_base_bias_nanos = current_bias_nanos
        previous_full_bias = current_full_bias

        gnss_id = CONSTELLATION_ANDROID.get(_safe_int(raw.get("ConstellationType")), "unknown")
        sv_id = _normalize_qzss_prn(gnss_id, raw.get("Svid"))
        gps_millis = _compute_gps_millis(raw)

        gps_week_nanos = (-current_full_bias // WEEK_NANOS) * WEEK_NANOS
        tx_rx_gnss_diff_ns = _safe_int(raw.get("TimeNanos")) - (
            _safe_int(segment_base_full_bias) + gps_week_nanos
        )
        tx_rx_gnss_diff_ns -= (tx_rx_gnss_diff_ns // WEEK_NANOS) * WEEK_NANOS

        t_rx_ns = float(tx_rx_gnss_diff_ns)
        if gnss_id == "beidou":
            t_rx_ns -= 14_000_000_000.0

        final_rx_tx_ns = (
            t_rx_ns
            - _safe_float(raw.get("ReceivedSvTimeNanos"), default=math.nan)
            + (_safe_float(raw.get("TimeOffsetNanos"), default=0.0) - segment_base_bias_nanos)
        )
        raw_pr_m = final_rx_tx_ns * (C * 1e-9)

        frequency_hz = _truncate_frequency_for_output(raw.get("CarrierFrequencyHz"))
        doppler_hz = -(
            _safe_float(raw.get("PseudorangeRateMetersPerSecond"), default=math.nan) / C
        ) * frequency_hz
        phase_cycle = (
            _safe_float(raw.get("AccumulatedDeltaRangeMeters"), default=math.nan) * frequency_hz / C
        )

        state = _safe_int(raw.get("State"))
        valid_pr = (
            ((state & STATE_CODE_LOCK) != 0)
            and (((state & STATE_TOW_DECODED) != 0) or ((state & STATE_TOW_KNOWN) != 0))
            and ((state & STATE_MSEC_AMBIGUOUS) == 0)
        )

        row = dict(raw)
        row.update(
            {
                "gnss_id": gnss_id,
                "sv_id": sv_id,
                "cn0_dbhz": _safe_float(raw.get("Cn0DbHz"), default=math.nan),
                "accumulated_delta_range_m": _safe_float(
                    raw.get("AccumulatedDeltaRangeMeters"), default=math.nan
                ),
                "accumulated_delta_range_sigma_m": _safe_float(
                    raw.get("AccumulatedDeltaRangeUncertaintyMeters"), default=math.nan
                ),
                "signal_type": CODE_TYPE_ANDROID.get(gnss_id, {}).get(str(raw.get("CodeType", "")), ""),
                "gps_millis": gps_millis,
                "gps_time": int(round(gps_millis / 1000.0)) if not math.isnan(gps_millis) else -1,
                "frequency_hz": frequency_hz,
                "doppler_hz": doppler_hz,
                "phase_cycle": phase_cycle,
                "raw_pr_m": raw_pr_m,
                "raw_pr_sigma_m": _safe_float(
                    raw.get("ReceivedSvTimeUncertaintyNanos"), default=math.nan
                )
                * (C * 1e-9),
                "loi": _safe_int(raw.get("AccumulatedDeltaRangeState")) & 4,
                "valid_pr": valid_pr,
                "segment_base_full_bias": segment_base_full_bias,
                "segment_base_bias_nanos": segment_base_bias_nanos,
            }
        )
        rows.append(row)

    return rows


class AndroidGnssTxtParser:
    def __init__(self, input_path: str | Path) -> None:
        self.input_path = Path(input_path).expanduser().resolve()
        if not self.input_path.exists():
            raise FileNotFoundError(self.input_path)

    def parse_rows(self) -> list[dict[str, object]]:
        return _compute_measurement_rows(_load_raw_rows(self.input_path))

    def parse_measurements(self, filter_valid_pr: bool = True) -> list[Measurement]:
        rows = self.parse_rows()
        if filter_valid_pr:
            rows = [row for row in rows if bool(row["valid_pr"])]
        return rows_to_measurements(rows)

    def convert_to_tsv(
        self,
        output_path: str | Path | None = None,
        filter_valid_pr: bool = True,
    ) -> ParseSummary:
        rows = self.parse_rows()
        target_rows = [row for row in rows if bool(row["valid_pr"])] if filter_valid_pr else rows
        output = Path(output_path) if output_path else self.input_path.with_suffix(".tsv")
        write_measurements_tsv(rows_to_measurements(target_rows), output)
        return ParseSummary(
            input_path=self.input_path,
            total_rows=len(rows),
            valid_rows=len(target_rows),
            output_path=output.resolve(),
        )


def parse_txt_to_rows(input_path: str | Path) -> list[dict[str, object]]:
    return AndroidGnssTxtParser(input_path).parse_rows()


def parse_txt_to_dataframe(input_path: str | Path):
    try:
        import pandas as pd
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "pandas is not installed. Use parse_txt_to_rows() or convert_txt_to_tsv() instead."
        ) from exc
    return pd.DataFrame(parse_txt_to_rows(input_path))


def rows_to_measurements(rows: Sequence[dict[str, object]]) -> list[Measurement]:
    measurements: list[Measurement] = []

    for row in rows:
        measurements.append(
            Measurement(
                time=TimeTag.from_sec(float(_safe_int(row.get("gps_time"), default=-1))),
                sat=SatelliteInfo(
                    sv_pos=[math.nan, math.nan, math.nan],
                    sv_vel=[math.nan, math.nan, math.nan],
                    sv_clock_bias=math.nan,
                    sv_clock_drift=math.nan,
                    pr_correction=math.nan,
                    iono_coff=[0.0] * 8,
                    constellation=str(row.get("gnss_id", "unknown")).upper(),
                    prn=_safe_int(row.get("sv_id"), default=-1),
                    frequency_hz=_safe_float(row.get("frequency_hz"), default=math.nan),
                    code_type=str(row.get("CodeType", "")),
                ),
                value=MeasurementValue(
                    pseudorange_m=_safe_float(row.get("raw_pr_m"), default=math.nan),
                    phase_cycle=_safe_float(row.get("phase_cycle"), default=math.nan),
                    doppler_hz=_safe_float(row.get("doppler_hz"), default=math.nan),
                    snr_dbhz=_safe_float(row.get("cn0_dbhz"), default=math.nan),
                    loi=bool(_safe_int(row.get("loi"))),
                ),
                ground_truth=[math.nan, math.nan, math.nan],
            )
        )

    return measurements


def dataframe_to_measurements(dataframe) -> list[Measurement]:
    if hasattr(dataframe, "to_dict"):
        return rows_to_measurements(dataframe.to_dict(orient="records"))
    return rows_to_measurements(dataframe)


def write_measurements_tsv(measurements: Sequence[Measurement], output_path: str | Path) -> Path:
    target = Path(output_path).expanduser()
    target.parent.mkdir(parents=True, exist_ok=True)

    with target.open("w", encoding="utf8", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        writer.writerow(Measurement.headers())
        for measurement in measurements:
            outfile.write(measurement.to_csv(sep="\t"))
            outfile.write("\n")

    return target.resolve()


def convert_txt_to_tsv(
    input_path: str | Path,
    output_path: str | Path | None = None,
    filter_valid_pr: bool = True,
) -> ParseSummary:
    parser = AndroidGnssTxtParser(input_path)
    return parser.convert_to_tsv(output_path=output_path, filter_valid_pr=filter_valid_pr)


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert Android GNSS txt logs to Measurement TSV.")
    parser.add_argument(
        "input_path",
        nargs="?",
        default="29740_gnss_log_2025_05_07_15_36_12.txt",
        help="Path to Android GNSS txt log.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="data/results/tsv/29740_gnss_log_2025_05_07_15_36_12.tsv",
        help="Output TSV path.",
    )
    parser.add_argument(
        "--keep-invalid",
        action="store_true",
        help="Keep rows that fail the current pseudorange validity filter.",
    )
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = _build_arg_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    summary = convert_txt_to_tsv(
        input_path=args.input_path,
        output_path=args.output,
        filter_valid_pr=not args.keep_invalid,
    )
    print(f"input: {summary.input_path}")
    print(f"rows: {summary.total_rows}")
    print(f"written: {summary.valid_rows}")
    print(f"output: {summary.output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

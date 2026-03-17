from __future__ import annotations

import argparse
import json
from pathlib import Path


CONSTELLATION_ENCODING = {
    0: "GPS",
    1: "GALILEO",
    2: "BEIDOU",
    3: "GLONASS",
    4: "QZSS",
    5: "SBAS",
    6: "IRNSS",
    7: "UNKNOWN",
}


def get_handoff_spec() -> dict:
    return {
        "file_format": {
            "physical_format": "TSV",
            "logical_format": "CSV-like tabular file",
            "delimiter": "\\t",
            "has_header": True,
            "nan_representation": "nan",
            "encoding": "utf-8",
        },
        "row_semantics": {
            "granularity": "one row per satellite measurement at one GNSS epoch",
            "primary_key_hint": [
                "t_sec",
                "constellation",
                "prn",
                "frequency_hz",
                "code_type",
            ],
            "time_reference": "GPS seconds since GPS epoch",
        },
        "constellation_encoding": CONSTELLATION_ENCODING,
        "filled_by_current_parser": [
            "t_sec",
            "gps_week",
            "tow_sec",
            "constellation",
            "prn",
            "iono_a0",
            "iono_a1",
            "iono_a2",
            "iono_a3",
            "iono_b0",
            "iono_b1",
            "iono_b2",
            "iono_b3",
            "frequency_hz",
            "code_type",
            "pseudorange_m",
            "phase_cycle",
            "doppler_hz",
            "snr_dbhz",
            "loi",
        ],
        "currently_empty_or_placeholder": [
            "sv_pos_x",
            "sv_pos_y",
            "sv_pos_z",
            "sv_vel_x",
            "sv_vel_y",
            "sv_vel_z",
            "sv_clock_bias",
            "sv_clock_drift",
            "pr_correction",
            "gt_pos_x",
            "gt_pos_y",
            "gt_pos_z",
        ],
        "columns": [
            {
                "name": "t_sec",
                "type": "int",
                "status": "filled",
                "description": "GPS absolute seconds at measurement epoch.",
            },
            {
                "name": "gps_week",
                "type": "int",
                "status": "filled",
                "description": "GPS week derived from t_sec.",
            },
            {
                "name": "tow_sec",
                "type": "int",
                "status": "filled",
                "description": "Time of week in GPS seconds.",
            },
            {
                "name": "constellation",
                "type": "int",
                "status": "filled",
                "description": "Constellation enum encoded as integer.",
                "encoding": CONSTELLATION_ENCODING,
            },
            {
                "name": "prn",
                "type": "int",
                "status": "filled",
                "description": "Satellite PRN or normalized QZSS SVN-like index.",
            },
            {
                "name": "sv_pos_x",
                "type": "float",
                "status": "empty",
                "description": "Satellite ECEF X position in meters.",
            },
            {
                "name": "sv_pos_y",
                "type": "float",
                "status": "empty",
                "description": "Satellite ECEF Y position in meters.",
            },
            {
                "name": "sv_pos_z",
                "type": "float",
                "status": "empty",
                "description": "Satellite ECEF Z position in meters.",
            },
            {
                "name": "sv_vel_x",
                "type": "float",
                "status": "empty",
                "description": "Satellite ECEF X velocity in m/s.",
            },
            {
                "name": "sv_vel_y",
                "type": "float",
                "status": "empty",
                "description": "Satellite ECEF Y velocity in m/s.",
            },
            {
                "name": "sv_vel_z",
                "type": "float",
                "status": "empty",
                "description": "Satellite ECEF Z velocity in m/s.",
            },
            {
                "name": "iono_a0",
                "type": "float",
                "status": "filled_constant",
                "description": "Ionospheric alpha coefficient 0. Currently placeholder 0.0.",
            },
            {
                "name": "iono_a1",
                "type": "float",
                "status": "filled_constant",
                "description": "Ionospheric alpha coefficient 1. Currently placeholder 0.0.",
            },
            {
                "name": "iono_a2",
                "type": "float",
                "status": "filled_constant",
                "description": "Ionospheric alpha coefficient 2. Currently placeholder 0.0.",
            },
            {
                "name": "iono_a3",
                "type": "float",
                "status": "filled_constant",
                "description": "Ionospheric alpha coefficient 3. Currently placeholder 0.0.",
            },
            {
                "name": "iono_b0",
                "type": "float",
                "status": "filled_constant",
                "description": "Ionospheric beta coefficient 0. Currently placeholder 0.0.",
            },
            {
                "name": "iono_b1",
                "type": "float",
                "status": "filled_constant",
                "description": "Ionospheric beta coefficient 1. Currently placeholder 0.0.",
            },
            {
                "name": "iono_b2",
                "type": "float",
                "status": "filled_constant",
                "description": "Ionospheric beta coefficient 2. Currently placeholder 0.0.",
            },
            {
                "name": "iono_b3",
                "type": "float",
                "status": "filled_constant",
                "description": "Ionospheric beta coefficient 3. Currently placeholder 0.0.",
            },
            {
                "name": "sv_clock_bias",
                "type": "float",
                "status": "empty",
                "description": "Satellite clock bias in meters.",
            },
            {
                "name": "sv_clock_drift",
                "type": "float",
                "status": "empty",
                "description": "Satellite clock drift in meters per second.",
            },
            {
                "name": "pr_correction",
                "type": "float",
                "status": "empty",
                "description": "Additional pseudorange correction in meters.",
            },
            {
                "name": "frequency_hz",
                "type": "float",
                "status": "filled",
                "description": "Carrier frequency in Hz.",
            },
            {
                "name": "code_type",
                "type": "str",
                "status": "filled",
                "description": "Android GnssMeasurement CodeType, such as C, Q, X, I, P.",
            },
            {
                "name": "pseudorange_m",
                "type": "float",
                "status": "filled",
                "description": "Measured pseudorange in meters.",
            },
            {
                "name": "phase_cycle",
                "type": "float",
                "status": "filled",
                "description": "Carrier phase in cycles.",
            },
            {
                "name": "doppler_hz",
                "type": "float",
                "status": "filled",
                "description": "Doppler in Hz.",
            },
            {
                "name": "snr_dbhz",
                "type": "float",
                "status": "filled",
                "description": "Carrier-to-noise density ratio in dB-Hz.",
            },
            {
                "name": "loi",
                "type": "int",
                "status": "filled",
                "description": "Loss-of-lock indicator encoded as 0 or 1-like bit result.",
            },
            {
                "name": "gt_pos_x",
                "type": "float",
                "status": "empty",
                "description": "Receiver ground-truth ECEF X in meters.",
            },
            {
                "name": "gt_pos_y",
                "type": "float",
                "status": "empty",
                "description": "Receiver ground-truth ECEF Y in meters.",
            },
            {
                "name": "gt_pos_z",
                "type": "float",
                "status": "empty",
                "description": "Receiver ground-truth ECEF Z in meters.",
            },
        ],
        "expected_fill_policy_for_other_model": {
            "must_preserve_existing_values": True,
            "should_only_fill_missing_fields": [
                "sv_pos_x",
                "sv_pos_y",
                "sv_pos_z",
                "sv_vel_x",
                "sv_vel_y",
                "sv_vel_z",
                "sv_clock_bias",
                "sv_clock_drift",
                "pr_correction",
                "gt_pos_x",
                "gt_pos_y",
                "gt_pos_z",
            ],
            "should_not_modify": [
                "t_sec",
                "gps_week",
                "tow_sec",
                "constellation",
                "prn",
                "frequency_hz",
                "code_type",
                "pseudorange_m",
                "phase_cycle",
                "doppler_hz",
                "snr_dbhz",
                "loi",
            ],
        },
    }


def get_handoff_prompt() -> str:
    return """You will receive a tab-separated GNSS measurement table.

File contract:
- It is physically TSV, even if we casually call it CSV.
- The first row is a header.
- Each row is one satellite measurement at one GNSS epoch.
- Missing numeric values are written as `nan`.

Column meanings:
- `t_sec`: GPS absolute seconds.
- `gps_week`, `tow_sec`: derived GPS time tags.
- `constellation`: integer enum where 0=GPS, 1=GALILEO, 2=BEIDOU, 3=GLONASS, 4=QZSS, 5=SBAS, 6=IRNSS, 7=UNKNOWN.
- `prn`: satellite PRN.
- `frequency_hz`, `code_type`: signal identifier.
- `pseudorange_m`, `phase_cycle`, `doppler_hz`, `snr_dbhz`, `loi`: already computed measurement values.

Columns that are already valid and must be preserved:
- `t_sec`, `gps_week`, `tow_sec`
- `constellation`, `prn`
- `frequency_hz`, `code_type`
- `pseudorange_m`, `phase_cycle`, `doppler_hz`, `snr_dbhz`, `loi`

Columns that are currently placeholders and should be filled if possible:
- `sv_pos_x`, `sv_pos_y`, `sv_pos_z`
- `sv_vel_x`, `sv_vel_y`, `sv_vel_z`
- `sv_clock_bias`, `sv_clock_drift`
- `pr_correction`
- `gt_pos_x`, `gt_pos_y`, `gt_pos_z`

Important rules:
- Do not alter existing non-missing measurement values.
- Keep row order and header names unchanged.
- Return the same tabular format with the same delimiter and header.
- If you cannot infer a missing field reliably, leave it as `nan`.
"""


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Print GNSS CSV handoff spec for another model.")
    parser.add_argument(
        "--format",
        choices=["json", "prompt"],
        default="json",
        help="Output format.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="",
        help="Optional output file path.",
    )
    return parser


def main() -> int:
    args = _build_arg_parser().parse_args()
    text = (
        json.dumps(get_handoff_spec(), indent=2, ensure_ascii=True)
        if args.format == "json"
        else get_handoff_prompt()
    )

    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(text + "\n", encoding="utf8")
    else:
        print(text)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

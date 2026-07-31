#!/usr/bin/env python3

import re
from pathlib import Path

import pandas as pd


INPUT_FILE = Path("data/raw/PXD043901/proteinGroups.txt")
OUTPUT_FILE = Path("data/processed/PXD043901_sample_metadata.csv")

LFQ_PREFIX = "LFQ intensity "


def parse_sample_name(sample_name: str) -> tuple[str, str]:
    """
    Parse patient and timepoint from MaxQuant sample names.

    Examples
    --------
    ..._BA_19T    -> patient 19, T0
    ..._BA_16_1T  -> patient 16, T1
    ..._BA_2_1T1  -> patient 2, T1
    ..._BA_11_6T  -> patient 11, T6
    """
    match = re.search(r"_BA_(\d+)(.*)$", sample_name)

    if match is None:
        raise ValueError(f"Could not parse sample name: {sample_name}")

    patient_id = match.group(1)
    suffix = match.group(2)

    if suffix == "T":
        timepoint = "T0"
    elif suffix in {"_1T", "_1T1"}:
        timepoint = "T1"
    elif suffix == "_6T":
        timepoint = "T6"
    else:
        raise ValueError(
            f"Unrecognised timepoint suffix '{suffix}' "
            f"in sample: {sample_name}"
        )

    return patient_id, timepoint


def main() -> None:
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

    columns = pd.read_csv(
        INPUT_FILE,
        sep="\t",
        nrows=0,
    ).columns

    lfq_columns = [
        column
        for column in columns
        if column.startswith(LFQ_PREFIX)
    ]

    if len(lfq_columns) != 30:
        raise ValueError(
            f"Expected 30 LFQ columns, found {len(lfq_columns)}"
        )

    records = []

    for lfq_column in lfq_columns:
        raw_sample_name = lfq_column.removeprefix(LFQ_PREFIX)
        patient_id, timepoint = parse_sample_name(raw_sample_name)

        condition = {
            "T0": "baseline",
            "T1": "post_treatment_1_month",
            "T6": "post_treatment_6_months",
        }[timepoint]

        records.append(
            {
                "sample_id": f"P{int(patient_id):02d}_{timepoint}",
                "raw_sample_name": raw_sample_name,
                "lfq_column": lfq_column,
                "patient_id": f"P{int(patient_id):02d}",
                "cohort": "PXD043901",
                "timepoint": timepoint,
                "condition": condition,
                "platform": "Q_Exactive_HF_X_DDA",
            }
        )

    metadata = pd.DataFrame(records)

    metadata["timepoint"] = pd.Categorical(
        metadata["timepoint"],
        categories=["T0", "T1", "T6"],
        ordered=True,
    )

    metadata = metadata.sort_values(
        ["patient_id", "timepoint"]
    ).reset_index(drop=True)

    duplicated = metadata["sample_id"].duplicated()

    if duplicated.any():
        raise ValueError(
            "Duplicated sample IDs detected: "
            + ", ".join(metadata.loc[duplicated, "sample_id"])
        )

    patient_counts = metadata.groupby(
        "patient_id",
        observed=True,
    )["timepoint"].nunique()

    incomplete_patients = patient_counts[patient_counts != 3]

    if not incomplete_patients.empty:
        raise ValueError(
            "Patients without exactly three timepoints:\n"
            + incomplete_patients.to_string()
        )

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(OUTPUT_FILE, index=False)

    print(f"Samples: {len(metadata)}")
    print(f"Patients: {metadata['patient_id'].nunique()}")
    print("\nTimepoints:")
    print(metadata["timepoint"].value_counts().sort_index())
    print(f"\nSaved: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()

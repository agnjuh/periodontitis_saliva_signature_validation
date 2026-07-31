#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd


PROTEIN_GROUPS = Path("data/raw/PXD043901/proteinGroups.txt")
METADATA = Path("data/processed/PXD043901_sample_metadata.csv")

OUTPUT_DIR = Path("data/processed")
ANNOTATION_OUT = OUTPUT_DIR / "PXD043901_protein_annotation.csv"
RAW_MATRIX_OUT = OUTPUT_DIR / "PXD043901_lfq_matrix_raw.csv"
LOG2_MATRIX_OUT = OUTPUT_DIR / "PXD043901_lfq_matrix_log2.csv"
SUMMARY_OUT = OUTPUT_DIR / "PXD043901_matrix_summary.txt"


def is_plus(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).eq("+")


def main() -> None:
    if not PROTEIN_GROUPS.exists():
        raise FileNotFoundError(f"Missing file: {PROTEIN_GROUPS}")

    if not METADATA.exists():
        raise FileNotFoundError(f"Missing file: {METADATA}")

    metadata = pd.read_csv(METADATA)
    protein_groups = pd.read_csv(PROTEIN_GROUPS, sep="\t", low_memory=False)

    lfq_columns = metadata["lfq_column"].tolist()
    missing_lfq = [col for col in lfq_columns if col not in protein_groups.columns]

    if missing_lfq:
        raise ValueError(
            "Missing LFQ columns in proteinGroups.txt:\n"
            + "\n".join(missing_lfq)
        )

    n_initial = len(protein_groups)

    filter_mask = pd.Series(False, index=protein_groups.index)

    for column in ["Reverse", "Potential contaminant", "Only identified by site"]:
        if column in protein_groups.columns:
            filter_mask |= is_plus(protein_groups[column])

    filtered = protein_groups.loc[~filter_mask].copy()
    n_after_flag_filter = len(filtered)

    if "Unique peptides" not in filtered.columns:
        raise ValueError("Column 'Unique peptides' not found")

    filtered["Unique peptides"] = pd.to_numeric(
        filtered["Unique peptides"],
        errors="coerce",
    )

    filtered = filtered.loc[filtered["Unique peptides"] >= 2].copy()
    n_after_peptide_filter = len(filtered)

    matrix = filtered[lfq_columns].apply(pd.to_numeric, errors="coerce")

    matrix = matrix.replace(0, np.nan)

    sample_names = metadata.set_index("lfq_column")["sample_id"].to_dict()
    matrix = matrix.rename(columns=sample_names)

    detected_count = matrix.notna().sum(axis=1)
    filtered["detected_samples"] = detected_count
    filtered["detected_fraction"] = detected_count / matrix.shape[1]

    keep_detected = filtered["detected_fraction"] >= 0.50

    filtered = filtered.loc[keep_detected].copy()
    matrix = matrix.loc[keep_detected].copy()

    n_after_detection_filter = len(filtered)

    preferred_annotation_columns = [
        "Protein IDs",
        "Majority protein IDs",
        "Protein names",
        "Gene names",
        "Fasta headers",
        "Number of proteins",
        "Peptides",
        "Razor + unique peptides",
        "Unique peptides",
        "Sequence coverage [%]",
        "Mol. weight [kDa]",
        "Q-value",
        "Score",
    ]

    annotation_columns = [
        column
        for column in preferred_annotation_columns
        if column in filtered.columns
    ]

    annotation = filtered[annotation_columns].copy()
    annotation.insert(0, "protein_group_id", filtered.index.astype(str))
    annotation["detected_samples"] = filtered["detected_samples"].values
    annotation["detected_fraction"] = filtered["detected_fraction"].values

    matrix.index = annotation["protein_group_id"]
    matrix.index.name = "protein_group_id"

    log2_matrix = np.log2(matrix)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    annotation.to_csv(ANNOTATION_OUT, index=False)
    matrix.to_csv(RAW_MATRIX_OUT)
    log2_matrix.to_csv(LOG2_MATRIX_OUT)

    summary = [
        f"Initial protein groups: {n_initial}",
        f"After reverse/contaminant/site filtering: {n_after_flag_filter}",
        f"After unique peptides >= 2: {n_after_peptide_filter}",
        f"After detection in >= 50% of samples: {n_after_detection_filter}",
        f"Samples: {matrix.shape[1]}",
        f"Missing values in retained matrix: {int(matrix.isna().sum().sum())}",
        f"Missing fraction: {matrix.isna().mean().mean():.4f}",
    ]

    SUMMARY_OUT.write_text("\n".join(summary) + "\n", encoding="utf-8")

    print("\n".join(summary))
    print(f"\nSaved annotation: {ANNOTATION_OUT}")
    print(f"Saved raw LFQ matrix: {RAW_MATRIX_OUT}")
    print(f"Saved log2 LFQ matrix: {LOG2_MATRIX_OUT}")
    print(f"Saved summary: {SUMMARY_OUT}")


if __name__ == "__main__":
    main()

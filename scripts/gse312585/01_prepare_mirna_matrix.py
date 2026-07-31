#!/usr/bin/env python3

from pathlib import Path
import re

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET = "GSE312585"


COUNT_FILE = PROJECT_ROOT / "data" / "processed" / DATASET / "raw_count_matrix.tsv"
METADATA_FILE = PROJECT_ROOT / "data" / "processed" / DATASET / "sample_metadata.tsv"
FEATURE_AUDIT_FILE = PROJECT_ROOT / "results" / DATASET / "tables" / "feature_audit.tsv"

OUTPUT_DIR = PROJECT_ROOT / "data" / "processed" / DATASET
TABLE_DIR = PROJECT_ROOT / "results" / DATASET / "tables"

MIRNA_RAW_OUTPUT = OUTPUT_DIR / "mirna_raw_counts.tsv"
MIRNA_FILTERED_OUTPUT = OUTPUT_DIR / "mirna_filtered_counts.tsv"
MIRNA_CPM_OUTPUT = OUTPUT_DIR / "mirna_log2_cpm.tsv"

MIRNA_AUDIT_OUTPUT = TABLE_DIR / "mirna_feature_audit.tsv"
SAMPLE_QC_OUTPUT = TABLE_DIR / "mirna_sample_qc.tsv"


def canonicalize_feature(feature: str) -> str:
    """
    Remove suffixes such as _1, _2 when they appear to represent
    duplicated annotation entries rather than distinct miRNA names.

    The original feature name is retained in the audit table.
    """
    return re.sub(r"_\d+$", "", feature)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(
        COUNT_FILE,
        sep="\t",
        index_col=0,
    )

    metadata = pd.read_csv(
        METADATA_FILE,
        sep="\t",
        dtype={"donor": str},
    )

    feature_audit = pd.read_csv(
        FEATURE_AUDIT_FILE,
        sep="\t",
    )

    expected_samples = metadata["sample_id"].tolist()

    missing_samples = [
        sample
        for sample in expected_samples
        if sample not in counts.columns
    ]

    if missing_samples:
        raise ValueError(
            f"Samples missing from count matrix: {missing_samples}"
        )

    counts = counts.loc[:, expected_samples]

    mirna_features = feature_audit.loc[
        feature_audit["category"] == "mirna_like",
        "feature",
    ].tolist()

    mirna_counts = counts.loc[
        counts.index.intersection(mirna_features)
    ].copy()

    mirna_counts.index.name = "original_feature"

    audit = pd.DataFrame(
        {
            "original_feature": mirna_counts.index,
            "canonical_feature": [
                canonicalize_feature(x)
                for x in mirna_counts.index
            ],
        }
    )

    canonical_sizes = (
        audit.groupby("canonical_feature")
        .size()
        .rename("n_annotation_rows")
    )

    audit = audit.merge(
        canonical_sizes,
        on="canonical_feature",
        how="left",
    )

    audit["annotation_duplicated"] = (
        audit["n_annotation_rows"] > 1
    )

    audit.to_csv(
        MIRNA_AUDIT_OUTPUT,
        sep="\t",
        index=False,
    )

    # Collapse duplicated annotation rows by summing counts.
    collapsed = mirna_counts.copy()
    collapsed["canonical_feature"] = [
        canonicalize_feature(x)
        for x in collapsed.index
    ]

    collapsed = (
        collapsed
        .groupby("canonical_feature")[expected_samples]
        .sum()
    )

    collapsed.index.name = "feature"

    collapsed.to_csv(
        MIRNA_RAW_OUTPUT,
        sep="\t",
    )

    # Retain miRNAs with at least 10 counts in at least 3 samples.
    keep = (
        collapsed.ge(10).sum(axis=1) >= 3
    )

    filtered = collapsed.loc[keep].copy()

    filtered.to_csv(
        MIRNA_FILTERED_OUTPUT,
        sep="\t",
    )

    library_sizes = filtered.sum(axis=0)

    cpm = filtered.divide(
        library_sizes.replace(0, np.nan),
        axis=1,
    ) * 1_000_000

    log2_cpm = np.log2(cpm + 1)

    log2_cpm.to_csv(
        MIRNA_CPM_OUTPUT,
        sep="\t",
    )

    sample_qc = metadata.copy()

    sample_qc["mirna_total_counts_before_filtering"] = (
        collapsed.sum(axis=0)
        .reindex(expected_samples)
        .values
    )

    sample_qc["mirna_detected_before_filtering"] = (
        collapsed.gt(0).sum(axis=0)
        .reindex(expected_samples)
        .values
    )

    sample_qc["mirna_total_counts_after_filtering"] = (
        filtered.sum(axis=0)
        .reindex(expected_samples)
        .values
    )

    sample_qc["mirna_detected_after_filtering"] = (
        filtered.gt(0).sum(axis=0)
        .reindex(expected_samples)
        .values
    )

    sample_qc.to_csv(
        SAMPLE_QC_OUTPUT,
        sep="\t",
        index=False,
    )

    duplicated_rows = int(
        audit["annotation_duplicated"].sum()
    )

    duplicated_groups = int(
        audit.loc[
            audit["annotation_duplicated"],
            "canonical_feature",
        ].nunique()
    )

    print(f"Original miRNA-like rows: {mirna_counts.shape[0]}")
    print(f"Canonical miRNA features: {collapsed.shape[0]}")
    print(
        "Rows belonging to duplicated annotation groups: "
        f"{duplicated_rows}"
    )
    print(
        "Duplicated canonical annotation groups: "
        f"{duplicated_groups}"
    )
    print(
        "Features retained after filtering: "
        f"{filtered.shape[0]}"
    )

    print("\nmiRNA sample QC:")
    print(
        sample_qc[
            [
                "sample_id",
                "donor",
                "method",
                "mirna_total_counts_before_filtering",
                "mirna_detected_before_filtering",
                "mirna_total_counts_after_filtering",
                "mirna_detected_after_filtering",
            ]
        ].to_string(index=False)
    )

    print("\nOutputs:")
    print(f"  {MIRNA_RAW_OUTPUT}")
    print(f"  {MIRNA_FILTERED_OUTPUT}")
    print(f"  {MIRNA_CPM_OUTPUT}")
    print(f"  {MIRNA_AUDIT_OUTPUT}")
    print(f"  {SAMPLE_QC_OUTPUT}")


if __name__ == "__main__":
    main()

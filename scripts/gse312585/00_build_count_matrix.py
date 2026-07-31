#!/usr/bin/env python3

from pathlib import Path
import re

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET = "GSE312585"


RAW_DIR = PROJECT_ROOT / "data" / "raw" / DATASET
METADATA_FILE = PROJECT_ROOT / "data" / "processed" / DATASET / "sample_metadata.tsv"
OUTPUT_DIR = PROJECT_ROOT / "data" / "processed" / DATASET
RESULTS_DIR = PROJECT_ROOT / "results" / DATASET / "tables"

COUNT_OUTPUT = OUTPUT_DIR / "raw_count_matrix.tsv"
TECHNICAL_OUTPUT = RESULTS_DIR / "technical_count_summary.tsv"
FEATURE_AUDIT_OUTPUT = RESULTS_DIR / "feature_audit.tsv"
SAMPLE_AUDIT_OUTPUT = RESULTS_DIR / "sample_count_audit.tsv"


def identify_sample(filename: str, metadata: pd.DataFrame) -> str:
    for row in metadata.itertuples(index=False):
        if row.gsm_id in filename:
            return row.sample_id

    raise ValueError(
        f"Could not match filename to metadata: {filename}"
    )


def feature_category(feature: str) -> str:
    upper = feature.upper()

    if feature.startswith("__"):
        return "technical"

    if (
        upper.startswith("MIR")
        or upper.startswith("HSA-MIR")
        or upper.startswith("LET-")
    ):
        return "mirna_like"

    if upper.endswith("-AS1") or "-AS" in upper:
        return "antisense_rna"

    if upper.startswith("LINC"):
        return "lincRNA"

    if upper.startswith("SNORD") or upper.startswith("SNORA"):
        return "snoRNA"

    if upper.startswith("RNU") or upper.startswith("RN7"):
        return "small_nuclear_rna"

    if upper.startswith("MT-"):
        return "mitochondrial"

    if re.fullmatch(r"[A-Z0-9._-]+", upper):
        return "gene_symbol_like"

    return "other"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(
        METADATA_FILE,
        sep="\t",
        dtype={"donor": str},
    )

    expected_columns = {
        "sample_id",
        "gsm_id",
        "donor",
        "method",
        "sample_type",
    }

    missing = expected_columns.difference(metadata.columns)
    if missing:
        raise ValueError(
            f"Metadata missing columns: {sorted(missing)}"
        )

    files = sorted(RAW_DIR.glob("GSM*.txt.gz"))

    if len(files) != len(metadata):
        raise ValueError(
            f"Found {len(files)} count files but "
            f"{len(metadata)} metadata rows"
        )

    biological_series = {}
    technical_records = []
    sample_records = []
    feature_sets = {}

    for path in files:
        sample_id = identify_sample(path.name, metadata)

        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["feature", "count"],
            compression="gzip",
        )

        if df["feature"].duplicated().any():
            duplicates = (
                df.loc[
                    df["feature"].duplicated(keep=False),
                    "feature",
                ]
                .unique()
                .tolist()
            )
            raise ValueError(
                f"Duplicated features in {path.name}: "
                f"{duplicates[:10]}"
            )

        technical = df.loc[
            df["feature"].str.startswith("__")
        ].copy()

        biological = df.loc[
            ~df["feature"].str.startswith("__")
        ].copy()

        biological_series[sample_id] = biological.set_index(
            "feature"
        )["count"]

        feature_sets[sample_id] = set(biological["feature"])

        technical_record = {"sample_id": sample_id}

        for row in technical.itertuples(index=False):
            technical_record[row.feature] = int(row.count)

        technical_records.append(technical_record)

        assigned_counts = int(biological["count"].sum())
        technical_counts = int(technical["count"].sum())
        total_counts = assigned_counts + technical_counts
        detected_features = int((biological["count"] > 0).sum())

        sample_records.append(
            {
                "sample_id": sample_id,
                "assigned_counts": assigned_counts,
                "technical_counts": technical_counts,
                "total_reported_counts": total_counts,
                "detected_features": detected_features,
                "n_biological_features": len(biological),
                "assigned_fraction": (
                    assigned_counts / total_counts
                    if total_counts > 0
                    else 0
                ),
            }
        )

    common_features = set.intersection(*feature_sets.values())
    union_features = set.union(*feature_sets.values())

    print(f"Count files: {len(files)}")
    print(f"Features in union: {len(union_features)}")
    print(f"Features common to all samples: {len(common_features)}")

    count_matrix = pd.concat(
        biological_series,
        axis=1,
    ).fillna(0)

    count_matrix = count_matrix.astype(int)
    count_matrix.index.name = "feature"
    count_matrix = count_matrix.sort_index()

    count_matrix.to_csv(
        COUNT_OUTPUT,
        sep="\t",
    )

    technical_df = pd.DataFrame(
        technical_records
    ).fillna(0)

    technical_df.to_csv(
        TECHNICAL_OUTPUT,
        sep="\t",
        index=False,
    )

    sample_audit = pd.DataFrame(sample_records)

    sample_audit = metadata.merge(
        sample_audit,
        on="sample_id",
        how="left",
        validate="one_to_one",
    )

    sample_audit.to_csv(
        SAMPLE_AUDIT_OUTPUT,
        sep="\t",
        index=False,
    )

    feature_audit = pd.DataFrame(
        {"feature": count_matrix.index}
    )

    feature_audit["category"] = feature_audit[
        "feature"
    ].map(feature_category)

    feature_audit["total_count"] = count_matrix.sum(
        axis=1
    ).values

    feature_audit["n_samples_detected"] = (
        count_matrix.gt(0).sum(axis=1).values
    )

    feature_audit.to_csv(
        FEATURE_AUDIT_OUTPUT,
        sep="\t",
        index=False,
    )

    print("\nFeature categories:")
    print(
        feature_audit["category"]
        .value_counts()
        .to_string()
    )

    print("\nSample audit:")
    print(
        sample_audit[
            [
                "sample_id",
                "method",
                "assigned_counts",
                "total_reported_counts",
                "assigned_fraction",
                "detected_features",
            ]
        ].to_string(index=False)
    )

    print("\nOutputs:")
    print(f"  {COUNT_OUTPUT}")
    print(f"  {TECHNICAL_OUTPUT}")
    print(f"  {FEATURE_AUDIT_OUTPUT}")
    print(f"  {SAMPLE_AUDIT_OUTPUT}")


if __name__ == "__main__":
    main()

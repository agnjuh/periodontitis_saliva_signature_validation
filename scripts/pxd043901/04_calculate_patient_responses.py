#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET = "PXD043901"

PROCESSED_DIR = PROJECT_ROOT / "data" / "processed" / DATASET
OUTPUT_DIR = PROJECT_ROOT / "results" / DATASET / "trajectories"

MATRIX_FILE = PROCESSED_DIR / "PXD043901_lfq_matrix_log2.csv"
METADATA_FILE = PROCESSED_DIR / "PXD043901_sample_metadata.csv"
ANNOTATION_FILE = PROCESSED_DIR / "PXD043901_protein_annotation_extended.csv"


def build_delta_matrix(
    matrix: pd.DataFrame,
    patients: list[str],
    later_timepoint: str,
) -> pd.DataFrame:
    deltas = {}

    for patient in patients:
        baseline_sample = f"{patient}_T0"
        later_sample = f"{patient}_{later_timepoint}"

        if baseline_sample not in matrix.columns:
            raise ValueError(f"Missing sample: {baseline_sample}")

        if later_sample not in matrix.columns:
            raise ValueError(f"Missing sample: {later_sample}")

        deltas[patient] = (
            matrix[later_sample] - matrix[baseline_sample]
        )

    delta_matrix = pd.DataFrame(deltas)
    delta_matrix.index.name = "protein_group_id"

    return delta_matrix


def summarise_response(
    delta_matrix: pd.DataFrame,
    comparison: str,
) -> pd.DataFrame:
    positive = delta_matrix.gt(0)
    negative = delta_matrix.lt(0)

    n_observed = delta_matrix.notna().sum(axis=1)
    n_positive = positive.sum(axis=1)
    n_negative = negative.sum(axis=1)

    direction_consistency = pd.concat(
        [n_positive, n_negative],
        axis=1,
    ).max(axis=1).div(n_observed.replace(0, np.nan))

    dominant_direction = np.where(
        n_positive > n_negative,
        "increase",
        np.where(
            n_negative > n_positive,
            "decrease",
            "tie",
        ),
    )

    summary = pd.DataFrame({
        "protein_group_id": delta_matrix.index.astype(str),
        "comparison": comparison,
        "n_observed_pairs": n_observed.values,
        "n_positive": n_positive.values,
        "n_negative": n_negative.values,
        "dominant_direction": dominant_direction,
        "direction_consistency": direction_consistency.values,
        "median_delta_log2": delta_matrix.median(axis=1).values,
        "mean_delta_log2": delta_matrix.mean(axis=1).values,
        "sd_delta_log2": delta_matrix.std(axis=1).values,
        "q25_delta_log2": delta_matrix.quantile(0.25, axis=1).values,
        "q75_delta_log2": delta_matrix.quantile(0.75, axis=1).values,
    })

    return summary


def main() -> None:
    for path in [
        MATRIX_FILE,
        METADATA_FILE,
        ANNOTATION_FILE,
    ]:
        if not path.exists():
            raise FileNotFoundError(f"Missing file: {path}")

    matrix = pd.read_csv(
        MATRIX_FILE,
        index_col=0,
    )

    matrix.index = matrix.index.astype(str)

    metadata = pd.read_csv(METADATA_FILE)

    annotation = pd.read_csv(
        ANNOTATION_FILE,
        dtype={"protein_group_id": str},
    )

    patients = sorted(metadata["patient_id"].unique())

    if len(patients) != 10:
        raise ValueError(
            f"Expected 10 patients, found {len(patients)}"
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    delta_t1 = build_delta_matrix(
        matrix=matrix,
        patients=patients,
        later_timepoint="T1",
    )

    delta_t6 = build_delta_matrix(
        matrix=matrix,
        patients=patients,
        later_timepoint="T6",
    )

    delta_t1.to_csv(
        PROCESSED_DIR / "PXD043901_delta_T1_vs_T0.csv"
    )

    delta_t6.to_csv(
        PROCESSED_DIR / "PXD043901_delta_T6_vs_T0.csv"
    )

    summary_t1 = summarise_response(
        delta_t1,
        comparison="T1_vs_T0",
    )

    summary_t6 = summarise_response(
        delta_t6,
        comparison="T6_vs_T0",
    )

    response_summary = pd.concat(
        [summary_t1, summary_t6],
        ignore_index=True,
    )

    annotation_columns = [
        "protein_group_id",
        "primary_accession",
        "organism_origin",
        "gene_symbol",
        "human_protein_name",
        "first_microbe_label",
        "detected_samples",
        "detected_fraction",
    ]

    response_summary = response_summary.merge(
        annotation[annotation_columns],
        on="protein_group_id",
        how="left",
        validate="many_to_one",
    )

    response_summary.to_csv(
        OUTPUT_DIR / "protein_response_summary.csv",
        index=False,
    )

    human_summary = response_summary.loc[
        response_summary["organism_origin"].eq("human")
    ].copy()

    robust_human = human_summary.loc[
        human_summary["n_observed_pairs"].ge(8)
        & human_summary["direction_consistency"].ge(0.75)
        & human_summary["median_delta_log2"].abs().ge(0.25)
    ].copy()

    robust_human = robust_human.sort_values(
        [
            "comparison",
            "direction_consistency",
            "median_delta_log2",
        ],
        ascending=[True, False, False],
    )

    robust_human.to_csv(
        OUTPUT_DIR / "robust_human_responses_preliminary.csv",
        index=False,
    )

    summary_lines = [
        f"Patients: {len(patients)}",
        f"Protein groups: {matrix.shape[0]}",
        "",
        "T1 versus T0:",
        (
            "  Median observed pairs per protein: "
            f"{summary_t1['n_observed_pairs'].median():.1f}"
        ),
        (
            "  Proteins observed in all 10 pairs: "
            f"{summary_t1['n_observed_pairs'].eq(10).sum()}"
        ),
        "",
        "T6 versus T0:",
        (
            "  Median observed pairs per protein: "
            f"{summary_t6['n_observed_pairs'].median():.1f}"
        ),
        (
            "  Proteins observed in all 10 pairs: "
            f"{summary_t6['n_observed_pairs'].eq(10).sum()}"
        ),
        "",
        (
            "Preliminary robust human responses "
            "(>=8 pairs, consistency >=0.75, "
            "|median delta| >=0.25): "
            f"{len(robust_human)}"
        ),
    ]

    summary_path = OUTPUT_DIR / "response_summary.txt"
    summary_path.write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    print("\n".join(summary_lines))
    print(f"\nSaved outputs to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET = "PXD043901"

OUTPUT_DIR = (
    PROJECT_ROOT
    / "results"
    / DATASET
    / "trajectories"
)

TEST_FILE = (
    OUTPUT_DIR
    / "longitudinal_statistical_tests.csv"
)


def classify_trajectory(row: pd.Series) -> str:
    d1 = row["median_delta_T1"]
    d6 = row["median_delta_T6"]

    c1 = row["consistency_T1"]
    c6 = row["consistency_T6"]

    reliable_t1 = (
        row["n_pairs_T1"] >= 8
        and c1 >= 0.75
        and abs(d1) >= 0.25
    )

    reliable_t6 = (
        row["n_pairs_T6"] >= 8
        and c6 >= 0.75
        and abs(d6) >= 0.25
    )

    if reliable_t1 and reliable_t6:
        if d1 > 0 and d6 > 0:
            if abs(d6) >= 0.75 * abs(d1):
                return "sustained_increase"
            return "early_increase_partial_return"

        if d1 < 0 and d6 < 0:
            if abs(d6) >= 0.75 * abs(d1):
                return "sustained_decrease"
            return "early_decrease_partial_return"

        if d1 > 0 and d6 < 0:
            return "increase_then_decrease"

        if d1 < 0 and d6 > 0:
            return "decrease_then_increase"

    if reliable_t1 and not reliable_t6:
        if d1 > 0:
            return "transient_early_increase"
        return "transient_early_decrease"

    if reliable_t6 and not reliable_t1:
        if d6 > 0:
            return "late_increase"
        return "late_decrease"

    return "weak_or_inconsistent"


def main() -> None:
    if not TEST_FILE.exists():
        raise FileNotFoundError(f"Missing file: {TEST_FILE}")

    tests = pd.read_csv(
        TEST_FILE,
        dtype={"protein_group_id": str},
    )

    human = tests.loc[
        tests["organism_origin"].eq("human")
    ].copy()

    value_columns = [
        "protein_group_id",
        "primary_accession",
        "gene_symbol",
        "human_protein_name",
        "comparison",
        "n_observed_pairs",
        "direction_consistency",
        "median_delta_log2",
        "wilcoxon_p",
        "wilcoxon_fdr",
        "sign_test_p",
        "reproducibility_score",
    ]

    human = human[value_columns]

    wide = human.pivot(
        index=[
            "protein_group_id",
            "primary_accession",
            "gene_symbol",
            "human_protein_name",
        ],
        columns="comparison",
        values=[
            "n_observed_pairs",
            "direction_consistency",
            "median_delta_log2",
            "wilcoxon_p",
            "wilcoxon_fdr",
            "sign_test_p",
            "reproducibility_score",
        ],
    )

    wide.columns = [
        f"{metric}_{comparison}"
        for metric, comparison in wide.columns
    ]

    wide = wide.reset_index()

    rename_map = {
        "n_observed_pairs_T1_vs_T0": "n_pairs_T1",
        "n_observed_pairs_T6_vs_T0": "n_pairs_T6",
        "direction_consistency_T1_vs_T0": "consistency_T1",
        "direction_consistency_T6_vs_T0": "consistency_T6",
        "median_delta_log2_T1_vs_T0": "median_delta_T1",
        "median_delta_log2_T6_vs_T0": "median_delta_T6",
        "wilcoxon_p_T1_vs_T0": "wilcoxon_p_T1",
        "wilcoxon_p_T6_vs_T0": "wilcoxon_p_T6",
        "wilcoxon_fdr_T1_vs_T0": "wilcoxon_fdr_T1",
        "wilcoxon_fdr_T6_vs_T0": "wilcoxon_fdr_T6",
        "sign_test_p_T1_vs_T0": "sign_test_p_T1",
        "sign_test_p_T6_vs_T0": "sign_test_p_T6",
        "reproducibility_score_T1_vs_T0": "score_T1",
        "reproducibility_score_T6_vs_T0": "score_T6",
    }

    wide = wide.rename(columns=rename_map)

    wide["trajectory_class"] = wide.apply(
        classify_trajectory,
        axis=1,
    )

    wide["combined_reproducibility_score"] = (
        wide[["score_T1", "score_T6"]]
        .fillna(0)
        .sum(axis=1)
    )

    wide["persistent_fdr_support"] = (
        wide["wilcoxon_fdr_T1"].le(0.10)
        & wide["wilcoxon_fdr_T6"].le(0.10)
    )

    wide = wide.sort_values(
        [
            "trajectory_class",
            "combined_reproducibility_score",
        ],
        ascending=[True, False],
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    wide.to_csv(
        OUTPUT_DIR / "human_protein_trajectory_classes.csv",
        index=False,
    )

    persistent = wide.loc[
        wide["trajectory_class"].isin(
            ["sustained_increase", "sustained_decrease"]
        )
        & wide["persistent_fdr_support"]
    ].copy()

    persistent = persistent.sort_values(
        "combined_reproducibility_score",
        ascending=False,
    )

    persistent.to_csv(
        OUTPUT_DIR / "persistent_human_responses.csv",
        index=False,
    )

    counts = (
        wide["trajectory_class"]
        .value_counts()
        .rename_axis("trajectory_class")
        .reset_index(name="protein_groups")
    )

    counts.to_csv(
        OUTPUT_DIR / "trajectory_class_counts.csv",
        index=False,
    )

    summary = [
        "Trajectory classes:",
        counts.to_string(index=False),
        "",
        (
            "Persistent responses supported at FDR <= 0.10 "
            f"at both timepoints: {len(persistent)}"
        ),
    ]

    summary_path = OUTPUT_DIR / "trajectory_summary.txt"
    summary_path.write_text(
        "\n".join(summary) + "\n",
        encoding="utf-8",
    )

    print("\n".join(summary))
    print(f"\nSaved: {OUTPUT_DIR / 'human_protein_trajectory_classes.csv'}")
    print(f"Saved: {OUTPUT_DIR / 'persistent_human_responses.csv'}")

if __name__ == "__main__":
    main()

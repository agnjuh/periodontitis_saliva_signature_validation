#!/usr/bin/env python3

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon
from statsmodels.stats.multitest import multipletests


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET = "PXD043901"

PROCESSED_DIR = (
    PROJECT_ROOT / "data" / "processed" / DATASET
)

RESULTS_DIR = (
    PROJECT_ROOT / "results" / DATASET
)

DELTA_FILES = {
    "T1_vs_T0": PROCESSED_DIR / "PXD043901_delta_T1_vs_T0.csv",
    "T6_vs_T0": PROCESSED_DIR / "PXD043901_delta_T6_vs_T0.csv",
}

ANNOTATION_FILE = (
    PROCESSED_DIR
    / "PXD043901_protein_annotation_extended.csv"
)

OUTPUT_DIR = RESULTS_DIR / "trajectories"


def test_protein(values: pd.Series) -> dict:
    values = pd.to_numeric(values, errors="coerce").dropna()
    values = values.loc[values != 0]

    n = len(values)

    if n < 6:
        return {
            "n_tested_pairs": n,
            "wilcoxon_statistic": np.nan,
            "wilcoxon_p": np.nan,
            "sign_test_p": np.nan,
        }

    n_positive = int((values > 0).sum())
    n_negative = int((values < 0).sum())

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            statistic, p_value = wilcoxon(
                values,
                alternative="two-sided",
                zero_method="wilcox",
                method="auto",
            )
    except ValueError:
        statistic = np.nan
        p_value = np.nan

    sign_p = binomtest(
        k=max(n_positive, n_negative),
        n=n,
        p=0.5,
        alternative="greater",
    ).pvalue

    return {
        "n_tested_pairs": n,
        "wilcoxon_statistic": statistic,
        "wilcoxon_p": p_value,
        "sign_test_p": sign_p,
    }


def add_fdr(
    table: pd.DataFrame,
    p_column: str,
    output_column: str,
) -> None:
    valid = table[p_column].notna()

    table[output_column] = np.nan

    if valid.any():
        table.loc[valid, output_column] = multipletests(
            table.loc[valid, p_column],
            method="fdr_bh",
        )[1]


def main() -> None:
    if not ANNOTATION_FILE.exists():
        raise FileNotFoundError(
            f"Missing annotation: {ANNOTATION_FILE}"
        )

    annotation = pd.read_csv(
        ANNOTATION_FILE,
        dtype={"protein_group_id": str},
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

    all_results = []

    for comparison, path in DELTA_FILES.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing delta file: {path}")

        delta = pd.read_csv(path, index_col=0)
        delta.index = delta.index.astype(str)

        records = []

        for protein_group_id, row in delta.iterrows():
            observed = row.dropna()

            n_positive = int((observed > 0).sum())
            n_negative = int((observed < 0).sum())
            n_observed = len(observed)

            consistency = (
                max(n_positive, n_negative) / n_observed
                if n_observed > 0
                else np.nan
            )

            if n_positive > n_negative:
                direction = "increase"
            elif n_negative > n_positive:
                direction = "decrease"
            else:
                direction = "tie"

            tests = test_protein(row)

            records.append({
                "protein_group_id": protein_group_id,
                "comparison": comparison,
                "n_observed_pairs": n_observed,
                "n_positive": n_positive,
                "n_negative": n_negative,
                "dominant_direction": direction,
                "direction_consistency": consistency,
                "median_delta_log2": observed.median(),
                "mean_delta_log2": observed.mean(),
                "sd_delta_log2": observed.std(),
                **tests,
            })

        result = pd.DataFrame(records)

        add_fdr(
            result,
            p_column="wilcoxon_p",
            output_column="wilcoxon_fdr",
        )

        add_fdr(
            result,
            p_column="sign_test_p",
            output_column="sign_test_fdr",
        )

        result = result.merge(
            annotation[annotation_columns],
            on="protein_group_id",
            how="left",
            validate="one_to_one",
        )

        all_results.append(result)

    combined = pd.concat(all_results, ignore_index=True)

    combined["absolute_median_delta"] = (
        combined["median_delta_log2"].abs()
    )

    combined["reproducibility_score"] = (
        combined["direction_consistency"]
        * combined["absolute_median_delta"]
        * np.sqrt(combined["n_observed_pairs"])
    )

    combined.to_csv(
        OUTPUT_DIR / "longitudinal_statistical_tests.csv",
        index=False,
    )

    human = combined.loc[
        combined["organism_origin"].eq("human")
    ].copy()

    stringent = human.loc[
        human["n_observed_pairs"].ge(8)
        & human["direction_consistency"].ge(0.80)
        & human["absolute_median_delta"].ge(0.50)
        & human["wilcoxon_fdr"].le(0.10)
    ].copy()

    stringent = stringent.sort_values(
        [
            "comparison",
            "wilcoxon_fdr",
            "reproducibility_score",
        ],
        ascending=[True, True, False],
    )

    stringent.to_csv(
        OUTPUT_DIR / "stringent_human_responses.csv",
        index=False,
    )

    ranked = human.loc[
        human["n_observed_pairs"].ge(8)
    ].sort_values(
        [
            "comparison",
            "reproducibility_score",
        ],
        ascending=[True, False],
    )

    ranked.to_csv(
        OUTPUT_DIR / "ranked_human_responses.csv",
        index=False,
    )

    summary_lines = []

    for comparison in DELTA_FILES:
        subset = human.loc[
            human["comparison"].eq(comparison)
        ]

        stringent_subset = stringent.loc[
            stringent["comparison"].eq(comparison)
        ]

        summary_lines.extend([
            comparison,
            (
                "  Human protein groups tested: "
                f"{len(subset)}"
            ),
            (
                "  Wilcoxon FDR < 0.05: "
                f"{subset['wilcoxon_fdr'].lt(0.05).sum()}"
            ),
            (
                "  Wilcoxon FDR < 0.10: "
                f"{subset['wilcoxon_fdr'].lt(0.10).sum()}"
            ),
            (
                "  Stringent reproducible responses: "
                f"{len(stringent_subset)}"
            ),
            "",
        ])

    summary_path = OUTPUT_DIR / "statistical_test_summary.txt"
    summary_path.write_text(
        "\n".join(summary_lines),
        encoding="utf-8",
    )

    print("\n".join(summary_lines))
    print(
        "Saved: "
        "results/trajectories/"
        "longitudinal_statistical_tests.csv"
    )
    print(
        "Saved: "
        "results/trajectories/"
        "stringent_human_responses.csv"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET = "PXD043901"

PROCESSED_DIR = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / DATASET
)

RESULTS_DIR = (
    PROJECT_ROOT
    / "results"
    / DATASET
)

OUTPUT_DIR = RESULTS_DIR / "modules"

MATRIX_FILE = (
    PROCESSED_DIR
    / "PXD043901_lfq_matrix_log2.csv"
)

ANNOTATION_FILE = (
    PROCESSED_DIR
    / "PXD043901_protein_annotation_extended.csv"
)

METADATA_FILE = (
    PROCESSED_DIR
    / "PXD043901_sample_metadata.csv"
)

MODULE_FILE = (
    PROJECT_ROOT
    / "config"
    / "biological_modules.csv"
)


def zscore_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    means = matrix.mean(axis=1, skipna=True)

    standard_deviations = matrix.std(
        axis=1,
        skipna=True,
        ddof=1,
    )

    standard_deviations = standard_deviations.replace(
        0,
        np.nan,
    )

    return matrix.sub(
        means,
        axis=0,
    ).div(
        standard_deviations,
        axis=0,
    )


def safe_wilcoxon(
    values: pd.Series,
) -> tuple[float, float]:
    values = pd.to_numeric(
        values,
        errors="coerce",
    ).dropna()

    values = values.loc[
        values != 0
    ]

    if len(values) < 6:
        return np.nan, np.nan

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            statistic, p_value = wilcoxon(
                values,
                alternative="two-sided",
                zero_method="wilcox",
                method="auto",
            )

        return statistic, p_value

    except ValueError:
        return np.nan, np.nan


def add_fdr(
    table: pd.DataFrame,
    p_column: str,
    output_column: str,
) -> None:
    table[output_column] = np.nan

    valid = table[p_column].notna()

    if valid.any():
        table.loc[
            valid,
            output_column,
        ] = multipletests(
            table.loc[
                valid,
                p_column,
            ],
            method="fdr_bh",
        )[1]


def main() -> None:
    required_files = [
        MATRIX_FILE,
        ANNOTATION_FILE,
        METADATA_FILE,
        MODULE_FILE,
    ]

    for path in required_files:
        if not path.exists():
            raise FileNotFoundError(
                f"Missing required file: {path}"
            )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    matrix = pd.read_csv(
        MATRIX_FILE,
        index_col=0,
    )

    matrix.index = matrix.index.astype(str)

    annotation = pd.read_csv(
        ANNOTATION_FILE,
        dtype={
            "protein_group_id": str,
        },
    )

    metadata = pd.read_csv(
        METADATA_FILE,
    )

    modules = pd.read_csv(
        MODULE_FILE,
    )

    required_module_columns = {
        "module",
        "gene_symbol",
        "expected_treatment_direction",
        "module_description",
    }

    missing_module_columns = (
        required_module_columns
        - set(modules.columns)
    )

    if missing_module_columns:
        raise ValueError(
            "Missing module columns: "
            f"{sorted(missing_module_columns)}"
        )

    annotation_subset = annotation[
        [
            "protein_group_id",
            "gene_symbol",
            "organism_origin",
        ]
    ].copy()

    human_annotation = annotation_subset.loc[
        annotation_subset[
            "organism_origin"
        ].eq("human")
        & annotation_subset[
            "gene_symbol"
        ].notna()
    ].copy()

    human_matrix = matrix.loc[
        matrix.index.intersection(
            human_annotation[
                "protein_group_id"
            ]
        )
    ].copy()

    mapped_annotation = human_annotation.loc[
        human_annotation[
            "protein_group_id"
        ].isin(
            human_matrix.index
        )
    ].copy()

    long_matrix = (
        human_matrix
        .reset_index()
        .rename(
            columns={
                (
                    human_matrix.index.name
                    if human_matrix.index.name
                    is not None
                    else "index"
                ): "protein_group_id"
            }
        )
    )

    if (
        "protein_group_id"
        not in long_matrix.columns
    ):
        long_matrix = long_matrix.rename(
            columns={
                long_matrix.columns[0]:
                    "protein_group_id"
            }
        )

    long_matrix = long_matrix.merge(
        mapped_annotation[
            [
                "protein_group_id",
                "gene_symbol",
            ]
        ],
        on="protein_group_id",
        how="inner",
        validate="one_to_one",
    )

    sample_columns = [
        column
        for column in matrix.columns
        if column in long_matrix.columns
    ]

    if not sample_columns:
        raise ValueError(
            "No sample columns from the LFQ matrix "
            "were retained after annotation mapping."
        )

    gene_matrix = (
        long_matrix
        .groupby(
            "gene_symbol",
            dropna=True,
        )[sample_columns]
        .median()
    )

    gene_matrix.index = (
        gene_matrix.index.astype(str)
    )

    standardized = zscore_rows(
        gene_matrix
    )

    module_records = []
    coverage_records = []

    for (
        module_name,
        module_table,
    ) in modules.groupby(
        "module",
        sort=True,
    ):
        requested_genes = (
            module_table[
                "gene_symbol"
            ]
            .dropna()
            .astype(str)
            .drop_duplicates()
            .tolist()
        )

        available_genes = [
            gene
            for gene in requested_genes
            if gene in standardized.index
        ]

        missing_genes = sorted(
            set(requested_genes)
            - set(available_genes)
        )

        direction_values = (
            module_table[
                "expected_treatment_direction"
            ]
            .dropna()
            .unique()
        )

        if len(direction_values) != 1:
            raise ValueError(
                f"Module {module_name} contains "
                "inconsistent expected directions."
            )

        expected_direction = (
            direction_values[0]
        )

        if expected_direction not in {
            "increase",
            "decrease",
        }:
            raise ValueError(
                "Invalid direction in module "
                f"{module_name}: "
                f"{expected_direction}"
            )

        direction_multiplier = (
            1.0
            if expected_direction == "increase"
            else -1.0
        )

        description_values = (
            module_table[
                "module_description"
            ]
            .dropna()
        )

        if description_values.empty:
            description = ""
        else:
            description = (
                description_values.iloc[0]
            )

        coverage_records.append(
            {
                "module": module_name,
                "module_description":
                    description,
                "expected_treatment_direction":
                    expected_direction,
                "requested_genes":
                    len(requested_genes),
                "available_genes":
                    len(available_genes),
                "coverage_fraction": (
                    len(available_genes)
                    / len(requested_genes)
                    if requested_genes
                    else np.nan
                ),
                "available_gene_symbols":
                    ";".join(
                        available_genes
                    ),
                "missing_gene_symbols":
                    ";".join(
                        missing_genes
                    ),
            }
        )

        if len(available_genes) < 2:
            print(
                f"Skipping {module_name}: "
                "fewer than two proteins "
                "were available."
            )
            continue

        module_abundance = (
            standardized.loc[
                available_genes
            ]
            .mean(
                axis=0,
                skipna=True,
            )
        )

        response_score = (
            module_abundance
            * direction_multiplier
        )

        n_available_per_sample = (
            standardized.loc[
                available_genes
            ]
            .notna()
            .sum(axis=0)
        )

        for sample_id in sample_columns:
            module_records.append(
                {
                    "sample_id":
                        sample_id,
                    "module":
                        module_name,
                    "module_description":
                        description,
                    "expected_treatment_direction":
                        expected_direction,
                    "module_abundance_score":
                        module_abundance.get(
                            sample_id,
                            np.nan,
                        ),
                    "response_oriented_score":
                        response_score.get(
                            sample_id,
                            np.nan,
                        ),
                    "available_proteins_in_module":
                        n_available_per_sample.get(
                            sample_id,
                            0,
                        ),
                    "module_size":
                        len(available_genes),
                }
            )

    coverage = pd.DataFrame(
        coverage_records
    )

    scores = pd.DataFrame(
        module_records
    )

    if scores.empty:
        raise ValueError(
            "No module scores were generated. "
            "Check module gene coverage."
        )

    metadata_columns = [
        "sample_id",
        "patient_id",
        "timepoint",
    ]

    missing_metadata_columns = [
        column
        for column in metadata_columns
        if column not in metadata.columns
    ]

    if missing_metadata_columns:
        raise ValueError(
            "Metadata is missing columns: "
            f"{missing_metadata_columns}"
        )

    scores = scores.merge(
        metadata[
            metadata_columns
        ],
        on="sample_id",
        how="left",
        validate="many_to_one",
    )

    if scores["patient_id"].isna().any():
        unmatched = (
            scores.loc[
                scores[
                    "patient_id"
                ].isna(),
                "sample_id",
            ]
            .drop_duplicates()
            .tolist()
        )

        raise ValueError(
            "Module scores could not be "
            "matched to metadata: "
            f"{unmatched}"
        )

    scores.to_csv(
        OUTPUT_DIR
        / "sample_module_scores.csv",
        index=False,
    )

    coverage.to_csv(
        OUTPUT_DIR
        / "module_coverage.csv",
        index=False,
    )

    patient_wide = scores.pivot_table(
        index=[
            "module",
            "patient_id",
        ],
        columns="timepoint",
        values="response_oriented_score",
        aggfunc="first",
    ).reset_index()

    required_timepoints = [
        "T0",
        "T1",
        "T6",
    ]

    for timepoint in required_timepoints:
        if timepoint not in patient_wide.columns:
            raise ValueError(
                "Missing timepoint in scores: "
                f"{timepoint}"
            )

    patient_wide[
        "delta_T1_vs_T0"
    ] = (
        patient_wide["T1"]
        - patient_wide["T0"]
    )

    patient_wide[
        "delta_T6_vs_T0"
    ] = (
        patient_wide["T6"]
        - patient_wide["T0"]
    )

    patient_wide.to_csv(
        OUTPUT_DIR
        / "patient_module_response_deltas.csv",
        index=False,
    )

    statistical_records = []

    for (
        module_name,
        module_data,
    ) in patient_wide.groupby(
        "module",
        sort=True,
    ):
        for (
            comparison,
            delta_column,
        ) in [
            (
                "T1_vs_T0",
                "delta_T1_vs_T0",
            ),
            (
                "T6_vs_T0",
                "delta_T6_vs_T0",
            ),
        ]:
            delta = (
                module_data[
                    delta_column
                ]
                .dropna()
            )

            n_pairs = len(delta)

            n_positive = int(
                (delta > 0).sum()
            )

            n_negative = int(
                (delta < 0).sum()
            )

            consistency = (
                max(
                    n_positive,
                    n_negative,
                )
                / n_pairs
                if n_pairs > 0
                else np.nan
            )

            statistic, p_value = (
                safe_wilcoxon(
                    delta
                )
            )

            statistical_records.append(
                {
                    "module":
                        module_name,
                    "comparison":
                        comparison,
                    "n_pairs":
                        n_pairs,
                    "n_positive_response":
                        n_positive,
                    "n_negative_response":
                        n_negative,
                    "direction_consistency":
                        consistency,
                    "median_response_delta":
                        delta.median(),
                    "mean_response_delta":
                        delta.mean(),
                    "sd_response_delta":
                        delta.std(ddof=1),
                    "wilcoxon_statistic":
                        statistic,
                    "wilcoxon_p":
                        p_value,
                }
            )

    statistics = pd.DataFrame(
        statistical_records
    )

    add_fdr(
        statistics,
        p_column="wilcoxon_p",
        output_column="wilcoxon_fdr",
    )

    statistics = statistics.sort_values(
        [
            "comparison",
            "wilcoxon_fdr",
            "median_response_delta",
        ],
        ascending=[
            True,
            True,
            False,
        ],
        na_position="last",
    )

    statistics.to_csv(
        OUTPUT_DIR
        / "module_longitudinal_statistics.csv",
        index=False,
    )

    summary_lines = [
        "Module coverage",
        coverage[
            [
                "module",
                "requested_genes",
                "available_genes",
                "coverage_fraction",
            ]
        ].to_string(
            index=False
        ),
        "",
        "Longitudinal module statistics",
        statistics[
            [
                "module",
                "comparison",
                "n_pairs",
                "n_positive_response",
                "direction_consistency",
                "median_response_delta",
                "wilcoxon_p",
                "wilcoxon_fdr",
            ]
        ].to_string(
            index=False
        ),
    ]

    summary_path = (
        OUTPUT_DIR
        / "module_analysis_summary.txt"
    )

    summary_path.write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    print(
        "\n".join(summary_lines)
    )

    print(
        f"\nSaved outputs to: {OUTPUT_DIR}"
    )


if __name__ == "__main__":
    main()

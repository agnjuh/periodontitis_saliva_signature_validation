#!/usr/bin/env python3

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon


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

OUTPUT_DIR = (
    RESULTS_DIR
    / "modules"
    / "sensitivity"
)


def zscore_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    means = matrix.mean(axis=1, skipna=True)
    sds = matrix.std(axis=1, skipna=True, ddof=1)
    sds = sds.replace(0, np.nan)

    return matrix.sub(means, axis=0).div(sds, axis=0)


def safe_wilcoxon(values: pd.Series) -> tuple[float, float]:
    values = pd.to_numeric(
        values,
        errors="coerce",
    ).dropna()

    values = values.loc[values != 0]

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

        return float(statistic), float(p_value)

    except ValueError:
        return np.nan, np.nan


def calculate_module_deltas(
    standardized: pd.DataFrame,
    genes: list[str],
    direction_multiplier: float,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    abundance = (
        standardized.loc[genes]
        .mean(axis=0, skipna=True)
    )

    oriented = abundance * direction_multiplier

    score_table = pd.DataFrame({
        "sample_id": oriented.index,
        "score": oriented.values,
    })

    score_table = score_table.merge(
        metadata[
            [
                "sample_id",
                "patient_id",
                "timepoint",
            ]
        ],
        on="sample_id",
        how="inner",
        validate="one_to_one",
    )

    wide = score_table.pivot(
        index="patient_id",
        columns="timepoint",
        values="score",
    )

    for required_timepoint in ["T0", "T1", "T6"]:
        if required_timepoint not in wide.columns:
            raise ValueError(
                f"Missing timepoint: {required_timepoint}"
            )

    return pd.DataFrame({
        "patient_id": wide.index,
        "delta_T1_vs_T0": wide["T1"] - wide["T0"],
        "delta_T6_vs_T0": wide["T6"] - wide["T0"],
    }).reset_index(drop=True)


def summarise_delta(
    values: pd.Series,
) -> dict:
    values = values.dropna()

    n_pairs = len(values)
    n_positive = int((values > 0).sum())
    n_negative = int((values < 0).sum())

    consistency = (
        max(n_positive, n_negative) / n_pairs
        if n_pairs > 0
        else np.nan
    )

    statistic, p_value = safe_wilcoxon(values)

    return {
        "n_pairs": n_pairs,
        "n_positive_response": n_positive,
        "n_negative_response": n_negative,
        "direction_consistency": consistency,
        "median_response_delta": values.median(),
        "mean_response_delta": values.mean(),
        "wilcoxon_statistic": statistic,
        "wilcoxon_p": p_value,
    }


def main() -> None:
    for path in [
        MATRIX_FILE,
        ANNOTATION_FILE,
        METADATA_FILE,
        MODULE_FILE,
    ]:
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
        dtype={"protein_group_id": str},
    )

    metadata = pd.read_csv(
        METADATA_FILE
    )

    modules = pd.read_csv(
        MODULE_FILE
    )

    mapping = annotation.loc[
        annotation["organism_origin"].eq("human")
        & annotation["gene_symbol"].notna(),
        [
            "protein_group_id",
            "gene_symbol",
        ],
    ].copy()

    retained_ids = matrix.index.intersection(
        mapping["protein_group_id"]
    )

    long_matrix = (
        matrix.loc[retained_ids]
        .reset_index()
    )

    long_matrix = long_matrix.rename(
        columns={
            long_matrix.columns[0]:
                "protein_group_id"
        }
    )

    long_matrix = long_matrix.merge(
        mapping,
        on="protein_group_id",
        how="inner",
        validate="one_to_one",
    )

    sample_columns = [
        sample_id
        for sample_id in metadata["sample_id"]
        if sample_id in long_matrix.columns
    ]

    gene_matrix = (
        long_matrix
        .groupby("gene_symbol")[sample_columns]
        .median()
    )

    standardized = zscore_rows(gene_matrix)

    leave_one_out_records = []
    full_module_deltas = []

    for module_name, module_table in modules.groupby(
        "module",
        sort=True,
    ):
        requested_genes = (
            module_table["gene_symbol"]
            .astype(str)
            .drop_duplicates()
            .tolist()
        )

        available_genes = [
            gene
            for gene in requested_genes
            if gene in standardized.index
        ]

        direction = (
            module_table[
                "expected_treatment_direction"
            ]
            .dropna()
            .unique()
        )

        if len(direction) != 1:
            raise ValueError(
                f"Inconsistent direction: {module_name}"
            )

        multiplier = (
            1.0
            if direction[0] == "increase"
            else -1.0
        )

        full_deltas = calculate_module_deltas(
            standardized=standardized,
            genes=available_genes,
            direction_multiplier=multiplier,
            metadata=metadata,
        )

        full_deltas["module"] = module_name

        full_module_deltas.append(
            full_deltas
        )

        gene_sets = [
            (
                "__none__",
                available_genes,
            )
        ]

        if len(available_genes) >= 3:
            gene_sets.extend([
                (
                    removed_gene,
                    [
                        gene
                        for gene in available_genes
                        if gene != removed_gene
                    ],
                )
                for removed_gene in available_genes
            ])

        for removed_gene, retained_genes in gene_sets:
            deltas = calculate_module_deltas(
                standardized=standardized,
                genes=retained_genes,
                direction_multiplier=multiplier,
                metadata=metadata,
            )

            for comparison in [
                "T1_vs_T0",
                "T6_vs_T0",
            ]:
                delta_column = (
                    f"delta_{comparison}"
                )

                summary = summarise_delta(
                    deltas[delta_column]
                )

                leave_one_out_records.append({
                    "module": module_name,
                    "removed_gene": removed_gene,
                    "retained_genes": len(
                        retained_genes
                    ),
                    "comparison": comparison,
                    **summary,
                })

    leave_one_out = pd.DataFrame(
        leave_one_out_records
    )

    leave_one_out.to_csv(
        OUTPUT_DIR
        / "leave_one_protein_out_results.csv",
        index=False,
    )

    sensitivity_summary = (
        leave_one_out.loc[
            leave_one_out["removed_gene"]
            .ne("__none__")
        ]
        .groupby(
            [
                "module",
                "comparison",
            ],
            as_index=False,
        )
        .agg(
            removals_tested=(
                "removed_gene",
                "nunique",
            ),
            minimum_positive_responses=(
                "n_positive_response",
                "min",
            ),
            minimum_direction_consistency=(
                "direction_consistency",
                "min",
            ),
            minimum_median_delta=(
                "median_response_delta",
                "min",
            ),
            maximum_median_delta=(
                "median_response_delta",
                "max",
            ),
            maximum_wilcoxon_p=(
                "wilcoxon_p",
                "max",
            ),
        )
    )

    sensitivity_summary[
        "robust_to_all_single_gene_removals"
    ] = (
        sensitivity_summary[
            "minimum_direction_consistency"
        ].ge(0.80)
        & sensitivity_summary[
            "minimum_median_delta"
        ].gt(0)
        & sensitivity_summary[
            "maximum_wilcoxon_p"
        ].lt(0.05)
    )

    sensitivity_summary.to_csv(
        OUTPUT_DIR
        / "leave_one_out_summary.csv",
        index=False,
    )

    combined_deltas = pd.concat(
        full_module_deltas,
        ignore_index=True,
    )

    correlation_records = []

    for comparison in [
        "T1_vs_T0",
        "T6_vs_T0",
    ]:
        delta_column = f"delta_{comparison}"

        wide = combined_deltas.pivot(
            index="patient_id",
            columns="module",
            values=delta_column,
        )

        module_names = list(wide.columns)

        for index_a, module_a in enumerate(
            module_names
        ):
            for module_b in module_names[
                index_a + 1:
            ]:
                valid = wide[
                    [
                        module_a,
                        module_b,
                    ]
                ].dropna()

                if len(valid) < 4:
                    rho = np.nan
                    p_value = np.nan
                else:
                    rho, p_value = spearmanr(
                        valid[module_a],
                        valid[module_b],
                    )

                correlation_records.append({
                    "comparison": comparison,
                    "module_a": module_a,
                    "module_b": module_b,
                    "n_patients": len(valid),
                    "spearman_rho": rho,
                    "spearman_p": p_value,
                })

    correlations = pd.DataFrame(
        correlation_records
    )

    correlations.to_csv(
        OUTPUT_DIR
        / "module_delta_correlations.csv",
        index=False,
    )

    high_correlations = correlations.loc[
        correlations["spearman_rho"]
        .abs()
        .ge(0.70)
    ].sort_values(
        [
            "comparison",
            "spearman_rho",
        ],
        ascending=[
            True,
            False,
        ],
    )

    high_correlations.to_csv(
        OUTPUT_DIR
        / "high_module_correlations.csv",
        index=False,
    )

    summary_lines = [
        "Leave-one-protein-out sensitivity",
        sensitivity_summary.to_string(
            index=False
        ),
        "",
        "Module-pair correlations with |rho| >= 0.70",
        (
            high_correlations.to_string(
                index=False
            )
            if not high_correlations.empty
            else "None"
        ),
        "",
        "Important interpretation:",
        (
            "These are exploratory sensitivity results. "
            "The modules were derived from the same cohort "
            "and require independent validation."
        ),
    ]

    summary_path = (
        OUTPUT_DIR
        / "module_sensitivity_summary.txt"
    )

    summary_path.write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    print("\n".join(summary_lines))

    print(f"\nSaved outputs to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()

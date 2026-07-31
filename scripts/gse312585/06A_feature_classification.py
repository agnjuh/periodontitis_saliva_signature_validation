#!/usr/bin/env python3

from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET = "GSE312585"


INPUT_FILE = PROJECT_ROOT / "results" / DATASET / "tables" / "consensus_method_sensitive_signature.tsv"

TABLE_DIR = PROJECT_ROOT / "results" / DATASET / "tables"
FIGURE_DIR = PROJECT_ROOT / "results" / DATASET / "figures"

TOP_FEATURES = 25


def classify_feature_name(feature: str) -> tuple[str, str]:
    """
    Classify a feature using conservative gene-symbol conventions.

    The classification is strictly name-based and does not replace
    authoritative annotation from GENCODE, Ensembl, or miRBase.

    Returns
    -------
    feature_class
        Broad feature category.

    classification_basis
        Rule used for classification.
    """
    symbol = str(feature).strip()
    upper = symbol.upper()

    # miRNA host genes must be detected before the generic MIR rule.
    if upper.endswith("HG") and upper.startswith("MIR"):
        return (
            "probable_miRNA_host_gene",
            "MIR-prefixed symbol ending in HG",
        )

    if upper.startswith("MIR"):
        return (
            "miRNA_associated_locus",
            "MIR-prefixed gene symbol",
        )

    if upper.startswith("LINC"):
        return (
            "long_intergenic_noncoding_RNA",
            "LINC-prefixed gene symbol",
        )

    if upper.startswith(("SNORD", "SNORA")):
        return (
            "small_nucleolar_RNA",
            "SNORD/SNORA-prefixed gene symbol",
        )

    if upper.startswith(
        (
            "RNU",
            "U1",
            "U2",
            "U3",
            "U4",
            "U5",
            "U6",
            "U7",
            "U8",
            "U11",
            "U12",
        )
    ):
        return (
            "small_nuclear_RNA",
            "RNU/U-prefixed gene symbol",
        )

    if upper.startswith(("MT-", "MT_")):
        return (
            "mitochondrial_feature",
            "MT-prefixed gene symbol",
        )

    if upper.startswith(("RPL", "RPS")):
        return (
            "ribosomal_protein_gene",
            "RPL/RPS-prefixed gene symbol",
        )

    if upper.startswith(("RNA5S", "RNA18S", "RNA28S", "RNA45S")):
        return (
            "ribosomal_RNA",
            "ribosomal RNA gene-symbol prefix",
        )

    if upper.startswith("Y_RNA") or re.match(r"^RNY\d+", upper):
        return (
            "Y_RNA",
            "Y_RNA/RNY-prefixed gene symbol",
        )

    if upper.startswith("VAULT") or upper.startswith("VTRNA"):
        return (
            "vault_RNA",
            "VAULT/VTRNA-prefixed gene symbol",
        )

    return (
        "other_or_unresolved_feature",
        "no conservative name-based rule matched",
    )


def validate_input(data: pd.DataFrame) -> None:
    required_columns = {
        "feature",
        "robust_contrasts",
        "cumulative_robust_effect",
        "maximum_absolute_effect",
        "consensus_rank",
    }

    missing = required_columns.difference(data.columns)

    if missing:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(sorted(missing))
        )

    if data["feature"].duplicated().any():
        duplicated = (
            data.loc[
                data["feature"].duplicated(keep=False),
                "feature",
            ]
            .astype(str)
            .tolist()
        )

        raise ValueError(
            "Duplicated feature names detected: "
            + ", ".join(sorted(set(duplicated)))
        )


def build_classified_table(
    data: pd.DataFrame,
) -> pd.DataFrame:
    classified = data.copy()

    assignments = classified["feature"].map(
        classify_feature_name
    )

    classified["feature_class"] = assignments.map(
        lambda value: value[0]
    )

    classified["classification_basis"] = assignments.map(
        lambda value: value[1]
    )

    classified["method_sensitivity_tier"] = pd.cut(
        classified["robust_contrasts"],
        bins=[-np.inf, 1, 3, 5, np.inf],
        labels=[
            "limited",
            "moderate",
            "high",
            "pan_contrast",
        ],
    ).astype(str)

    preferred_columns = [
        "consensus_rank",
        "feature",
        "feature_class",
        "classification_basis",
        "method_sensitivity_tier",
        "robust_contrasts",
        "robust_contrast_names",
        "cumulative_robust_effect",
        "maximum_absolute_effect",
        "median_absolute_effect",
        "mean_absolute_effect",
        "mean_direction_consistency",
        "fully_direction_consistent_contrasts",
        "leave_one_out_stable_contrasts",
        "contrasts_evaluated",
    ]

    existing_columns = [
        column
        for column in preferred_columns
        if column in classified.columns
    ]

    remaining_columns = [
        column
        for column in classified.columns
        if column not in existing_columns
    ]

    classified = classified[
        existing_columns + remaining_columns
    ]

    classified = classified.sort_values(
        "consensus_rank"
    ).reset_index(drop=True)

    return classified


def build_class_summary(
    classified: pd.DataFrame,
) -> pd.DataFrame:
    summary = (
        classified.groupby(
            "feature_class",
            dropna=False,
        )
        .agg(
            feature_count=("feature", "size"),
            median_robust_contrasts=(
                "robust_contrasts",
                "median",
            ),
            maximum_robust_contrasts=(
                "robust_contrasts",
                "max",
            ),
            median_cumulative_robust_effect=(
                "cumulative_robust_effect",
                "median",
            ),
            maximum_cumulative_robust_effect=(
                "cumulative_robust_effect",
                "max",
            ),
            median_maximum_absolute_effect=(
                "maximum_absolute_effect",
                "median",
            ),
        )
        .reset_index()
    )

    summary["percentage_of_signature"] = (
        100.0
        * summary["feature_count"]
        / summary["feature_count"].sum()
    )

    summary = summary.sort_values(
        [
            "feature_count",
            "median_cumulative_robust_effect",
        ],
        ascending=[False, False],
    ).reset_index(drop=True)

    return summary


def build_tier_summary(
    classified: pd.DataFrame,
) -> pd.DataFrame:
    tier_order = [
        "pan_contrast",
        "high",
        "moderate",
        "limited",
    ]

    summary = (
        classified.groupby(
            "method_sensitivity_tier",
            observed=False,
        )
        .agg(
            feature_count=("feature", "size"),
            median_robust_contrasts=(
                "robust_contrasts",
                "median",
            ),
            median_cumulative_robust_effect=(
                "cumulative_robust_effect",
                "median",
            ),
            maximum_cumulative_robust_effect=(
                "cumulative_robust_effect",
                "max",
            ),
        )
        .reset_index()
    )

    summary["method_sensitivity_tier"] = pd.Categorical(
        summary["method_sensitivity_tier"],
        categories=tier_order,
        ordered=True,
    )

    summary = summary.sort_values(
        "method_sensitivity_tier"
    ).reset_index(drop=True)

    summary["method_sensitivity_tier"] = (
        summary["method_sensitivity_tier"].astype(str)
    )

    return summary


def save_class_distribution_plot(
    class_summary: pd.DataFrame,
) -> None:
    plot_data = class_summary.sort_values(
        "feature_count",
        ascending=True,
    )

    figure_height = max(
        4.5,
        0.55 * len(plot_data),
    )

    fig, ax = plt.subplots(
        figsize=(9.0, figure_height)
    )

    positions = np.arange(len(plot_data))

    bars = ax.barh(
        positions,
        plot_data["feature_count"],
    )

    ax.set_yticks(positions)
    ax.set_yticklabels(
        plot_data["feature_class"],
        fontsize=9,
    )

    ax.set_xlabel(
        "Number of consensus-signature features"
    )

    ax.set_ylabel("Name-based feature class")

    ax.set_title(
        "Composition of the consensus method-sensitive signature"
    )

    for bar, count, percentage in zip(
        bars,
        plot_data["feature_count"],
        plot_data["percentage_of_signature"],
    ):
        ax.text(
            bar.get_width(),
            bar.get_y() + bar.get_height() / 2,
            f"  {int(count)} ({percentage:.1f}%)",
            va="center",
            fontsize=8,
        )

    fig.tight_layout()

    fig.savefig(
        FIGURE_DIR
        / "feature_class_distribution.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def save_top_features_plot(
    classified: pd.DataFrame,
) -> None:
    top = (
        classified.head(TOP_FEATURES)
        .sort_values(
            "cumulative_robust_effect",
            ascending=True,
        )
    )

    figure_height = max(
        7.0,
        0.34 * len(top),
    )

    fig, ax = plt.subplots(
        figsize=(10.0, figure_height)
    )

    positions = np.arange(len(top))

    bars = ax.barh(
        positions,
        top["cumulative_robust_effect"],
    )

    ax.set_yticks(positions)

    labels = [
        f"{feature} [{feature_class}]"
        for feature, feature_class in zip(
            top["feature"],
            top["feature_class"],
        )
    ]

    ax.set_yticklabels(
        labels,
        fontsize=8,
    )

    ax.set_xlabel(
        "Cumulative absolute effect across robust contrasts"
    )

    ax.set_ylabel("Feature")

    ax.set_title(
        "Top consensus method-sensitive features by name-based class"
    )

    for bar, robust_count in zip(
        bars,
        top["robust_contrasts"],
    ):
        ax.text(
            bar.get_width(),
            bar.get_y() + bar.get_height() / 2,
            f"  n={int(robust_count)}",
            va="center",
            fontsize=8,
        )

    fig.tight_layout()

    fig.savefig(
        FIGURE_DIR
        / "top_consensus_features_by_class.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def main() -> None:
    TABLE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    FIGURE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    data = pd.read_csv(
        INPUT_FILE,
        sep="\t",
    )

    validate_input(data)

    classified = build_classified_table(data)
    class_summary = build_class_summary(classified)
    tier_summary = build_tier_summary(classified)

    classified.to_csv(
        TABLE_DIR / "feature_classification.tsv",
        sep="\t",
        index=False,
    )

    class_summary.to_csv(
        TABLE_DIR / "feature_class_summary.tsv",
        sep="\t",
        index=False,
    )

    tier_summary.to_csv(
        TABLE_DIR / "method_sensitivity_tier_summary.tsv",
        sep="\t",
        index=False,
    )

    save_class_distribution_plot(class_summary)
    save_top_features_plot(classified)

    print(
        f"Consensus-signature features classified: "
        f"{classified.shape[0]}"
    )

    print("\nFeature-class summary:")
    print(
        class_summary.to_string(
            index=False,
            float_format=lambda value: f"{value:.3f}",
        )
    )

    print("\nMethod-sensitivity tiers:")
    print(
        tier_summary.to_string(
            index=False,
            float_format=lambda value: f"{value:.3f}",
        )
    )

    print("\nTop classified consensus features:")
    display_columns = [
        "consensus_rank",
        "feature",
        "feature_class",
        "method_sensitivity_tier",
        "robust_contrasts",
        "cumulative_robust_effect",
        "maximum_absolute_effect",
    ]

    print(
        classified[
            display_columns
        ]
        .head(25)
        .to_string(
            index=False,
            float_format=lambda value: f"{value:.3f}",
        )
    )

    print("\nTables:")
    print(
        "  results/GSE312585/tables/feature_classification.tsv"
    )
    print(
        "  results/GSE312585/tables/feature_class_summary.tsv"
    )
    print(
        "  results/GSE312585/tables/method_sensitivity_tier_summary.tsv"
    )

    print("\nFigures:")
    print(
        "  results/GSE312585/figures/feature_class_distribution.png"
    )
    print(
        "  results/GSE312585/figures/top_consensus_features_by_class.png"
    )

    print(
        "\nInterpretation note: feature classes are based only "
        "on conservative gene-symbol conventions. Genomic "
        "coordinates, biotypes, host-gene relationships, and "
        "mature miRNA products require external annotation."
    )


if __name__ == "__main__":
    main()

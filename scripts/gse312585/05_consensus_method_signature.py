#!/usr/bin/env python3

from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


from matplotlib.colors import TwoSlopeNorm


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET = "GSE312585"


INPUT_FILE = PROJECT_ROOT / "results" / DATASET / "tables" / "all_paired_feature_effects.tsv"

TABLE_DIR = PROJECT_ROOT / "results" / DATASET / "tables"
FIGURE_DIR = PROJECT_ROOT / "results" / DATASET / "figures"

MIN_ABSOLUTE_EFFECT = 1.0
TOP_HEATMAP_FEATURES = 30
TOP_RANK_FEATURES = 25

CONTRAST_ORDER = [
    "UC_minus_WS",
    "Q_minus_WS",
    "M_minus_WS",
    "Q_minus_UC",
    "M_minus_UC",
    "M_minus_Q",
]


def classify_feature(feature: str) -> str:
    """
    Provisional name-based classification only.

    This does not replace an authoritative annotation database.
    """
    feature_upper = feature.upper()

    if (
        feature_upper.endswith("HG")
        or re.search(r"-\d+HG$", feature_upper)
    ):
        return "probable_miRNA_host_gene"

    if feature_upper.startswith("MIR"):
        return "miRNA_associated_locus"

    return "other_feature"


def validate_input(data: pd.DataFrame) -> None:
    required_columns = {
        "feature",
        "contrast",
        "mean_log2_cpm_difference",
        "absolute_mean_difference",
        "direction_consistency",
        "loo_direction_stable",
    }

    missing = required_columns.difference(data.columns)

    if missing:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(sorted(missing))
        )

    duplicated_pairs = data.duplicated(
        subset=["feature", "contrast"]
    )

    if duplicated_pairs.any():
        duplicates = data.loc[
            duplicated_pairs,
            ["feature", "contrast"],
        ]

        raise ValueError(
            "Duplicated feature-contrast rows detected:\n"
            + duplicates.to_string(index=False)
        )


def build_feature_summary(
    data: pd.DataFrame,
) -> pd.DataFrame:
    data = data.copy()

    data["robust_contrast"] = (
        (data["direction_consistency"] == 1.0)
        & data["loo_direction_stable"].astype(bool)
        & (
            data["absolute_mean_difference"]
            >= MIN_ABSOLUTE_EFFECT
        )
    )

    data["robust_absolute_effect"] = np.where(
        data["robust_contrast"],
        data["absolute_mean_difference"],
        0.0,
    )

    data["robust_signed_effect"] = np.where(
        data["robust_contrast"],
        data["mean_log2_cpm_difference"],
        0.0,
    )

    summary = (
        data.groupby("feature")
        .agg(
            contrasts_evaluated=(
                "contrast",
                "nunique",
            ),
            fully_direction_consistent_contrasts=(
                "direction_consistency",
                lambda x: int((x == 1.0).sum()),
            ),
            leave_one_out_stable_contrasts=(
                "loo_direction_stable",
                lambda x: int(
                    pd.Series(x).astype(bool).sum()
                ),
            ),
            robust_contrasts=(
                "robust_contrast",
                "sum",
            ),
            maximum_absolute_effect=(
                "absolute_mean_difference",
                "max",
            ),
            median_absolute_effect=(
                "absolute_mean_difference",
                "median",
            ),
            mean_absolute_effect=(
                "absolute_mean_difference",
                "mean",
            ),
            cumulative_robust_effect=(
                "robust_absolute_effect",
                "sum",
            ),
            mean_direction_consistency=(
                "direction_consistency",
                "mean",
            ),
        )
        .reset_index()
    )

    summary["robust_contrasts"] = (
        summary["robust_contrasts"].astype(int)
    )

    summary["provisional_feature_class"] = (
        summary["feature"].map(classify_feature)
    )

    robust_contrast_names = (
        data.loc[data["robust_contrast"]]
        .groupby("feature")["contrast"]
        .apply(lambda x: ";".join(
            contrast
            for contrast in CONTRAST_ORDER
            if contrast in set(x)
        ))
        .rename("robust_contrast_names")
    )

    summary = summary.merge(
        robust_contrast_names,
        on="feature",
        how="left",
    )

    summary["robust_contrast_names"] = (
        summary["robust_contrast_names"]
        .fillna("")
    )

    summary = summary.sort_values(
        [
            "robust_contrasts",
            "cumulative_robust_effect",
            "maximum_absolute_effect",
            "mean_direction_consistency",
            "feature",
        ],
        ascending=[
            False,
            False,
            False,
            False,
            True,
        ],
    ).reset_index(drop=True)

    summary["consensus_rank"] = (
        np.arange(1, len(summary) + 1)
    )

    return data, summary


def build_effect_matrices(
    data: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    effect_matrix = data.pivot(
        index="feature",
        columns="contrast",
        values="mean_log2_cpm_difference",
    )

    robust_matrix = data.pivot(
        index="feature",
        columns="contrast",
        values="robust_contrast",
    )

    existing_order = [
        contrast
        for contrast in CONTRAST_ORDER
        if contrast in effect_matrix.columns
    ]

    effect_matrix = effect_matrix.loc[
        :,
        existing_order,
    ]

    robust_matrix = robust_matrix.loc[
        :,
        existing_order,
    ].astype(bool)

    return effect_matrix, robust_matrix


def save_consensus_heatmap(
    summary: pd.DataFrame,
    effect_matrix: pd.DataFrame,
    robust_matrix: pd.DataFrame,
) -> None:
    selected_features = (
        summary.loc[
            summary["robust_contrasts"] > 0,
            "feature",
        ]
        .head(TOP_HEATMAP_FEATURES)
        .tolist()
    )

    if not selected_features:
        selected_features = (
            summary["feature"]
            .head(TOP_HEATMAP_FEATURES)
            .tolist()
        )

    plot_matrix = effect_matrix.loc[
        selected_features
    ]

    plot_robust = robust_matrix.loc[
        selected_features
    ]

    maximum = float(
        np.nanmax(
            np.abs(plot_matrix.to_numpy())
        )
    )

    if maximum == 0:
        maximum = 1.0

    figure_height = max(
        7.0,
        0.30 * len(selected_features),
    )

    fig, ax = plt.subplots(
        figsize=(9.5, figure_height)
    )

    norm = TwoSlopeNorm(
        vmin=-maximum,
        vcenter=0,
        vmax=maximum,
    )

    image = ax.imshow(
        plot_matrix.to_numpy(),
        aspect="auto",
        cmap="coolwarm",
        norm=norm,
    )

    ax.set_xticks(
        np.arange(plot_matrix.shape[1])
    )

    ax.set_xticklabels(
        plot_matrix.columns,
        rotation=45,
        ha="right",
        fontsize=9,
    )

    ax.set_yticks(
        np.arange(plot_matrix.shape[0])
    )

    ax.set_yticklabels(
        plot_matrix.index,
        fontsize=8,
    )

    # Mark robust feature-contrast combinations.
    for row_index in range(
        plot_robust.shape[0]
    ):
        for column_index in range(
            plot_robust.shape[1]
        ):
            if plot_robust.iloc[
                row_index,
                column_index,
            ]:
                ax.text(
                    column_index,
                    row_index,
                    "●",
                    ha="center",
                    va="center",
                    fontsize=7,
                )

    colorbar = fig.colorbar(
        image,
        ax=ax,
        fraction=0.035,
        pad=0.03,
    )

    colorbar.set_label(
        "Mean paired log₂ CPM difference"
    )

    ax.set_title(
        "Consensus isolation-method-sensitive signature\n"
        "Dots indicate robust feature–contrast effects"
    )

    ax.set_xlabel("Method contrast")
    ax.set_ylabel("miRNA-associated feature")

    fig.tight_layout()

    fig.savefig(
        FIGURE_DIR
        / "consensus_method_signature_heatmap.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def save_rank_plot(
    summary: pd.DataFrame,
) -> None:
    top = (
        summary.loc[
            summary["robust_contrasts"] > 0
        ]
        .head(TOP_RANK_FEATURES)
        .sort_values(
            [
                "robust_contrasts",
                "cumulative_robust_effect",
            ],
            ascending=True,
        )
    )

    if top.empty:
        top = (
            summary.head(TOP_RANK_FEATURES)
            .sort_values(
                "maximum_absolute_effect",
                ascending=True,
            )
        )

    positions = np.arange(len(top))

    fig, ax = plt.subplots(
        figsize=(9.0, 7.5)
    )

    bars = ax.barh(
        positions,
        top["cumulative_robust_effect"],
    )

    ax.set_yticks(positions)
    ax.set_yticklabels(
        top["feature"],
        fontsize=8,
    )

    ax.set_xlabel(
        "Cumulative absolute effect across robust contrasts"
    )

    ax.set_ylabel(
        "miRNA-associated feature"
    )

    ax.set_title(
        "Ranked consensus method-sensitive features"
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
        / "consensus_signature_rank.png",
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

    # The feature name may have been written as an
    # unnamed index column.
    if "feature" not in data.columns:
        unnamed_columns = [
            column
            for column in data.columns
            if column.startswith("Unnamed:")
        ]

        if len(unnamed_columns) == 1:
            data = data.rename(
                columns={
                    unnamed_columns[0]: "feature"
                }
            )

    validate_input(data)

    data, summary = build_feature_summary(
        data
    )

    effect_matrix, robust_matrix = (
        build_effect_matrices(data)
    )

    summary.to_csv(
        TABLE_DIR
        / "feature_method_sensitivity_summary.tsv",
        sep="\t",
        index=False,
    )

    consensus_signature = summary.loc[
        summary["robust_contrasts"] > 0
    ].copy()

    consensus_signature.to_csv(
        TABLE_DIR
        / "consensus_method_sensitive_signature.tsv",
        sep="\t",
        index=False,
    )

    effect_matrix.to_csv(
        TABLE_DIR
        / "feature_contrast_effect_matrix.tsv",
        sep="\t",
    )

    robust_matrix.astype(int).to_csv(
        TABLE_DIR
        / "feature_contrast_robustness_matrix.tsv",
        sep="\t",
    )

    data.to_csv(
        TABLE_DIR
        / "all_paired_feature_effects_with_robustness.tsv",
        sep="\t",
        index=False,
    )

    save_consensus_heatmap(
        summary=summary,
        effect_matrix=effect_matrix,
        robust_matrix=robust_matrix,
    )

    save_rank_plot(summary)

    robust_distribution = (
        summary["robust_contrasts"]
        .value_counts()
        .sort_index()
        .rename_axis("robust_contrasts")
        .reset_index(name="features")
    )

    robust_distribution.to_csv(
        TABLE_DIR
        / "consensus_robust_contrast_distribution.tsv",
        sep="\t",
        index=False,
    )

    print(
        f"Input rows: {data.shape[0]}"
    )

    print(
        f"Unique features: "
        f"{summary.shape[0]}"
    )

    print(
        "Features with at least one robust contrast: "
        f"{consensus_signature.shape[0]}"
    )

    print(
        "\nRobust-contrast distribution:"
    )

    print(
        robust_distribution.to_string(
            index=False
        )
    )

    print(
        "\nTop consensus features:"
    )

    display_columns = [
        "consensus_rank",
        "feature",
        "provisional_feature_class",
        "robust_contrasts",
        "robust_contrast_names",
        "cumulative_robust_effect",
        "maximum_absolute_effect",
        "mean_direction_consistency",
    ]

    print(
        consensus_signature[
            display_columns
        ]
        .head(25)
        .to_string(
            index=False,
            float_format=lambda value: (
                f"{value:.3f}"
            ),
        )
    )

    print("\nTables:")
    print(
        "  results/GSE312585/tables/"
        "consensus_method_sensitive_signature.tsv"
    )
    print(
        "  results/GSE312585/tables/"
        "feature_method_sensitivity_summary.tsv"
    )
    print(
        "  results/GSE312585/tables/"
        "feature_contrast_effect_matrix.tsv"
    )
    print(
        "  results/GSE312585/tables/"
        "feature_contrast_robustness_matrix.tsv"
    )

    print("\nFigures:")
    print(
        "  results/GSE312585/figures/"
        "consensus_method_signature_heatmap.png"
    )
    print(
        "  results/GSE312585/figures/"
        "consensus_signature_rank.png"
    )


if __name__ == "__main__":
    main()

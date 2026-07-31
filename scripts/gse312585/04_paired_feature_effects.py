#!/usr/bin/env python3

from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET = "GSE312585"


EXPRESSION_FILE = PROJECT_ROOT / "data" / "processed" / DATASET / "mirna_log2_cpm.tsv"
METADATA_FILE = PROJECT_ROOT / "data" / "processed" / DATASET / "sample_metadata.tsv"

TABLE_DIR = PROJECT_ROOT / "results" / DATASET / "tables"
FIGURE_DIR = PROJECT_ROOT / "results" / DATASET / "figures"

METHOD_ORDER = ["WS", "UC", "Q", "M"]

COMPARISONS = [
    ("WS", "UC"),
    ("WS", "Q"),
    ("WS", "M"),
    ("UC", "Q"),
    ("UC", "M"),
    ("Q", "M"),
]


def calculate_contrast(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    reference_method: str,
    comparison_method: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate comparison_method - reference_method separately
    for each donor.
    """
    donor_effects = {}

    for donor, donor_meta in metadata.groupby("donor"):
        donor_meta = donor_meta.set_index("method")

        reference_sample = donor_meta.loc[
            reference_method,
            "sample_id",
        ]

        comparison_sample = donor_meta.loc[
            comparison_method,
            "sample_id",
        ]

        donor_effects[donor] = (
            expression[comparison_sample]
            - expression[reference_sample]
        )

    effects = pd.DataFrame(donor_effects)
    effects.index.name = "feature"

    values = effects.to_numpy()

    mean_effect = effects.mean(axis=1)
    sd_effect = effects.std(axis=1, ddof=1)
    median_effect = effects.median(axis=1)

    positive_count = (effects > 0).sum(axis=1)
    negative_count = (effects < 0).sum(axis=1)
    zero_count = (effects == 0).sum(axis=1)

    direction_consistency = np.maximum(
        positive_count,
        negative_count,
    ) / effects.shape[1]

    consistent_direction = np.where(
        positive_count == effects.shape[1],
        "positive",
        np.where(
            negative_count == effects.shape[1],
            "negative",
            "mixed",
        ),
    )

    # Leave-one-donor-out mean effects.
    loo_means = []

    donor_columns = list(effects.columns)

    for omitted_donor in donor_columns:
        retained = [
            donor
            for donor in donor_columns
            if donor != omitted_donor
        ]

        loo_means.append(
            effects[retained].mean(axis=1).rename(
                f"mean_without_{omitted_donor}"
            )
        )

    loo_table = pd.concat(loo_means, axis=1)

    loo_min = loo_table.min(axis=1)
    loo_max = loo_table.max(axis=1)

    loo_direction_stable = (
        ((loo_table > 0).all(axis=1))
        | ((loo_table < 0).all(axis=1))
    )

    # Ranking prioritises large, directionally consistent effects.
    stability_score = (
        mean_effect.abs()
        * direction_consistency
        * (1 + loo_direction_stable.astype(float))
    )

    summary = effects.copy()

    summary["reference_method"] = reference_method
    summary["comparison_method"] = comparison_method
    summary["contrast"] = (
        f"{comparison_method}_minus_{reference_method}"
    )

    summary["mean_log2_cpm_difference"] = mean_effect
    summary["median_log2_cpm_difference"] = median_effect
    summary["sd_log2_cpm_difference"] = sd_effect
    summary["absolute_mean_difference"] = mean_effect.abs()

    summary["positive_donors"] = positive_count
    summary["negative_donors"] = negative_count
    summary["zero_donors"] = zero_count

    summary["direction_consistency"] = direction_consistency
    summary["consistent_direction"] = consistent_direction

    summary["loo_minimum_mean"] = loo_min
    summary["loo_maximum_mean"] = loo_max
    summary["loo_direction_stable"] = loo_direction_stable

    summary["stability_score"] = stability_score

    summary = summary.sort_values(
        [
            "loo_direction_stable",
            "direction_consistency",
            "absolute_mean_difference",
        ],
        ascending=[False, False, False],
    )

    return effects, summary


def save_top_effect_plot(
    summary: pd.DataFrame,
    contrast_name: str,
    output_path: Path,
    top_n: int = 20,
) -> None:
    top = (
        summary.loc[
            summary["direction_consistency"] == 1.0
        ]
        .nlargest(
            top_n,
            "absolute_mean_difference",
        )
        .sort_values(
            "mean_log2_cpm_difference"
        )
    )

    if top.empty:
        top = (
            summary.nlargest(
                top_n,
                "absolute_mean_difference",
            )
            .sort_values(
                "mean_log2_cpm_difference"
            )
        )

    fig, ax = plt.subplots(figsize=(8.5, 7.0))

    positions = np.arange(len(top))

    ax.barh(
        positions,
        top["mean_log2_cpm_difference"],
    )

    ax.set_yticks(positions)
    ax.set_yticklabels(
        top.index,
        fontsize=8,
    )

    ax.axvline(0, linewidth=0.8)

    ax.set_xlabel(
        "Mean paired log₂ CPM difference"
    )
    ax.set_ylabel("miRNA-associated feature")
    ax.set_title(
        f"Consistent feature effects: {contrast_name}"
    )

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def save_effect_heatmap(
    combined_summary: pd.DataFrame,
    output_path: Path,
    features_per_contrast: int = 8,
) -> None:
    selected_features = []

    for contrast, contrast_table in combined_summary.groupby(
        "contrast",
        sort=False,
    ):
        selected = (
            contrast_table.loc[
                contrast_table["direction_consistency"] == 1.0
            ]
            .nlargest(
                features_per_contrast,
                "absolute_mean_difference",
            )
            .index
            .tolist()
        )

        selected_features.extend(selected)

    selected_features = list(
        dict.fromkeys(selected_features)
    )

    effect_matrix = (
        combined_summary.loc[
            combined_summary.index.isin(selected_features)
        ]
        .reset_index()
        .pivot_table(
            index="feature",
            columns="contrast",
            values="mean_log2_cpm_difference",
            aggfunc="first",
        )
    )

    desired_columns = [
        f"{comparison}_minus_{reference}"
        for reference, comparison in COMPARISONS
    ]

    desired_columns = [
        column
        for column in desired_columns
        if column in effect_matrix.columns
    ]

    effect_matrix = effect_matrix.loc[
        :,
        desired_columns,
    ]

    effect_matrix.to_csv(
        TABLE_DIR / "selected_feature_effect_matrix.tsv",
        sep="\t",
    )

    fig_height = max(
        6.0,
        0.28 * effect_matrix.shape[0],
    )

    fig, ax = plt.subplots(
        figsize=(9.0, fig_height)
    )

    maximum = np.nanmax(
        np.abs(effect_matrix.to_numpy())
    )

    image = ax.imshow(
        effect_matrix.to_numpy(),
        aspect="auto",
        vmin=-maximum,
        vmax=maximum,
    )

    ax.set_xticks(
        np.arange(effect_matrix.shape[1])
    )

    ax.set_xticklabels(
        effect_matrix.columns,
        rotation=45,
        ha="right",
        fontsize=8,
    )

    ax.set_yticks(
        np.arange(effect_matrix.shape[0])
    )

    ax.set_yticklabels(
        effect_matrix.index,
        fontsize=7,
    )

    colorbar = fig.colorbar(
        image,
        ax=ax,
    )

    colorbar.set_label(
        "Mean paired log₂ CPM difference"
    )

    ax.set_title(
        "Isolation-method effect profiles"
    )

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    expression = pd.read_csv(
        EXPRESSION_FILE,
        sep="\t",
        index_col=0,
    )

    metadata = pd.read_csv(
        METADATA_FILE,
        sep="\t",
        dtype={"donor": str},
    )

    sample_ids = metadata["sample_id"].tolist()
    expression = expression.loc[:, sample_ids]

    all_summaries = []
    contrast_overview = []

    for reference_method, comparison_method in COMPARISONS:
        contrast_name = (
            f"{comparison_method}_minus_{reference_method}"
        )

        effects, summary = calculate_contrast(
            expression=expression,
            metadata=metadata,
            reference_method=reference_method,
            comparison_method=comparison_method,
        )

        effects.to_csv(
            TABLE_DIR
            / f"paired_effects_{contrast_name}.tsv",
            sep="\t",
        )

        summary.to_csv(
            TABLE_DIR
            / f"paired_effect_summary_{contrast_name}.tsv",
            sep="\t",
        )

        all_summaries.append(summary)

        n_fully_consistent = int(
            (summary["direction_consistency"] == 1.0).sum()
        )

        n_loo_stable = int(
            summary["loo_direction_stable"].sum()
        )

        n_large_consistent = int(
            (
                (summary["direction_consistency"] == 1.0)
                & (
                    summary[
                        "absolute_mean_difference"
                    ] >= 1.0
                )
            ).sum()
        )

        contrast_overview.append(
            {
                "contrast": contrast_name,
                "reference_method": reference_method,
                "comparison_method": comparison_method,
                "features": summary.shape[0],
                "fully_direction_consistent": n_fully_consistent,
                "leave_one_out_direction_stable": n_loo_stable,
                "consistent_with_abs_mean_ge_1": (
                    n_large_consistent
                ),
                "median_absolute_effect": summary[
                    "absolute_mean_difference"
                ].median(),
                "maximum_absolute_effect": summary[
                    "absolute_mean_difference"
                ].max(),
            }
        )

        save_top_effect_plot(
            summary=summary,
            contrast_name=contrast_name,
            output_path=(
                FIGURE_DIR
                / f"top_paired_effects_{contrast_name}.png"
            ),
        )

    combined = pd.concat(
        all_summaries,
        axis=0,
    )

    combined.index.name = "feature"

    combined.to_csv(
        TABLE_DIR / "all_paired_feature_effects.tsv",
        sep="\t",
    )

    overview = pd.DataFrame(
        contrast_overview
    ).sort_values(
        "median_absolute_effect",
        ascending=False,
    )

    overview.to_csv(
        TABLE_DIR / "paired_effect_contrast_overview.tsv",
        sep="\t",
        index=False,
    )

    save_effect_heatmap(
        combined_summary=combined,
        output_path=(
            FIGURE_DIR
            / "paired_method_effect_heatmap.png"
        ),
    )

    print(
        f"Expression matrix: {expression.shape[0]} features "
        f"x {expression.shape[1]} samples"
    )

    print("\nContrast overview:")
    print(
        overview.to_string(
            index=False,
            float_format=lambda value: f"{value:.3f}",
        )
    )

    for contrast in overview["contrast"]:
        subset = combined.loc[
            combined["contrast"] == contrast
        ]

        top = subset.loc[
            subset["direction_consistency"] == 1.0
        ].nlargest(
            5,
            "absolute_mean_difference",
        )

        print(f"\nTop consistent features: {contrast}")

        print(
            top[
                [
                    "mean_log2_cpm_difference",
                    "sd_log2_cpm_difference",
                    "direction_consistency",
                    "loo_direction_stable",
                ]
            ].to_string(
                float_format=lambda value: f"{value:.3f}",
            )
        )

    print("\nOutputs:")
    print(
        f"  {TABLE_DIR / 'paired_effect_contrast_overview.tsv'}"
    )
    print(
        f"  {TABLE_DIR / 'all_paired_feature_effects.tsv'}"
    )
    print(
        f"  {FIGURE_DIR / 'paired_method_effect_heatmap.png'}"
    )


if __name__ == "__main__":
    main()

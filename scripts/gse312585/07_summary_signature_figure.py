#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET = "GSE312585"



# Configuration
# 

INPUT_FILE = PROJECT_ROOT / "results" / DATASET / "tables" / "feature_classification.tsv"

FIGURE_DIR = PROJECT_ROOT / "results" / DATASET / "figures"
TABLE_DIR = PROJECT_ROOT / "results" / DATASET / "tables"
TEXT_DIR = PROJECT_ROOT / "results" / DATASET / "text"

OUTPUT_BASENAME = "summary_signature_scatter"

DPI = 300

FIGURE_WIDTH = 10.5
FIGURE_HEIGHT = 7.2

MIN_POINT_SIZE = 55
MAX_POINT_SIZE = 360

LABEL_FEATURES = [
    "MIR1249",
    "MIR32",
    "MIR136",
    "MIR885",
    "MIR193B",
    "MIR29B1",
]

LABEL_OFFSETS = {
    "MIR1249": (-0.45, -0.08),
    "MIR32": (-0.45, -0.08),
    "MIR136": (-0.45, -0.08),
    "MIR885": (-0.45, -0.08),
    "MIR193B": (0.18, -0.08),
    "MIR29B1": (0.18, 0.08),
}

CLASS_MARKERS = {
    "miRNA_associated_locus": "o",
    "probable_miRNA_host_gene": "s",
}

CLASS_LABELS = {
    "miRNA_associated_locus":
        "miRNA-associated locus",
    "probable_miRNA_host_gene":
        "probable miRNA host gene",
}




# Validation
# 

def validate_input(
    data: pd.DataFrame,
) -> None:
    """
    Validate required input columns and values.
    """
    required_columns = {
        "feature",
        "feature_class",
        "consensus_rank",
        "robust_contrasts",
        "cumulative_robust_effect",
        "maximum_absolute_effect",
    }

    missing = required_columns.difference(
        data.columns
    )

    if missing:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(sorted(missing))
        )

    if data.empty:
        raise ValueError(
            "The input table contains no rows."
        )

    if data["feature"].duplicated().any():
        duplicated = (
            data.loc[
                data["feature"].duplicated(
                    keep=False
                ),
                "feature",
            ]
            .astype(str)
            .tolist()
        )

        raise ValueError(
            "Duplicated feature names detected: "
            + ", ".join(
                sorted(set(duplicated))
            )
        )

    numeric_columns = [
        "consensus_rank",
        "robust_contrasts",
        "cumulative_robust_effect",
        "maximum_absolute_effect",
    ]

    for column in numeric_columns:
        converted = pd.to_numeric(
            data[column],
            errors="coerce",
        )

        invalid = (
            converted.isna()
            | ~np.isfinite(converted)
        )

        if invalid.any():
            invalid_features = (
                data.loc[
                    invalid,
                    "feature",
                ]
                .astype(str)
                .tolist()
            )

            raise ValueError(
                f"Invalid values in '{column}' for: "
                + ", ".join(
                    invalid_features
                )
            )



# Data preparation
# 

def prepare_data(
    data: pd.DataFrame,
) -> pd.DataFrame:
    """
    Prepare plotting columns.
    """
    prepared = data.copy()

    numeric_columns = [
        "consensus_rank",
        "robust_contrasts",
        "cumulative_robust_effect",
        "maximum_absolute_effect",
    ]

    for column in numeric_columns:
        prepared[column] = pd.to_numeric(
            prepared[column]
        )

    prepared["feature"] = (
        prepared["feature"]
        .astype(str)
        .str.strip()
    )

    prepared["feature_class"] = (
        prepared["feature_class"]
        .fillna("unclassified")
        .astype(str)
        .str.strip()
    )

    prepared["label_feature"] = (
        prepared["feature"].isin(
            LABEL_FEATURES
        )
    )

    prepared = prepared.sort_values(
        [
            "robust_contrasts",
            "cumulative_robust_effect",
            "consensus_rank",
        ],
        ascending=[
            True,
            True,
            False,
        ],
    ).reset_index(drop=True)

    return prepared



# Point-size scaling
# 

def scale_point_sizes(
    values: pd.Series,
) -> np.ndarray:
    """
    Scale maximum absolute effect values to point areas.
    """
    numeric = pd.to_numeric(
        values,
        errors="coerce",
    ).to_numpy(dtype=float)

    minimum = float(
        np.min(numeric)
    )

    maximum = float(
        np.max(numeric)
    )

    if np.isclose(
        minimum,
        maximum,
    ):
        return np.full(
            numeric.shape,
            (
                MIN_POINT_SIZE
                + MAX_POINT_SIZE
            )
            / 2,
        )

    scaled = (
        numeric - minimum
    ) / (
        maximum - minimum
    )

    return (
        MIN_POINT_SIZE
        + scaled
        * (
            MAX_POINT_SIZE
            - MIN_POINT_SIZE
        )
    )


def scale_reference_size(
    value: float,
    minimum: float,
    maximum: float,
) -> float:
    """
    Scale one representative effect value to point area.
    """
    if np.isclose(
        minimum,
        maximum,
    ):
        return (
            MIN_POINT_SIZE
            + MAX_POINT_SIZE
        ) / 2

    scaled = (
        value - minimum
    ) / (
        maximum - minimum
    )

    scaled = float(
        np.clip(
            scaled,
            0.0,
            1.0,
        )
    )

    return (
        MIN_POINT_SIZE
        + scaled
        * (
            MAX_POINT_SIZE
            - MIN_POINT_SIZE
        )
    )



# Coordinates table
# 

def build_coordinates_table(
    data: pd.DataFrame,
    point_sizes: np.ndarray,
) -> pd.DataFrame:
    """
    Save plotting coordinates and aesthetics.
    """
    output = data.copy()

    output["x_coordinate"] = (
        output[
            "cumulative_robust_effect"
        ]
    )

    output["y_coordinate"] = (
        output[
            "robust_contrasts"
        ]
    )

    output["point_size"] = (
        point_sizes
    )

    preferred_columns = [
        "feature",
        "consensus_rank",
        "feature_class",
        "method_sensitivity_tier",
        "robust_contrasts",
        "cumulative_robust_effect",
        "maximum_absolute_effect",
        "x_coordinate",
        "y_coordinate",
        "point_size",
        "label_feature",
        "robust_contrast_names",
    ]

    existing_columns = [
        column
        for column in preferred_columns
        if column in output.columns
    ]

    remaining_columns = [
        column
        for column in output.columns
        if column not in existing_columns
    ]

    return (
        output[
            existing_columns
            + remaining_columns
        ]
        .sort_values(
            "consensus_rank"
        )
        .reset_index(drop=True)
    )



# Figure
# 

def create_scatter_plot(
    data: pd.DataFrame,
    point_sizes: np.ndarray,
) -> None:
    """
    Create a clean summary scatter plot.
    """
    fig, ax = plt.subplots(
        figsize=(
            FIGURE_WIDTH,
            FIGURE_HEIGHT,
        )
    )

    feature_classes = sorted(
        data[
            "feature_class"
        ].unique()
    )

    for feature_class in feature_classes:
        mask = (
            data[
                "feature_class"
            ]
            == feature_class
        )

        subset = data.loc[
            mask
        ]

        subset_sizes = (
            point_sizes[
                mask.to_numpy()
            ]
        )

        marker = CLASS_MARKERS.get(
            feature_class,
            "D",
        )

        label = CLASS_LABELS.get(
            feature_class,
            feature_class.replace(
                "_",
                " ",
            ),
        )

        ax.scatter(
            subset[
                "cumulative_robust_effect"
            ],
            subset[
                "robust_contrasts"
            ],
            s=subset_sizes,
            marker=marker,
            alpha=0.78,
            edgecolors="black",
            linewidths=0.55,
            label=label,
            zorder=3,
        )

    ax.set_title(
        "Consensus method-sensitive salivary RNA signature",
        fontsize=15,
        pad=14,
    )

    ax.set_xlabel(
        "Cumulative absolute effect across robust paired contrasts",
        fontsize=11,
        labelpad=8,
    )

    ax.set_ylabel(
        "Number of robust paired contrasts",
        fontsize=11,
        labelpad=8,
    )

    maximum_robust = int(
        data[
            "robust_contrasts"
        ].max()
    )

    ax.set_yticks(
        np.arange(
            1,
            maximum_robust + 1,
        )
    )

    ax.set_ylim(
        0.65,
        maximum_robust + 0.45,
    )

    x_maximum = float(
        data[
            "cumulative_robust_effect"
        ].max()
    )

    ax.set_xlim(
        0,
        x_maximum + 1.8,
    )

    ax.grid(
        axis="y",
        linestyle="--",
        linewidth=0.55,
        alpha=0.35,
        zorder=0,
    )

    ax.grid(
        axis="x",
        linestyle=":",
        linewidth=0.5,
        alpha=0.22,
        zorder=0,
    )

    for spine in [
        "top",
        "right",
    ]:
        ax.spines[
            spine
        ].set_visible(
            False
        )

    # Stable manual labels for six key features.
    labelled = data.loc[
        data[
            "label_feature"
        ]
    ].copy()

    for _, row in labelled.iterrows():
        feature = str(
            row["feature"]
        )

        dx, dy = LABEL_OFFSETS.get(
            feature,
            (0.15, 0.08),
        )

        horizontal_alignment = (
            "right"
            if dx < 0
            else "left"
        )

        ax.text(
            row[
                "cumulative_robust_effect"
            ]
            + dx,
            row[
                "robust_contrasts"
            ]
            + dy,
            feature,
            fontsize=8.5,
            fontweight=(
                "bold"
                if row[
                    "consensus_rank"
                ]
                <= 3
                else "normal"
            ),
            ha=horizontal_alignment,
            va="center",
            zorder=5,
        )

    # Compact feature-class legend.
    class_legend = ax.legend(
        title="Feature class",
        loc="upper left",
        frameon=True,
        fontsize=8,
        title_fontsize=9,
        borderpad=0.7,
        labelspacing=0.5,
        handletextpad=0.6,
    )

    ax.add_artist(
        class_legend
    )

    # Compact point-size legend.
    effects = (
        data[
            "maximum_absolute_effect"
        ]
        .astype(float)
        .to_numpy()
    )

    effect_minimum = float(
        np.min(effects)
    )

    effect_maximum = float(
        np.max(effects)
    )

    representative_effects = [
        2.0,
        4.0,
        6.0,
    ]

    size_handles = []

    for value in representative_effects:
        marker_size = scale_reference_size(
            value,
            effect_minimum,
            effect_maximum,
        )

        handle = ax.scatter(
            [],
            [],
            s=marker_size,
            facecolors="none",
            edgecolors="black",
            linewidths=0.8,
            label=f"{value:.0f}",
        )

        size_handles.append(
            handle
        )

    size_legend = ax.legend(
        handles=size_handles,
        title="Maximum effect",
        loc="lower right",
        frameon=True,
        fontsize=8,
        title_fontsize=9,
        borderpad=0.7,
        labelspacing=0.7,
    )

    ax.add_artist(
        size_legend
    )

    fig.tight_layout()

    png_path = (
        FIGURE_DIR
        / f"{OUTPUT_BASENAME}.png"
    )

    pdf_path = (
        FIGURE_DIR
        / f"{OUTPUT_BASENAME}.pdf"
    )

    fig.savefig(
        png_path,
        dpi=DPI,
        bbox_inches="tight",
    )

    fig.savefig(
        pdf_path,
        bbox_inches="tight",
    )

    plt.close(
        fig
    )



# Figure legend
# 

def write_figure_legend(
    data: pd.DataFrame,
) -> None:
    """
    Write a Markdown figure legend.
    """
    top_feature = (
        data.sort_values(
            "consensus_rank"
        )
        .iloc[0]
    )

    labelled_features = (
        data.loc[
            data[
                "label_feature"
            ],
            "feature",
        ]
        .astype(str)
        .tolist()
    )

    legend = (
        "**Figure 7. Consensus method-sensitive salivary RNA "
        "signature across extracellular-vesicle isolation "
        "procedures.**\n\n"
        f"Each point represents one of {len(data)} consensus "
        "features supported by at least one robust paired contrast. "
        "The x-axis shows the cumulative absolute effect across "
        "robust contrasts, while the y-axis shows the number of "
        "robust paired contrasts. Point size represents the maximum "
        "absolute paired effect. Marker shape denotes the "
        "conservative name-based feature class. Labels identify "
        "selected high-ranking or otherwise informative features: "
        f"{', '.join(labelled_features)}. "
        f"{top_feature['feature']} was the highest-ranked feature, "
        f"supported by {int(top_feature['robust_contrasts'])} "
        "robust contrasts and a cumulative robust effect of "
        f"{float(top_feature['cumulative_robust_effect']):.2f}. "
        "The three-donor analysis is exploratory and describes "
        "isolation-method sensitivity rather than formal "
        "differential expression.\n"
    )

    output_path = (
        TEXT_DIR
        / "Figure7_legend.md"
    )

    output_path.write_text(
        legend,
        encoding="utf-8",
    )



# Console summary
# 

def print_summary(
    data: pd.DataFrame,
) -> None:
    """
    Print a concise execution summary.
    """
    labelled = (
        data.loc[
            data[
                "label_feature"
            ]
        ]
        .sort_values(
            "consensus_rank"
        )
    )

    print(
        "Summary-signature features plotted: "
        f"{data.shape[0]}"
    )

    print(
        "Features labelled: "
        f"{labelled.shape[0]}"
    )

    print(
        "Maximum number of robust contrasts: "
        f"{int(data['robust_contrasts'].max())}"
    )

    print(
        "\nLabelled features:"
    )

    columns = [
        "consensus_rank",
        "feature",
        "robust_contrasts",
        "cumulative_robust_effect",
        "maximum_absolute_effect",
    ]

    print(
        labelled[
            columns
        ].to_string(
            index=False,
            float_format=lambda value: f"{value:.3f}",
        )
    )

    print(
        "\nOutputs:"
    )

    print(
        "  results/GSE312585/figures/"
        "summary_signature_scatter.png"
    )

    print(
        "  results/GSE312585/figures/"
        "summary_signature_scatter.pdf"
    )

    print(
        "  results/GSE312585/tables/"
        "summary_signature_coordinates.tsv"
    )

    print(
        "  results/GSE312585/text/"
        "Figure7_legend.md"
    )



# Main
# 

def main() -> None:
    """
    Run the corrected summary-figure workflow.
    """
    FIGURE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    TABLE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    TEXT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    data = pd.read_csv(
        INPUT_FILE,
        sep="\t",
    )

    validate_input(
        data
    )

    prepared = prepare_data(
        data
    )

    point_sizes = scale_point_sizes(
        prepared[
            "maximum_absolute_effect"
        ]
    )

    coordinates = build_coordinates_table(
        prepared,
        point_sizes,
    )

    coordinates.to_csv(
        TABLE_DIR
        / "summary_signature_coordinates.tsv",
        sep="\t",
        index=False,
    )

    create_scatter_plot(
        prepared,
        point_sizes,
    )

    write_figure_legend(
        prepared
    )

    print_summary(
        prepared
    )


if __name__ == "__main__":
    main()

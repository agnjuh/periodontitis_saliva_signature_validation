#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET = "GSE312585"


EXPRESSION_FILE = PROJECT_ROOT / "data" / "processed" / DATASET / "mirna_log2_cpm.tsv"
METADATA_FILE = PROJECT_ROOT / "data" / "processed" / DATASET / "sample_metadata.tsv"
SAMPLE_QC_FILE = PROJECT_ROOT / "results" / DATASET / "tables" / "mirna_sample_qc.tsv"

FIGURE_DIR = PROJECT_ROOT / "results" / DATASET / "figures"
TABLE_DIR = PROJECT_ROOT / "results" / DATASET / "tables"

METHOD_ORDER = ["WS", "UC", "Q", "M"]


def validate_inputs(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_qc: pd.DataFrame,
) -> list[str]:
    required_metadata = {
        "sample_id",
        "gsm_id",
        "donor",
        "method",
        "sample_type",
    }

    missing_metadata = required_metadata.difference(metadata.columns)

    if missing_metadata:
        raise ValueError(
            f"Metadata columns missing: {sorted(missing_metadata)}"
        )

    required_qc = {
        "sample_id",
        "mirna_total_counts_after_filtering",
        "mirna_detected_after_filtering",
    }

    missing_qc = required_qc.difference(sample_qc.columns)

    if missing_qc:
        raise ValueError(
            f"QC columns missing: {sorted(missing_qc)}"
        )

    sample_ids = metadata["sample_id"].tolist()

    missing_expression = [
        sample for sample in sample_ids
        if sample not in expression.columns
    ]

    if missing_expression:
        raise ValueError(
            "Samples missing from expression matrix: "
            f"{missing_expression}"
        )

    duplicated_samples = metadata.loc[
        metadata["sample_id"].duplicated(),
        "sample_id",
    ].tolist()

    if duplicated_samples:
        raise ValueError(
            f"Duplicated metadata sample IDs: {duplicated_samples}"
        )

    if expression.isna().any().any():
        raise ValueError(
            "Expression matrix contains missing values."
        )

    return sample_ids


def assign_categories(values: pd.Series) -> tuple[dict, dict]:
    unique_values = list(dict.fromkeys(values.astype(str)))

    default_colors = plt.rcParams[
        "axes.prop_cycle"
    ].by_key()["color"]

    markers = ["o", "s", "^", "D", "P", "X", "v", "<", ">"]

    color_map = {
        value: default_colors[index % len(default_colors)]
        for index, value in enumerate(unique_values)
    }

    marker_map = {
        value: markers[index % len(markers)]
        for index, value in enumerate(unique_values)
    }

    return color_map, marker_map


def save_pca_by_method(
    pca_table: pd.DataFrame,
    explained_variance: np.ndarray,
) -> None:
    method_colors, _ = assign_categories(
        pd.Series(METHOD_ORDER)
    )
    _, donor_markers = assign_categories(
        pca_table["donor"]
    )

    fig, ax = plt.subplots(figsize=(8.2, 6.5))

    for row in pca_table.itertuples(index=False):
        ax.scatter(
            row.PC1,
            row.PC2,
            color=method_colors[row.method],
            marker=donor_markers[str(row.donor)],
            s=95,
            edgecolor="black",
            linewidth=0.6,
        )

        ax.annotate(
            row.sample_id,
            (row.PC1, row.PC2),
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=8,
        )

    method_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=method_colors[method],
            markeredgecolor="black",
            label=method,
            markersize=8,
        )
        for method in METHOD_ORDER
        if method in set(pca_table["method"])
    ]

    donor_handles = [
        plt.Line2D(
            [0],
            [0],
            marker=donor_markers[str(donor)],
            linestyle="",
            markerfacecolor="white",
            markeredgecolor="black",
            label=str(donor),
            markersize=8,
        )
        for donor in pca_table["donor"].astype(str).unique()
    ]

    method_legend = ax.legend(
        handles=method_handles,
        title="Isolation method",
        loc="best",
    )

    ax.add_artist(method_legend)

    ax.legend(
        handles=donor_handles,
        title="Donor",
        loc="lower right",
    )

    ax.axhline(0, linewidth=0.7)
    ax.axvline(0, linewidth=0.7)

    ax.set_xlabel(
        f"PC1 ({explained_variance[0]:.1f}% variance)"
    )
    ax.set_ylabel(
        f"PC2 ({explained_variance[1]:.1f}% variance)"
    )
    ax.set_title(
        "Salivary EV miRNA-associated profiles"
    )

    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "pca_method_and_donor.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def save_pca_by_donor(
    pca_table: pd.DataFrame,
    explained_variance: np.ndarray,
) -> None:
    donor_colors, _ = assign_categories(
        pca_table["donor"].astype(str)
    )

    method_markers = {
        "WS": "o",
        "UC": "s",
        "Q": "^",
        "M": "D",
    }

    fig, ax = plt.subplots(figsize=(8.2, 6.5))

    for row in pca_table.itertuples(index=False):
        ax.scatter(
            row.PC1,
            row.PC2,
            color=donor_colors[str(row.donor)],
            marker=method_markers[row.method],
            s=95,
            edgecolor="black",
            linewidth=0.6,
        )

        ax.annotate(
            row.sample_id,
            (row.PC1, row.PC2),
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=8,
        )

    ax.axhline(0, linewidth=0.7)
    ax.axvline(0, linewidth=0.7)

    ax.set_xlabel(
        f"PC1 ({explained_variance[0]:.1f}% variance)"
    )
    ax.set_ylabel(
        f"PC2 ({explained_variance[1]:.1f}% variance)"
    )
    ax.set_title(
        "PCA coloured by donor"
    )

    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "pca_by_donor.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def save_correlation_heatmap(
    correlation: pd.DataFrame,
    ordered_samples: list[str],
) -> None:
    ordered = correlation.loc[
        ordered_samples,
        ordered_samples,
    ]

    fig, ax = plt.subplots(figsize=(9.2, 8.0))

    image = ax.imshow(
        ordered.values,
        vmin=0,
        vmax=1,
        aspect="auto",
    )

    positions = np.arange(len(ordered_samples))

    ax.set_xticks(positions)
    ax.set_yticks(positions)

    ax.set_xticklabels(
        ordered_samples,
        rotation=90,
        fontsize=8,
    )
    ax.set_yticklabels(
        ordered_samples,
        fontsize=8,
    )

    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("Pearson correlation")

    ax.set_title(
        "Sample correlation: miRNA-associated log₂ CPM"
    )

    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "sample_correlation_heatmap.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def save_dendrogram(
    clustering: np.ndarray,
    sample_ids: list[str],
) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.8))

    dendrogram(
        clustering,
        labels=sample_ids,
        leaf_rotation=90,
        leaf_font_size=9,
        ax=ax,
    )

    ax.set_ylabel("Correlation distance")
    ax.set_title(
        "Hierarchical clustering of samples"
    )

    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "sample_hierarchical_clustering.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def save_library_size_plot(
    sample_qc: pd.DataFrame,
) -> None:
    plot_data = sample_qc.copy()

    method_rank = {
        method: index
        for index, method in enumerate(METHOD_ORDER)
    }

    plot_data["method_rank"] = plot_data[
        "method"
    ].map(method_rank)

    plot_data = plot_data.sort_values(
        ["donor", "method_rank"]
    )

    colors, _ = assign_categories(
        pd.Series(METHOD_ORDER)
    )

    bar_colors = [
        colors[method]
        for method in plot_data["method"]
    ]

    fig, ax = plt.subplots(figsize=(10.0, 5.8))

    ax.bar(
        plot_data["sample_id"],
        plot_data[
            "mirna_total_counts_after_filtering"
        ],
        color=bar_colors,
    )

    ax.set_ylabel(
        "Counts assigned to retained miRNA-associated loci"
    )
    ax.set_xlabel("Sample")
    ax.set_title(
        "miRNA-associated library size after filtering"
    )

    ax.tick_params(
        axis="x",
        rotation=60,
        labelsize=8,
    )

    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "mirna_library_size_qc.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def save_detected_feature_plot(
    sample_qc: pd.DataFrame,
) -> None:
    plot_data = sample_qc.copy()

    method_rank = {
        method: index
        for index, method in enumerate(METHOD_ORDER)
    }

    plot_data["method_rank"] = plot_data[
        "method"
    ].map(method_rank)

    plot_data = plot_data.sort_values(
        ["donor", "method_rank"]
    )

    colors, _ = assign_categories(
        pd.Series(METHOD_ORDER)
    )

    bar_colors = [
        colors[method]
        for method in plot_data["method"]
    ]

    fig, ax = plt.subplots(figsize=(10.0, 5.8))

    ax.bar(
        plot_data["sample_id"],
        plot_data[
            "mirna_detected_after_filtering"
        ],
        color=bar_colors,
    )

    ax.set_ylim(
        0,
        max(
            plot_data[
                "mirna_detected_after_filtering"
            ]
        ) * 1.08,
    )

    ax.set_ylabel(
        "Detected retained miRNA-associated loci"
    )
    ax.set_xlabel("Sample")
    ax.set_title(
        "Feature detection after abundance filtering"
    )

    ax.tick_params(
        axis="x",
        rotation=60,
        labelsize=8,
    )

    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "mirna_detected_features_qc.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

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

    sample_qc = pd.read_csv(
        SAMPLE_QC_FILE,
        sep="\t",
        dtype={"donor": str},
    )

    sample_ids = validate_inputs(
        expression,
        metadata,
        sample_qc,
    )

    expression = expression.loc[:, sample_ids]

    # Samples are observations and miRNA-associated loci are variables.
    sample_matrix = expression.T

    # PCA is performed on centred features without variance scaling.
    scaled_matrix = StandardScaler(
        with_mean=True,
        with_std=False,
    ).fit_transform(sample_matrix)

    pca = PCA()
    pca_scores = pca.fit_transform(scaled_matrix)

    explained_variance = (
        pca.explained_variance_ratio_ * 100
    )

    pca_table = metadata.copy()
    pca_table["PC1"] = pca_scores[:, 0]
    pca_table["PC2"] = pca_scores[:, 1]
    pca_table["PC3"] = pca_scores[:, 2]

    pca_table.to_csv(
        TABLE_DIR / "pca_sample_coordinates.tsv",
        sep="\t",
        index=False,
    )

    variance_table = pd.DataFrame(
        {
            "component": [
                f"PC{index + 1}"
                for index in range(len(explained_variance))
            ],
            "explained_variance_percent": explained_variance,
            "cumulative_variance_percent": np.cumsum(
                explained_variance
            ),
        }
    )

    variance_table.to_csv(
        TABLE_DIR / "pca_explained_variance.tsv",
        sep="\t",
        index=False,
    )

    correlation = sample_matrix.T.corr(
        method="pearson"
    )

    correlation.to_csv(
        TABLE_DIR / "sample_pearson_correlation.tsv",
        sep="\t",
    )

    distance = 1 - correlation

    # Make a writable NumPy copy for NumPy 2.x compatibility.
    distance_matrix = distance.to_numpy(copy=True)

    np.fill_diagonal(
        distance_matrix,
        0,
    )

    condensed_distance = squareform(
        distance_matrix,
        checks=False,
    )

    clustering = linkage(
        condensed_distance,
        method="complete",
    )

    leaf_indices = leaves_list(clustering)

    ordered_samples = [
        sample_ids[index]
        for index in leaf_indices
    ]

    pd.DataFrame(
        {
            "cluster_order": np.arange(
                1,
                len(ordered_samples) + 1,
            ),
            "sample_id": ordered_samples,
        }
    ).to_csv(
        TABLE_DIR / "sample_clustering_order.tsv",
        sep="\t",
        index=False,
    )

    save_pca_by_method(
        pca_table,
        explained_variance,
    )

    save_pca_by_donor(
        pca_table,
        explained_variance,
    )

    save_correlation_heatmap(
        correlation,
        ordered_samples,
    )

    save_dendrogram(
        clustering,
        sample_ids,
    )

    save_library_size_plot(sample_qc)
    save_detected_feature_plot(sample_qc)

    print(
        f"Input matrix: {expression.shape[0]} features "
        f"x {expression.shape[1]} samples"
    )

    print("\nExplained variance:")
    print(
        variance_table.head(5).to_string(
            index=False
        )
    )

    print("\nHierarchical clustering order:")
    print(" -> ".join(ordered_samples))

    print("\nPCA coordinates:")
    print(
        pca_table[
            [
                "sample_id",
                "donor",
                "method",
                "PC1",
                "PC2",
            ]
        ].to_string(
            index=False,
            float_format=lambda value: f"{value:.3f}",
        )
    )

    print("\nFigures written to:")
    print(f"  {FIGURE_DIR}")

    print("\nTables written to:")
    print(f"  {TABLE_DIR}")


if __name__ == "__main__":
    main()

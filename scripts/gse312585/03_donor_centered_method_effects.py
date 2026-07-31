#!/usr/bin/env python3

from pathlib import Path
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


from sklearn.decomposition import PCA


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATASET = "GSE312585"


EXPRESSION_FILE = PROJECT_ROOT / "data" / "processed" / DATASET / "mirna_log2_cpm.tsv"
METADATA_FILE = PROJECT_ROOT / "data" / "processed" / DATASET / "sample_metadata.tsv"

FIGURE_DIR = PROJECT_ROOT / "results" / DATASET / "figures"
TABLE_DIR = PROJECT_ROOT / "results" / DATASET / "tables"

METHOD_ORDER = ["WS", "UC", "Q", "M"]


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

    sample_ids = metadata["sample_id"].tolist()
    expression = expression.loc[:, sample_ids]

    # Samples x features.
    sample_matrix = expression.T.copy()

    # Remove each donor's feature-wise mean.
    donor_centered = sample_matrix.copy()

    for donor, donor_metadata in metadata.groupby("donor"):
        donor_samples = donor_metadata["sample_id"].tolist()

        donor_centered.loc[donor_samples] = (
            sample_matrix.loc[donor_samples]
            - sample_matrix.loc[donor_samples].mean(axis=0)
        )

    donor_centered.T.to_csv(
        TABLE_DIR / "donor_centered_log2_cpm.tsv",
        sep="\t",
    )

    # PCA on donor-centred profiles.
    pca = PCA()
    scores = pca.fit_transform(donor_centered)

    explained = pca.explained_variance_ratio_ * 100

    pca_table = metadata.copy()
    pca_table["PC1"] = scores[:, 0]
    pca_table["PC2"] = scores[:, 1]
    pca_table["PC3"] = scores[:, 2]

    pca_table.to_csv(
        TABLE_DIR / "donor_centered_pca_coordinates.tsv",
        sep="\t",
        index=False,
    )

    variance_table = pd.DataFrame(
        {
            "component": [
                f"PC{i + 1}"
                for i in range(len(explained))
            ],
            "explained_variance_percent": explained,
            "cumulative_variance_percent": np.cumsum(explained),
        }
    )

    variance_table.to_csv(
        TABLE_DIR / "donor_centered_pca_variance.tsv",
        sep="\t",
        index=False,
    )

    # PCA loadings.
    loading_table = pd.DataFrame(
        {
            "feature": donor_centered.columns,
            "PC1_loading": pca.components_[0],
            "PC2_loading": pca.components_[1],
            "PC1_absolute_loading": np.abs(pca.components_[0]),
            "PC2_absolute_loading": np.abs(pca.components_[1]),
        }
    )

    loading_table["maximum_PC1_PC2_loading"] = loading_table[
        ["PC1_absolute_loading", "PC2_absolute_loading"]
    ].max(axis=1)

    loading_table = loading_table.sort_values(
        "maximum_PC1_PC2_loading",
        ascending=False,
    )

    loading_table.to_csv(
        TABLE_DIR / "donor_centered_pca_loadings.tsv",
        sep="\t",
        index=False,
    )

    # Euclidean distances between methods within each donor.
    distance_rows = []

    for donor, donor_metadata in metadata.groupby("donor"):
        donor_metadata = donor_metadata.set_index("method")

        for method_a, method_b in combinations(METHOD_ORDER, 2):
            sample_a = donor_metadata.loc[method_a, "sample_id"]
            sample_b = donor_metadata.loc[method_b, "sample_id"]

            vector_a = donor_centered.loc[sample_a].to_numpy()
            vector_b = donor_centered.loc[sample_b].to_numpy()

            distance = np.linalg.norm(vector_a - vector_b)

            distance_rows.append(
                {
                    "donor": donor,
                    "method_a": method_a,
                    "method_b": method_b,
                    "comparison": f"{method_a}_vs_{method_b}",
                    "euclidean_distance": distance,
                }
            )

    distance_table = pd.DataFrame(distance_rows)

    distance_table.to_csv(
        TABLE_DIR / "within_donor_method_distances.tsv",
        sep="\t",
        index=False,
    )

    distance_summary = (
        distance_table
        .groupby("comparison", as_index=False)
        .agg(
            mean_distance=("euclidean_distance", "mean"),
            sd_distance=("euclidean_distance", "std"),
            minimum_distance=("euclidean_distance", "min"),
            maximum_distance=("euclidean_distance", "max"),
        )
        .sort_values("mean_distance", ascending=False)
    )

    distance_summary.to_csv(
        TABLE_DIR / "within_donor_method_distance_summary.tsv",
        sep="\t",
        index=False,
    )

    # Method centroids after donor centring.
    centroid_rows = []

    for method in METHOD_ORDER:
        method_samples = metadata.loc[
            metadata["method"] == method,
            "sample_id",
        ].tolist()

        centroid = donor_centered.loc[method_samples].mean(axis=0)

        centroid_rows.append(
            pd.Series(
                centroid,
                name=method,
            )
        )

    method_centroids = pd.DataFrame(centroid_rows)
    method_centroids.index.name = "method"

    method_centroids.to_csv(
        TABLE_DIR / "donor_centered_method_centroids.tsv",
        sep="\t",
    )

    # Plot donor-centred PCA.
    default_colors = plt.rcParams[
        "axes.prop_cycle"
    ].by_key()["color"]

    method_colors = {
        method: default_colors[i]
        for i, method in enumerate(METHOD_ORDER)
    }

    donor_markers = {
        donor: marker
        for donor, marker in zip(
            metadata["donor"].unique(),
            ["o", "s", "^"],
        )
    }

    fig, ax = plt.subplots(figsize=(8.5, 6.8))

    # Connect methods belonging to the same donor.
    for donor, donor_metadata in pca_table.groupby("donor"):
        ordered = (
            donor_metadata
            .assign(
                method_rank=lambda x: x["method"].map(
                    {
                        method: i
                        for i, method in enumerate(METHOD_ORDER)
                    }
                )
            )
            .sort_values("method_rank")
        )

        ax.plot(
            ordered["PC1"],
            ordered["PC2"],
            linewidth=0.8,
            alpha=0.6,
        )

    for row in pca_table.itertuples(index=False):
        ax.scatter(
            row.PC1,
            row.PC2,
            color=method_colors[row.method],
            marker=donor_markers[row.donor],
            s=105,
            edgecolor="black",
            linewidth=0.6,
            zorder=3,
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
    ]

    donor_handles = [
        plt.Line2D(
            [0],
            [0],
            marker=donor_markers[donor],
            linestyle="",
            markerfacecolor="white",
            markeredgecolor="black",
            label=donor,
            markersize=8,
        )
        for donor in donor_markers
    ]

    first_legend = ax.legend(
        handles=method_handles,
        title="Method",
        loc="upper left",
    )

    ax.add_artist(first_legend)

    ax.legend(
        handles=donor_handles,
        title="Donor",
        loc="lower right",
    )

    ax.axhline(0, linewidth=0.7)
    ax.axvline(0, linewidth=0.7)

    ax.set_xlabel(
        f"PC1 ({explained[0]:.1f}% variance)"
    )
    ax.set_ylabel(
        f"PC2 ({explained[1]:.1f}% variance)"
    )
    ax.set_title(
        "Donor-centred salivary miRNA-associated profiles"
    )

    fig.tight_layout()

    fig.savefig(
        FIGURE_DIR / "donor_centered_pca.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    # Plot within-donor method distances.
    comparison_order = distance_summary["comparison"].tolist()

    fig, ax = plt.subplots(figsize=(9.5, 5.8))

    positions = {
        comparison: i
        for i, comparison in enumerate(comparison_order)
    }

    for row in distance_table.itertuples(index=False):
        ax.scatter(
            positions[row.comparison],
            row.euclidean_distance,
            s=70,
        )

    means = distance_summary.set_index(
        "comparison"
    ).loc[comparison_order, "mean_distance"]

    ax.scatter(
        range(len(comparison_order)),
        means,
        marker="_",
        s=450,
        linewidth=2.5,
        label="Mean",
    )

    ax.set_xticks(range(len(comparison_order)))
    ax.set_xticklabels(
        comparison_order,
        rotation=45,
        ha="right",
    )

    ax.set_ylabel(
        "Within-donor Euclidean distance"
    )
    ax.set_title(
        "Magnitude of isolation-method effects"
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        FIGURE_DIR / "within_donor_method_distances.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(
        f"Input: {donor_centered.shape[0]} samples x "
        f"{donor_centered.shape[1]} features"
    )

    print("\nDonor-centred PCA variance:")
    print(
        variance_table.head(5).to_string(
            index=False,
            float_format=lambda x: f"{x:.3f}",
        )
    )

    print("\nMethod-distance summary:")
    print(
        distance_summary.to_string(
            index=False,
            float_format=lambda x: f"{x:.3f}",
        )
    )

    print("\nTop PCA-associated features:")
    print(
        loading_table[
            [
                "feature",
                "PC1_loading",
                "PC2_loading",
                "maximum_PC1_PC2_loading",
            ]
        ]
        .head(20)
        .to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}",
        )
    )

    print("\nFigures:")
    print(f"  {FIGURE_DIR / 'donor_centered_pca.png'}")
    print(
        f"  {FIGURE_DIR / 'within_donor_method_distances.png'}"
    )


if __name__ == "__main__":
    main()

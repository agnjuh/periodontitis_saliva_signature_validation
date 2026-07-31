#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


MATRIX_FILE = Path("data/processed/PXD043901_lfq_matrix_log2.csv")
METADATA_FILE = Path("data/processed/PXD043901_sample_metadata.csv")

QC_DIR = Path("results/qc")
FIGURE_DIR = Path("results/figures/qc")


def main() -> None:
    if not MATRIX_FILE.exists():
        raise FileNotFoundError(f"Missing matrix: {MATRIX_FILE}")

    if not METADATA_FILE.exists():
        raise FileNotFoundError(f"Missing metadata: {METADATA_FILE}")

    matrix = pd.read_csv(MATRIX_FILE, index_col=0)
    metadata = pd.read_csv(METADATA_FILE)

    QC_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    expected_samples = metadata["sample_id"].tolist()

    missing_samples = [
        sample for sample in expected_samples
        if sample not in matrix.columns
    ]

    if missing_samples:
        raise ValueError(
            "Samples missing from matrix: "
            + ", ".join(missing_samples)
        )

    matrix = matrix[expected_samples]

    sample_metrics = pd.DataFrame({
        "sample_id": matrix.columns,
        "detected_proteins": matrix.notna().sum(axis=0).values,
        "missing_values": matrix.isna().sum(axis=0).values,
        "missing_fraction": matrix.isna().mean(axis=0).values,
        "median_log2_lfq": matrix.median(axis=0).values,
        "mean_log2_lfq": matrix.mean(axis=0).values,
    })

    sample_metrics = sample_metrics.merge(
        metadata[
            [
                "sample_id",
                "patient_id",
                "timepoint",
                "condition",
            ]
        ],
        on="sample_id",
        how="left",
        validate="one_to_one",
    )

    sample_metrics.to_csv(
        QC_DIR / "sample_qc_metrics.csv",
        index=False,
    )

    protein_metrics = pd.DataFrame({
        "protein_group_id": matrix.index,
        "detected_samples": matrix.notna().sum(axis=1).values,
        "missing_fraction": matrix.isna().mean(axis=1).values,
        "median_log2_lfq": matrix.median(axis=1).values,
        "variance_log2_lfq": matrix.var(axis=1).values,
    })

    protein_metrics.to_csv(
        QC_DIR / "protein_qc_metrics.csv",
        index=False,
    )

    plt.figure(figsize=(11, 6))
    plt.bar(
        sample_metrics["sample_id"],
        sample_metrics["detected_proteins"],
    )
    plt.xticks(rotation=90)
    plt.ylabel("Detected protein groups")
    plt.xlabel("Sample")
    plt.title("Protein detection per saliva sample")
    plt.tight_layout()
    plt.savefig(
        FIGURE_DIR / "detected_proteins_per_sample.png",
        dpi=300,
    )
    plt.close()

    plt.figure(figsize=(11, 6))
    plt.bar(
        sample_metrics["sample_id"],
        sample_metrics["missing_fraction"],
    )
    plt.xticks(rotation=90)
    plt.ylabel("Missing fraction")
    plt.xlabel("Sample")
    plt.title("Missing LFQ values per saliva sample")
    plt.tight_layout()
    plt.savefig(
        FIGURE_DIR / "missing_fraction_per_sample.png",
        dpi=300,
    )
    plt.close()

    plt.figure(figsize=(12, 7))
    matrix.boxplot(
        column=list(matrix.columns),
        rot=90,
        grid=False,
    )
    plt.ylabel("log2 LFQ intensity")
    plt.xlabel("Sample")
    plt.title("Sample-wise LFQ intensity distributions")
    plt.tight_layout()
    plt.savefig(
        FIGURE_DIR / "lfq_intensity_boxplots.png",
        dpi=300,
    )
    plt.close()

    correlation = matrix.corr(
        method="spearman",
        min_periods=100,
    )

    correlation.to_csv(
        QC_DIR / "sample_spearman_correlation.csv"
    )

    plt.figure(figsize=(10, 9))
    image = plt.imshow(
        correlation,
        aspect="auto",
        vmin=-1,
        vmax=1,
    )
    plt.colorbar(image, label="Spearman correlation")
    plt.xticks(
        range(len(correlation.columns)),
        correlation.columns,
        rotation=90,
        fontsize=7,
    )
    plt.yticks(
        range(len(correlation.index)),
        correlation.index,
        fontsize=7,
    )
    plt.title("Sample-to-sample Spearman correlation")
    plt.tight_layout()
    plt.savefig(
        FIGURE_DIR / "sample_correlation_heatmap.png",
        dpi=300,
    )
    plt.close()

    pca_input = matrix.T.copy()

    protein_medians = pca_input.median(axis=0)
    pca_input = pca_input.fillna(protein_medians)

    zero_variance = pca_input.var(axis=0).eq(0)
    pca_input = pca_input.loc[:, ~zero_variance]

    scaled = StandardScaler().fit_transform(pca_input)

    pca = PCA(n_components=5)
    coordinates = pca.fit_transform(scaled)

    pca_table = metadata[
        [
            "sample_id",
            "patient_id",
            "timepoint",
            "condition",
        ]
    ].copy()

    for component in range(5):
        pca_table[f"PC{component + 1}"] = coordinates[:, component]

    pca_table.to_csv(
        QC_DIR / "sample_pca_coordinates.csv",
        index=False,
    )

    variance_table = pd.DataFrame({
        "component": [
            f"PC{i + 1}"
            for i in range(len(pca.explained_variance_ratio_))
        ],
        "explained_variance_ratio": pca.explained_variance_ratio_,
        "explained_variance_percent":
            pca.explained_variance_ratio_ * 100,
    })

    variance_table.to_csv(
        QC_DIR / "pca_explained_variance.csv",
        index=False,
    )

    plt.figure(figsize=(8, 7))

    timepoint_markers = {
        "T0": "o",
        "T1": "s",
        "T6": "^",
    }

    for timepoint, group in pca_table.groupby(
        "timepoint",
        sort=False,
    ):
        plt.scatter(
            group["PC1"],
            group["PC2"],
            marker=timepoint_markers[timepoint],
            s=70,
            label=timepoint,
        )

        for _, row in group.iterrows():
            plt.annotate(
                row["sample_id"],
                (row["PC1"], row["PC2"]),
                fontsize=7,
                xytext=(3, 3),
                textcoords="offset points",
            )

    pc1_variance = variance_table.loc[
        variance_table["component"] == "PC1",
        "explained_variance_percent",
    ].iloc[0]

    pc2_variance = variance_table.loc[
        variance_table["component"] == "PC2",
        "explained_variance_percent",
    ].iloc[0]

    plt.xlabel(f"PC1 ({pc1_variance:.1f}%)")
    plt.ylabel(f"PC2 ({pc2_variance:.1f}%)")
    plt.title("PCA of salivary LFQ proteomes")
    plt.legend(title="Timepoint")
    plt.tight_layout()
    plt.savefig(
        FIGURE_DIR / "pca_samples.png",
        dpi=300,
    )
    plt.close()

    summary_lines = [
        f"Samples: {matrix.shape[1]}",
        f"Retained protein groups: {matrix.shape[0]}",
        (
            "Detected proteins per sample: "
            f"{sample_metrics['detected_proteins'].min()}-"
            f"{sample_metrics['detected_proteins'].max()}"
        ),
        (
            "Median detected proteins per sample: "
            f"{sample_metrics['detected_proteins'].median():.1f}"
        ),
        (
            "Missing fraction range: "
            f"{sample_metrics['missing_fraction'].min():.4f}-"
            f"{sample_metrics['missing_fraction'].max():.4f}"
        ),
        (
            "Median sample missing fraction: "
            f"{sample_metrics['missing_fraction'].median():.4f}"
        ),
        f"PC1 explained variance: {pc1_variance:.2f}%",
        f"PC2 explained variance: {pc2_variance:.2f}%",
    ]

    summary_path = QC_DIR / "qc_summary.txt"
    summary_path.write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    print("\n".join(summary_lines))
    print(f"\nQC tables saved to: {QC_DIR}")
    print(f"QC figures saved to: {FIGURE_DIR}")


if __name__ == "__main__":
    main()

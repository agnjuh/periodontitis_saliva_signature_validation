#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


INPUT = Path(
    "results/tables/"
    "known_mirna_collapsed_unique_sequences.tsv"
)

TABLE_DIR = Path("results/tables")
FIGURE_DIR = Path("results/figures")

PRIORITY_OUTPUT = TABLE_DIR / "exploratory_mirna_priorities.tsv"

SCATTER_OUTPUT = FIGURE_DIR / "case_vs_control_tpm.png"
EFFECT_OUTPUT = FIGURE_DIR / "effect_abundance_plot.png"
TOP_OUTPUT = FIGURE_DIR / "top20_mirna_abundance.png"


def save_figure(path: Path) -> None:
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(f"Missing input file: {INPUT}")

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(INPUT, sep="\t")

    required = {
        "Known miRNA",
        "tpm_case",
        "tpm_control",
        "log2FC_case_vs_control",
        "mean_tpm",
    }

    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns: {sorted(missing)}"
        )

    # Exploratory thresholds only.
    # These are not statistical significance thresholds.
    df["passes_abundance_filter"] = df["mean_tpm"] >= 10
    df["passes_effect_filter"] = (
        df["log2FC_case_vs_control"].abs() >= 1
    )
    df["exploratory_priority"] = (
        df["passes_abundance_filter"]
        & df["passes_effect_filter"]
    )

    df["direction"] = np.select(
        [
            df["log2FC_case_vs_control"] >= 1,
            df["log2FC_case_vs_control"] <= -1,
        ],
        [
            "higher_in_case",
            "higher_in_control",
        ],
        default="limited_change",
    )

    priorities = (
        df.loc[df["exploratory_priority"]]
        .sort_values(
            ["mean_tpm", "log2FC_case_vs_control"],
            ascending=[False, False],
        )
        .copy()
    )

    priorities.to_csv(
        PRIORITY_OUTPUT,
        sep="\t",
        index=False,
    )

 
    # Figure 1: case versus control TPM
   
    plt.figure(figsize=(7, 7))

    x = np.log10(df["tpm_control"] + 1)
    y = np.log10(df["tpm_case"] + 1)

    plt.scatter(x, y, alpha=0.75)

    lower = min(x.min(), y.min())
    upper = max(x.max(), y.max())

    plt.plot(
        [lower, upper],
        [lower, upper],
        linestyle="--",
        linewidth=1,
    )

    plt.xlabel("log10(TPM control + 1)")
    plt.ylabel("log10(TPM case + 1)")
    plt.title(
        "GSE121910: salivary miRNA abundance\n"
        "periodontitis case versus control"
    )

    labels = (
        df.assign(abs_effect=df["log2FC_case_vs_control"].abs())
        .sort_values(
            ["mean_tpm", "abs_effect"],
            ascending=False,
        )
        .head(10)
    )

    for _, row in labels.iterrows():
        plt.annotate(
            row["Known miRNA"],
            (
                np.log10(row["tpm_control"] + 1),
                np.log10(row["tpm_case"] + 1),
            ),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )

    save_figure(SCATTER_OUTPUT)

    
    # Figure 2: effect-abundance plot
   
    plt.figure(figsize=(8, 6))

    x = np.log10(df["mean_tpm"] + 1)
    y = df["log2FC_case_vs_control"]

    plt.scatter(x, y, alpha=0.75)

    plt.axhline(0, linewidth=1)
    plt.axhline(1, linestyle="--", linewidth=1)
    plt.axhline(-1, linestyle="--", linewidth=1)
    plt.axvline(
        np.log10(10 + 1),
        linestyle="--",
        linewidth=1,
    )

    plt.xlabel("log10(mean TPM + 1)")
    plt.ylabel("Descriptive log2 fold change: case/control")
    plt.title(
        "GSE121910: exploratory effect–abundance profile"
    )

    priority_labels = (
        priorities.assign(
            abs_effect=priorities[
                "log2FC_case_vs_control"
            ].abs()
        )
        .sort_values(
            ["mean_tpm", "abs_effect"],
            ascending=False,
        )
        .head(12)
    )

    for _, row in priority_labels.iterrows():
        plt.annotate(
            row["Known miRNA"],
            (
                np.log10(row["mean_tpm"] + 1),
                row["log2FC_case_vs_control"],
            ),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )

    save_figure(EFFECT_OUTPUT)

    
    # Figure 3: top 20 by mean abundance
   
    top20 = (
        df.sort_values("mean_tpm", ascending=False)
        .head(20)
        .sort_values("mean_tpm", ascending=True)
    )

    plt.figure(figsize=(8, 7))
    plt.barh(
        top20["Known miRNA"],
        np.log10(top20["mean_tpm"] + 1),
    )

    plt.xlabel("log10(mean TPM + 1)")
    plt.ylabel("Known miRNA")
    plt.title(
        "GSE121910: top 20 salivary miRNAs by mean abundance"
    )

    save_figure(TOP_OUTPUT)

    print(f"Input miRNAs: {len(df)}")
    print(
        "Exploratory priority miRNAs "
        "(mean TPM >= 10 and |log2FC| >= 1): "
        f"{len(priorities)}"
    )

    print("\nPriority miRNAs:")
    if priorities.empty:
        print("None")
    else:
        print(
            priorities[
                [
                    "Known miRNA",
                    "tpm_case",
                    "tpm_control",
                    "mean_tpm",
                    "log2FC_case_vs_control",
                    "direction",
                ]
            ].to_string(index=False)
        )

    print("\nOutputs:")
    print(f"  {PRIORITY_OUTPUT}")
    print(f"  {SCATTER_OUTPUT}")
    print(f"  {EFFECT_OUTPUT}")
    print(f"  {TOP_OUTPUT}")


if __name__ == "__main__":
    main()

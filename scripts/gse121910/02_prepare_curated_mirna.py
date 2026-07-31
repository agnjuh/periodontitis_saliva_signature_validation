#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd


INPUT = Path("data/processed/known_miRNA_case_vs_control.tsv")
OUTPUT_DIR = Path("results/tables")

SEQUENCE_AUDIT = OUTPUT_DIR / "known_mirna_sequence_audit.tsv"
ALL_SUMMARY = OUTPUT_DIR / "known_mirna_collapsed_all.tsv"
UNIQUE_SUMMARY = OUTPUT_DIR / "known_mirna_collapsed_unique_sequences.tsv"


def collapse_by_mirna(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse isomiR/sequence rows to the reported miRNA identifier."""
    summary = (
        df.groupby("Known miRNA", as_index=False)
        .agg(
            n_sequence_rows=("miRNA Sequence", "size"),
            n_unique_sequences=("miRNA Sequence", "nunique"),
            case_raw=("Case Raw expression", "sum"),
            control_raw=("Control Raw expression", "sum"),
            tpm_case=("TPM Case", "sum"),
            tpm_control=("TPM Control", "sum"),
            max_genome_matches=("Genome Matches", "max"),
        )
    )

    # Descriptive effect size only: pooled case versus pooled control.
    pseudocount = 0.5
    summary["log2FC_case_vs_control"] = np.log2(
        (summary["tpm_case"] + pseudocount)
        / (summary["tpm_control"] + pseudocount)
    )

    summary["mean_tpm"] = (
        summary["tpm_case"] + summary["tpm_control"]
    ) / 2

    return summary.sort_values(
        ["mean_tpm", "Known miRNA"],
        ascending=[False, True],
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(INPUT, sep="\t")

    required_columns = {
        "miRNA Sequence",
        "Known miRNA",
        "Genome Matches",
        "Case Raw expression",
        "Control Raw expression",
        "TPM Case",
        "TPM Control",
    }

    missing = required_columns.difference(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns: {sorted(missing)}"
        )

    # Number of different miRNA labels assigned to each sequence.
    annotation_count = (
        df.groupby("miRNA Sequence")["Known miRNA"]
        .nunique()
        .rename("n_mirna_annotations")
    )

    df = df.merge(
        annotation_count,
        on="miRNA Sequence",
        how="left",
        validate="many_to_one",
    )

    df["ambiguous_mirna_mapping"] = (
        df["n_mirna_annotations"] > 1
    )

    df.to_csv(SEQUENCE_AUDIT, sep="\t", index=False)

    all_summary = collapse_by_mirna(df)
    all_summary.to_csv(ALL_SUMMARY, sep="\t", index=False)

    unique_df = df.loc[
        ~df["ambiguous_mirna_mapping"]
    ].copy()

    unique_summary = collapse_by_mirna(unique_df)
    unique_summary.to_csv(
        UNIQUE_SUMMARY,
        sep="\t",
        index=False,
    )

    print(f"Input rows: {len(df)}")
    print(
        "Unique sequences: "
        f"{df['miRNA Sequence'].nunique()}"
    )
    print(
        "Unique miRNA identifiers: "
        f"{df['Known miRNA'].nunique()}"
    )
    print(
        "Rows with ambiguous miRNA assignment: "
        f"{df['ambiguous_mirna_mapping'].sum()}"
    )
    print(
        "Sequences with ambiguous miRNA assignment: "
        f"{df.loc[df['ambiguous_mirna_mapping'], 'miRNA Sequence'].nunique()}"
    )
    print(
        "miRNAs in all-sequence summary: "
        f"{len(all_summary)}"
    )
    print(
        "miRNAs in unique-sequence summary: "
        f"{len(unique_summary)}"
    )

    print("\nTop 10 miRNAs by mean TPM — unique sequences only:")
    print(
        unique_summary[
            [
                "Known miRNA",
                "tpm_case",
                "tpm_control",
                "log2FC_case_vs_control",
                "mean_tpm",
            ]
        ]
        .head(10)
        .to_string(index=False)
    )

    print("\nOutputs:")
    print(f"  {SEQUENCE_AUDIT}")
    print(f"  {ALL_SUMMARY}")
    print(f"  {UNIQUE_SUMMARY}")


if __name__ == "__main__":
    main()

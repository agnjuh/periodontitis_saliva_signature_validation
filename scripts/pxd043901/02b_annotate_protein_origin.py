#!/usr/bin/env python3

import re
from pathlib import Path

import pandas as pd


PROTEIN_GROUPS_FILE = Path(
    "data/raw/PXD043901/proteinGroups.txt"
)
CURRENT_ANNOTATION_FILE = Path(
    "data/processed/PXD043901_protein_annotation.csv"
)
OUTPUT_FILE = Path(
    "data/processed/PXD043901_protein_annotation_extended.csv"
)
SUMMARY_FILE = Path(
    "data/processed/PXD043901_organism_summary.txt"
)


def first_entry(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).split(";")[0].strip()


def extract_primary_accession(protein_ids: object) -> str:
    entry = first_entry(protein_ids)

    match = re.match(r"^(?:sp|tr)\|([^|]+)\|", entry)
    if match:
        return match.group(1)

    return entry


def extract_gene_symbol(fasta_headers: object) -> str:
    if pd.isna(fasta_headers):
        return ""

    header = str(fasta_headers)

    match = re.search(
        r"OS=Homo sapiens\s+GN=([^\s;]+)",
        header,
    )

    if match:
        return match.group(1)

    return ""


def extract_human_protein_name(fasta_headers: object) -> str:
    if pd.isna(fasta_headers):
        return ""

    entries = [
        entry.strip()
        for entry in str(fasta_headers).split(";")
        if entry.strip()
    ]

    for entry in entries:
        if "OS=Homo sapiens" not in entry:
            continue

        match = re.match(
            r"^(?:sp|tr)\|[^|]+\|[^\s]+\s+(.+?)\s+OS=Homo sapiens",
            entry,
        )

        if match:
            return match.group(1).strip()

    return ""


def classify_origin(fasta_headers: object) -> str:
    if pd.isna(fasta_headers) or not str(fasta_headers).strip():
        return "unknown"

    entries = [
        entry.strip()
        for entry in str(fasta_headers).split(";")
        if entry.strip()
    ]

    human_flags = [
        "OS=Homo sapiens" in entry
        for entry in entries
    ]

    if all(human_flags):
        return "human"

    if any(human_flags):
        return "mixed_human_microbiome"

    return "microbiome"


def extract_first_microbe_label(fasta_headers: object) -> str:
    if pd.isna(fasta_headers):
        return ""

    entries = [
        entry.strip()
        for entry in str(fasta_headers).split(";")
        if entry.strip()
    ]

    for entry in entries:
        if "OS=Homo sapiens" in entry:
            continue

        entry = re.sub(r"^\S+\s+", "", entry)
        entry = re.sub(r"\s+\[[^\]]+\].*$", "", entry)
        entry = re.sub(r"\s+\(REVERSE SENSE\).*$", "", entry)
        entry = re.sub(r"\s+ORF$", "", entry)

        return entry.strip()

    return ""


def main() -> None:
    if not PROTEIN_GROUPS_FILE.exists():
        raise FileNotFoundError(
            f"Missing file: {PROTEIN_GROUPS_FILE}"
        )

    if not CURRENT_ANNOTATION_FILE.exists():
        raise FileNotFoundError(
            f"Missing file: {CURRENT_ANNOTATION_FILE}"
        )

    protein_groups = pd.read_csv(
        PROTEIN_GROUPS_FILE,
        sep="\t",
        low_memory=False,
    )

    annotation = pd.read_csv(
        CURRENT_ANNOTATION_FILE,
        dtype={"protein_group_id": str},
    )

    protein_groups = protein_groups.reset_index().rename(
        columns={"index": "protein_group_id"}
    )

    protein_groups["protein_group_id"] = (
        protein_groups["protein_group_id"].astype(str)
    )

    required = [
        "protein_group_id",
        "Protein IDs",
        "Majority protein IDs",
        "Fasta headers",
    ]

    source = protein_groups[required].copy()

    extended = annotation.merge(
        source,
        on="protein_group_id",
        how="left",
        validate="one_to_one",
        suffixes=("", "_source"),
    )

    protein_id_column = (
        "Majority protein IDs_source"
        if "Majority protein IDs_source" in extended.columns
        else "Majority protein IDs"
    )

    fasta_column = (
        "Fasta headers_source"
        if "Fasta headers_source" in extended.columns
        else "Fasta headers"
    )

    extended["primary_accession"] = extended[
        protein_id_column
    ].apply(extract_primary_accession)

    extended["organism_origin"] = extended[
        fasta_column
    ].apply(classify_origin)

    extended["gene_symbol"] = extended[
        fasta_column
    ].apply(extract_gene_symbol)

    extended["human_protein_name"] = extended[
        fasta_column
    ].apply(extract_human_protein_name)

    extended["first_microbe_label"] = extended[
        fasta_column
    ].apply(extract_first_microbe_label)

    duplicate_columns = [
        column
        for column in extended.columns
        if column.endswith("_source")
    ]

    extended = extended.drop(columns=duplicate_columns)

    extended.to_csv(OUTPUT_FILE, index=False)

    origin_counts = (
        extended["organism_origin"]
        .value_counts(dropna=False)
        .rename_axis("organism_origin")
        .reset_index(name="protein_groups")
    )

    summary_lines = [
        f"Retained protein groups: {len(extended)}",
        "",
        "Origin classification:",
        origin_counts.to_string(index=False),
        "",
        (
            "Human protein groups with gene symbol: "
            f"{extended['gene_symbol'].ne('').sum()}"
        ),
        (
            "Human protein groups with parsed protein name: "
            f"{extended['human_protein_name'].ne('').sum()}"
        ),
    ]

    SUMMARY_FILE.write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    print("\n".join(summary_lines))
    print(f"\nSaved: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()

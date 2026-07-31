#!/usr/bin/env python3

"""
Build experimentally validated miRNA–protein regulatory edges.

Inputs
------
1. Aggregated validated miRNA-target interactions from multiMiR.
2. Candidate network node registry produced by the network-preparation step.

Outputs
-------
- regulatory_edges_all.tsv
- regulatory_edges_high_confidence.tsv
- regulatory_target_overlap.tsv
- unmatched_candidate_proteins.tsv
- regulatory_edge_summary.tsv
- regulatory_edge_report.txt

The script creates only experimentally validated edges already present in
multiMiR. It does not add predicted interactions.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd



# Constants


DEFAULT_MULTIMIR_PATH = Path(
    "data/reference/multimir/multimir_validated_targets.tsv"
)

DEFAULT_OUTPUT_DIR = Path(
    "results/integration/regulatory_network"
)

NODE_REGISTRY_CANDIDATES = [
    Path("results/integration/network_node_registry.tsv"),
    Path("results/integration/biological_network/network_node_registry.tsv"),
    Path("results/network/network_node_registry.tsv"),
    Path("results/network_node_registry.tsv"),
    Path("data/processed/network_node_registry.tsv"),
]

MIRNA_ALIASES = {
    "MIR-143-3P": "hsa-miR-143-3p",
    "MIR143-3P": "hsa-miR-143-3p",
    "HSA-MIR-143-3P": "hsa-miR-143-3p",
    "HSA-MIR143-3P": "hsa-miR-143-3p",
    "HSA-MIR-143-3P-": "hsa-miR-143-3p",
}

# Curated aliases encountered frequently in miRNA-target resources.
# This is deliberately conservative rather than exhaustive.
GENE_ALIASES = {
    "ERK5": "MAPK7",
    "BCL-2": "BCL2",
    "BCL 2": "BCL2",
    "P53": "TP53",
    "ERBB1": "EGFR",
    "HER1": "EGFR",
    "HER2": "ERBB2",
    "C-MYC": "MYC",
    "CMYC": "MYC",
    "C-JUN": "JUN",
    "CJUN": "JUN",
    "C-FOS": "FOS",
    "CFOS": "FOS",
    "PKB": "AKT1",
    "PKBΑ": "AKT1",
    "PKBALPHA": "AKT1",
    "P21": "CDKN1A",
    "P27": "CDKN1B",
    "P16": "CDKN2A",
    "VEGF": "VEGFA",
    "TGF-BETA1": "TGFB1",
    "TGFΒ1": "TGFB1",
    "COX2": "PTGS2",
    "COX-2": "PTGS2",
    "MMP-9": "MMP9",
    "MMP-2": "MMP2",
    "IL-6": "IL6",
    "IL-8": "CXCL8",
    "TNF-ALPHA": "TNF",
    "TNFΑ": "TNF",
    "NFKB1A": "NFKBIA",
}



# Utility functions


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Intersect experimentally validated multiMiR targets with "
            "candidate protein nodes."
        )
    )

    parser.add_argument(
        "--multimir",
        type=Path,
        default=DEFAULT_MULTIMIR_PATH,
        help=(
            "Path to multimir_validated_targets.tsv "
            f"(default: {DEFAULT_MULTIMIR_PATH})"
        ),
    )

    parser.add_argument(
        "--node-registry",
        type=Path,
        default=None,
        help=(
            "Path to network_node_registry.tsv. If omitted, the script "
            "searches common repository locations."
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )

    parser.add_argument(
        "--minimum-databases",
        type=int,
        default=1,
        help=(
            "Minimum number of supporting databases for the "
            "high-confidence edge list."
        ),
    )

    parser.add_argument(
        "--minimum-publications",
        type=int,
        default=1,
        help=(
            "Minimum number of supporting publications for the "
            "high-confidence edge list."
        ),
    )

    parser.add_argument(
        "--require-strong-evidence",
        action="store_true",
        help=(
            "Require maximum_evidence_strength == strong for the "
            "high-confidence edge list."
        ),
    )

    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    frame = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    if frame.empty:
        raise ValueError(f"Input file contains no rows: {path}")

    return frame


def detect_column(
    frame: pd.DataFrame,
    candidates: Iterable[str],
    *,
    required: bool = True,
) -> str | None:
    exact_lookup = {
        column.lower(): column
        for column in frame.columns
    }

    for candidate in candidates:
        if candidate.lower() in exact_lookup:
            return exact_lookup[candidate.lower()]

    normalised_columns = {
        re.sub(r"[^a-z0-9]+", "", column.lower()): column
        for column in frame.columns
    }

    for candidate in candidates:
        key = re.sub(r"[^a-z0-9]+", "", candidate.lower())

        if key in normalised_columns:
            return normalised_columns[key]

    if required:
        raise KeyError(
            "Could not identify required column. "
            f"Candidates: {list(candidates)}. "
            f"Available columns: {list(frame.columns)}"
        )

    return None


def find_node_registry(explicit_path: Path | None) -> Path:
    if explicit_path is not None:
        if not explicit_path.exists():
            raise FileNotFoundError(
                f"Node registry not found: {explicit_path}"
            )
        return explicit_path

    for candidate in NODE_REGISTRY_CANDIDATES:
        if candidate.exists():
            return candidate

    matches = sorted(
        Path(".").glob("**/network_node_registry.tsv")
    )

    if len(matches) == 1:
        return matches[0]

    if len(matches) > 1:
        formatted = "\n".join(
            f"  - {match}"
            for match in matches
        )
        raise FileNotFoundError(
            "Multiple network_node_registry.tsv files were found. "
            "Specify one with --node-registry:\n"
            f"{formatted}"
        )

    raise FileNotFoundError(
        "Could not locate network_node_registry.tsv. "
        "Provide its path with --node-registry."
    )


def clean_text(value: object) -> str:
    if value is None:
        return ""

    text = str(value).strip()

    if text.lower() in {"nan", "none", "na", "n/a"}:
        return ""

    return text


def normalise_gene_symbol(value: object) -> str:
    text = clean_text(value).upper()

    if not text:
        return ""

    text = (
        text.replace("–", "-")
        .replace("—", "-")
        .replace("−", "-")
        .replace("_", "-")
    )

    text = re.sub(r"\s+", " ", text).strip()

    if text in GENE_ALIASES:
        return GENE_ALIASES[text]

    compact = re.sub(r"[^A-Z0-9]", "", text)

    if compact in GENE_ALIASES:
        return GENE_ALIASES[compact]

    # Standard HGNC symbols generally contain letters, numbers, and
    # occasionally hyphens. For matching proteomics symbols, remove
    # punctuation left over from legacy aliases.
    return compact


def normalise_mirna(value: object) -> str:
    text = clean_text(value)

    if not text:
        return ""

    text = (
        text.replace("_", "-")
        .replace("–", "-")
        .replace("—", "-")
    )

    text = re.sub(r"\s+", "", text)

    alias_key = text.upper()

    if alias_key in MIRNA_ALIASES:
        return MIRNA_ALIASES[alias_key]

    lower = text.lower()

    if lower.startswith("hsa-mir-"):
        suffix = text[len("hsa-miR-"):]
        return f"hsa-miR-{suffix}"

    if lower.startswith("mir-"):
        suffix = text[len("miR-"):]
        return f"hsa-miR-{suffix}"

    return text


def parse_integer_series(series: pd.Series) -> pd.Series:
    return pd.to_numeric(
        series.replace("", "0"),
        errors="coerce",
    ).fillna(0).astype(int)


def parse_boolean_series(series: pd.Series) -> pd.Series:
    true_values = {
        "true",
        "t",
        "yes",
        "y",
        "1",
    }

    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin(true_values)
    )


def collapse_unique(values: Iterable[object]) -> str:
    unique_values: set[str] = set()

    for value in values:
        text = clean_text(value)

        if not text:
            continue

        parts = re.split(r"[;|]+", text)

        for part in parts:
            part = part.strip()

            if part:
                unique_values.add(part)

    return ";".join(sorted(unique_values))


def select_candidate_proteins(
    node_registry: pd.DataFrame,
) -> pd.DataFrame:
    feature_column = detect_column(
        node_registry,
        [
            "gene_symbol",
            "feature_name",
            "node_name",
            "feature",
            "target_gene",
            "feature_id",
            "node_id",
        ],
    )

    feature_type_column = detect_column(
        node_registry,
        [
            "feature_type",
            "node_type",
            "entity_type",
            "molecule_type",
        ],
    )

    network_eligible_column = detect_column(
        node_registry,
        [
            "network_eligible",
            "biological_network_eligible",
            "is_network_eligible",
        ],
        required=False,
    )

    biological_node_column = detect_column(
        node_registry,
        [
            "biological_node",
            "is_biological_node",
            "biological_evidence",
        ],
        required=False,
    )

    candidate_class_column = detect_column(
        node_registry,
        [
            "candidate_class",
            "confidence_class",
            "priority_class",
            "evidence_class",
        ],
        required=False,
    )

    output = node_registry.copy()

    output["candidate_feature_id"] = (
        output[feature_column]
        .astype(str)
        .map(clean_text)
    )

    output["candidate_gene_symbol"] = (
        output["candidate_feature_id"]
        .map(normalise_gene_symbol)
    )

    feature_types = (
        output[feature_type_column]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    protein_mask = feature_types.isin(
        {
            "protein",
            "proteins",
            "protein_gene",
            "protein-coding gene",
            "protein_coding_gene",
        }
    )

    if not protein_mask.any():
        protein_mask = feature_types.str.contains(
            "protein",
            regex=False,
            na=False,
        ) & ~feature_types.str.contains(
            "module",
            regex=False,
            na=False,
        )

    output = output.loc[
        protein_mask
    ].copy()

    if network_eligible_column is not None:
        output = output.loc[
            parse_boolean_series(
                output[network_eligible_column]
            )
        ].copy()

    if biological_node_column is not None:
        biological_values = parse_boolean_series(
            output[biological_node_column]
        )

        if biological_values.any():
            output = output.loc[
                biological_values
            ].copy()

    output = output.loc[
        output["candidate_gene_symbol"].ne("")
    ].copy()

    selected_columns = [
        "candidate_feature_id",
        "candidate_gene_symbol",
    ]

    if candidate_class_column is not None:
        output["candidate_class"] = (
            output[candidate_class_column]
            .astype(str)
            .map(clean_text)
        )
        selected_columns.append("candidate_class")

    for optional_column in [
        "reproducibility_score",
        "aggregate_score",
        "priority_score",
        "confidence_score",
        "n_datasets",
        "n_cohorts",
        "direction_consistency",
    ]:
        detected = detect_column(
            output,
            [optional_column],
            required=False,
        )

        if (
            detected is not None
            and detected not in selected_columns
        ):
            output[optional_column] = output[detected]
            selected_columns.append(optional_column)

    output = (
        output[selected_columns]
        .drop_duplicates(
            subset=["candidate_gene_symbol"]
        )
        .sort_values(
            ["candidate_gene_symbol"]
        )
        .reset_index(drop=True)
    )

    if output.empty:
        raise ValueError(
            "No eligible protein nodes were identified in the node registry."
        )

    return output


def standardise_multimir(
    multimir: pd.DataFrame,
) -> pd.DataFrame:
    mirna_column = detect_column(
        multimir,
        ["mirna", "mature_mirna_id"],
    )

    target_column = detect_column(
        multimir,
        ["target_gene", "target_symbol", "gene_symbol"],
    )

    strength_column = detect_column(
        multimir,
        [
            "maximum_evidence_strength",
            "evidence_strength",
        ],
        required=False,
    )

    databases_column = detect_column(
        multimir,
        [
            "evidence_databases",
            "evidence_database",
            "database",
        ],
        required=False,
    )

    publications_column = detect_column(
        multimir,
        [
            "publication_ids",
            "publication_id",
            "pubmed_id",
            "pmid",
        ],
        required=False,
    )

    methods_column = detect_column(
        multimir,
        [
            "experimental_methods",
            "experimental_method",
        ],
        required=False,
    )

    support_column = detect_column(
        multimir,
        [
            "support_types",
            "support_type",
        ],
        required=False,
    )

    n_database_column = detect_column(
        multimir,
        ["n_databases"],
        required=False,
    )

    n_publication_column = detect_column(
        multimir,
        ["n_publications"],
        required=False,
    )

    n_method_column = detect_column(
        multimir,
        ["n_experimental_methods"],
        required=False,
    )

    n_source_column = detect_column(
        multimir,
        ["n_source_records"],
        required=False,
    )

    output = pd.DataFrame(
        {
            "source_mirna": (
                multimir[mirna_column]
                .map(normalise_mirna)
            ),
            "source_target_gene": (
                multimir[target_column]
                .astype(str)
                .map(clean_text)
            ),
            "target_gene_symbol": (
                multimir[target_column]
                .map(normalise_gene_symbol)
            ),
        }
    )

    output["maximum_evidence_strength"] = (
        multimir[strength_column].astype(str).map(clean_text)
        if strength_column is not None
        else ""
    )

    output["evidence_databases"] = (
        multimir[databases_column].astype(str).map(clean_text)
        if databases_column is not None
        else ""
    )

    output["publication_ids"] = (
        multimir[publications_column].astype(str).map(clean_text)
        if publications_column is not None
        else ""
    )

    output["experimental_methods"] = (
        multimir[methods_column].astype(str).map(clean_text)
        if methods_column is not None
        else ""
    )

    output["support_types"] = (
        multimir[support_column].astype(str).map(clean_text)
        if support_column is not None
        else ""
    )

    output["n_databases"] = (
        parse_integer_series(multimir[n_database_column])
        if n_database_column is not None
        else output["evidence_databases"].map(
            lambda value: len(
                {
                    part.strip()
                    for part in value.split(";")
                    if part.strip()
                }
            )
        )
    )

    output["n_publications"] = (
        parse_integer_series(multimir[n_publication_column])
        if n_publication_column is not None
        else output["publication_ids"].map(
            lambda value: len(
                {
                    part.strip()
                    for part in value.split(";")
                    if part.strip()
                }
            )
        )
    )

    output["n_experimental_methods"] = (
        parse_integer_series(multimir[n_method_column])
        if n_method_column is not None
        else output["experimental_methods"].map(
            lambda value: len(
                {
                    part.strip()
                    for part in value.split(";")
                    if part.strip()
                }
            )
        )
    )

    output["n_source_records"] = (
        parse_integer_series(multimir[n_source_column])
        if n_source_column is not None
        else 1
    )

    output = output.loc[
        output["source_mirna"].ne("")
        & output["target_gene_symbol"].ne("")
    ].copy()

    return output


def aggregate_alias_collisions(
    frame: pd.DataFrame,
) -> pd.DataFrame:
    """
    Re-aggregate records after gene alias normalisation.

    Example:
        ERK5 and MAPK7 collapse to MAPK7.
        BCL-2 and BCL2 collapse to BCL2.
    """

    grouped_rows: list[dict[str, object]] = []

    grouping_columns = [
        "source_mirna",
        "target_gene_symbol",
    ]

    for (
        source_mirna,
        target_gene_symbol,
    ), group in frame.groupby(
        grouping_columns,
        sort=True,
        dropna=False,
    ):
        strength_order = {
            "strong": 4,
            "experimentally_supported_unspecified": 3,
            "supportive": 2,
            "unclassified": 1,
            "": 0,
        }

        strengths = [
            clean_text(value).lower()
            for value in group["maximum_evidence_strength"]
        ]

        strongest = max(
            strengths,
            key=lambda value: strength_order.get(value, 0),
            default="",
        )

        source_symbols = sorted(
            {
                clean_text(value)
                for value in group["source_target_gene"]
                if clean_text(value)
            }
        )

        grouped_rows.append(
            {
                "source_mirna": source_mirna,
                "target_gene_symbol": target_gene_symbol,
                "source_target_symbols": ";".join(
                    source_symbols
                ),
                "maximum_evidence_strength": strongest,
                "evidence_databases": collapse_unique(
                    group["evidence_databases"]
                ),
                "publication_ids": collapse_unique(
                    group["publication_ids"]
                ),
                "experimental_methods": collapse_unique(
                    group["experimental_methods"]
                ),
                "support_types": collapse_unique(
                    group["support_types"]
                ),
                "n_source_records": int(
                    group["n_source_records"].sum()
                ),
            }
        )

    output = pd.DataFrame(grouped_rows)

    output["n_databases"] = output[
        "evidence_databases"
    ].map(
        lambda value: len(
            {
                part.strip()
                for part in value.split(";")
                if part.strip()
            }
        )
    )

    output["n_publications"] = output[
        "publication_ids"
    ].map(
        lambda value: len(
            {
                part.strip()
                for part in value.split(";")
                if part.strip()
            }
        )
    )

    output["n_experimental_methods"] = output[
        "experimental_methods"
    ].map(
        lambda value: len(
            {
                part.strip()
                for part in value.split(";")
                if part.strip()
            }
        )
    )

    return output


def build_edges(
    interactions: pd.DataFrame,
    candidate_proteins: pd.DataFrame,
) -> pd.DataFrame:
    edges = candidate_proteins.merge(
        interactions,
        how="inner",
        left_on="candidate_gene_symbol",
        right_on="target_gene_symbol",
        validate="one_to_many",
    )

    if edges.empty:
        return edges

    edges["source"] = edges["source_mirna"]
    edges["target"] = edges["candidate_gene_symbol"]
    edges["source_type"] = "miRNA"
    edges["target_type"] = "protein"
    edges["edge_type"] = "experimentally_validated_regulation"
    edges["interaction_direction"] = "miRNA_to_target"
    edges["evidence_source"] = "multiMiR"
    edges["is_experimentally_validated"] = True

    edges["edge_id"] = (
        edges["source"]
        + "::"
        + edges["target"]
        + "::validated_regulation"
    )

    preferred_columns = [
        "edge_id",
        "source",
        "target",
        "source_type",
        "target_type",
        "edge_type",
        "interaction_direction",
        "is_experimentally_validated",
        "evidence_source",
        "maximum_evidence_strength",
        "evidence_databases",
        "n_databases",
        "publication_ids",
        "n_publications",
        "experimental_methods",
        "n_experimental_methods",
        "support_types",
        "n_source_records",
        "source_target_symbols",
        "candidate_feature_id",
        "candidate_gene_symbol",
    ]

    optional_columns = [
        "candidate_class",
        "reproducibility_score",
        "aggregate_score",
        "priority_score",
        "confidence_score",
        "n_datasets",
        "n_cohorts",
        "direction_consistency",
    ]

    selected_columns = [
        column
        for column in preferred_columns + optional_columns
        if column in edges.columns
    ]

    edges = edges[selected_columns].copy()

    strength_rank = {
        "strong": 4,
        "experimentally_supported_unspecified": 3,
        "supportive": 2,
        "unclassified": 1,
        "": 0,
    }

    edges["_strength_rank"] = (
        edges["maximum_evidence_strength"]
        .astype(str)
        .str.lower()
        .map(strength_rank)
        .fillna(0)
    )

    edges = (
        edges.sort_values(
            [
                "_strength_rank",
                "n_databases",
                "n_publications",
                "n_experimental_methods",
                "target",
            ],
            ascending=[
                False,
                False,
                False,
                False,
                True,
            ],
        )
        .drop(columns="_strength_rank")
        .reset_index(drop=True)
    )

    return edges


def build_overlap_table(
    candidate_proteins: pd.DataFrame,
    edges: pd.DataFrame,
) -> pd.DataFrame:
    validated_targets = set(
        edges["target"]
    ) if not edges.empty else set()

    output = candidate_proteins.copy()

    output["has_validated_mirna_edge"] = output[
        "candidate_gene_symbol"
    ].isin(validated_targets)

    if not edges.empty:
        edge_summary = (
            edges.groupby("target", as_index=False)
            .agg(
                n_regulatory_mirnas=(
                    "source",
                    "nunique",
                ),
                maximum_n_databases=(
                    "n_databases",
                    "max",
                ),
                maximum_n_publications=(
                    "n_publications",
                    "max",
                ),
                regulatory_mirnas=(
                    "source",
                    collapse_unique,
                ),
                strongest_evidence=(
                    "maximum_evidence_strength",
                    lambda values: (
                        "strong"
                        if any(
                            clean_text(value).lower() == "strong"
                            for value in values
                        )
                        else collapse_unique(values)
                    ),
                ),
            )
            .rename(
                columns={
                    "target": "candidate_gene_symbol",
                }
            )
        )

        output = output.merge(
            edge_summary,
            on="candidate_gene_symbol",
            how="left",
        )

    output = output.sort_values(
        [
            "has_validated_mirna_edge",
            "candidate_gene_symbol",
        ],
        ascending=[
            False,
            True,
        ],
    ).reset_index(drop=True)

    return output


def write_report(
    output_path: Path,
    *,
    multimir_path: Path,
    node_registry_path: Path,
    candidate_proteins: pd.DataFrame,
    interactions: pd.DataFrame,
    edges: pd.DataFrame,
    high_confidence_edges: pd.DataFrame,
    unmatched: pd.DataFrame,
    minimum_databases: int,
    minimum_publications: int,
    require_strong_evidence: bool,
) -> None:
    candidate_count = len(candidate_proteins)
    matched_count = (
        edges["target"].nunique()
        if not edges.empty
        else 0
    )

    overlap_fraction = (
        matched_count / candidate_count
        if candidate_count
        else 0
    )

    lines = [
        "Experimentally validated miRNA–protein regulatory edges",
        "=======================================================",
        "",
        f"multiMiR input: {multimir_path}",
        f"Node registry input: {node_registry_path}",
        "",
        f"Candidate protein nodes: {candidate_count}",
        (
            "Alias-normalised validated miRNA-target interactions: "
            f"{len(interactions)}"
        ),
        f"Matched candidate proteins: {matched_count}",
        f"Unmatched candidate proteins: {len(unmatched)}",
        f"Candidate protein overlap fraction: {overlap_fraction:.4f}",
        f"All validated regulatory edges: {len(edges)}",
        (
            "High-confidence regulatory edges: "
            f"{len(high_confidence_edges)}"
        ),
        "",
        "High-confidence criteria:",
        f"  Minimum supporting databases: {minimum_databases}",
        f"  Minimum supporting publications: {minimum_publications}",
        (
            "  Strong evidence required: "
            f"{require_strong_evidence}"
        ),
        "",
        "Evidence policy:",
        (
            "  Only experimentally validated multiMiR records were "
            "used."
        ),
        "  Predicted interactions were not included.",
        (
            "  Gene aliases were conservatively normalised before "
            "matching."
        ),
        (
            "  Multiple source records were aggregated after alias "
            "normalisation."
        ),
        "",
    ]

    if not edges.empty:
        lines.extend(
            [
                "Matched targets:",
                *[
                    (
                        f"  {row.target}: "
                        f"{row.maximum_evidence_strength}; "
                        f"{row.n_databases} database(s); "
                        f"{row.n_publications} publication(s)"
                    )
                    for row in edges.itertuples()
                ],
                "",
            ]
        )

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()

    try:
        node_registry_path = find_node_registry(
            args.node_registry
        )

        multimir = read_tsv(
            args.multimir
        )

        node_registry = read_tsv(
            node_registry_path
        )

        candidate_proteins = select_candidate_proteins(
            node_registry
        )

        standardised_interactions = standardise_multimir(
            multimir
        )

        aggregated_interactions = aggregate_alias_collisions(
            standardised_interactions
        )

        edges = build_edges(
            aggregated_interactions,
            candidate_proteins,
        )

        high_confidence_mask = pd.Series(
            True,
            index=edges.index,
        )

        if not edges.empty:
            high_confidence_mask &= (
                edges["n_databases"]
                >= args.minimum_databases
            )

            high_confidence_mask &= (
                edges["n_publications"]
                >= args.minimum_publications
            )

            if args.require_strong_evidence:
                high_confidence_mask &= (
                    edges["maximum_evidence_strength"]
                    .astype(str)
                    .str.lower()
                    .eq("strong")
                )

        high_confidence_edges = edges.loc[
            high_confidence_mask
        ].copy()

        overlap = build_overlap_table(
            candidate_proteins,
            edges,
        )

        unmatched = overlap.loc[
            ~overlap["has_validated_mirna_edge"]
        ].copy()

        args.output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        all_edges_path = (
            args.output_dir
            / "regulatory_edges_all.tsv"
        )

        high_confidence_path = (
            args.output_dir
            / "regulatory_edges_high_confidence.tsv"
        )

        overlap_path = (
            args.output_dir
            / "regulatory_target_overlap.tsv"
        )

        unmatched_path = (
            args.output_dir
            / "unmatched_candidate_proteins.tsv"
        )

        summary_path = (
            args.output_dir
            / "regulatory_edge_summary.tsv"
        )

        report_path = (
            args.output_dir
            / "regulatory_edge_report.txt"
        )

        edges.to_csv(
            all_edges_path,
            sep="\t",
            index=False,
        )

        high_confidence_edges.to_csv(
            high_confidence_path,
            sep="\t",
            index=False,
        )

        overlap.to_csv(
            overlap_path,
            sep="\t",
            index=False,
        )

        unmatched.to_csv(
            unmatched_path,
            sep="\t",
            index=False,
        )

        candidate_count = len(candidate_proteins)

        matched_count = (
            edges["target"].nunique()
            if not edges.empty
            else 0
        )

        summary = pd.DataFrame(
            {
                "metric": [
                    "candidate_protein_nodes",
                    "validated_interactions_after_alias_normalisation",
                    "matched_candidate_proteins",
                    "unmatched_candidate_proteins",
                    "candidate_protein_overlap_fraction",
                    "all_validated_regulatory_edges",
                    "high_confidence_regulatory_edges",
                    "unique_regulatory_mirnas",
                ],
                "value": [
                    candidate_count,
                    len(aggregated_interactions),
                    matched_count,
                    len(unmatched),
                    (
                        matched_count / candidate_count
                        if candidate_count
                        else 0
                    ),
                    len(edges),
                    len(high_confidence_edges),
                    (
                        edges["source"].nunique()
                        if not edges.empty
                        else 0
                    ),
                ],
            }
        )

        summary.to_csv(
            summary_path,
            sep="\t",
            index=False,
        )

        write_report(
            report_path,
            multimir_path=args.multimir,
            node_registry_path=node_registry_path,
            candidate_proteins=candidate_proteins,
            interactions=aggregated_interactions,
            edges=edges,
            high_confidence_edges=high_confidence_edges,
            unmatched=unmatched,
            minimum_databases=args.minimum_databases,
            minimum_publications=args.minimum_publications,
            require_strong_evidence=args.require_strong_evidence,
        )

        print(f"Node registry: {node_registry_path}")
        print(
            "Candidate protein nodes: "
            f"{candidate_count}"
        )
        print(
            "Validated interactions after alias normalisation: "
            f"{len(aggregated_interactions)}"
        )
        print(
            "Matched candidate proteins: "
            f"{matched_count}"
        )
        print(
            "All validated regulatory edges: "
            f"{len(edges)}"
        )
        print(
            "High-confidence regulatory edges: "
            f"{len(high_confidence_edges)}"
        )
        print(
            f"Wrote all edges: {all_edges_path}"
        )
        print(
            "Wrote high-confidence edges: "
            f"{high_confidence_path}"
        )
        print(
            f"Wrote target overlap: {overlap_path}"
        )
        print(
            f"Wrote unmatched proteins: {unmatched_path}"
        )
        print(
            f"Wrote summary: {summary_path}"
        )
        print(
            f"Wrote report: {report_path}"
        )

        return 0

    except Exception as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

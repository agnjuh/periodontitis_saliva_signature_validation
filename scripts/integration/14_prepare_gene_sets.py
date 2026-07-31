#!/usr/bin/env python3

"""
Prepare standardized gene sets for downstream pathway enrichment.

Inputs
------
1. Biological network node registry.
2. Experimentally supported miRNA–protein regulatory edges.

Outputs
-------
- candidate_panel_genes.txt
- mir143_linked_genes.txt
- nonlinked_candidate_genes.txt
- candidate_panel_mapping.tsv
- gene_set_membership.tsv
- gene_set_summary.tsv
- gene_set_report.txt

No pathway enrichment is performed in this step.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


DEFAULT_NODE_REGISTRY = Path(
    "results/integration/biological_network/network_node_registry.tsv"
)

DEFAULT_EDGE_FILE = Path(
    "results/integration/regulatory_network/regulatory_edges_all.tsv"
)

DEFAULT_OUTPUT_DIR = Path(
    "results/integration/gene_sets"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare standardized gene sets for enrichment analysis."
    )

    parser.add_argument(
        "--node-registry",
        type=Path,
        default=DEFAULT_NODE_REGISTRY,
        help=f"Node registry path (default: {DEFAULT_NODE_REGISTRY})",
    )

    parser.add_argument(
        "--edge-file",
        type=Path,
        default=DEFAULT_EDGE_FILE,
        help=f"Regulatory edge path (default: {DEFAULT_EDGE_FILE})",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
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

    return frame


def parse_boolean(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "t", "yes", "y", "1"})
    )


def clean_gene_symbol(value: object) -> str:
    if value is None:
        return ""

    symbol = str(value).strip().upper()

    if symbol.lower() in {"", "nan", "none", "na", "n/a"}:
        return ""

    return symbol


def write_gene_list(path: Path, genes: list[str]) -> None:
    text = "\n".join(genes)

    if text:
        text += "\n"

    path.write_text(
        text,
        encoding="utf-8",
    )


def build_candidate_panel(
    node_registry: pd.DataFrame,
) -> pd.DataFrame:
    required_columns = {
        "feature_id",
        "feature_name",
        "feature_type",
        "network_eligible",
    }

    missing = required_columns - set(node_registry.columns)

    if missing:
        raise KeyError(
            "Node registry is missing required columns: "
            + ", ".join(sorted(missing))
        )

    feature_types = (
        node_registry["feature_type"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    eligible = parse_boolean(
        node_registry["network_eligible"]
    )

    proteins = node_registry.loc[
        feature_types.eq("protein") & eligible
    ].copy()

    proteins["gene_symbol"] = (
        proteins["feature_name"]
        .map(clean_gene_symbol)
    )

    proteins["uniprot_accession"] = (
        proteins["feature_id"]
        .astype(str)
        .str.strip()
    )

    proteins = proteins.loc[
        proteins["gene_symbol"].ne("")
    ].copy()

    output_columns = [
        "gene_symbol",
        "uniprot_accession",
    ]

    optional_columns = [
        "candidate_tier",
        "aggregated_reproducibility_score",
        "within_feature_type_percentile",
        "within_feature_type_rank",
        "n_evidence_records",
        "n_comparisons",
        "comparisons",
        "effect_directions",
        "direction_agreement",
        "maximum_absolute_effect",
        "median_absolute_effect",
        "mean_direction_consistency",
        "mean_robustness_value",
        "dataset",
        "study_role",
        "evidence_role",
    ]

    output_columns.extend(
        column
        for column in optional_columns
        if column in proteins.columns
    )

    proteins = (
        proteins[output_columns]
        .drop_duplicates(
            subset=["gene_symbol"]
        )
        .sort_values("gene_symbol")
        .reset_index(drop=True)
    )

    if proteins.empty:
        raise ValueError(
            "No eligible protein nodes were found in the node registry."
        )

    return proteins


def build_linked_gene_set(
    edges: pd.DataFrame,
) -> tuple[list[str], dict[str, dict[str, object]]]:
    if edges.empty:
        return [], {}

    required_columns = {
        "source",
        "target",
    }

    missing = required_columns - set(edges.columns)

    if missing:
        raise KeyError(
            "Regulatory edge file is missing required columns: "
            + ", ".join(sorted(missing))
        )

    edges = edges.copy()

    edges["target_gene"] = (
        edges["target"]
        .map(clean_gene_symbol)
    )

    edges = edges.loc[
        edges["target_gene"].ne("")
    ].copy()

    linked_genes = sorted(
        edges["target_gene"]
        .drop_duplicates()
        .tolist()
    )

    evidence_lookup: dict[str, dict[str, object]] = {}

    for gene, group in edges.groupby(
        "target_gene",
        sort=True,
    ):
        record: dict[str, object] = {
            "regulatory_mirnas": ";".join(
                sorted(
                    {
                        str(value).strip()
                        for value in group["source"]
                        if str(value).strip()
                    }
                )
            ),
            "n_regulatory_edges": len(group),
        }

        for column in [
            "maximum_evidence_strength",
            "evidence_databases",
            "n_databases",
            "publication_ids",
            "n_publications",
            "experimental_methods",
            "n_experimental_methods",
            "support_types",
            "n_source_records",
        ]:
            if column not in group.columns:
                continue

            if column.startswith("n_"):
                numeric = pd.to_numeric(
                    group[column],
                    errors="coerce",
                )

                record[column] = (
                    numeric.max()
                    if numeric.notna().any()
                    else 0
                )

            else:
                values: set[str] = set()

                for value in group[column]:
                    text = str(value).strip()

                    if not text:
                        continue

                    for part in text.split(";"):
                        part = part.strip()

                        if part:
                            values.add(part)

                record[column] = ";".join(
                    sorted(values)
                )

        evidence_lookup[gene] = record

    return linked_genes, evidence_lookup


def build_membership_table(
    candidate_panel: pd.DataFrame,
    linked_genes: list[str],
    evidence_lookup: dict[str, dict[str, object]],
) -> pd.DataFrame:
    linked_set = set(linked_genes)

    membership = candidate_panel.copy()

    membership["in_candidate_panel"] = True

    membership["mir143_linked"] = (
        membership["gene_symbol"]
        .isin(linked_set)
    )

    membership["nonlinked_candidate"] = (
        ~membership["mir143_linked"]
    )

    membership["gene_sets"] = membership.apply(
        lambda row: ";".join(
            [
                "candidate_panel",
                *(
                    ["mir143_linked"]
                    if row["mir143_linked"]
                    else ["nonlinked_candidate"]
                ),
            ]
        ),
        axis=1,
    )

    evidence_columns = [
        "regulatory_mirnas",
        "n_regulatory_edges",
        "maximum_evidence_strength",
        "evidence_databases",
        "n_databases",
        "publication_ids",
        "n_publications",
        "experimental_methods",
        "n_experimental_methods",
        "support_types",
        "n_source_records",
    ]

    for column in evidence_columns:
        membership[column] = membership[
            "gene_symbol"
        ].map(
            lambda gene: evidence_lookup.get(
                gene,
                {},
            ).get(
                column,
                "",
            )
        )

    membership = membership.sort_values(
        [
            "mir143_linked",
            "gene_symbol",
        ],
        ascending=[
            False,
            True,
        ],
    ).reset_index(drop=True)

    return membership


def build_summary(
    candidate_genes: list[str],
    linked_genes: list[str],
    nonlinked_genes: list[str],
) -> pd.DataFrame:
    n_candidate = len(candidate_genes)
    n_linked = len(linked_genes)
    n_nonlinked = len(nonlinked_genes)

    return pd.DataFrame(
        {
            "gene_set": [
                "candidate_panel",
                "mir143_linked",
                "nonlinked_candidate",
            ],
            "n_genes": [
                n_candidate,
                n_linked,
                n_nonlinked,
            ],
            "fraction_of_candidate_panel": [
                1.0 if n_candidate else 0.0,
                (
                    n_linked / n_candidate
                    if n_candidate
                    else 0.0
                ),
                (
                    n_nonlinked / n_candidate
                    if n_candidate
                    else 0.0
                ),
            ],
            "description": [
                (
                    "All network-eligible candidate proteins from "
                    "the integrated evidence registry."
                ),
                (
                    "Candidate proteins with experimentally supported "
                    "hsa-miR-143-3p interactions in multiMiR."
                ),
                (
                    "Candidate proteins without an experimentally "
                    "supported hsa-miR-143-3p interaction in the "
                    "current multiMiR query."
                ),
            ],
        }
    )


def write_report(
    path: Path,
    *,
    node_registry_path: Path,
    edge_file_path: Path,
    candidate_genes: list[str],
    linked_genes: list[str],
    nonlinked_genes: list[str],
) -> None:
    lines = [
        "Gene-set preparation report",
        "===========================",
        "",
        f"Node registry: {node_registry_path}",
        f"Regulatory edge file: {edge_file_path}",
        "",
        f"Candidate panel genes: {len(candidate_genes)}",
        f"miR-143-3p-linked genes: {len(linked_genes)}",
        f"Nonlinked candidate genes: {len(nonlinked_genes)}",
        (
            "Linked fraction of candidate panel: "
            f"{len(linked_genes) / len(candidate_genes):.4f}"
            if candidate_genes
            else "Linked fraction of candidate panel: 0.0000"
        ),
        "",
        "miR-143-3p-linked genes:",
    ]

    if linked_genes:
        lines.extend(
            f"  {gene}"
            for gene in linked_genes
        )
    else:
        lines.append("  None")

    lines.extend(
        [
            "",
            "Interpretation policy:",
            (
                "  The linked gene set contains experimentally "
                "supported multiMiR interactions."
            ),
            (
                "  Supportive CLIP-based evidence is retained but "
                "must not be described as definitive functional "
                "regulation."
            ),
            (
                "  Absence from the linked set does not prove absence "
                "of biological regulation."
            ),
            (
                "  The candidate panel is the primary gene set for "
                "pathway enrichment."
            ),
            (
                "  The linked subset should be interpreted as an "
                "exploratory mechanistic subset because of its small "
                "size."
            ),
            "",
        ]
    )

    path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()

    try:
        node_registry = read_tsv(
            args.node_registry
        )

        edges = read_tsv(
            args.edge_file
        )

        candidate_panel = build_candidate_panel(
            node_registry
        )

        linked_genes, evidence_lookup = build_linked_gene_set(
            edges
        )

        candidate_genes = sorted(
            candidate_panel["gene_symbol"]
            .drop_duplicates()
            .tolist()
        )

        candidate_set = set(candidate_genes)

        linked_genes = sorted(
            set(linked_genes) & candidate_set
        )

        nonlinked_genes = sorted(
            candidate_set - set(linked_genes)
        )

        membership = build_membership_table(
            candidate_panel,
            linked_genes,
            evidence_lookup,
        )

        summary = build_summary(
            candidate_genes,
            linked_genes,
            nonlinked_genes,
        )

        args.output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        candidate_path = (
            args.output_dir
            / "candidate_panel_genes.txt"
        )

        linked_path = (
            args.output_dir
            / "mir143_linked_genes.txt"
        )

        nonlinked_path = (
            args.output_dir
            / "nonlinked_candidate_genes.txt"
        )

        mapping_path = (
            args.output_dir
            / "candidate_panel_mapping.tsv"
        )

        membership_path = (
            args.output_dir
            / "gene_set_membership.tsv"
        )

        summary_path = (
            args.output_dir
            / "gene_set_summary.tsv"
        )

        report_path = (
            args.output_dir
            / "gene_set_report.txt"
        )

        write_gene_list(
            candidate_path,
            candidate_genes,
        )

        write_gene_list(
            linked_path,
            linked_genes,
        )

        write_gene_list(
            nonlinked_path,
            nonlinked_genes,
        )

        candidate_panel.to_csv(
            mapping_path,
            sep="\t",
            index=False,
        )

        membership.to_csv(
            membership_path,
            sep="\t",
            index=False,
        )

        summary.to_csv(
            summary_path,
            sep="\t",
            index=False,
        )

        write_report(
            report_path,
            node_registry_path=args.node_registry,
            edge_file_path=args.edge_file,
            candidate_genes=candidate_genes,
            linked_genes=linked_genes,
            nonlinked_genes=nonlinked_genes,
        )

        print(
            f"Candidate panel genes: {len(candidate_genes)}"
        )
        print(
            f"miR-143-3p-linked genes: {len(linked_genes)}"
        )
        print(
            f"Nonlinked candidate genes: {len(nonlinked_genes)}"
        )
        print(
            f"Wrote candidate genes: {candidate_path}"
        )
        print(
            f"Wrote linked genes: {linked_path}"
        )
        print(
            f"Wrote nonlinked genes: {nonlinked_path}"
        )
        print(
            f"Wrote mapping table: {mapping_path}"
        )
        print(
            f"Wrote membership table: {membership_path}"
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

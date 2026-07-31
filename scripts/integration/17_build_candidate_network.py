#!/usr/bin/env python3

"""
Build a multilayer candidate–regulatory–pathway network.

Inputs
------
1. Validated miRNA–target regulatory edges
2. Significant pathway enrichment results
3. Optional candidate panel / gene-set files

Outputs
-------
- candidate_network_nodes.tsv
- candidate_network_edges.tsv
- candidate_network.graphml
- candidate_network.gml
- candidate_network_summary.tsv
- candidate_network_report.txt

The network contains:
- miRNA nodes
- protein/gene nodes
- pathway nodes
- validated miRNA-target edges
- protein-pathway membership edges

Interpretation policy
---------------------
- miRNA-target edges represent validated database support, not proof of
  functional regulation in saliva or periodontitis.
- pathway enrichment represents over-representation, not pathway activation.
- very small gene sets, especially the miR-143-linked subset, are exploratory.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable

import networkx as nx
import pandas as pd


DEFAULT_REGULATORY_EDGE_CANDIDATES = [
    Path("results/integration/biological_network/regulatory_edges.tsv"),
    Path("results/integration/regulatory_edges.tsv"),
    Path("results/integration/network/regulatory_edges.tsv"),
    Path("results/integration/multimir/regulatory_edges.tsv"),
    Path("results/integration/mirna_target_edges.tsv"),
]

DEFAULT_ENRICHMENT_PATH = Path(
    "results/integration/pathway_enrichment/tables/"
    "significant_pathway_enrichment_results.tsv"
)

DEFAULT_OUTPUT_DIR = Path(
    "results/integration/candidate_network"
)

PRIMARY_PATHWAYS = {
    "Antimicrobial peptides",
    "Signaling by Rho GTPases",
    "Signaling by Rho GTPases, Miro GTPases and RHOBTB3",
    "Cell-Cell communication",
    "RHO GTPases activate PAKs",
    "RHO GTPase cycle",
}

LOW_SPECIFICITY_OR_CONTEXT_MISMATCH_PATHWAYS = {
    "Sensory processing of sound by outer hair cells of the cochlea",
    "Sensory processing of sound by inner hair cells of the cochlea",
    "Sensory processing of sound",
    "Sensory Perception",
    "NCAM signaling for neurite out-growth",
    "Nephrin family interactions",
    "RHOU GTPase cycle",
    "RHOV GTPase cycle",
    "L1CAM interactions",
    "Signal Transduction",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build candidate miRNA–protein–pathway network."
    )
    parser.add_argument(
        "--regulatory-edges",
        type=Path,
        default=None,
        help=(
            "Validated miRNA-target edge table. If omitted, the script "
            "searches known output locations."
        ),
    )
    parser.add_argument(
        "--enrichment",
        type=Path,
        default=DEFAULT_ENRICHMENT_PATH,
        help="Significant pathway enrichment table.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Output directory.",
    )
    parser.add_argument(
        "--include-all-significant-pathways",
        action="store_true",
        help=(
            "Include every significant pathway. By default, all significant "
            "pathways are retained but context-mismatched pathways are marked "
            "as low interpretation priority."
        ),
    )
    return parser.parse_args()


def normalize_column_name(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def normalize_dataframe_columns(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()
    result.columns = [normalize_column_name(str(col)) for col in result.columns]
    return result


def first_existing_column(
    columns: Iterable[str],
    candidates: Iterable[str],
) -> str | None:
    column_set = set(columns)
    for candidate in candidates:
        if candidate in column_set:
            return candidate
    return None


def clean_string(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def clean_gene_symbol(value: object) -> str:
    symbol = clean_string(value)
    return symbol.upper()


def split_gene_symbols(value: object) -> list[str]:
    text = clean_string(value)
    if not text:
        return []

    genes = re.split(r"[;,/|]\s*", text)
    cleaned = []

    for gene in genes:
        gene = clean_gene_symbol(gene)
        if gene and gene not in cleaned:
            cleaned.append(gene)

    return cleaned


def safe_float(value: object) -> float | None:
    try:
        if pd.isna(value):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def find_regulatory_edge_file(
    explicitly_requested: Path | None,
) -> Path:
    if explicitly_requested is not None:
        if not explicitly_requested.exists():
            raise FileNotFoundError(
                f"Regulatory edge file not found: {explicitly_requested}"
            )
        return explicitly_requested

    for candidate in DEFAULT_REGULATORY_EDGE_CANDIDATES:
        if candidate.exists():
            return candidate

    search_patterns = [
        "results/**/*regulatory*edge*.tsv",
        "results/**/*mirna*target*.tsv",
        "results/**/*validated*target*.tsv",
    ]

    matches: list[Path] = []
    for pattern in search_patterns:
        matches.extend(Path(".").glob(pattern))

    matches = sorted(
        {
            path
            for path in matches
            if path.is_file()
            and "candidate_network" not in str(path)
        }
    )

    if not matches:
        searched = "\n".join(
            f"  - {path}" for path in DEFAULT_REGULATORY_EDGE_CANDIDATES
        )
        raise FileNotFoundError(
            "Could not locate the validated regulatory edge table.\n"
            "Checked the following standard paths:\n"
            f"{searched}\n\n"
            "Run with an explicit path, for example:\n"
            "python scripts/integration/17_build_candidate_network.py "
            "--regulatory-edges path/to/regulatory_edges.tsv"
        )

    print(
        "Automatically selected regulatory edge table:",
        matches[0],
    )
    return matches[0]


def read_regulatory_edges(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = normalize_dataframe_columns(df)

    mirna_col = first_existing_column(
        df.columns,
        [
            "mirna",
            "mirna_name",
            "mature_mirna",
            "mirna_id",
            "source_mirna",
            "source",
        ],
    )

    target_col = first_existing_column(
        df.columns,
        [
            "target_gene",
            "target_symbol",
            "gene_symbol",
            "target",
            "protein",
            "target_gene_symbol",
        ],
    )

    evidence_col = first_existing_column(
        df.columns,
        [
            "evidence",
            "evidence_type",
            "support",
            "validation",
            "interaction_support",
        ],
    )

    database_col = first_existing_column(
        df.columns,
        [
            "database",
            "source_database",
            "resource",
            "databases",
        ],
    )

    if mirna_col is None or target_col is None:
        raise ValueError(
            "Could not identify miRNA and target-gene columns in "
            f"{path}.\nAvailable columns:\n  "
            + "\n  ".join(df.columns)
        )

    result = pd.DataFrame(
        {
            "mirna": df[mirna_col].map(clean_string),
            "target_gene": df[target_col].map(clean_gene_symbol),
        }
    )

    if evidence_col is not None:
        result["evidence"] = df[evidence_col].map(clean_string)
    else:
        result["evidence"] = "validated_database_support"

    if database_col is not None:
        result["database"] = df[database_col].map(clean_string)
    else:
        result["database"] = "multiMiR"

    result = result.loc[
        result["mirna"].ne("")
        & result["target_gene"].ne("")
    ].copy()

    result["interaction"] = "validated_target"
    result["support_level"] = "supportive"
    result["functional_validation_in_study"] = False

    result = result.drop_duplicates(
        subset=["mirna", "target_gene", "database"]
    ).reset_index(drop=True)

    return result


def classify_pathway_priority(
    description: str,
    gene_set: str,
) -> str:
    if description in PRIMARY_PATHWAYS:
        if gene_set == "candidate_panel":
            return "primary"
        if gene_set == "nonlinked_candidate":
            return "primary_supportive"
        return "exploratory_relevant"

    if description in LOW_SPECIFICITY_OR_CONTEXT_MISMATCH_PATHWAYS:
        return "low_interpretation_priority"

    if gene_set == "mir143_linked":
        return "exploratory"

    return "supportive"


def classify_biological_module(
    description: str,
) -> str:
    lower = description.lower()

    if "antimicrobial" in lower:
        return "innate_antimicrobial_defence"

    if (
        "rho gtpase" in lower
        or "cytoskeleton" in lower
        or "cell-cell communication" in lower
        or "ncam" in lower
        or "l1cam" in lower
        or "nephrin" in lower
    ):
        return "cytoskeleton_cell_architecture"

    if "sensory" in lower or "cochlea" in lower:
        return "shared_structural_annotation"

    return "other_pathway"


def read_enrichment_results(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Enrichment table not found: {path}"
        )

    df = pd.read_csv(path, sep="\t")
    df = normalize_dataframe_columns(df)

    required = {
        "gene_set",
        "database",
        "description",
        "gene_symbols",
    }

    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(
            "Enrichment table is missing required columns: "
            + ", ".join(missing)
            + "\nAvailable columns:\n  "
            + "\n  ".join(df.columns)
        )

    p_adjust_col = first_existing_column(
        df.columns,
        ["p_adjust", "p_adj", "adjusted_p_value", "fdr"],
    )

    fold_col = first_existing_column(
        df.columns,
        ["fold_enrichment", "foldenrichment"],
    )

    id_col = first_existing_column(
        df.columns,
        ["id", "pathway_id", "term_id"],
    )

    result = pd.DataFrame(
        {
            "gene_set": df["gene_set"].map(clean_string),
            "database": df["database"].map(clean_string),
            "description": df["description"].map(clean_string),
            "gene_symbols": df["gene_symbols"].map(clean_string),
        }
    )

    result["pathway_id"] = (
        df[id_col].map(clean_string)
        if id_col is not None
        else ""
    )

    result["adjusted_p_value"] = (
        df[p_adjust_col].map(safe_float)
        if p_adjust_col is not None
        else None
    )

    result["fold_enrichment"] = (
        df[fold_col].map(safe_float)
        if fold_col is not None
        else None
    )

    result["interpretation_priority"] = result.apply(
        lambda row: classify_pathway_priority(
            row["description"],
            row["gene_set"],
        ),
        axis=1,
    )

    result["biological_module"] = result["description"].map(
        classify_biological_module
    )

    return result


def make_node_id(node_type: str, label: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.:-]+", "_", label.strip())
    return f"{node_type}:{cleaned}"


def add_node(
    graph: nx.Graph,
    node_registry: dict[str, dict[str, object]],
    node_id: str,
    **attributes: object,
) -> None:
    cleaned_attributes = {}

    for key, value in attributes.items():
        if value is None:
            cleaned_attributes[key] = ""
        elif isinstance(value, bool):
            cleaned_attributes[key] = bool(value)
        elif isinstance(value, (int, float)):
            cleaned_attributes[key] = value
        else:
            cleaned_attributes[key] = str(value)

    if node_id not in node_registry:
        node_registry[node_id] = {
            "node_id": node_id,
            **cleaned_attributes,
        }
        graph.add_node(node_id, **cleaned_attributes)
        return

    for key, value in cleaned_attributes.items():
        old_value = node_registry[node_id].get(key, "")

        if old_value in ("", None, False) and value not in ("", None):
            node_registry[node_id][key] = value
            graph.nodes[node_id][key] = value


def add_edge(
    graph: nx.MultiDiGraph,
    edge_records: list[dict[str, object]],
    source: str,
    target: str,
    interaction: str,
    **attributes: object,
) -> None:
    record = {
        "source": source,
        "target": target,
        "interaction": interaction,
        **attributes,
    }

    edge_key = (
        source,
        target,
        interaction,
        str(attributes.get("database", "")),
        str(attributes.get("gene_set", "")),
    )

    existing_keys = {
        (
            str(row["source"]),
            str(row["target"]),
            str(row["interaction"]),
            str(row.get("database", "")),
            str(row.get("gene_set", "")),
        )
        for row in edge_records
    }

    if edge_key in existing_keys:
        return

    graph_attributes = {}

    for key, value in {
        "interaction": interaction,
        **attributes,
    }.items():
        if value is None:
            graph_attributes[key] = ""
        elif isinstance(value, bool):
            graph_attributes[key] = bool(value)
        elif isinstance(value, (int, float)):
            graph_attributes[key] = value
        else:
            graph_attributes[key] = str(value)

    graph.add_edge(
        source,
        target,
        **graph_attributes,
    )

    edge_records.append(record)


def build_network(
    regulatory_edges: pd.DataFrame,
    enrichment_results: pd.DataFrame,
) -> tuple[nx.MultiDiGraph, pd.DataFrame, pd.DataFrame]:
    graph = nx.MultiDiGraph()
    node_registry: dict[str, dict[str, object]] = {}
    edge_records: list[dict[str, object]] = []

    regulatory_targets = set(
        regulatory_edges["target_gene"]
    )

    pathway_genes = set()
    for value in enrichment_results["gene_symbols"]:
        pathway_genes.update(split_gene_symbols(value))

    all_proteins = sorted(regulatory_targets | pathway_genes)

    for protein in all_proteins:
        node_id = make_node_id("protein", protein)

        add_node(
            graph,
            node_registry,
            node_id,
            label=protein,
            node_type="protein",
            molecular_entity="gene_or_protein",
            is_regulatory_target=protein in regulatory_targets,
            is_pathway_member=protein in pathway_genes,
            biological_module="unassigned",
        )

    for _, row in regulatory_edges.iterrows():
        mirna = row["mirna"]
        target = row["target_gene"]

        mirna_id = make_node_id("mirna", mirna)
        target_id = make_node_id("protein", target)

        add_node(
            graph,
            node_registry,
            mirna_id,
            label=mirna,
            node_type="miRNA",
            molecular_entity="mature_miRNA",
            is_regulatory_target=False,
            is_pathway_member=False,
            biological_module="regulatory",
        )

        add_edge(
            graph,
            edge_records,
            mirna_id,
            target_id,
            "validated_target",
            database=row["database"],
            evidence=row["evidence"],
            support_level=row["support_level"],
            functional_validation_in_study=False,
            gene_set="mir143_linked",
        )

    for _, row in enrichment_results.iterrows():
        pathway_label = row["description"]
        pathway_db = row["database"]
        pathway_term_id = row["pathway_id"]
        gene_set = row["gene_set"]
        priority = row["interpretation_priority"]
        module = row["biological_module"]

        pathway_key = (
            pathway_term_id
            if pathway_term_id
            else pathway_label
        )

        pathway_node_id = make_node_id(
            "pathway",
            f"{pathway_db}_{pathway_key}",
        )

        add_node(
            graph,
            node_registry,
            pathway_node_id,
            label=pathway_label,
            node_type="pathway",
            molecular_entity="functional_pathway",
            database=pathway_db,
            pathway_id=pathway_term_id,
            gene_set=gene_set,
            adjusted_p_value=row["adjusted_p_value"],
            fold_enrichment=row["fold_enrichment"],
            interpretation_priority=priority,
            biological_module=module,
            is_regulatory_target=False,
            is_pathway_member=False,
        )

        genes = split_gene_symbols(row["gene_symbols"])

        for gene in genes:
            protein_id = make_node_id("protein", gene)

            add_node(
                graph,
                node_registry,
                protein_id,
                label=gene,
                node_type="protein",
                molecular_entity="gene_or_protein",
                is_regulatory_target=gene in regulatory_targets,
                is_pathway_member=True,
                biological_module=module,
            )

            add_edge(
                graph,
                edge_records,
                protein_id,
                pathway_node_id,
                "member_of_pathway",
                database=pathway_db,
                pathway_id=pathway_term_id,
                gene_set=gene_set,
                adjusted_p_value=row["adjusted_p_value"],
                fold_enrichment=row["fold_enrichment"],
                interpretation_priority=priority,
                biological_module=module,
            )

    nodes_df = pd.DataFrame(node_registry.values())

    if not nodes_df.empty:
        node_degree = dict(graph.degree())
        in_degree = dict(graph.in_degree())
        out_degree = dict(graph.out_degree())

        nodes_df["degree"] = nodes_df["node_id"].map(node_degree).fillna(0)
        nodes_df["in_degree"] = nodes_df["node_id"].map(in_degree).fillna(0)
        nodes_df["out_degree"] = nodes_df["node_id"].map(out_degree).fillna(0)

        nodes_df = nodes_df.sort_values(
            ["node_type", "label"],
            kind="stable",
        ).reset_index(drop=True)

    edges_df = pd.DataFrame(edge_records)

    if not edges_df.empty:
        edges_df = edges_df.sort_values(
            ["interaction", "source", "target"],
            kind="stable",
        ).reset_index(drop=True)

        edges_df.insert(
            0,
            "edge_id",
            [
                f"edge_{index:04d}"
                for index in range(1, len(edges_df) + 1)
            ],
        )

    return graph, nodes_df, edges_df


def sanitize_graph_for_export(
    graph: nx.MultiDiGraph,
) -> nx.MultiDiGraph:
    clean_graph = nx.MultiDiGraph()

    for node_id, attrs in graph.nodes(data=True):
        clean_attrs = {}
        for key, value in attrs.items():
            if value is None:
                clean_attrs[key] = ""
            elif isinstance(value, bool):
                clean_attrs[key] = int(value)
            elif isinstance(value, (int, float, str)):
                clean_attrs[key] = value
            else:
                clean_attrs[key] = str(value)

        clean_graph.add_node(node_id, **clean_attrs)

    for source, target, key, attrs in graph.edges(
        keys=True,
        data=True,
    ):
        clean_attrs = {}
        for attr_name, value in attrs.items():
            if value is None:
                clean_attrs[attr_name] = ""
            elif isinstance(value, bool):
                clean_attrs[attr_name] = int(value)
            elif isinstance(value, (int, float, str)):
                clean_attrs[attr_name] = value
            else:
                clean_attrs[attr_name] = str(value)

        clean_graph.add_edge(
            source,
            target,
            key=key,
            **clean_attrs,
        )

    return clean_graph


def write_summary(
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
    output_path: Path,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    rows.append(
        {
            "category": "network",
            "metric": "total_nodes",
            "value": len(nodes_df),
        }
    )

    rows.append(
        {
            "category": "network",
            "metric": "total_edges",
            "value": len(edges_df),
        }
    )

    if not nodes_df.empty:
        for node_type, count in (
            nodes_df["node_type"]
            .value_counts()
            .sort_index()
            .items()
        ):
            rows.append(
                {
                    "category": "node_type",
                    "metric": node_type,
                    "value": int(count),
                }
            )

    if not edges_df.empty:
        for interaction, count in (
            edges_df["interaction"]
            .value_counts()
            .sort_index()
            .items()
        ):
            rows.append(
                {
                    "category": "edge_type",
                    "metric": interaction,
                    "value": int(count),
                }
            )

    if (
        not nodes_df.empty
        and "interpretation_priority" in nodes_df.columns
    ):
        pathway_nodes = nodes_df.loc[
            nodes_df["node_type"].eq("pathway")
        ]

        for priority, count in (
            pathway_nodes["interpretation_priority"]
            .fillna("")
            .value_counts()
            .sort_index()
            .items()
        ):
            if priority:
                rows.append(
                    {
                        "category": "pathway_priority",
                        "metric": priority,
                        "value": int(count),
                    }
                )

    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(output_path, sep="\t", index=False)

    return summary_df


def write_report(
    output_path: Path,
    regulatory_path: Path,
    enrichment_path: Path,
    graph: nx.MultiDiGraph,
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
) -> None:
    node_counts = (
        nodes_df["node_type"].value_counts().to_dict()
        if not nodes_df.empty
        else {}
    )

    edge_counts = (
        edges_df["interaction"].value_counts().to_dict()
        if not edges_df.empty
        else {}
    )

    pathway_nodes = (
        nodes_df.loc[nodes_df["node_type"].eq("pathway")].copy()
        if not nodes_df.empty
        else pd.DataFrame()
    )

    lines = [
        "Candidate network report",
        "========================",
        "",
        f"Regulatory edge input: {regulatory_path}",
        f"Pathway enrichment input: {enrichment_path}",
        "",
        "Network dimensions:",
        f"  Nodes: {graph.number_of_nodes()}",
        f"  Edges: {graph.number_of_edges()}",
        "",
        "Node composition:",
    ]

    for node_type, count in sorted(node_counts.items()):
        lines.append(f"  {node_type}: {count}")

    lines.extend(
        [
            "",
            "Edge composition:",
        ]
    )

    for edge_type, count in sorted(edge_counts.items()):
        lines.append(f"  {edge_type}: {count}")

    lines.extend(
        [
            "",
            "Biological interpretation:",
            "  The primary proteomic signal is the Reactome "
            "'Antimicrobial peptides' pathway.",
            "  The miR-143-linked subset is enriched for overlapping "
            "cytoskeletal and Rho-GTPase-related annotations.",
            "  Cochlear, sensory, neurite, nephrin and related pathway "
            "labels are driven by shared structural proteins and should "
            "not be interpreted as tissue-specific findings.",
            "",
            "Interpretation constraints:",
            "  miRNA-target edges indicate validated database support.",
            "  They do not establish functional miRNA regulation in the "
            "analysed saliva samples.",
            "  Pathway enrichment indicates over-representation, not "
            "pathway activation, direction or causality.",
            "  The miR-143-linked set is small and exploratory.",
            "  The nonlinked candidate set is descriptive and is not an "
            "independent control group.",
            "",
            "Pathways included:",
        ]
    )

    if pathway_nodes.empty:
        lines.append("  None")
    else:
        pathway_nodes = pathway_nodes.sort_values(
            [
                "interpretation_priority",
                "adjusted_p_value",
                "label",
            ],
            kind="stable",
        )

        for _, row in pathway_nodes.iterrows():
            fdr = row.get("adjusted_p_value", "")
            fold = row.get("fold_enrichment", "")
            priority = row.get("interpretation_priority", "")
            gene_set = row.get("gene_set", "")
            database = row.get("database", "")
            label = row.get("label", "")

            lines.append(
                f"  {database} | {gene_set} | {label} | "
                f"FDR={fdr} | fold={fold} | priority={priority}"
            )

    lines.extend(
        [
            "",
            "Main outputs:",
            "  candidate_network_nodes.tsv",
            "  candidate_network_edges.tsv",
            "  candidate_network.graphml",
            "  candidate_network.gml",
            "  candidate_network_summary.tsv",
        ]
    )

    output_path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()

    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    regulatory_path = find_regulatory_edge_file(
        args.regulatory_edges
    )

    regulatory_edges = read_regulatory_edges(
        regulatory_path
    )

    enrichment_results = read_enrichment_results(
        args.enrichment
    )

    graph, nodes_df, edges_df = build_network(
        regulatory_edges,
        enrichment_results,
    )

    nodes_path = output_dir / "candidate_network_nodes.tsv"
    edges_path = output_dir / "candidate_network_edges.tsv"
    graphml_path = output_dir / "candidate_network.graphml"
    gml_path = output_dir / "candidate_network.gml"
    summary_path = output_dir / "candidate_network_summary.tsv"
    report_path = output_dir / "candidate_network_report.txt"

    nodes_df.to_csv(
        nodes_path,
        sep="\t",
        index=False,
    )

    edges_df.to_csv(
        edges_path,
        sep="\t",
        index=False,
    )

    clean_graph = sanitize_graph_for_export(graph)

    nx.write_graphml(
        clean_graph,
        graphml_path,
    )

    nx.write_gml(
        clean_graph,
        gml_path,
    )

    write_summary(
        nodes_df,
        edges_df,
        summary_path,
    )

    write_report(
        report_path,
        regulatory_path,
        args.enrichment,
        graph,
        nodes_df,
        edges_df,
    )

    print("Candidate network completed.")
    print(f"Regulatory edge input: {regulatory_path}")
    print(f"Pathway enrichment input: {args.enrichment}")
    print(f"Nodes: {graph.number_of_nodes()}")
    print(f"Edges: {graph.number_of_edges()}")
    print(f"Wrote nodes: {nodes_path}")
    print(f"Wrote edges: {edges_path}")
    print(f"Wrote GraphML: {graphml_path}")
    print(f"Wrote GML: {gml_path}")
    print(f"Wrote summary: {summary_path}")
    print(f"Wrote report: {report_path}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(
            f"ERROR: {exc}",
            file=sys.stderr,
        )
        raise

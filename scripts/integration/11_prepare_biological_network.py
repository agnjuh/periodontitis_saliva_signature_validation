#!/usr/bin/env python3
"""
Prepare evidence-aware biological-network inputs from the final candidate panel.

This script does not invent miRNA-target or protein-module relationships.
Instead, it:

1. creates a harmonised node registry;
2. distinguishes biological evidence from technical method-sensitivity evidence;
3. identifies features eligible for regulatory network integration;
4. creates auditable mapping templates for validated external relationships;
5. reports network-readiness and evidence coverage.

The next network-building step should consume manually reviewed or
database-derived mapping tables.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = [
    "feature_id",
    "feature_name",
    "feature_type",
    "modality",
    "dataset",
    "study_role",
    "n_evidence_records",
    "n_comparisons",
    "comparisons",
    "effect_directions",
    "direction_agreement",
    "maximum_absolute_effect",
    "median_absolute_effect",
    "mean_direction_consistency",
    "mean_robustness_value",
    "aggregated_reproducibility_score",
    "within_feature_type_percentile",
    "within_feature_type_rank",
    "candidate_tier",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare biological-network nodes and auditable "
            "relationship templates."
        )
    )

    parser.add_argument(
        "--panel",
        type=Path,
        default=Path(
            "results/integration/cross_cohort_validation/"
            "final_candidate_panel.tsv"
        ),
        help="Feature-level final candidate panel.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "results/integration/biological_network"
        ),
        help="Directory for biological-network preparation outputs.",
    )

    return parser.parse_args()


def require_columns(frame: pd.DataFrame) -> None:
    missing = sorted(
        set(REQUIRED_COLUMNS) - set(frame.columns)
    )

    if missing:
        raise ValueError(
            "Final candidate panel is missing required columns: "
            + ", ".join(missing)
        )


def normalise_text(value: object) -> str:
    if pd.isna(value):
        return ""

    return str(value).strip()


def make_node_id(row: pd.Series) -> str:
    feature_type = normalise_text(row["feature_type"])
    feature_id = normalise_text(row["feature_id"])

    safe_type = (
        feature_type.lower()
        .replace(" ", "_")
        .replace("/", "_")
    )

    return f"{safe_type}:{feature_id}"


def assign_node_class(feature_type: str) -> str:
    mapping = {
        "protein": "molecular_feature",
        "protein_module": "biological_module",
        "mature_miRNA": "regulatory_rna",
        "miRNA_associated_locus": "technical_feature",
        "probable_miRNA_host_gene": (
            "technical_or_regulatory_locus"
        ),
    }

    return mapping.get(
        feature_type,
        "other_feature",
    )


def assign_evidence_role(row: pd.Series) -> str:
    dataset = normalise_text(row["dataset"])
    study_role = normalise_text(row["study_role"])
    feature_type = normalise_text(row["feature_type"])

    if (
        dataset == "PXD043901"
        and feature_type == "protein"
    ):
        return "longitudinal_biological_response"

    if (
        dataset == "PXD043901"
        and feature_type == "protein_module"
    ):
        return "longitudinal_functional_response"

    if (
        dataset == "GSE121910"
        and feature_type == "mature_miRNA"
    ):
        return "case_control_biological_signal"

    if dataset == "GSE312585":
        return "technical_method_sensitivity"

    if "technical" in study_role.lower():
        return "technical_method_sensitivity"

    return "other_evidence"


def assign_network_eligibility(
    row: pd.Series,
) -> bool:
    evidence_role = normalise_text(
        row["evidence_role"]
    )
    feature_type = normalise_text(
        row["feature_type"]
    )

    biologically_eligible_roles = {
        "longitudinal_biological_response",
        "longitudinal_functional_response",
        "case_control_biological_signal",
    }

    eligible_feature_types = {
        "protein",
        "protein_module",
        "mature_miRNA",
    }

    return bool(
        evidence_role in biologically_eligible_roles
        and feature_type in eligible_feature_types
    )


def assign_regulatory_source_eligibility(
    row: pd.Series,
) -> bool:
    return bool(
        row["feature_type"] == "mature_miRNA"
        and row["evidence_role"]
        == "case_control_biological_signal"
    )


def assign_regulatory_target_eligibility(
    row: pd.Series,
) -> bool:
    return bool(
        row["feature_type"] == "protein"
        and row["evidence_role"]
        == "longitudinal_biological_response"
    )


def assign_module_mapping_eligibility(
    row: pd.Series,
) -> bool:
    return bool(
        row["feature_type"]
        in {"protein", "protein_module"}
        and row["dataset"] == "PXD043901"
    )


def build_node_registry(
    panel: pd.DataFrame,
) -> pd.DataFrame:
    nodes = panel.copy()

    nodes.insert(
        0,
        "node_id",
        nodes.apply(
            make_node_id,
            axis=1,
        ),
    )

    nodes["node_class"] = nodes[
        "feature_type"
    ].map(assign_node_class)

    nodes["evidence_role"] = nodes.apply(
        assign_evidence_role,
        axis=1,
    )

    nodes["network_eligible"] = nodes.apply(
        assign_network_eligibility,
        axis=1,
    )

    nodes[
        "regulatory_source_eligible"
    ] = nodes.apply(
        assign_regulatory_source_eligibility,
        axis=1,
    )

    nodes[
        "regulatory_target_eligible"
    ] = nodes.apply(
        assign_regulatory_target_eligibility,
        axis=1,
    )

    nodes[
        "module_mapping_eligible"
    ] = nodes.apply(
        assign_module_mapping_eligibility,
        axis=1,
    )

    selected_columns = [
        "node_id",
        "feature_id",
        "feature_name",
        "feature_type",
        "node_class",
        "modality",
        "dataset",
        "study_role",
        "evidence_role",
        "network_eligible",
        "regulatory_source_eligible",
        "regulatory_target_eligible",
        "module_mapping_eligible",
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
    ]

    return (
        nodes[selected_columns]
        .sort_values(
            [
                "node_class",
                "feature_type",
                "within_feature_type_rank",
                "feature_name",
            ]
        )
        .reset_index(drop=True)
    )


def build_mirna_target_template(
    nodes: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    mirnas = nodes.loc[
        nodes["regulatory_source_eligible"]
    ].copy()

    proteins = nodes.loc[
        nodes["regulatory_target_eligible"]
    ].copy()

    rows: list[dict[str, object]] = []

    for _, mirna in mirnas.iterrows():
        rows.append(
            {
                "mirna_node_id": mirna["node_id"],
                "mirna_id": mirna["feature_id"],
                "mirna_name": mirna["feature_name"],
                "mirna_dataset": mirna["dataset"],
                "mirna_candidate_tier": (
                    mirna["candidate_tier"]
                ),
                "mirna_reproducibility_score": (
                    mirna[
                        "aggregated_reproducibility_score"
                    ]
                ),
                "target_gene_symbol": "",
                "target_protein_node_id": "",
                "target_protein_detected_in_panel": False,
                "relationship_type": "",
                "evidence_database": "",
                "database_version": "",
                "evidence_category": "",
                "experimental_method": "",
                "publication_id": "",
                "interaction_direction": "",
                "curation_status": "not_reviewed",
                "review_notes": "",
            }
        )

    columns = [
        "mirna_node_id",
        "mirna_id",
        "mirna_name",
        "mirna_dataset",
        "mirna_candidate_tier",
        "mirna_reproducibility_score",
        "target_gene_symbol",
        "target_protein_node_id",
        "target_protein_detected_in_panel",
        "relationship_type",
        "evidence_database",
        "database_version",
        "evidence_category",
        "experimental_method",
        "publication_id",
        "interaction_direction",
        "curation_status",
        "review_notes",
    ]

    template = pd.DataFrame(
        rows,
        columns=columns,
    )

    protein_reference = proteins[
        [
            "node_id",
            "feature_id",
            "feature_name",
            "candidate_tier",
            "aggregated_reproducibility_score",
            "effect_directions",
            "n_comparisons",
        ]
    ].copy()

    protein_reference = protein_reference.rename(
        columns={
            "node_id": "target_protein_node_id",
            "feature_id": "protein_feature_id",
            "feature_name": "target_gene_symbol",
            "candidate_tier": (
                "protein_candidate_tier"
            ),
            "aggregated_reproducibility_score": (
                "protein_reproducibility_score"
            ),
            "effect_directions": (
                "protein_effect_directions"
            ),
            "n_comparisons": (
                "protein_n_comparisons"
            ),
        }
    )

    protein_reference = protein_reference.sort_values(
        [
            "protein_reproducibility_score",
            "target_gene_symbol",
        ],
        ascending=[False, True],
    ).reset_index(drop=True)

    return template, protein_reference


def build_protein_module_template(
    nodes: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    proteins = nodes.loc[
        (
            nodes["feature_type"].eq("protein")
            & nodes["module_mapping_eligible"]
        )
    ].copy()

    modules = nodes.loc[
        (
            nodes["feature_type"].eq(
                "protein_module"
            )
            & nodes["module_mapping_eligible"]
        )
    ].copy()

    rows: list[dict[str, object]] = []

    for _, module in modules.iterrows():
        rows.append(
            {
                "module_node_id": (
                    module["node_id"]
                ),
                "module_name": (
                    module["feature_name"]
                ),
                "module_candidate_tier": (
                    module["candidate_tier"]
                ),
                "module_reproducibility_score": (
                    module[
                        "aggregated_reproducibility_score"
                    ]
                ),
                "protein_node_id": "",
                "protein_gene_symbol": "",
                "membership_source": "",
                "membership_definition": "",
                "membership_weight": "",
                "curation_status": "not_reviewed",
                "review_notes": "",
            }
        )

    columns = [
        "module_node_id",
        "module_name",
        "module_candidate_tier",
        "module_reproducibility_score",
        "protein_node_id",
        "protein_gene_symbol",
        "membership_source",
        "membership_definition",
        "membership_weight",
        "curation_status",
        "review_notes",
    ]

    template = pd.DataFrame(
        rows,
        columns=columns,
    )

    protein_reference = proteins[
        [
            "node_id",
            "feature_id",
            "feature_name",
            "candidate_tier",
            "aggregated_reproducibility_score",
            "effect_directions",
            "n_comparisons",
        ]
    ].copy()

    protein_reference = protein_reference.rename(
        columns={
            "node_id": "protein_node_id",
            "feature_id": "protein_feature_id",
            "feature_name": "protein_gene_symbol",
            "candidate_tier": (
                "protein_candidate_tier"
            ),
            "aggregated_reproducibility_score": (
                "protein_reproducibility_score"
            ),
            "effect_directions": (
                "protein_effect_directions"
            ),
            "n_comparisons": (
                "protein_n_comparisons"
            ),
        }
    )

    protein_reference = protein_reference.sort_values(
        [
            "protein_reproducibility_score",
            "protein_gene_symbol",
        ],
        ascending=[False, True],
    ).reset_index(drop=True)

    return template, protein_reference


def build_network_readiness_summary(
    nodes: pd.DataFrame,
) -> pd.DataFrame:
    summary = (
        nodes.groupby(
            [
                "feature_type",
                "node_class",
                "evidence_role",
            ],
            dropna=False,
        )
        .agg(
            n_nodes=(
                "node_id",
                "size",
            ),
            n_network_eligible=(
                "network_eligible",
                "sum",
            ),
            n_regulatory_sources=(
                "regulatory_source_eligible",
                "sum",
            ),
            n_regulatory_targets=(
                "regulatory_target_eligible",
                "sum",
            ),
            n_module_mapping_eligible=(
                "module_mapping_eligible",
                "sum",
            ),
            mean_reproducibility_score=(
                "aggregated_reproducibility_score",
                "mean",
            ),
            median_reproducibility_score=(
                "aggregated_reproducibility_score",
                "median",
            ),
            minimum_reproducibility_score=(
                "aggregated_reproducibility_score",
                "min",
            ),
            maximum_reproducibility_score=(
                "aggregated_reproducibility_score",
                "max",
            ),
        )
        .reset_index()
    )

    score_columns = [
        "mean_reproducibility_score",
        "median_reproducibility_score",
        "minimum_reproducibility_score",
        "maximum_reproducibility_score",
    ]

    summary[score_columns] = (
        summary[score_columns]
        .round(3)
    )

    return (
        summary.sort_values(
            [
                "n_network_eligible",
                "mean_reproducibility_score",
                "feature_type",
            ],
            ascending=[
                False,
                False,
                True,
            ],
        )
        .reset_index(drop=True)
    )


def write_report(
    nodes: pd.DataFrame,
    mirna_template: pd.DataFrame,
    mirna_protein_reference: pd.DataFrame,
    module_template: pd.DataFrame,
    module_protein_reference: pd.DataFrame,
    readiness: pd.DataFrame,
    output_path: Path,
) -> None:
    total_nodes = len(nodes)

    network_eligible = int(
        nodes["network_eligible"].sum()
    )

    regulatory_sources = int(
        nodes[
            "regulatory_source_eligible"
        ].sum()
    )

    regulatory_targets = int(
        nodes[
            "regulatory_target_eligible"
        ].sum()
    )

    module_mapping_eligible = int(
        nodes[
            "module_mapping_eligible"
        ].sum()
    )

    technical_nodes = int(
        nodes["evidence_role"]
        .eq("technical_method_sensitivity")
        .sum()
    )

    biological_nodes = int(
        nodes["evidence_role"]
        .isin(
            [
                "longitudinal_biological_response",
                "longitudinal_functional_response",
                "case_control_biological_signal",
            ]
        )
        .sum()
    )

    lines = [
        "Biological-network preparation",
        "==============================",
        "",
        f"Total candidate nodes: {total_nodes}",
        (
            "Biological evidence nodes: "
            f"{biological_nodes}"
        ),
        (
            "Biological-network eligible nodes: "
            f"{network_eligible}"
        ),
        (
            "Regulatory miRNA source nodes: "
            f"{regulatory_sources}"
        ),
        (
            "Eligible protein target nodes: "
            f"{regulatory_targets}"
        ),
        (
            "Module-mapping eligible nodes: "
            f"{module_mapping_eligible}"
        ),
        (
            "Technical method-sensitivity nodes: "
            f"{technical_nodes}"
        ),
        "",
        "Mapping templates:",
        (
            "  miRNA-target template rows: "
            f"{len(mirna_template)}"
        ),
        (
            "  eligible protein target reference rows: "
            f"{len(mirna_protein_reference)}"
        ),
        (
            "  protein-module template rows: "
            f"{len(module_template)}"
        ),
        (
            "  protein-module candidate reference rows: "
            f"{len(module_protein_reference)}"
        ),
        (
            "  network-readiness summary rows: "
            f"{len(readiness)}"
        ),
        "",
        "Interpretation:",
        (
            "  - GSE121910 mature miRNAs are eligible as biological "
            "regulatory-source nodes."
        ),
        (
            "  - PXD043901 proteins are eligible as experimentally "
            "observed biological target nodes."
        ),
        (
            "  - PXD043901 protein modules are retained as longitudinal "
            "functional-response nodes."
        ),
        (
            "  - GSE312585 loci are retained as technical robustness "
            "nodes but are excluded from disease-regulatory edges."
        ),
        (
            "  - No miRNA-target or protein-module edge is generated "
            "without an explicit reviewed mapping."
        ),
        (
            "  - Validated interaction databases, database versions, "
            "publication identifiers, and curation status must be "
            "recorded before network construction."
        ),
        "",
    ]

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()

    if not args.panel.is_file():
        raise FileNotFoundError(
            "Final candidate panel not found: "
            f"{args.panel}"
        )

    panel = pd.read_csv(
        args.panel,
        sep="\t",
    )

    require_columns(panel)

    if panel.empty:
        raise ValueError(
            "Final candidate panel is empty."
        )

    output_dir = args.output_dir
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    nodes = build_node_registry(panel)

    (
        mirna_template,
        mirna_protein_reference,
    ) = build_mirna_target_template(nodes)

    (
        module_template,
        module_protein_reference,
    ) = build_protein_module_template(nodes)

    readiness = build_network_readiness_summary(
        nodes
    )

    paths = {
        "nodes": (
            output_dir
            / "network_node_registry.tsv"
        ),
        "readiness": (
            output_dir
            / "network_readiness_summary.tsv"
        ),
        "mirna_template": (
            output_dir
            / "mirna_target_mapping_template.tsv"
        ),
        "mirna_protein_reference": (
            output_dir
            / "eligible_protein_target_reference.tsv"
        ),
        "module_template": (
            output_dir
            / "protein_module_mapping_template.tsv"
        ),
        "module_protein_reference": (
            output_dir
            / "protein_module_candidate_reference.tsv"
        ),
        "report": (
            output_dir
            / "network_preparation_report.txt"
        ),
    }

    nodes.to_csv(
        paths["nodes"],
        sep="\t",
        index=False,
    )

    readiness.to_csv(
        paths["readiness"],
        sep="\t",
        index=False,
    )

    mirna_template.to_csv(
        paths["mirna_template"],
        sep="\t",
        index=False,
    )

    mirna_protein_reference.to_csv(
        paths["mirna_protein_reference"],
        sep="\t",
        index=False,
    )

    module_template.to_csv(
        paths["module_template"],
        sep="\t",
        index=False,
    )

    module_protein_reference.to_csv(
        paths["module_protein_reference"],
        sep="\t",
        index=False,
    )

    write_report(
        nodes=nodes,
        mirna_template=mirna_template,
        mirna_protein_reference=(
            mirna_protein_reference
        ),
        module_template=module_template,
        module_protein_reference=(
            module_protein_reference
        ),
        readiness=readiness,
        output_path=paths["report"],
    )

    print(
        f"Wrote node registry: "
        f"{paths['nodes']}"
    )

    print(
        f"Wrote readiness summary: "
        f"{paths['readiness']}"
    )

    print(
        f"Wrote miRNA-target template: "
        f"{paths['mirna_template']}"
    )

    print(
        "Wrote eligible protein reference: "
        f"{paths['mirna_protein_reference']}"
    )

    print(
        "Wrote protein-module template: "
        f"{paths['module_template']}"
    )

    print(
        "Wrote module protein reference: "
        f"{paths['module_protein_reference']}"
    )

    print(
        f"Wrote report: "
        f"{paths['report']}"
    )

    print(
        f"Total nodes: {len(nodes)}"
    )

    print(
        "Biological-network eligible nodes: "
        f"{int(nodes['network_eligible'].sum())}"
    )

    print(
        "Regulatory source nodes: "
        f"{int(nodes['regulatory_source_eligible'].sum())}"
    )

    print(
        "Regulatory target nodes: "
        f"{int(nodes['regulatory_target_eligible'].sum())}"
    )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        FileNotFoundError,
        ValueError,
        pd.errors.EmptyDataError,
        pd.errors.ParserError,
    ) as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        raise SystemExit(1)

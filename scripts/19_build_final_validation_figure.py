#!/usr/bin/env python3
"""
Build the final validation framework figure for the
periodontitis_saliva_signature_validation repository.

The script integrates repository-derived evidence from:

1. Longitudinal protein trajectories
2. Biological module responses
3. Leave-one-protein-out robustness analysis
4. Cross-cohort candidate prioritisation
5. Regulatory-network support
6. Pathway-enrichment support
7. Candidate-network integration

Outputs
-------
results/final_figure/final_validation_framework.svg
results/final_figure/final_validation_framework.pdf
results/final_figure/final_validation_framework.png
results/final_figure/final_validation_framework_caption.txt
results/final_figure/final_validation_framework_report.txt

Usage
-----
Run from the repository root:

    python scripts/19_build_final_validation_figure.py

The script performs strict input and schema validation before generating
any output. It does not reconstruct upstream analyses and does not insert
placeholder evidence.
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.ticker import MaxNLocator


SCRIPT_NAME = "19_build_final_validation_figure.py"
DEFAULT_OUTPUT_DPI = 450



# Repository paths


@dataclass(frozen=True)
class RepositoryPaths:
    """Absolute paths used by the final validation-figure workflow."""

    repository_root: Path
    results_root: Path
    processed_data_root: Path
    output_directory: Path

    protein_response_summary: Path
    robust_human_responses: Path
    persistent_human_responses: Path
    ranked_human_responses: Path
    stringent_human_responses: Path
    human_trajectory_classes: Path
    trajectory_class_counts: Path
    longitudinal_statistical_tests: Path

    module_longitudinal_statistics: Path
    module_coverage: Path
    patient_module_response_deltas: Path
    module_delta_correlations: Path
    leave_one_out_summary: Path
    leave_one_protein_out_results: Path

    aggregated_candidate_evidence: Path
    cross_cohort_evidence_registry: Path
    final_candidate_panel: Path
    prioritised_candidate_panel: Path
    scored_evidence_registry: Path

    regulatory_edges_high_confidence: Path
    regulatory_edge_summary: Path

    pathway_enrichment_summary: Path
    significant_pathway_enrichment_results: Path

    gene_set_summary: Path
    candidate_panel_mapping: Path

    candidate_network_nodes: Path
    candidate_network_edges: Path
    candidate_network_summary: Path

    protein_annotation: Path
    delta_t1_vs_t0: Path
    delta_t6_vs_t0: Path

    output_svg: Path
    output_pdf: Path
    output_png: Path
    output_caption: Path
    output_report: Path
    output_manifest: Path


def build_repository_paths(repository_root: Path) -> RepositoryPaths:
    """Construct all expected repository input and output paths."""

    root = repository_root.expanduser().resolve()
    results = root / "results"
    processed = root / "data" / "processed"
    output = results / "final_figure"

    trajectories = results / "PXD043901" / "trajectories"
    modules = results / "PXD043901" / "modules"
    sensitivity = modules / "sensitivity"

    integration = results / "integration"
    cross_cohort = integration / "cross_cohort_validation"
    regulatory = integration / "regulatory_network"
    pathway = integration / "pathway_enrichment"
    gene_sets = integration / "gene_sets"
    candidate_network = integration / "candidate_network"

    pxd_processed = processed / "PXD043901"

    return RepositoryPaths(
        repository_root=root,
        results_root=results,
        processed_data_root=processed,
        output_directory=output,

        protein_response_summary=(
            trajectories / "protein_response_summary.csv"
        ),
        robust_human_responses=(
            trajectories / "robust_human_responses_preliminary.csv"
        ),
        persistent_human_responses=(
            trajectories / "persistent_human_responses.csv"
        ),
        ranked_human_responses=(
            trajectories / "ranked_human_responses.csv"
        ),
        stringent_human_responses=(
            trajectories / "stringent_human_responses.csv"
        ),
        human_trajectory_classes=(
            trajectories / "human_protein_trajectory_classes.csv"
        ),
        trajectory_class_counts=(
            trajectories / "trajectory_class_counts.csv"
        ),
        longitudinal_statistical_tests=(
            trajectories / "longitudinal_statistical_tests.csv"
        ),

        module_longitudinal_statistics=(
            modules / "module_longitudinal_statistics.csv"
        ),
        module_coverage=(
            modules / "module_coverage.csv"
        ),
        patient_module_response_deltas=(
            modules / "patient_module_response_deltas.csv"
        ),
        module_delta_correlations=(
            sensitivity / "module_delta_correlations.csv"
        ),
        leave_one_out_summary=(
            sensitivity / "leave_one_out_summary.csv"
        ),
        leave_one_protein_out_results=(
            sensitivity / "leave_one_protein_out_results.csv"
        ),

        aggregated_candidate_evidence=(
            cross_cohort / "aggregated_candidate_evidence.tsv"
        ),
        cross_cohort_evidence_registry=(
            cross_cohort / "cross_cohort_evidence_registry.tsv"
        ),
        final_candidate_panel=(
            cross_cohort / "final_candidate_panel.tsv"
        ),
        prioritised_candidate_panel=(
            cross_cohort / "prioritised_candidate_panel.tsv"
        ),
        scored_evidence_registry=(
            cross_cohort / "scored_evidence_registry.tsv"
        ),

        regulatory_edges_high_confidence=(
            regulatory / "regulatory_edges_high_confidence.tsv"
        ),
        regulatory_edge_summary=(
            regulatory / "regulatory_edge_summary.tsv"
        ),

        pathway_enrichment_summary=(
            pathway / "pathway_enrichment_summary.tsv"
        ),
        significant_pathway_enrichment_results=(
            pathway
            / "tables"
            / "significant_pathway_enrichment_results.tsv"
        ),

        gene_set_summary=(
            gene_sets / "gene_set_summary.tsv"
        ),
        candidate_panel_mapping=(
            gene_sets / "candidate_panel_mapping.tsv"
        ),

        candidate_network_nodes=(
            candidate_network / "candidate_network_nodes.tsv"
        ),
        candidate_network_edges=(
            candidate_network / "candidate_network_edges.tsv"
        ),
        candidate_network_summary=(
            candidate_network / "candidate_network_summary.tsv"
        ),

        protein_annotation=(
            pxd_processed / "PXD043901_protein_annotation_extended.csv"
        ),
        delta_t1_vs_t0=(
            pxd_processed / "PXD043901_delta_T1_vs_T0.csv"
        ),
        delta_t6_vs_t0=(
            pxd_processed / "PXD043901_delta_T6_vs_T0.csv"
        ),

        output_svg=(
            output / "final_validation_framework.svg"
        ),
        output_pdf=(
            output / "final_validation_framework.pdf"
        ),
        output_png=(
            output / "final_validation_framework.png"
        ),
        output_caption=(
            output / "final_validation_framework_caption.txt"
        ),
        output_report=(
            output / "final_validation_framework_report.txt"
        ),
        output_manifest=(
            output / "final_validation_framework_manifest.json"
        ),
    )



# Required input schemas


REQUIRED_COLUMNS: dict[str, tuple[str, ...]] = {
    "protein_response_summary": (
        "protein_group_id",
        "comparison",
        "n_observed_pairs",
        "n_positive",
        "n_negative",
        "dominant_direction",
        "direction_consistency",
        "median_delta_log2",
        "mean_delta_log2",
        "sd_delta_log2",
        "q25_delta_log2",
        "q75_delta_log2",
        "primary_accession",
        "organism_origin",
        "gene_symbol",
        "human_protein_name",
        "first_microbe_label",
        "detected_samples",
        "detected_fraction",
    ),
    "robust_human_responses": (
        "protein_group_id",
        "comparison",
        "n_observed_pairs",
        "n_positive",
        "n_negative",
        "dominant_direction",
        "direction_consistency",
        "median_delta_log2",
        "mean_delta_log2",
        "sd_delta_log2",
        "q25_delta_log2",
        "q75_delta_log2",
        "primary_accession",
        "organism_origin",
        "gene_symbol",
        "human_protein_name",
        "first_microbe_label",
        "detected_samples",
        "detected_fraction",
    ),
    "persistent_human_responses": (
        "protein_group_id",
        "primary_accession",
        "gene_symbol",
        "human_protein_name",
        "n_pairs_T1",
        "n_pairs_T6",
        "consistency_T1",
        "consistency_T6",
        "median_delta_T1",
        "median_delta_T6",
        "wilcoxon_p_T1",
        "wilcoxon_p_T6",
        "wilcoxon_fdr_T1",
        "wilcoxon_fdr_T6",
        "sign_test_p_T1",
        "sign_test_p_T6",
        "score_T1",
        "score_T6",
        "trajectory_class",
        "combined_reproducibility_score",
        "persistent_fdr_support",
    ),
    "ranked_human_responses": (
        "protein_group_id",
        "comparison",
        "n_observed_pairs",
        "n_positive",
        "n_negative",
        "dominant_direction",
        "direction_consistency",
        "median_delta_log2",
        "mean_delta_log2",
        "sd_delta_log2",
        "n_tested_pairs",
        "wilcoxon_statistic",
        "wilcoxon_p",
        "sign_test_p",
        "wilcoxon_fdr",
        "sign_test_fdr",
        "primary_accession",
        "organism_origin",
        "gene_symbol",
        "human_protein_name",
        "first_microbe_label",
        "detected_samples",
        "detected_fraction",
        "absolute_median_delta",
        "reproducibility_score",
    ),
    "stringent_human_responses": (
        "protein_group_id",
        "comparison",
        "n_observed_pairs",
        "n_positive",
        "n_negative",
        "dominant_direction",
        "direction_consistency",
        "median_delta_log2",
        "mean_delta_log2",
        "sd_delta_log2",
        "n_tested_pairs",
        "wilcoxon_statistic",
        "wilcoxon_p",
        "sign_test_p",
        "wilcoxon_fdr",
        "sign_test_fdr",
        "primary_accession",
        "organism_origin",
        "gene_symbol",
        "human_protein_name",
        "first_microbe_label",
        "detected_samples",
        "detected_fraction",
        "absolute_median_delta",
        "reproducibility_score",
    ),
    "human_trajectory_classes": (
        "protein_group_id",
        "primary_accession",
        "gene_symbol",
        "human_protein_name",
        "n_pairs_T1",
        "n_pairs_T6",
        "consistency_T1",
        "consistency_T6",
        "median_delta_T1",
        "median_delta_T6",
        "wilcoxon_p_T1",
        "wilcoxon_p_T6",
        "wilcoxon_fdr_T1",
        "wilcoxon_fdr_T6",
        "sign_test_p_T1",
        "sign_test_p_T6",
        "score_T1",
        "score_T6",
        "trajectory_class",
        "combined_reproducibility_score",
        "persistent_fdr_support",
    ),
    "trajectory_class_counts": (
        "trajectory_class",
        "protein_groups",
    ),
    "longitudinal_statistical_tests": (
        "protein_group_id",
        "comparison",
        "n_observed_pairs",
        "n_positive",
        "n_negative",
        "dominant_direction",
        "direction_consistency",
        "median_delta_log2",
        "mean_delta_log2",
        "sd_delta_log2",
        "n_tested_pairs",
        "wilcoxon_statistic",
        "wilcoxon_p",
        "sign_test_p",
        "wilcoxon_fdr",
        "sign_test_fdr",
        "primary_accession",
        "organism_origin",
        "gene_symbol",
        "human_protein_name",
        "first_microbe_label",
        "detected_samples",
        "detected_fraction",
        "absolute_median_delta",
        "reproducibility_score",
    ),
    "module_longitudinal_statistics": (
        "module",
        "comparison",
        "n_pairs",
        "n_positive_response",
        "n_negative_response",
        "direction_consistency",
        "median_response_delta",
        "mean_response_delta",
        "sd_response_delta",
        "wilcoxon_statistic",
        "wilcoxon_p",
        "wilcoxon_fdr",
    ),
    "module_coverage": (
        "module",
        "module_description",
        "expected_treatment_direction",
        "requested_genes",
        "available_genes",
        "coverage_fraction",
        "available_gene_symbols",
        "missing_gene_symbols",
    ),
    "patient_module_response_deltas": (
        "module",
        "patient_id",
        "T0",
        "T1",
        "T6",
        "delta_T1_vs_T0",
        "delta_T6_vs_T0",
    ),
    "module_delta_correlations": (
        "comparison",
        "module_a",
        "module_b",
        "n_patients",
        "spearman_rho",
        "spearman_p",
    ),
    "leave_one_out_summary": (
        "module",
        "comparison",
        "removals_tested",
        "minimum_positive_responses",
        "minimum_direction_consistency",
        "minimum_median_delta",
        "maximum_median_delta",
        "maximum_wilcoxon_p",
        "robust_to_all_single_gene_removals",
    ),
    "leave_one_protein_out_results": (
        "module",
        "removed_gene",
        "retained_genes",
        "comparison",
        "n_pairs",
        "n_positive_response",
        "n_negative_response",
        "direction_consistency",
        "median_response_delta",
        "mean_response_delta",
        "wilcoxon_statistic",
        "wilcoxon_p",
    ),
    "aggregated_candidate_evidence": (
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
        "maximum_record_score",
        "mean_record_score",
        "minimum_record_score",
        "strongest_evidence_level",
        "any_prioritised",
        "base_feature_score",
        "replication_adjustment",
        "aggregated_reproducibility_score",
        "within_feature_type_percentile",
        "within_feature_type_rank",
        "candidate_tier",
    ),
    "cross_cohort_evidence_registry": (
        "evidence_id",
        "feature_id",
        "feature_name",
        "feature_type",
        "modality",
        "dataset",
        "study_role",
        "comparison",
        "effect_size",
        "effect_metric",
        "effect_direction",
        "direction_consistency",
        "robustness_value",
        "robustness_metric",
        "statistic",
        "p_value",
        "adjusted_p_value",
        "n_observations",
        "abundance_or_detection",
        "evidence_level",
        "is_prioritised",
        "source_file",
    ),
    "final_candidate_panel": (
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
        "maximum_record_score",
        "mean_record_score",
        "minimum_record_score",
        "strongest_evidence_level",
        "any_prioritised",
        "base_feature_score",
        "replication_adjustment",
        "aggregated_reproducibility_score",
        "within_feature_type_percentile",
        "within_feature_type_rank",
        "candidate_tier",
    ),

    "prioritised_candidate_panel": (
        "within_feature_type_rank",
        "within_dataset_rank",
        "feature_id",
        "feature_name",
        "feature_type",
        "modality",
        "dataset",
        "study_role",
        "comparison",
        "effect_size",
        "effect_metric",
        "effect_direction",
        "effect_score",
        "consistency_score",
        "robustness_score",
        "support_quality_score",
        "reproducibility_score",
        "score_band",
        "evidence_level",
        "is_prioritised",
        "source_file",
    ),
    "scored_evidence_registry": (
        "evidence_id",
        "feature_id",
        "feature_name",
        "feature_type",
        "modality",
        "dataset",
        "study_role",
        "comparison",
        "effect_size",
        "effect_metric",
        "effect_direction",
        "direction_consistency",
        "robustness_value",
        "robustness_metric",
        "statistic",
        "p_value",
        "adjusted_p_value",
        "n_observations",
        "abundance_or_detection",
        "evidence_level",
        "is_prioritised",
        "source_file",
        "effect_score",
        "consistency_score",
        "robustness_score",
        "support_quality_score",
        "reproducibility_score",
        "score_band",
        "within_feature_type_rank",
        "within_dataset_rank",
    ),
    "regulatory_edges_high_confidence": (
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
    ),
    "regulatory_edge_summary": (
        "metric",
        "value",
    ),
    "pathway_enrichment_summary": (
        "gene_set",
        "database",
        "input_symbols",
        "mapped_entrez_ids",
        "universe_symbols",
        "universe_entrez_ids",
        "tested_terms",
        "significant_terms_fdr_0_05",
        "minimum_adjusted_p_value",
    ),
    "significant_pathway_enrichment_results": (
        "gene_set",
        "database",
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "RichFactor",
        "FoldEnrichment",
        "zScore",
        "pvalue",
        "p.adjust",
        "qvalue",
        "geneID",
        "Count",
        "foreground_count",
        "foreground_size",
        "background_count",
        "background_size",
        "foreground_fraction",
        "background_fraction",
        "fold_enrichment",
        "gene_symbols",
        "significant_fdr_0_05",
        "category",
        "subcategory",
    ),
    "gene_set_summary": (
        "gene_set",
        "n_genes",
        "fraction_of_candidate_panel",
        "description",
    ),
    "candidate_panel_mapping": (
        "gene_symbol",
        "uniprot_accession",
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
    ),
    "candidate_network_nodes": (
        "node_id",
        "label",
        "node_type",
        "molecular_entity",
        "is_regulatory_target",
        "is_pathway_member",
        "biological_module",
        "database",
        "pathway_id",
        "gene_set",
        "adjusted_p_value",
        "fold_enrichment",
        "interpretation_priority",
        "degree",
        "in_degree",
        "out_degree",
    ),
    "candidate_network_edges": (
        "edge_id",
        "source",
        "target",
        "interaction",
        "database",
        "pathway_id",
        "gene_set",
        "adjusted_p_value",
        "fold_enrichment",
        "interpretation_priority",
        "biological_module",
    ),
    "candidate_network_summary": (
        "metric",
        "value",
    ),

}



# Input-table registry


@dataclass
class InputTables:
    """Container holding all validated input tables."""

    protein_response_summary: pd.DataFrame
    robust_human_responses: pd.DataFrame
    persistent_human_responses: pd.DataFrame
    ranked_human_responses: pd.DataFrame
    stringent_human_responses: pd.DataFrame
    human_trajectory_classes: pd.DataFrame
    trajectory_class_counts: pd.DataFrame
    longitudinal_statistical_tests: pd.DataFrame

    module_longitudinal_statistics: pd.DataFrame
    module_coverage: pd.DataFrame
    patient_module_response_deltas: pd.DataFrame
    module_delta_correlations: pd.DataFrame
    leave_one_out_summary: pd.DataFrame
    leave_one_protein_out_results: pd.DataFrame

    aggregated_candidate_evidence: pd.DataFrame
    cross_cohort_evidence_registry: pd.DataFrame
    final_candidate_panel: pd.DataFrame
    prioritised_candidate_panel: pd.DataFrame
    scored_evidence_registry: pd.DataFrame

    regulatory_edges_high_confidence: pd.DataFrame
    regulatory_edge_summary: pd.DataFrame

    pathway_enrichment_summary: pd.DataFrame
    significant_pathway_enrichment_results: pd.DataFrame

    gene_set_summary: pd.DataFrame
    candidate_panel_mapping: pd.DataFrame

    candidate_network_nodes: pd.DataFrame
    candidate_network_edges: pd.DataFrame
    candidate_network_summary: pd.DataFrame


def configure_logging(verbose: bool = False) -> None:
    """Configure console logging for the script."""

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )


def read_delimited_table(path: Path) -> pd.DataFrame:
    """
    Read a CSV or TSV table based on its filename suffix.

    Empty tables are permitted when their header is present, which is
    relevant for regulatory_edges_high_confidence.tsv.
    """

    suffix = path.suffix.lower()

    if suffix == ".csv":
        separator = ","
    elif suffix in {".tsv", ".txt"}:
        separator = "\t"
    else:
        raise ValueError(
            f"Unsupported table extension for {path}: {suffix}"
        )

    try:
        frame = pd.read_csv(
            path,
            sep=separator,
            low_memory=False,
        )
    except pd.errors.EmptyDataError as exc:
        raise ValueError(
            f"Input table contains no header or data: {path}"
        ) from exc
    except Exception as exc:
        raise RuntimeError(
            f"Failed to read input table: {path}"
        ) from exc

    return frame


def validate_required_columns(
    frame: pd.DataFrame,
    table_name: str,
) -> None:
    """Require every documented column for a named input table."""

    required = REQUIRED_COLUMNS[table_name]
    missing = [
        column
        for column in required
        if column not in frame.columns
    ]

    if missing:
        raise ValueError(
            f"{table_name} is missing required columns: "
            + ", ".join(missing)
        )


def validate_no_duplicate_columns(
    frame: pd.DataFrame,
    table_name: str,
) -> None:
    """Reject duplicated column labels."""

    duplicated = frame.columns[
        frame.columns.duplicated()
    ].tolist()

    if duplicated:
        raise ValueError(
            f"{table_name} contains duplicate column names: "
            + ", ".join(map(str, duplicated))
        )


def validate_nonempty_table(
    frame: pd.DataFrame,
    table_name: str,
) -> None:
    """Require at least one row for tables expected to contain results."""

    allowed_empty = {
        "regulatory_edges_high_confidence",
    }

    if frame.empty and table_name not in allowed_empty:
        raise ValueError(
            f"{table_name} contains zero rows."
        )


def coerce_boolean_series(series: pd.Series) -> pd.Series:
    """
    Convert common serialized boolean representations to nullable BooleanDtype.
    """

    if pd.api.types.is_bool_dtype(series):
        return series.astype("boolean")

    normalized = (
        series.astype("string")
        .str.strip()
        .str.lower()
    )

    mapping = {
        "true": True,
        "false": False,
        "1": True,
        "0": False,
        "yes": True,
        "no": False,
    }

    converted = normalized.map(mapping)
    converted[series.isna()] = pd.NA

    invalid = (
        converted.isna()
        & series.notna()
        & normalized.ne("")
        & normalized.ne("nan")
        & normalized.ne("<na>")
    )

    if invalid.any():
        invalid_values = sorted(
            normalized.loc[invalid]
            .dropna()
            .unique()
            .tolist()
        )
        raise ValueError(
            "Unrecognized boolean values: "
            + ", ".join(invalid_values)
        )

    return converted.astype("boolean")


def normalize_boolean_columns(
    tables: InputTables,
) -> None:
    """Normalize boolean columns used by downstream summaries."""

    boolean_columns: dict[str, tuple[str, ...]] = {
        "persistent_human_responses": (
            "persistent_fdr_support",
        ),
        "human_trajectory_classes": (
            "persistent_fdr_support",
        ),
        "leave_one_out_summary": (
            "robust_to_all_single_gene_removals",
        ),
        "aggregated_candidate_evidence": (
            "any_prioritised",
        ),
        "cross_cohort_evidence_registry": (
            "is_prioritised",
        ),
        "final_candidate_panel": (
            "any_prioritised",
        ),
        "prioritised_candidate_panel": (
            "is_prioritised",
        ),
        "scored_evidence_registry": (
            "is_prioritised",
        ),
        "regulatory_edges_high_confidence": (
            "is_experimentally_validated",
        ),
        "significant_pathway_enrichment_results": (
            "significant_fdr_0_05",
        ),
        "candidate_network_nodes": (
            "is_regulatory_target",
            "is_pathway_member",
        ),
    }

    for table_name, columns in boolean_columns.items():
        frame = getattr(tables, table_name)

        for column in columns:
            if column in frame.columns:
                frame[column] = coerce_boolean_series(
                    frame[column]
                )


def table_path_registry(
    paths: RepositoryPaths,
) -> dict[str, Path]:
    """Map each table name to its repository path."""

    return {
        "protein_response_summary":
            paths.protein_response_summary,
        "robust_human_responses":
            paths.robust_human_responses,
        "persistent_human_responses":
            paths.persistent_human_responses,
        "ranked_human_responses":
            paths.ranked_human_responses,
        "stringent_human_responses":
            paths.stringent_human_responses,
        "human_trajectory_classes":
            paths.human_trajectory_classes,
        "trajectory_class_counts":
            paths.trajectory_class_counts,
        "longitudinal_statistical_tests":
            paths.longitudinal_statistical_tests,

        "module_longitudinal_statistics":
            paths.module_longitudinal_statistics,
        "module_coverage":
            paths.module_coverage,
        "patient_module_response_deltas":
            paths.patient_module_response_deltas,
        "module_delta_correlations":
            paths.module_delta_correlations,
        "leave_one_out_summary":
            paths.leave_one_out_summary,
        "leave_one_protein_out_results":
            paths.leave_one_protein_out_results,

        "aggregated_candidate_evidence":
            paths.aggregated_candidate_evidence,
        "cross_cohort_evidence_registry":
            paths.cross_cohort_evidence_registry,
        "final_candidate_panel":
            paths.final_candidate_panel,
        "prioritised_candidate_panel":
            paths.prioritised_candidate_panel,
        "scored_evidence_registry":
            paths.scored_evidence_registry,

        "regulatory_edges_high_confidence":
            paths.regulatory_edges_high_confidence,
        "regulatory_edge_summary":
            paths.regulatory_edge_summary,

        "pathway_enrichment_summary":
            paths.pathway_enrichment_summary,
        "significant_pathway_enrichment_results":
            paths.significant_pathway_enrichment_results,

        "gene_set_summary":
            paths.gene_set_summary,
        "candidate_panel_mapping":
            paths.candidate_panel_mapping,

        "candidate_network_nodes":
            paths.candidate_network_nodes,
        "candidate_network_edges":
            paths.candidate_network_edges,
        "candidate_network_summary":
            paths.candidate_network_summary,
    }


def validate_input_files_exist(
    registry: Mapping[str, Path],
) -> None:
    """Require every repository-derived input file to exist."""

    missing = [
        (table_name, path)
        for table_name, path in registry.items()
        if not path.is_file()
    ]

    if missing:
        details = "\n".join(
            f"  - {table_name}: {path}"
            for table_name, path in missing
        )
        raise FileNotFoundError(
            "Required final-figure inputs are missing:\n"
            + details
        )


def load_input_tables(
    paths: RepositoryPaths,
) -> InputTables:
    """Load and validate every final-figure input table."""

    registry = table_path_registry(paths)
    validate_input_files_exist(registry)

    loaded: dict[str, pd.DataFrame] = {}

    for table_name, path in registry.items():
        logging.info("Reading %s", path.relative_to(
            paths.repository_root
        ))

        frame = read_delimited_table(path)

        validate_no_duplicate_columns(
            frame=frame,
            table_name=table_name,
        )
        validate_required_columns(
            frame=frame,
            table_name=table_name,
        )
        validate_nonempty_table(
            frame=frame,
            table_name=table_name,
        )

        loaded[table_name] = frame

    tables = InputTables(**loaded)
    normalize_boolean_columns(tables)

    return tables


# Cross-table consistency validation


NUMERIC_COLUMNS: dict[str, tuple[str, ...]] = {
    "protein_response_summary": (
        "n_observed_pairs",
        "n_positive",
        "n_negative",
        "direction_consistency",
        "median_delta_log2",
        "mean_delta_log2",
        "sd_delta_log2",
        "q25_delta_log2",
        "q75_delta_log2",
        "detected_samples",
        "detected_fraction",
    ),
    "robust_human_responses": (
        "n_observed_pairs",
        "n_positive",
        "n_negative",
        "direction_consistency",
        "median_delta_log2",
        "mean_delta_log2",
        "sd_delta_log2",
        "q25_delta_log2",
        "q75_delta_log2",
        "detected_samples",
        "detected_fraction",
    ),
    "persistent_human_responses": (
        "n_pairs_T1",
        "n_pairs_T6",
        "consistency_T1",
        "consistency_T6",
        "median_delta_T1",
        "median_delta_T6",
        "wilcoxon_p_T1",
        "wilcoxon_p_T6",
        "wilcoxon_fdr_T1",
        "wilcoxon_fdr_T6",
        "sign_test_p_T1",
        "sign_test_p_T6",
        "score_T1",
        "score_T6",
        "combined_reproducibility_score",
    ),
    "ranked_human_responses": (
        "n_observed_pairs",
        "n_positive",
        "n_negative",
        "direction_consistency",
        "median_delta_log2",
        "mean_delta_log2",
        "sd_delta_log2",
        "n_tested_pairs",
        "wilcoxon_statistic",
        "wilcoxon_p",
        "sign_test_p",
        "wilcoxon_fdr",
        "sign_test_fdr",
        "detected_samples",
        "detected_fraction",
        "absolute_median_delta",
        "reproducibility_score",
    ),
    "stringent_human_responses": (
        "n_observed_pairs",
        "n_positive",
        "n_negative",
        "direction_consistency",
        "median_delta_log2",
        "mean_delta_log2",
        "sd_delta_log2",
        "n_tested_pairs",
        "wilcoxon_statistic",
        "wilcoxon_p",
        "sign_test_p",
        "wilcoxon_fdr",
        "sign_test_fdr",
        "detected_samples",
        "detected_fraction",
        "absolute_median_delta",
        "reproducibility_score",
    ),
    "human_trajectory_classes": (
        "n_pairs_T1",
        "n_pairs_T6",
        "consistency_T1",
        "consistency_T6",
        "median_delta_T1",
        "median_delta_T6",
        "wilcoxon_p_T1",
        "wilcoxon_p_T6",
        "wilcoxon_fdr_T1",
        "wilcoxon_fdr_T6",
        "sign_test_p_T1",
        "sign_test_p_T6",
        "score_T1",
        "score_T6",
        "combined_reproducibility_score",
    ),
    "trajectory_class_counts": (
        "protein_groups",
    ),
    "longitudinal_statistical_tests": (
        "n_observed_pairs",
        "n_positive",
        "n_negative",
        "direction_consistency",
        "median_delta_log2",
        "mean_delta_log2",
        "sd_delta_log2",
        "n_tested_pairs",
        "wilcoxon_statistic",
        "wilcoxon_p",
        "sign_test_p",
        "wilcoxon_fdr",
        "sign_test_fdr",
        "detected_samples",
        "detected_fraction",
        "absolute_median_delta",
        "reproducibility_score",
    ),
    "module_longitudinal_statistics": (
        "n_pairs",
        "n_positive_response",
        "n_negative_response",
        "direction_consistency",
        "median_response_delta",
        "mean_response_delta",
        "sd_response_delta",
        "wilcoxon_statistic",
        "wilcoxon_p",
        "wilcoxon_fdr",
    ),
    "module_coverage": (
        "requested_genes",
        "available_genes",
        "coverage_fraction",
    ),
    "patient_module_response_deltas": (
        "T0",
        "T1",
        "T6",
        "delta_T1_vs_T0",
        "delta_T6_vs_T0",
    ),
    "module_delta_correlations": (
        "n_patients",
        "spearman_rho",
        "spearman_p",
    ),
    "leave_one_out_summary": (
        "removals_tested",
        "minimum_positive_responses",
        "minimum_direction_consistency",
        "minimum_median_delta",
        "maximum_median_delta",
        "maximum_wilcoxon_p",
    ),
    "leave_one_protein_out_results": (
        "retained_genes",
        "n_pairs",
        "n_positive_response",
        "n_negative_response",
        "direction_consistency",
        "median_response_delta",
        "mean_response_delta",
        "wilcoxon_statistic",
        "wilcoxon_p",
    ),
    "aggregated_candidate_evidence": (
        "n_evidence_records",
        "n_comparisons",
        "direction_agreement",
        "maximum_absolute_effect",
        "median_absolute_effect",
        "mean_direction_consistency",
        "mean_robustness_value",
        "maximum_record_score",
        "mean_record_score",
        "minimum_record_score",
        "base_feature_score",
        "replication_adjustment",
        "aggregated_reproducibility_score",
        "within_feature_type_percentile",
        "within_feature_type_rank",
    ),
    "cross_cohort_evidence_registry": (
        "effect_size",
        "direction_consistency",
        "robustness_value",
        "statistic",
        "p_value",
        "adjusted_p_value",
        "n_observations",
        "abundance_or_detection",
    ),
    "final_candidate_panel": (
        "n_evidence_records",
        "n_comparisons",
        "direction_agreement",
        "maximum_absolute_effect",
        "median_absolute_effect",
        "mean_direction_consistency",
        "mean_robustness_value",
        "maximum_record_score",
        "mean_record_score",
        "minimum_record_score",
        "base_feature_score",
        "replication_adjustment",
        "aggregated_reproducibility_score",
        "within_feature_type_percentile",
        "within_feature_type_rank",
    ),
    "prioritised_candidate_panel": (
        "within_feature_type_rank",
        "within_dataset_rank",
        "effect_size",
        "effect_score",
        "consistency_score",
        "robustness_score",
        "support_quality_score",
        "reproducibility_score",
    ),
    "scored_evidence_registry": (
        "effect_size",
        "direction_consistency",
        "robustness_value",
        "statistic",
        "p_value",
        "adjusted_p_value",
        "n_observations",
        "abundance_or_detection",
        "effect_score",
        "consistency_score",
        "robustness_score",
        "support_quality_score",
        "reproducibility_score",
        "within_feature_type_rank",
        "within_dataset_rank",
    ),
    "regulatory_edges_high_confidence": (
        "n_databases",
        "n_publications",
        "n_experimental_methods",
        "n_source_records",
    ),
    "regulatory_edge_summary": (
        "value",
    ),
    "pathway_enrichment_summary": (
        "input_symbols",
        "mapped_entrez_ids",
        "universe_symbols",
        "universe_entrez_ids",
        "tested_terms",
        "significant_terms_fdr_0_05",
        "minimum_adjusted_p_value",
    ),
    "significant_pathway_enrichment_results": (
        "RichFactor",
        "FoldEnrichment",
        "zScore",
        "pvalue",
        "p.adjust",
        "qvalue",
        "Count",
        "foreground_count",
        "foreground_size",
        "background_count",
        "background_size",
        "foreground_fraction",
        "background_fraction",
        "fold_enrichment",
    ),
    "gene_set_summary": (
        "n_genes",
        "fraction_of_candidate_panel",
    ),
    "candidate_panel_mapping": (
        "aggregated_reproducibility_score",
        "within_feature_type_percentile",
        "within_feature_type_rank",
        "n_evidence_records",
        "n_comparisons",
        "direction_agreement",
        "maximum_absolute_effect",
        "median_absolute_effect",
        "mean_direction_consistency",
        "mean_robustness_value",
    ),
    "candidate_network_nodes": (
        "adjusted_p_value",
        "fold_enrichment",
        "degree",
        "in_degree",
        "out_degree",
    ),
    "candidate_network_edges": (
        "adjusted_p_value",
        "fold_enrichment",
    ),
    "candidate_network_summary": (
        "value",
    ),
}


PROPORTION_COLUMNS: dict[str, tuple[str, ...]] = {
    "protein_response_summary": (
        "direction_consistency",
        "detected_fraction",
    ),
    "robust_human_responses": (
        "direction_consistency",
        "detected_fraction",
    ),
    "persistent_human_responses": (
        "consistency_T1",
        "consistency_T6",
    ),
    "ranked_human_responses": (
        "direction_consistency",
        "detected_fraction",
    ),
    "stringent_human_responses": (
        "direction_consistency",
        "detected_fraction",
    ),
    "human_trajectory_classes": (
        "consistency_T1",
        "consistency_T6",
    ),
    "module_longitudinal_statistics": (
        "direction_consistency",
    ),
    "module_coverage": (
        "coverage_fraction",
    ),
    "leave_one_out_summary": (
        "minimum_direction_consistency",
    ),
    "leave_one_protein_out_results": (
        "direction_consistency",
    ),
    "aggregated_candidate_evidence": (
        "direction_agreement",
        "mean_direction_consistency",
    ),
    "cross_cohort_evidence_registry": (
        "direction_consistency",
    ),
    "final_candidate_panel": (
        "direction_agreement",
        "mean_direction_consistency",
    ),
    "scored_evidence_registry": (
        "direction_consistency",
    ),
    "gene_set_summary": (
        "fraction_of_candidate_panel",
    ),
    "candidate_panel_mapping": (
        "direction_agreement",
        "mean_direction_consistency",
    ),
}


P_VALUE_COLUMNS: dict[str, tuple[str, ...]] = {
    "persistent_human_responses": (
        "wilcoxon_p_T1",
        "wilcoxon_p_T6",
        "wilcoxon_fdr_T1",
        "wilcoxon_fdr_T6",
        "sign_test_p_T1",
        "sign_test_p_T6",
    ),
    "ranked_human_responses": (
        "wilcoxon_p",
        "sign_test_p",
        "wilcoxon_fdr",
        "sign_test_fdr",
    ),
    "stringent_human_responses": (
        "wilcoxon_p",
        "sign_test_p",
        "wilcoxon_fdr",
        "sign_test_fdr",
    ),
    "human_trajectory_classes": (
        "wilcoxon_p_T1",
        "wilcoxon_p_T6",
        "wilcoxon_fdr_T1",
        "wilcoxon_fdr_T6",
        "sign_test_p_T1",
        "sign_test_p_T6",
    ),
    "longitudinal_statistical_tests": (
        "wilcoxon_p",
        "sign_test_p",
        "wilcoxon_fdr",
        "sign_test_fdr",
    ),
    "module_longitudinal_statistics": (
        "wilcoxon_p",
        "wilcoxon_fdr",
    ),
    "module_delta_correlations": (
        "spearman_p",
    ),
    "leave_one_out_summary": (
        "maximum_wilcoxon_p",
    ),
    "leave_one_protein_out_results": (
        "wilcoxon_p",
    ),
    "cross_cohort_evidence_registry": (
        "p_value",
        "adjusted_p_value",
    ),
    "scored_evidence_registry": (
        "p_value",
        "adjusted_p_value",
    ),
    "pathway_enrichment_summary": (
        "minimum_adjusted_p_value",
    ),
    "significant_pathway_enrichment_results": (
        "pvalue",
        "p.adjust",
        "qvalue",
    ),
    "candidate_network_nodes": (
        "adjusted_p_value",
    ),
    "candidate_network_edges": (
        "adjusted_p_value",
    ),
}


def coerce_numeric_columns(
    tables: InputTables,
) -> None:
    """Convert documented numeric fields to numeric dtypes."""

    for table_name, columns in NUMERIC_COLUMNS.items():
        frame = getattr(tables, table_name)

        for column in columns:
            if column not in frame.columns:
                continue

            original_nonmissing = frame[column].notna()

            converted = pd.to_numeric(
                frame[column],
                errors="coerce",
            )

            newly_missing = (
                original_nonmissing
                & converted.isna()
            )

            if newly_missing.any():
                bad_values = (
                    frame.loc[newly_missing, column]
                    .astype("string")
                    .dropna()
                    .unique()
                    .tolist()
                )

                raise ValueError(
                    f"{table_name}.{column} contains "
                    f"non-numeric values: {bad_values[:10]}"
                )

            frame[column] = converted


def validate_proportion_columns(
    tables: InputTables,
    tolerance: float = 1e-9,
) -> None:
    """Require proportions to remain within [0, 1]."""

    for table_name, columns in PROPORTION_COLUMNS.items():
        frame = getattr(tables, table_name)

        for column in columns:
            values = frame[column].dropna()

            invalid = (
                (values < -tolerance)
                | (values > 1.0 + tolerance)
            )

            if invalid.any():
                examples = values.loc[invalid].head(10).tolist()

                raise ValueError(
                    f"{table_name}.{column} contains values "
                    f"outside [0, 1]: {examples}"
                )


def validate_p_value_columns(
    tables: InputTables,
    tolerance: float = 1e-9,
) -> None:
    """Require raw and adjusted p-values to remain within [0, 1]."""

    for table_name, columns in P_VALUE_COLUMNS.items():
        frame = getattr(tables, table_name)

        for column in columns:
            values = frame[column].dropna()

            invalid = (
                (values < -tolerance)
                | (values > 1.0 + tolerance)
            )

            if invalid.any():
                examples = values.loc[invalid].head(10).tolist()

                raise ValueError(
                    f"{table_name}.{column} contains invalid "
                    f"probabilities: {examples}"
                )


def validate_count_balance(
    frame: pd.DataFrame,
    total_column: str,
    positive_column: str,
    negative_column: str,
    table_name: str,
) -> None:
    """
    Check that positive and negative response counts do not exceed totals.
    """

    available = frame[
        [
            total_column,
            positive_column,
            negative_column,
        ]
    ].dropna()

    invalid = (
        available[positive_column]
        + available[negative_column]
        > available[total_column]
    )

    if invalid.any():
        examples = available.loc[invalid].head(10)

        raise ValueError(
            f"{table_name} contains positive/negative counts "
            f"that exceed {total_column}:\n"
            f"{examples.to_string(index=False)}"
        )


def validate_response_count_balances(
    tables: InputTables,
) -> None:
    """Validate count arithmetic in longitudinal response tables."""

    for table_name in (
        "protein_response_summary",
        "robust_human_responses",
        "ranked_human_responses",
        "stringent_human_responses",
        "longitudinal_statistical_tests",
    ):
        validate_count_balance(
            frame=getattr(tables, table_name),
            total_column="n_observed_pairs",
            positive_column="n_positive",
            negative_column="n_negative",
            table_name=table_name,
        )

    validate_count_balance(
        frame=tables.module_longitudinal_statistics,
        total_column="n_pairs",
        positive_column="n_positive_response",
        negative_column="n_negative_response",
        table_name="module_longitudinal_statistics",
    )

    validate_count_balance(
        frame=tables.leave_one_protein_out_results,
        total_column="n_pairs",
        positive_column="n_positive_response",
        negative_column="n_negative_response",
        table_name="leave_one_protein_out_results",
    )


def validate_unique_keys(
    frame: pd.DataFrame,
    columns: Sequence[str],
    table_name: str,
) -> None:
    """Reject duplicate records for a documented logical key."""

    duplicated = frame.duplicated(
        subset=list(columns),
        keep=False,
    )

    if duplicated.any():
        examples = (
            frame.loc[duplicated, list(columns)]
            .head(20)
            .to_string(index=False)
        )

        raise ValueError(
            f"{table_name} contains duplicate logical keys "
            f"{tuple(columns)}:\n{examples}"
        )


def validate_logical_keys(
    tables: InputTables,
) -> None:
    """Validate uniqueness for tables with stable record identifiers."""

    validate_unique_keys(
        tables.protein_response_summary,
        ("protein_group_id", "comparison"),
        "protein_response_summary",
    )

    validate_unique_keys(
        tables.robust_human_responses,
        ("protein_group_id", "comparison"),
        "robust_human_responses",
    )

    validate_unique_keys(
        tables.persistent_human_responses,
        ("protein_group_id",),
        "persistent_human_responses",
    )

    validate_unique_keys(
        tables.human_trajectory_classes,
        ("protein_group_id",),
        "human_trajectory_classes",
    )

    validate_unique_keys(
        tables.trajectory_class_counts,
        ("trajectory_class",),
        "trajectory_class_counts",
    )

    validate_unique_keys(
        tables.module_longitudinal_statistics,
        ("module", "comparison"),
        "module_longitudinal_statistics",
    )

    validate_unique_keys(
        tables.module_coverage,
        ("module",),
        "module_coverage",
    )

    validate_unique_keys(
        tables.patient_module_response_deltas,
        ("module", "patient_id"),
        "patient_module_response_deltas",
    )

    validate_unique_keys(
        tables.module_delta_correlations,
        ("comparison", "module_a", "module_b"),
        "module_delta_correlations",
    )

    validate_unique_keys(
        tables.leave_one_out_summary,
        ("module", "comparison"),
        "leave_one_out_summary",
    )

    validate_unique_keys(
        tables.leave_one_protein_out_results,
        ("module", "removed_gene", "comparison"),
        "leave_one_protein_out_results",
    )

    validate_unique_keys(
        tables.cross_cohort_evidence_registry,
        ("evidence_id",),
        "cross_cohort_evidence_registry",
    )

    validate_unique_keys(
        tables.scored_evidence_registry,
        ("evidence_id",),
        "scored_evidence_registry",
    )

    validate_unique_keys(
        tables.aggregated_candidate_evidence,
        ("feature_id", "dataset", "study_role"),
        "aggregated_candidate_evidence",
    )

    validate_unique_keys(
        tables.final_candidate_panel,
        ("feature_id", "dataset", "study_role"),
        "final_candidate_panel",
    )

    validate_unique_keys(
        tables.candidate_panel_mapping,
        ("gene_symbol",),
        "candidate_panel_mapping",
    )

    validate_unique_keys(
        tables.candidate_network_nodes,
        ("node_id",),
        "candidate_network_nodes",
    )

    validate_unique_keys(
        tables.candidate_network_edges,
        ("edge_id",),
        "candidate_network_edges",
    )


def validate_expected_comparisons(
    tables: InputTables,
) -> None:
    """Require the expected paired longitudinal comparisons."""

    expected = {
        "T1_vs_T0",
        "T6_vs_T0",
    }

    comparison_tables = {
        "protein_response_summary":
            tables.protein_response_summary,
        "robust_human_responses":
            tables.robust_human_responses,
        "ranked_human_responses":
            tables.ranked_human_responses,
        "stringent_human_responses":
            tables.stringent_human_responses,
        "longitudinal_statistical_tests":
            tables.longitudinal_statistical_tests,
        "module_longitudinal_statistics":
            tables.module_longitudinal_statistics,
        "leave_one_out_summary":
            tables.leave_one_out_summary,
        "leave_one_protein_out_results":
            tables.leave_one_protein_out_results,
    }

    for table_name, frame in comparison_tables.items():
        observed = set(
            frame["comparison"]
            .dropna()
            .astype(str)
            .unique()
            .tolist()
        )

        missing = expected - observed

        if missing:
            raise ValueError(
                f"{table_name} lacks expected comparisons: "
                + ", ".join(sorted(missing))
            )


def validate_trajectory_class_counts(
    tables: InputTables,
) -> None:
    """
    Confirm that the supplied trajectory count table reproduces the classes
    present in human_trajectory_classes.csv.
    """

    observed = (
        tables.human_trajectory_classes
        .groupby(
            "trajectory_class",
            dropna=False,
        )
        .size()
        .rename("observed_protein_groups")
        .reset_index()
    )

    reported = (
        tables.trajectory_class_counts
        .rename(
            columns={
                "protein_groups":
                    "reported_protein_groups",
            }
        )
        .copy()
    )

    comparison = reported.merge(
        observed,
        on="trajectory_class",
        how="outer",
    ).fillna(0)

    comparison[
        "reported_protein_groups"
    ] = comparison[
        "reported_protein_groups"
    ].astype(int)

    comparison[
        "observed_protein_groups"
    ] = comparison[
        "observed_protein_groups"
    ].astype(int)

    mismatch = (
        comparison["reported_protein_groups"]
        != comparison["observed_protein_groups"]
    )

    if mismatch.any():
        raise ValueError(
            "trajectory_class_counts.csv does not match "
            "human_trajectory_classes.csv:\n"
            + comparison.loc[mismatch].to_string(index=False)
        )


def validate_robust_subset_relationship(
    tables: InputTables,
) -> None:
    """
    Confirm that robust response rows are present in the full response table.
    """

    full_keys = set(
        zip(
            tables.protein_response_summary[
                "protein_group_id"
            ].astype(str),
            tables.protein_response_summary[
                "comparison"
            ].astype(str),
        )
    )

    robust_keys = set(
        zip(
            tables.robust_human_responses[
                "protein_group_id"
            ].astype(str),
            tables.robust_human_responses[
                "comparison"
            ].astype(str),
        )
    )

    missing = robust_keys - full_keys

    if missing:
        examples = sorted(missing)[:20]

        raise ValueError(
            "robust_human_responses_preliminary.csv contains "
            "records absent from protein_response_summary.csv: "
            f"{examples}"
        )


def validate_stringent_subset_relationship(
    tables: InputTables,
) -> None:
    """
    Confirm that stringent response rows are present in ranked responses.
    """

    ranked_keys = set(
        zip(
            tables.ranked_human_responses[
                "protein_group_id"
            ].astype(str),
            tables.ranked_human_responses[
                "comparison"
            ].astype(str),
        )
    )

    stringent_keys = set(
        zip(
            tables.stringent_human_responses[
                "protein_group_id"
            ].astype(str),
            tables.stringent_human_responses[
                "comparison"
            ].astype(str),
        )
    )

    missing = stringent_keys - ranked_keys

    if missing:
        examples = sorted(missing)[:20]

        raise ValueError(
            "stringent_human_responses.csv contains records "
            "absent from ranked_human_responses.csv: "
            f"{examples}"
        )


def validate_persistent_subset_relationship(
    tables: InputTables,
) -> None:
    """
    Confirm that persistent responses are represented in trajectory classes.
    """

    trajectory_ids = set(
        tables.human_trajectory_classes[
            "protein_group_id"
        ].astype(str)
    )

    persistent_ids = set(
        tables.persistent_human_responses[
            "protein_group_id"
        ].astype(str)
    )

    missing = persistent_ids - trajectory_ids

    if missing:
        raise ValueError(
            "persistent_human_responses.csv contains proteins "
            "absent from human_protein_trajectory_classes.csv: "
            f"{sorted(missing)[:20]}"
        )


def validate_final_candidate_subset(
    tables: InputTables,
) -> None:
    """
    Require every final-panel feature to occur in aggregated evidence.
    """

    aggregated_keys = set(
        zip(
            tables.aggregated_candidate_evidence[
                "feature_id"
            ].astype(str),
            tables.aggregated_candidate_evidence[
                "dataset"
            ].astype(str),
            tables.aggregated_candidate_evidence[
                "study_role"
            ].astype(str),
        )
    )

    final_keys = set(
        zip(
            tables.final_candidate_panel[
                "feature_id"
            ].astype(str),
            tables.final_candidate_panel[
                "dataset"
            ].astype(str),
            tables.final_candidate_panel[
                "study_role"
            ].astype(str),
        )
    )

    missing = final_keys - aggregated_keys

    if missing:
        raise ValueError(
            "final_candidate_panel.tsv contains records absent "
            "from aggregated_candidate_evidence.tsv: "
            f"{sorted(missing)[:20]}"
        )


def validate_candidate_mapping_subset(
    tables: InputTables,
) -> None:
    """
    Require mapped candidate proteins to occur in the final candidate panel.
    """

    final_proteins = set(
        tables.final_candidate_panel.loc[
            tables.final_candidate_panel[
                "feature_type"
            ].eq("protein"),
            "feature_name",
        ]
        .dropna()
        .astype(str)
    )

    mapped_symbols = set(
        tables.candidate_panel_mapping[
            "gene_symbol"
        ]
        .dropna()
        .astype(str)
    )

    missing = mapped_symbols - final_proteins

    if missing:
        raise ValueError(
            "candidate_panel_mapping.tsv contains gene symbols "
            "absent from protein entries in final_candidate_panel.tsv: "
            f"{sorted(missing)[:20]}"
        )


def validate_pathway_significance(
    tables: InputTables,
    alpha: float = 0.05,
) -> None:
    """
    Confirm that the significant-pathway table contains only FDR-significant
    records and that its Boolean significance flag agrees with p.adjust.
    """

    frame = tables.significant_pathway_enrichment_results

    if frame.empty:
        raise ValueError(
            "significant_pathway_enrichment_results.tsv "
            "contains no significant pathways."
        )

    adjusted = frame["p.adjust"]
    flag = frame["significant_fdr_0_05"].fillna(False)

    invalid_adjusted = adjusted > alpha + 1e-12

    if invalid_adjusted.any():
        examples = frame.loc[
            invalid_adjusted,
            [
                "gene_set",
                "database",
                "ID",
                "Description",
                "p.adjust",
            ],
        ].head(20)

        raise ValueError(
            "The significant-pathway table contains entries "
            f"with p.adjust > {alpha}:\n"
            f"{examples.to_string(index=False)}"
        )

    disagreement = flag.ne(adjusted.le(alpha))

    if disagreement.any():
        examples = frame.loc[
            disagreement,
            [
                "gene_set",
                "database",
                "ID",
                "p.adjust",
                "significant_fdr_0_05",
            ],
        ].head(20)

        raise ValueError(
            "Pathway significance flags disagree with p.adjust:\n"
            f"{examples.to_string(index=False)}"
        )


def validate_regulatory_summary(
    tables: InputTables,
) -> None:
    """
    Check consistency between the regulatory edge summary and edge table.

    The high-confidence edge table may legitimately contain zero rows.
    """

    summary = tables.regulatory_edge_summary.copy()

    metric_values = {
        str(row.metric): float(row.value)
        for row in summary.itertuples(index=False)
    }

    required_metrics = {
        "candidate_protein_nodes",
        "validated_interactions_after_alias_normalisation",
        "matched_candidate_proteins",
        "unmatched_candidate_proteins",
        "candidate_protein_overlap_fraction",
        "all_validated_regulatory_edges",
        "high_confidence_regulatory_edges",
        "unique_regulatory_mirnas",
    }

    missing_metrics = required_metrics - set(metric_values)

    if missing_metrics:
        raise ValueError(
            "regulatory_edge_summary.tsv lacks metrics: "
            + ", ".join(sorted(missing_metrics))
        )

    actual_high_confidence = len(
        tables.regulatory_edges_high_confidence
    )

    reported_high_confidence = int(
        round(
            metric_values[
                "high_confidence_regulatory_edges"
            ]
        )
    )

    if actual_high_confidence != reported_high_confidence:
        raise ValueError(
            "High-confidence regulatory-edge count mismatch: "
            f"table={actual_high_confidence}, "
            f"summary={reported_high_confidence}"
        )

    candidate_nodes = metric_values[
        "candidate_protein_nodes"
    ]
    matched = metric_values[
        "matched_candidate_proteins"
    ]
    reported_fraction = metric_values[
        "candidate_protein_overlap_fraction"
    ]

    calculated_fraction = (
        matched / candidate_nodes
        if candidate_nodes > 0
        else math.nan
    )

    if (
        math.isfinite(calculated_fraction)
        and not math.isclose(
            reported_fraction,
            calculated_fraction,
            rel_tol=1e-9,
            abs_tol=1e-9,
        )
    ):
        raise ValueError(
            "candidate_protein_overlap_fraction does not "
            "match matched/candidate counts."
        )


def validate_module_coverage(
    tables: InputTables,
) -> None:
    """Check arithmetic in module_coverage.csv."""

    frame = tables.module_coverage.copy()

    invalid_counts = (
        frame["available_genes"]
        > frame["requested_genes"]
    )

    if invalid_counts.any():
        raise ValueError(
            "module_coverage.csv has available_genes greater "
            "than requested_genes."
        )

    expected_fraction = np.where(
        frame["requested_genes"] > 0,
        frame["available_genes"]
        / frame["requested_genes"],
        np.nan,
    )

    discrepancy = (
        frame["coverage_fraction"].notna()
        & ~np.isclose(
            frame["coverage_fraction"],
            expected_fraction,
            rtol=1e-8,
            atol=1e-8,
            equal_nan=True,
        )
    )

    if discrepancy.any():
        examples = frame.loc[
            discrepancy,
            [
                "module",
                "requested_genes",
                "available_genes",
                "coverage_fraction",
            ],
        ]

        raise ValueError(
            "module_coverage.csv contains inconsistent "
            "coverage fractions:\n"
            f"{examples.to_string(index=False)}"
        )


def validate_all_input_tables(
    tables: InputTables,
) -> None:
    """Run all numerical and cross-table validation checks."""

    coerce_numeric_columns(tables)
    validate_proportion_columns(tables)
    validate_p_value_columns(tables)
    validate_response_count_balances(tables)
    validate_logical_keys(tables)
    validate_expected_comparisons(tables)

    validate_trajectory_class_counts(tables)
    validate_robust_subset_relationship(tables)
    validate_stringent_subset_relationship(tables)
    validate_persistent_subset_relationship(tables)

    validate_final_candidate_subset(tables)
    validate_candidate_mapping_subset(tables)
    validate_pathway_significance(tables)
    validate_regulatory_summary(tables)
    validate_module_coverage(tables)



# Derived summary objects


@dataclass(frozen=True)
class LongitudinalSummary:
    """Repository-derived longitudinal response summary."""

    total_human_proteins: int
    total_response_records: int
    robust_response_records: int
    robust_unique_proteins: int
    stringent_response_records: int
    stringent_unique_proteins: int
    persistent_unique_proteins: int
    persistent_fdr_supported: int
    trajectory_class_counts: Mapping[str, int]
    median_consistency_t1: float
    median_consistency_t6: float
    median_absolute_delta_t1: float
    median_absolute_delta_t6: float


@dataclass(frozen=True)
class ModuleSummary:
    """Repository-derived biological module validation summary."""

    total_modules: int
    total_module_comparisons: int
    robust_module_comparisons: int
    modules_robust_both_timepoints: tuple[str, ...]
    median_module_coverage: float
    minimum_module_coverage: float
    maximum_module_coverage: float
    median_direction_consistency_t1: float
    median_direction_consistency_t6: float
    median_response_delta_t1: float
    median_response_delta_t6: float


@dataclass(frozen=True)
class CrossCohortSummary:
    """Integrated evidence and candidate prioritisation summary."""

    evidence_records: int
    aggregated_features: int
    final_candidates: int
    final_candidate_types: Mapping[str, int]
    final_candidate_datasets: Mapping[str, int]
    high_confidence_candidates: int
    prioritised_candidates: int
    protein_candidates: int
    mirna_candidates: int
    replicated_candidates: int
    median_reproducibility_score: float
    maximum_reproducibility_score: float


@dataclass(frozen=True)
class RegulatorySummary:
    """Regulatory-network evidence summary."""

    candidate_protein_nodes: int
    validated_interactions_after_alias_normalisation: int
    matched_candidate_proteins: int
    unmatched_candidate_proteins: int
    candidate_protein_overlap_fraction: float
    all_validated_regulatory_edges: int
    high_confidence_regulatory_edges: int
    unique_regulatory_mirnas: int


@dataclass(frozen=True)
class PathwaySummary:
    """Pathway enrichment and integrated network summary."""

    tested_gene_set_database_pairs: int
    significant_pathway_terms: int
    significant_candidate_panel_terms: int
    significant_mir143_linked_terms: int
    significant_nonlinked_candidate_terms: int
    minimum_adjusted_p_value: float
    maximum_fold_enrichment: float
    candidate_network_nodes: int
    candidate_network_edges: int
    pathway_nodes: int
    protein_nodes: int
    primary_pathways: tuple[str, ...]


@dataclass(frozen=True)
class ValidationSummary:
    """All evidence summaries required for the final figure."""

    longitudinal: LongitudinalSummary
    modules: ModuleSummary
    cross_cohort: CrossCohortSummary
    regulatory: RegulatorySummary
    pathways: PathwaySummary


def safe_median(
    series: pd.Series,
) -> float:
    """Return a finite median or NaN for an empty numeric series."""

    values = pd.to_numeric(
        series,
        errors="coerce",
    ).dropna()

    if values.empty:
        return math.nan

    return float(values.median())


def safe_minimum(
    series: pd.Series,
) -> float:
    """Return a finite minimum or NaN for an empty numeric series."""

    values = pd.to_numeric(
        series,
        errors="coerce",
    ).dropna()

    if values.empty:
        return math.nan

    return float(values.min())


def safe_maximum(
    series: pd.Series,
) -> float:
    """Return a finite maximum or NaN for an empty numeric series."""

    values = pd.to_numeric(
        series,
        errors="coerce",
    ).dropna()

    if values.empty:
        return math.nan

    return float(values.max())


def value_counts_mapping(
    series: pd.Series,
) -> dict[str, int]:
    """Convert value counts into a deterministic string-to-int mapping."""

    counts = (
        series.dropna()
        .astype(str)
        .value_counts()
        .sort_index()
    )

    return {
        str(key): int(value)
        for key, value in counts.items()
    }


def count_unique_feature_ids(
    frame: pd.DataFrame,
) -> int:
    """Count unique non-null feature identifiers."""

    return int(
        frame["feature_id"]
        .dropna()
        .astype(str)
        .nunique()
    )


def build_longitudinal_summary(
    tables: InputTables,
) -> LongitudinalSummary:
    """Create the longitudinal evidence summary."""

    trajectories = tables.human_trajectory_classes
    robust = tables.robust_human_responses
    stringent = tables.stringent_human_responses
    persistent = tables.persistent_human_responses

    class_counts = {
        str(row.trajectory_class):
            int(row.protein_groups)
        for row in (
            tables.trajectory_class_counts
            .sort_values(
                "protein_groups",
                ascending=False,
            )
            .itertuples(index=False)
        )
    }

    return LongitudinalSummary(
        total_human_proteins=int(
            trajectories["protein_group_id"].nunique()
        ),
        total_response_records=int(
            len(tables.protein_response_summary)
        ),
        robust_response_records=int(len(robust)),
        robust_unique_proteins=int(
            robust["protein_group_id"].nunique()
        ),
        stringent_response_records=int(len(stringent)),
        stringent_unique_proteins=int(
            stringent["protein_group_id"].nunique()
        ),
        persistent_unique_proteins=int(
            persistent["protein_group_id"].nunique()
        ),
        persistent_fdr_supported=int(
            persistent[
                "persistent_fdr_support"
            ]
            .fillna(False)
            .sum()
        ),
        trajectory_class_counts=class_counts,
        median_consistency_t1=safe_median(
            trajectories["consistency_T1"]
        ),
        median_consistency_t6=safe_median(
            trajectories["consistency_T6"]
        ),
        median_absolute_delta_t1=safe_median(
            trajectories["median_delta_T1"].abs()
        ),
        median_absolute_delta_t6=safe_median(
            trajectories["median_delta_T6"].abs()
        ),
    )


def build_module_summary(
    tables: InputTables,
) -> ModuleSummary:
    """Create the module-level robustness summary."""

    statistics = tables.module_longitudinal_statistics
    coverage = tables.module_coverage
    leave_one_out = tables.leave_one_out_summary

    robust = leave_one_out[
        leave_one_out[
            "robust_to_all_single_gene_removals"
        ].fillna(False)
    ].copy()

    robust_counts = (
        robust.groupby("module")["comparison"]
        .nunique()
    )

    robust_both = tuple(
        sorted(
            robust_counts[
                robust_counts >= 2
            ].index.astype(str).tolist()
        )
    )

    t1 = statistics[
        statistics["comparison"].eq("T1_vs_T0")
    ]
    t6 = statistics[
        statistics["comparison"].eq("T6_vs_T0")
    ]

    return ModuleSummary(
        total_modules=int(
            statistics["module"].nunique()
        ),
        total_module_comparisons=int(
            len(statistics)
        ),
        robust_module_comparisons=int(
            len(robust)
        ),
        modules_robust_both_timepoints=robust_both,
        median_module_coverage=safe_median(
            coverage["coverage_fraction"]
        ),
        minimum_module_coverage=safe_minimum(
            coverage["coverage_fraction"]
        ),
        maximum_module_coverage=safe_maximum(
            coverage["coverage_fraction"]
        ),
        median_direction_consistency_t1=safe_median(
            t1["direction_consistency"]
        ),
        median_direction_consistency_t6=safe_median(
            t6["direction_consistency"]
        ),
        median_response_delta_t1=safe_median(
            t1["median_response_delta"]
        ),
        median_response_delta_t6=safe_median(
            t6["median_response_delta"]
        ),
    )


def build_cross_cohort_summary(
    tables: InputTables,
) -> CrossCohortSummary:
    """Create integrated cross-cohort candidate summary."""

    evidence = tables.cross_cohort_evidence_registry
    aggregated = tables.aggregated_candidate_evidence
    final = tables.final_candidate_panel

    feature_types = value_counts_mapping(
        final["feature_type"]
    )
    datasets = value_counts_mapping(
        final["dataset"]
    )

    high_confidence = int(
        final["candidate_tier"]
        .astype(str)
        .eq("high_confidence")
        .sum()
    )

    prioritised = int(
        final["candidate_tier"]
        .astype(str)
        .eq("prioritised")
        .sum()
    )

    protein_candidates = int(
        final["feature_type"]
        .astype(str)
        .eq("protein")
        .sum()
    )

    mirna_candidates = int(
        final["feature_type"]
        .astype(str)
        .isin(
            {
                "mature_miRNA",
                "miRNA_associated_locus",
                "probable_miRNA_host_gene",
            }
        )
        .sum()
    )

    replicated_candidates = int(
        final["n_evidence_records"]
        .fillna(0)
        .ge(2)
        .sum()
    )

    return CrossCohortSummary(
        evidence_records=int(len(evidence)),
        aggregated_features=count_unique_feature_ids(
            aggregated
        ),
        final_candidates=count_unique_feature_ids(
            final
        ),
        final_candidate_types=feature_types,
        final_candidate_datasets=datasets,
        high_confidence_candidates=high_confidence,
        prioritised_candidates=prioritised,
        protein_candidates=protein_candidates,
        mirna_candidates=mirna_candidates,
        replicated_candidates=replicated_candidates,
        median_reproducibility_score=safe_median(
            final[
                "aggregated_reproducibility_score"
            ]
        ),
        maximum_reproducibility_score=safe_maximum(
            final[
                "aggregated_reproducibility_score"
            ]
        ),
    )


def summary_metric_mapping(
    frame: pd.DataFrame,
) -> dict[str, float]:
    """Convert metric/value summary tables to a mapping."""

    return {
        str(row.metric): float(row.value)
        for row in frame.itertuples(index=False)
    }


def build_regulatory_summary(
    tables: InputTables,
) -> RegulatorySummary:
    """Create the regulatory-network evidence summary."""

    metrics = summary_metric_mapping(
        tables.regulatory_edge_summary
    )

    return RegulatorySummary(
        candidate_protein_nodes=int(
            round(metrics["candidate_protein_nodes"])
        ),
        validated_interactions_after_alias_normalisation=int(
            round(
                metrics[
                    "validated_interactions_after_alias_normalisation"
                ]
            )
        ),
        matched_candidate_proteins=int(
            round(metrics["matched_candidate_proteins"])
        ),
        unmatched_candidate_proteins=int(
            round(metrics["unmatched_candidate_proteins"])
        ),
        candidate_protein_overlap_fraction=float(
            metrics[
                "candidate_protein_overlap_fraction"
            ]
        ),
        all_validated_regulatory_edges=int(
            round(
                metrics[
                    "all_validated_regulatory_edges"
                ]
            )
        ),
        high_confidence_regulatory_edges=int(
            round(
                metrics[
                    "high_confidence_regulatory_edges"
                ]
            )
        ),
        unique_regulatory_mirnas=int(
            round(metrics["unique_regulatory_mirnas"])
        ),
    )


def build_pathway_summary(
    tables: InputTables,
) -> PathwaySummary:
    """Create pathway and candidate-network summary."""

    significant = (
        tables.significant_pathway_enrichment_results
        .loc[
            lambda frame:
                frame[
                    "significant_fdr_0_05"
                ].fillna(False)
        ]
        .copy()
    )

    nodes = tables.candidate_network_nodes
    edges = tables.candidate_network_edges

    primary_nodes = nodes[
        nodes["interpretation_priority"]
        .astype("string")
        .isin(
            {
                "primary",
                "primary_supportive",
            }
        )
        & nodes["node_type"].eq("pathway")
    ]

    primary_pathways = tuple(
        primary_nodes.sort_values(
            [
                "adjusted_p_value",
                "fold_enrichment",
            ],
            ascending=[
                True,
                False,
            ],
            na_position="last",
        )["label"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    return PathwaySummary(
        tested_gene_set_database_pairs=int(
            len(tables.pathway_enrichment_summary)
        ),
        significant_pathway_terms=int(
            len(significant)
        ),
        significant_candidate_panel_terms=int(
            significant[
                "gene_set"
            ].eq("candidate_panel").sum()
        ),
        significant_mir143_linked_terms=int(
            significant[
                "gene_set"
            ].eq("mir143_linked").sum()
        ),
        significant_nonlinked_candidate_terms=int(
            significant[
                "gene_set"
            ].eq("nonlinked_candidate").sum()
        ),
        minimum_adjusted_p_value=safe_minimum(
            significant["p.adjust"]
        ),
        maximum_fold_enrichment=safe_maximum(
            significant["fold_enrichment"]
        ),
        candidate_network_nodes=int(len(nodes)),
        candidate_network_edges=int(len(edges)),
        pathway_nodes=int(
            nodes["node_type"].eq("pathway").sum()
        ),
        protein_nodes=int(
            nodes["node_type"].eq("protein").sum()
        ),
        primary_pathways=primary_pathways,
    )


def build_validation_summary(
    tables: InputTables,
) -> ValidationSummary:
    """Build the complete final-figure evidence summary."""

    return ValidationSummary(
        longitudinal=build_longitudinal_summary(
            tables
        ),
        modules=build_module_summary(
            tables
        ),
        cross_cohort=build_cross_cohort_summary(
            tables
        ),
        regulatory=build_regulatory_summary(
            tables
        ),
        pathways=build_pathway_summary(
            tables
        ),
    )


# Figure configuration


@dataclass(frozen=True)
class FigureConfig:
    """Visual configuration for the final validation framework figure."""

    width_inches: float = 14.5
    height_inches: float = 11.0
    dpi: int = DEFAULT_OUTPUT_DPI

    title_fontsize: float = 15.5
    subtitle_fontsize: float = 9.5
    panel_label_fontsize: float = 12.0
    panel_title_fontsize: float = 10.5
    axis_label_fontsize: float = 8.5
    tick_fontsize: float = 7.5
    annotation_fontsize: float = 7.5
    small_fontsize: float = 6.8

    line_width: float = 0.9
    grid_line_width: float = 0.5
    marker_size: float = 32.0
    bar_edge_width: float = 0.5


@dataclass(frozen=True)
class FigurePalette:
    """Consistent figure colours."""

    navy: str = "#17365D"
    blue: str = "#3B75AF"
    light_blue: str = "#B9D3EA"
    teal: str = "#2A7F7F"
    light_teal: str = "#BFE1DF"
    green: str = "#4D8B57"
    light_green: str = "#CBE2CD"
    orange: str = "#C7761D"
    light_orange: str = "#F1D2A8"
    red: str = "#B34A4A"
    light_red: str = "#E8BABA"
    purple: str = "#7158A5"
    light_purple: str = "#D5CCE8"
    gold: str = "#B68A1E"
    light_gold: str = "#E7D9A8"
    grey: str = "#6E6E6E"
    mid_grey: str = "#A6A6A6"
    light_grey: str = "#E8E8E8"
    very_light_grey: str = "#F5F5F5"
    dark: str = "#222222"
    white: str = "#FFFFFF"


FIGURE_CONFIG = FigureConfig()
PALETTE = FigurePalette()


def configure_matplotlib() -> None:
    """Apply publication-oriented Matplotlib defaults."""

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": FIGURE_CONFIG.tick_fontsize,
            "axes.titlesize": FIGURE_CONFIG.panel_title_fontsize,
            "axes.labelsize": FIGURE_CONFIG.axis_label_fontsize,
            "xtick.labelsize": FIGURE_CONFIG.tick_fontsize,
            "ytick.labelsize": FIGURE_CONFIG.tick_fontsize,
            "axes.linewidth": 0.7,
            "axes.edgecolor": PALETTE.dark,
            "axes.labelcolor": PALETTE.dark,
            "xtick.color": PALETTE.dark,
            "ytick.color": PALETTE.dark,
            "text.color": PALETTE.dark,
            "figure.facecolor": PALETTE.white,
            "axes.facecolor": PALETTE.white,
            "savefig.facecolor": PALETTE.white,
            "savefig.edgecolor": PALETTE.white,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def remove_top_right_spines(ax: Axes) -> None:
    """Remove non-essential top and right axis spines."""

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def add_panel_label(
    ax: Axes,
    label: str,
    x: float = -0.12,
    y: float = 1.08,
) -> None:
    """Add a bold panel label in axes-relative coordinates."""

    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=FIGURE_CONFIG.panel_label_fontsize,
        fontweight="bold",
        va="top",
        ha="left",
    )


def format_module_name(module: str) -> str:
    """Convert repository module identifiers to readable labels."""

    replacements = {
        "blood_erythrocyte": "Blood / erythrocyte",
        "epithelial_barrier": "Epithelial barrier",
        "mucosal_defence": "Mucosal defence",
        "neutrophil_inflammation": "Neutrophil inflammation",
        "oxidative_stress": "Oxidative stress",
        "tissue_damage_cytoskeleton": "Tissue damage / cytoskeleton",
    }

    return replacements.get(
        str(module),
        str(module).replace("_", " ").title(),
    )


def format_trajectory_class(value: str) -> str:
    """Convert trajectory-class identifiers to readable labels."""

    return (
        str(value)
        .replace("_", " ")
        .replace("T1", "T1")
        .replace("T6", "T6")
        .title()
    )


def format_feature_type(value: str) -> str:
    """Convert feature-type identifiers to compact labels."""

    replacements = {
        "mature_miRNA": "Mature miRNA",
        "miRNA_associated_locus": "miRNA-associated locus",
        "probable_miRNA_host_gene": "Probable miRNA host gene",
        "protein": "Protein",
        "protein_module": "Protein module",
    }

    return replacements.get(
        str(value),
        str(value).replace("_", " ").title(),
    )


def format_dataset(value: str) -> str:
    """Return a compact dataset label."""

    return str(value)


def finite_or_zero(value: float) -> float:
    """Return zero for non-finite values."""

    if value is None:
        return 0.0

    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return 0.0

    return numeric if math.isfinite(numeric) else 0.0


def format_float(
    value: float,
    digits: int = 2,
) -> str:
    """Format finite floats while preserving unavailable values."""

    if value is None:
        return "NA"

    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "NA"

    if not math.isfinite(numeric):
        return "NA"

    return f"{numeric:.{digits}f}"


def format_percent(
    value: float,
    digits: int = 0,
) -> str:
    """Format a fraction as a percentage."""

    if value is None:
        return "NA"

    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "NA"

    if not math.isfinite(numeric):
        return "NA"

    return f"{100.0 * numeric:.{digits}f}%"


def wrap_label(
    value: str,
    width: int,
) -> str:
    """Wrap a label without breaking words."""

    return "\n".join(
        textwrap.wrap(
            str(value),
            width=width,
            break_long_words=False,
            break_on_hyphens=False,
        )
    )


def comparison_order(series: pd.Series) -> pd.Categorical:
    """Return ordered comparison labels."""

    return pd.Categorical(
        series,
        categories=[
            "T1_vs_T0",
            "T6_vs_T0",
        ],
        ordered=True,
    )


def comparison_label(value: str) -> str:
    """Convert comparison identifiers to compact labels."""

    mapping = {
        "T1_vs_T0": "T1 − T0",
        "T6_vs_T0": "T6 − T0",
    }

    return mapping.get(str(value), str(value))


def add_zero_line(ax: Axes) -> None:
    """Add a subtle horizontal zero reference line."""

    ax.axhline(
        0.0,
        color=PALETTE.grey,
        linewidth=0.7,
        linestyle="--",
        zorder=0,
    )


def annotate_bars(
    ax: Axes,
    bars: Iterable[Any],
    fmt: str = "{:.0f}",
    offset_fraction: float = 0.02,
    fontsize: float | None = None,
) -> None:
    """Annotate vertical bars with their heights."""

    font_size = (
        FIGURE_CONFIG.annotation_fontsize
        if fontsize is None
        else fontsize
    )

    y_min, y_max = ax.get_ylim()
    span = max(y_max - y_min, 1e-9)
    offset = span * offset_fraction

    for bar in bars:
        height = float(bar.get_height())

        if not math.isfinite(height):
            continue

        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            height + offset,
            fmt.format(height),
            ha="center",
            va="bottom",
            fontsize=font_size,
        )


def add_metric_box(
    ax: Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    value: str,
    subtitle: str,
    facecolor: str,
    edgecolor: str,
) -> None:
    """Add one rounded metric box to an axes."""

    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        transform=ax.transAxes,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=0.8,
        edgecolor=edgecolor,
        facecolor=facecolor,
        clip_on=False,
    )

    ax.add_patch(patch)

    ax.text(
        x + width / 2.0,
        y + height * 0.73,
        title,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=FIGURE_CONFIG.small_fontsize,
        color=PALETTE.dark,
    )

    ax.text(
        x + width / 2.0,
        y + height * 0.44,
        value,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=12.0,
        fontweight="bold",
        color=edgecolor,
    )

    ax.text(
        x + width / 2.0,
        y + height * 0.16,
        subtitle,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=FIGURE_CONFIG.small_fontsize,
        color=PALETTE.grey,
    )



# Panel A: longitudinal trajectory classes


def plot_trajectory_class_panel(
    ax: Axes,
    tables: InputTables,
    summary: ValidationSummary,
) -> None:
    """Plot counts of repository-derived longitudinal trajectory classes."""

    frame = (
        tables.trajectory_class_counts
        .sort_values(
            "protein_groups",
            ascending=True,
        )
        .copy()
    )

    labels = [
        format_trajectory_class(value)
        for value in frame["trajectory_class"]
    ]

    values = frame["protein_groups"].astype(float)

    colours = []
    for trajectory in frame["trajectory_class"].astype(str):
        lower = trajectory.lower()

        if "increase" in lower or "up" in lower:
            colours.append(PALETTE.blue)
        elif "decrease" in lower or "down" in lower:
            colours.append(PALETTE.orange)
        elif "persistent" in lower:
            colours.append(PALETTE.teal)
        elif "transient" in lower:
            colours.append(PALETTE.purple)
        else:
            colours.append(PALETTE.mid_grey)

    bars = ax.barh(
        np.arange(len(frame)),
        values,
        color=colours,
        edgecolor=PALETTE.white,
        linewidth=FIGURE_CONFIG.bar_edge_width,
    )

    ax.set_yticks(np.arange(len(frame)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("Protein groups")
    ax.set_title(
        "Longitudinal protein-response trajectories",
        loc="left",
        fontweight="bold",
    )

    remove_top_right_spines(ax)
    ax.xaxis.set_major_locator(
        MaxNLocator(
            integer=True,
            nbins=5,
        )
    )

    x_max = max(
        float(values.max()) if len(values) else 1.0,
        1.0,
    )
    ax.set_xlim(0, x_max * 1.18)

    for bar, value in zip(bars, values):
        ax.text(
            float(value) + x_max * 0.02,
            bar.get_y() + bar.get_height() / 2.0,
            f"{int(value)}",
            va="center",
            ha="left",
            fontsize=FIGURE_CONFIG.annotation_fontsize,
        )

    longitudinal = summary.longitudinal

    annotation = (
        f"Robust responses: "
        f"{longitudinal.robust_unique_proteins} proteins\n"
        f"Stringent responses: "
        f"{longitudinal.stringent_unique_proteins} proteins\n"
        f"Persistent responses: "
        f"{longitudinal.persistent_unique_proteins} proteins"
    )

    ax.text(
        0.98,
        0.03,
        annotation,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=FIGURE_CONFIG.small_fontsize,
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": PALETTE.very_light_grey,
            "edgecolor": PALETTE.light_grey,
            "linewidth": 0.7,
        },
    )

    add_panel_label(ax, "A")



# Panel B: module-level longitudinal responses




def plot_module_response_panel(
    ax: Axes,
    tables: InputTables,
    summary: ValidationSummary,
) -> None:
    """Plot module response medians for T1 and T6 relative to baseline."""

    frame = tables.module_longitudinal_statistics.copy()

    frame["comparison"] = comparison_order(
        frame["comparison"]
    )

    pivot = frame.pivot(
        index="module",
        columns="comparison",
        values="median_response_delta",
    )

    order = (
        pivot.abs()
        .max(axis=1)
        .sort_values(ascending=True)
        .index
    )

    pivot = pivot.loc[order]

    y_positions = np.arange(len(pivot))
    bar_height = 0.34

    t1_values = (
        pivot["T1_vs_T0"]
        if "T1_vs_T0" in pivot.columns
        else pd.Series(np.nan, index=pivot.index)
    )

    t6_values = (
        pivot["T6_vs_T0"]
        if "T6_vs_T0" in pivot.columns
        else pd.Series(np.nan, index=pivot.index)
    )

    ax.barh(
        y_positions - bar_height / 2.0,
        t1_values,
        height=bar_height,
        color=PALETTE.blue,
        edgecolor=PALETTE.white,
        linewidth=FIGURE_CONFIG.bar_edge_width,
        label="T1 − T0",
    )

    ax.barh(
        y_positions + bar_height / 2.0,
        t6_values,
        height=bar_height,
        color=PALETTE.teal,
        edgecolor=PALETTE.white,
        linewidth=FIGURE_CONFIG.bar_edge_width,
        label="T6 − T0",
    )

    ax.set_yticks(y_positions)

    ax.set_yticklabels(
        [
            format_module_name(module)
            for module in pivot.index
        ]
    )

    ax.tick_params(
        axis="y",
        pad=5,
    )

    ax.set_xlabel(
        "Median paired module-score change"
    )

    ax.set_title(
        "Biological module responses",
        loc="left",
        fontweight="bold",
        pad=7,
    )

    # Dedicated upper band for the legend.
    ax.set_ylim(
        -0.65,
        len(pivot) + 0.95,
    )

    finite_values = np.concatenate(
        [
            np.asarray(t1_values, dtype=float),
            np.asarray(t6_values, dtype=float),
        ]
    )

    finite_values = finite_values[
        np.isfinite(finite_values)
    ]

    if finite_values.size:
        minimum_value = min(
            0.0,
            float(finite_values.min()),
        )

        maximum_value = max(
            0.0,
            float(finite_values.max()),
        )

        data_range = max(
            maximum_value - minimum_value,
            0.1,
        )

        ax.set_xlim(
            minimum_value - data_range * 0.025,
            maximum_value + data_range * 0.10,
        )

    add_zero_line(ax)
    remove_top_right_spines(ax)

    legend = ax.legend(
        loc="upper right",
        bbox_to_anchor=(0.995, 0.955),
        frameon=True,
        fontsize=FIGURE_CONFIG.small_fontsize,
        handlelength=1.35,
        handletextpad=0.42,
        ncol=2,
        columnspacing=1.05,
        borderaxespad=0.0,
        borderpad=0.25,
    )

    legend.get_frame().set_facecolor(
        PALETTE.white
    )

    legend.get_frame().set_edgecolor(
        PALETTE.very_light_grey
    )

    legend.get_frame().set_linewidth(0.5)
    legend.get_frame().set_alpha(0.96)

    add_panel_label(ax, "B")


# Panel C: leave-one-protein-out robustness



def plot_module_robustness_panel(
    ax: Axes,
    tables: InputTables,
) -> None:
    """Plot minimum direction consistency after single-protein removal."""

    frame = (
        tables.leave_one_out_summary
        .copy()
    )

    frame["comparison"] = comparison_order(
        frame["comparison"]
    )

    pivot = frame.pivot(
        index="module",
        columns="comparison",
        values="minimum_direction_consistency",
    )

    order = (
        pivot.mean(axis=1)
        .sort_values(ascending=True)
        .index
    )

    pivot = pivot.loc[order]

    y = np.arange(len(pivot))
    offset = 0.10

    t1 = (
        pivot["T1_vs_T0"]
        if "T1_vs_T0" in pivot.columns
        else pd.Series(np.nan, index=pivot.index)
    )

    t6 = (
        pivot["T6_vs_T0"]
        if "T6_vs_T0" in pivot.columns
        else pd.Series(np.nan, index=pivot.index)
    )

    ax.scatter(
        t1,
        y - offset,
        s=FIGURE_CONFIG.marker_size,
        color=PALETTE.blue,
        edgecolor=PALETTE.white,
        linewidth=0.5,
        label="T1 − T0",
        zorder=3,
    )

    ax.scatter(
        t6,
        y + offset,
        s=FIGURE_CONFIG.marker_size,
        color=PALETTE.teal,
        edgecolor=PALETTE.white,
        linewidth=0.5,
        label="T6 − T0",
        zorder=3,
    )

    for index, module in enumerate(pivot.index):
        values = [
            finite_or_zero(t1.loc[module]),
            finite_or_zero(t6.loc[module]),
        ]

        ax.plot(
            values,
            [index - offset, index + offset],
            color=PALETTE.light_grey,
            linewidth=1.2,
            zorder=1,
        )

    ax.axvline(
        0.8,
        color=PALETTE.orange,
        linestyle="--",
        linewidth=0.8,
        label="0.80 consistency",
    )

    ax.set_xlim(0.74, 1.015)
    ax.set_xticks(
        [
            0.75,
            0.80,
            0.85,
            0.90,
            0.95,
            1.00,
        ]
    )
    ax.set_xticklabels(
        [
            "75%",
            "80%",
            "85%",
            "90%",
            "95%",
            "100%",
        ]
    )
    ax.set_yticks(y)
    ax.set_yticklabels(
        [
            format_module_name(module)
            for module in pivot.index
        ]
    )
    ax.set_xlabel(
        "Minimum direction consistency after removal"
    )
    ax.set_title(
        "Leave-one-protein-out sensitivity",
        loc="left",
        fontweight="bold",
    )

    remove_top_right_spines(ax)

    ax.legend(
        loc="lower right",
        frameon=False,
        fontsize=FIGURE_CONFIG.small_fontsize,
        handlelength=1.5,
    )

    add_panel_label(ax, "C")



# Panel D: integrated evidence composition





def plot_evidence_composition_panel(
    ax: Axes,
    tables: InputTables,
    summary: ValidationSummary,
) -> None:
    """Visualise the progressive integration of candidate evidence."""

    final = tables.final_candidate_panel.copy()

    dataset_counts = (
        final.groupby("dataset")["feature_id"]
        .nunique()
        .sort_values(ascending=False)
    )

    feature_counts = (
        final.groupby("feature_type")["feature_id"]
        .nunique()
        .sort_values(ascending=False)
    )

    ax.axis("off")

    ax.set_title(
        "Integrated candidate evidence",
        loc="left",
        fontweight="bold",
        pad=8,
    )

    cross = summary.cross_cohort

    box_width = 0.205
    box_height = 0.22
    box_y = 0.69

    x_positions = [
        0.00,
        0.265,
        0.530,
        0.795,
    ]

    metric_specs = [
        (
            "Evidence records",
            str(cross.evidence_records),
            "all observations",
            PALETTE.light_blue,
            PALETTE.blue,
        ),
        (
            "Aggregated features",
            str(cross.aggregated_features),
            "unique feature IDs",
            PALETTE.light_teal,
            PALETTE.teal,
        ),
        (
            "Final candidates",
            str(cross.final_candidates),
            "integrated panel",
            PALETTE.light_green,
            PALETTE.green,
        ),
        (
            "High confidence",
            str(cross.high_confidence_candidates),
            "candidate tier",
            PALETTE.light_gold,
            PALETTE.gold,
        ),
    ]

    for x_position, specification in zip(
        x_positions,
        metric_specs,
    ):
        add_metric_box(
            ax=ax,
            x=x_position,
            y=box_y,
            width=box_width,
            height=box_height,
            title=specification[0],
            value=specification[1],
            subtitle=specification[2],
            facecolor=specification[3],
            edgecolor=specification[4],
        )

    for index in range(len(x_positions) - 1):
        arrow_start = (
            x_positions[index]
            + box_width
            + 0.008
        )
        arrow_end = (
            x_positions[index + 1]
            - 0.008
        )

        ax.annotate(
            "",
            xy=(arrow_end, box_y + box_height / 2.0),
            xytext=(arrow_start, box_y + box_height / 2.0),
            xycoords=ax.transAxes,
            textcoords=ax.transAxes,
            arrowprops={
                "arrowstyle": "-|>",
                "color": PALETTE.grey,
                "linewidth": 1.0,
                "shrinkA": 0,
                "shrinkB": 0,
            },
        )

    process_labels = [
        "harmonise",
        "aggregate",
        "prioritise",
    ]

    for index, label in enumerate(process_labels):
        midpoint = (
            x_positions[index]
            + box_width
            + x_positions[index + 1]
        ) / 2.0

        ax.text(
            midpoint,
            box_y - 0.035,
            label,
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=FIGURE_CONFIG.small_fontsize,
            color=PALETTE.grey,
        )

    ax.plot(
        [0.02, 0.98],
        [0.59, 0.59],
        transform=ax.transAxes,
        color=PALETTE.light_grey,
        linewidth=0.8,
    )

    dataset_text = "\n".join(
        f"{format_dataset(dataset)}: {int(count)}"
        for dataset, count in dataset_counts.items()
    )

    feature_text = "\n".join(
        f"{format_feature_type(feature_type)}: {int(count)}"
        for feature_type, count in feature_counts.items()
    )

    ax.text(
        0.02,
        0.52,
        "Candidate representation by dataset",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=FIGURE_CONFIG.annotation_fontsize,
        fontweight="bold",
    )

    ax.text(
        0.02,
        0.45,
        dataset_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=FIGURE_CONFIG.small_fontsize,
        linespacing=1.40,
    )

    ax.text(
        0.53,
        0.52,
        "Candidate representation by feature type",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=FIGURE_CONFIG.annotation_fontsize,
        fontweight="bold",
    )

    ax.text(
        0.53,
        0.45,
        feature_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=FIGURE_CONFIG.small_fontsize,
        linespacing=1.40,
    )

    footer = (
        f"Median reproducibility score "
        f"{format_float(cross.median_reproducibility_score, 1)}"
        f"  •  Maximum "
        f"{format_float(cross.maximum_reproducibility_score, 1)}"
    )

    ax.text(
        0.02,
        0.04,
        footer,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=FIGURE_CONFIG.small_fontsize,
        color=PALETTE.grey,
        bbox={
            "boxstyle": "round,pad=0.28",
            "facecolor": PALETTE.very_light_grey,
            "edgecolor": PALETTE.light_grey,
            "linewidth": 0.6,
        },
    )

    add_panel_label(
        ax,
        "D",
        x=-0.05,
        y=1.08,
    )


# Panel E: top integrated candidates



def select_top_integrated_candidates(
    tables: InputTables,
    top_n: int = 12,
) -> pd.DataFrame:
    """Select top final candidates by integrated reproducibility score."""

    frame = (
        tables.final_candidate_panel
        .copy()
    )

    frame = frame.sort_values(
        [
            "aggregated_reproducibility_score",
            "within_feature_type_percentile",
            "maximum_absolute_effect",
        ],
        ascending=[
            False,
            False,
            False,
        ],
        na_position="last",
    )

    frame = frame.drop_duplicates(
        subset=["feature_id"],
        keep="first",
    )

    return frame.head(top_n).copy()







def plot_top_candidate_panel(
    ax: Axes,
    tables: InputTables,
    top_n: int = 12,
) -> None:
    """Plot top integrated candidate reproducibility scores."""

    frame = select_top_integrated_candidates(
        tables,
        top_n=top_n,
    )

    frame = frame.sort_values(
        "aggregated_reproducibility_score",
        ascending=True,
    )

    labels = [
        str(name)
        for name in frame["feature_name"]
    ]

    colours_by_type = {
        "mature_miRNA": PALETTE.purple,
        "miRNA_associated_locus": PALETTE.light_purple,
        "probable_miRNA_host_gene": PALETTE.gold,
        "protein": PALETTE.blue,
        "protein_module": PALETTE.teal,
    }

    colours = [
        colours_by_type.get(
            str(feature_type),
            PALETTE.grey,
        )
        for feature_type in frame["feature_type"]
    ]

    values = frame[
        "aggregated_reproducibility_score"
    ].astype(float)

    y_positions = np.arange(len(frame))

    bars = ax.barh(
        y_positions,
        values,
        height=0.74,
        color=colours,
        edgecolor=PALETTE.white,
        linewidth=FIGURE_CONFIG.bar_edge_width,
    )

    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels)

    ax.tick_params(
        axis="y",
        pad=4,
    )

    ax.set_xlabel(
        "Aggregated reproducibility score"
    )

    ax.set_title(
        "Highest-ranked integrated candidates",
        loc="left",
        fontweight="bold",
        pad=7,
    )

    remove_top_right_spines(ax)

    # A distinct upper legend band above the highest candidate.
    ax.set_ylim(
        -0.75,
        len(frame) + 1.65,
    )

    x_min = 0.0

    x_max = max(
        112.0,
        float(values.max()) + 11.0,
    )

    ax.set_xlim(
        x_min,
        x_max,
    )

    for bar, score, tier in zip(
        bars,
        values,
        frame["candidate_tier"],
    ):
        ax.text(
            float(score) + 1.45,
            bar.get_y() + bar.get_height() / 2.0,
            f"{float(score):.1f}",
            va="center",
            ha="left",
            fontsize=FIGURE_CONFIG.small_fontsize,
            color=PALETTE.dark,
        )

        if str(tier) == "high_confidence":
            ax.text(
                1.35,
                bar.get_y() + bar.get_height() / 2.0,
                "●",
                va="center",
                ha="center",
                fontsize=5.7,
                color=PALETTE.gold,
                zorder=4,
            )

    observed_types = (
        frame["feature_type"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    preferred_type_order = [
        "protein",
        "protein_module",
        "miRNA_associated_locus",
        "mature_miRNA",
        "probable_miRNA_host_gene",
    ]

    observed_types = [
        feature_type
        for feature_type in preferred_type_order
        if feature_type in observed_types
    ]

    legend_handles = []
    legend_labels = []

    for feature_type in observed_types:
        legend_handles.append(
            mpl.lines.Line2D(
                [0],
                [0],
                marker="s",
                linestyle="",
                markerfacecolor=colours_by_type.get(
                    feature_type,
                    PALETTE.grey,
                ),
                markeredgecolor="none",
                markersize=5.5,
            )
        )

        legend_labels.append(
            format_feature_type(feature_type)
        )

    high_confidence_handle = mpl.lines.Line2D(
        [0],
        [0],
        marker="o",
        linestyle="",
        markerfacecolor=PALETTE.gold,
        markeredgecolor="none",
        markersize=5.5,
    )

    legend_handles.append(high_confidence_handle)
    legend_labels.append("High-confidence candidate")

    legend = ax.legend(
        legend_handles,
        legend_labels,
        loc="upper right",
        bbox_to_anchor=(0.995, 0.975),
        frameon=True,
        fontsize=FIGURE_CONFIG.small_fontsize,
        ncol=2,
        columnspacing=1.15,
        handletextpad=0.35,
        borderaxespad=0.0,
        borderpad=0.28,
    )

    legend.get_frame().set_facecolor(
        PALETTE.white
    )

    legend.get_frame().set_edgecolor(
        PALETTE.very_light_grey
    )

    legend.get_frame().set_linewidth(0.5)
    legend.get_frame().set_alpha(0.96)

    add_panel_label(ax, "E")


# Panel F: pathway enrichment and regulatory integration








def plot_pathway_summary_panel(
    ax: Axes,
    summary: ValidationSummary,
) -> None:
    """Summarise pathway enrichment and regulatory support."""

    ax.axis("off")

    pathway = summary.pathways
    regulatory = summary.regulatory

    ax.set_title(
        "Biological interpretation",
        loc="left",
        fontweight="bold",
        pad=8,
    )

    metrics = [
        (
            "Significant\npathways",
            str(pathway.significant_pathway_terms),
            PALETTE.light_blue,
            PALETTE.blue,
        ),
        (
            "Candidate\nnetwork nodes",
            str(pathway.candidate_network_nodes),
            PALETTE.light_teal,
            PALETTE.teal,
        ),
        (
            "Matched\ncandidate proteins",
            str(regulatory.matched_candidate_proteins),
            PALETTE.light_green,
            PALETTE.green,
        ),
        (
            "Maximum\nenrichment",
            f"{format_float(pathway.maximum_fold_enrichment,1)}×",
            PALETTE.light_gold,
            PALETTE.gold,
        ),
    ]

    box_width = 0.215
    box_height = 0.225
    box_y = 0.715

    x_positions = [
        0.00,
        0.262,
        0.524,
        0.786,
    ]

    for x_position, metric in zip(
        x_positions,
        metrics,
    ):
        title, value, facecolor, edgecolor = metric

        patch = FancyBboxPatch(
            (x_position, box_y),
            box_width,
            box_height,
            transform=ax.transAxes,
            boxstyle=(
                "round,pad=0.012,"
                "rounding_size=0.018"
            ),
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=1.0,
            clip_on=False,
        )

        ax.add_patch(patch)

        # Title kept in the upper half.
        ax.text(
            x_position + box_width / 2.0,
            box_y + box_height * 0.77,
            title,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=FIGURE_CONFIG.small_fontsize,
            color=PALETTE.dark,
            linespacing=1.02,
        )

        # Value placed clearly in the lower half.
        ax.text(
            x_position + box_width / 2.0,
            box_y + box_height * 0.29,
            value,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=FIGURE_CONFIG.annotation_fontsize,
            fontweight="bold",
            color=edgecolor,
        )

    ax.plot(
        [0.02, 0.98],
        [0.655, 0.655],
        transform=ax.transAxes,
        color=PALETTE.light_grey,
        linewidth=0.8,
    )

    ax.text(
        0.02,
        0.598,
        "Principal enriched biological process",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=FIGURE_CONFIG.annotation_fontsize,
        fontweight="bold",
    )

    if pathway.primary_pathways:
        primary_names = [
            str(name)
            for name in pathway.primary_pathways[:5]
        ]

        if len(primary_names) == 1:
            pathway_name = wrap_label(
                primary_names[0],
                34,
            )

            callout = FancyBboxPatch(
                (0.02, 0.305),
                0.53,
                0.205,
                transform=ax.transAxes,
                boxstyle=(
                    "round,pad=0.018,"
                    "rounding_size=0.018"
                ),
                facecolor=PALETTE.very_light_grey,
                edgecolor=PALETTE.light_blue,
                linewidth=0.9,
            )

            ax.add_patch(callout)

            circle = mpl.patches.Circle(
                (0.075, 0.407),
                radius=0.027,
                transform=ax.transAxes,
                facecolor=PALETTE.light_blue,
                edgecolor=PALETTE.blue,
                linewidth=0.9,
            )

            ax.add_patch(circle)

            ax.text(
                0.075,
                0.407,
                "1",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=FIGURE_CONFIG.annotation_fontsize,
                fontweight="bold",
                color=PALETTE.blue,
            )

            ax.text(
                0.125,
                0.415,
                pathway_name,
                transform=ax.transAxes,
                ha="left",
                va="center",
                fontsize=FIGURE_CONFIG.annotation_fontsize,
                fontweight="bold",
                color=PALETTE.dark,
            )


        else:
            row_y_positions = np.linspace(
                0.50,
                0.16,
                len(primary_names),
            )

            for index, (name, row_y) in enumerate(
                zip(primary_names, row_y_positions),
                start=1,
            ):
                circle = mpl.patches.Circle(
                    (0.045, row_y),
                    radius=0.016,
                    transform=ax.transAxes,
                    facecolor=PALETTE.light_blue,
                    edgecolor=PALETTE.blue,
                    linewidth=0.8,
                )

                ax.add_patch(circle)

                ax.text(
                    0.045,
                    row_y,
                    str(index),
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=FIGURE_CONFIG.small_fontsize,
                    fontweight="bold",
                    color=PALETTE.blue,
                )

                ax.text(
                    0.080,
                    row_y,
                    wrap_label(name, 34),
                    transform=ax.transAxes,
                    ha="left",
                    va="center",
                    fontsize=FIGURE_CONFIG.small_fontsize,
                )

    else:
        ax.text(
            0.02,
            0.49,
            "No significant pathways available.",
            transform=ax.transAxes,
            fontsize=FIGURE_CONFIG.small_fontsize,
            va="top",
            color=PALETTE.grey,
        )

    ax.plot(
        [0.605, 0.605],
        [0.10, 0.595],
        transform=ax.transAxes,
        color=PALETTE.light_grey,
        linewidth=0.8,
    )

    ax.text(
        0.645,
        0.598,
        "Biological evidence summary",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=FIGURE_CONFIG.annotation_fontsize,
        fontweight="bold",
    )

    evidence_rows = [
        (
            "Minimum adjusted P value",
            format_float(
                pathway.minimum_adjusted_p_value,
                4,
            ),
        ),
        (
            "Candidate-protein overlap",
            format_percent(
                regulatory.candidate_protein_overlap_fraction
            ),
        ),
        (
            "Matched candidate proteins",
            str(
                regulatory.matched_candidate_proteins
            ),
        ),
        (
            "Significant pathways",
            str(
                pathway.significant_pathway_terms
            ),
        ),
    ]

    evidence_y_positions = [
        0.485,
        0.380,
        0.275,
        0.170,
    ]

    for (label, value), row_y in zip(
        evidence_rows,
        evidence_y_positions,
    ):
        ax.text(
            0.645,
            row_y,
            label,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=FIGURE_CONFIG.small_fontsize,
            color=PALETTE.grey,
        )

        ax.text(
            0.965,
            row_y,
            value,
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=FIGURE_CONFIG.annotation_fontsize,
            fontweight="bold",
            color=PALETTE.dark,
        )

        ax.plot(
            [0.645, 0.965],
            [row_y - 0.052, row_y - 0.052],
            transform=ax.transAxes,
            color=PALETTE.very_light_grey,
            linewidth=0.9,
        )

    add_panel_label(ax, "F")


# Figure assembly



def build_final_validation_figure(
    tables: InputTables,
    summary: ValidationSummary,
) -> Figure:
    """
    Assemble the complete publication figure.
    """

    configure_matplotlib()

    figure = plt.figure(
        figsize=(
            FIGURE_CONFIG.width_inches,
            FIGURE_CONFIG.height_inches,
        ),
        constrained_layout=False,
    )

    grid = figure.add_gridspec(
        nrows=3,
        ncols=2,
        left=0.055,
        right=0.985,
        top=0.90,
        bottom=0.075,
        hspace=0.42,
        wspace=0.24,
        height_ratios=[
            1.0,
            1.0,
            0.92,
        ],
    )

    ax_a = figure.add_subplot(grid[0, 0])
    ax_b = figure.add_subplot(grid[0, 1])

    ax_c = figure.add_subplot(grid[1, 0])
    ax_d = figure.add_subplot(grid[1, 1])

    ax_e = figure.add_subplot(grid[2, 0])
    ax_f = figure.add_subplot(grid[2, 1])

    plot_trajectory_class_panel(
        ax_a,
        tables,
        summary,
    )

    plot_module_response_panel(
        ax_b,
        tables,
        summary,
    )

    plot_module_robustness_panel(
        ax_c,
        tables,
    )

    plot_evidence_composition_panel(
        ax_d,
        tables,
        summary,
    )

    plot_top_candidate_panel(
        ax_e,
        tables,
    )

    plot_pathway_summary_panel(
        ax_f,
        summary,
    )

    figure.suptitle(
        (
            "Integrated validation framework for longitudinal "
            "salivary biomarker prioritisation"
        ),
        fontsize=FIGURE_CONFIG.title_fontsize,
        fontweight="bold",
        y=0.965,
    )

    figure.text(
        0.5,
        0.935,
        (
            "Longitudinal protein responses were integrated with "
            "robustness, cross-cohort, regulatory, and pathway-level evidence."
        ),
        ha="center",
        fontsize=FIGURE_CONFIG.subtitle_fontsize,
        color=PALETTE.grey,
    )

    return figure



# Caption generation



def build_caption(
    summary: ValidationSummary,
) -> str:
    """
    Generate the manuscript figure legend.
    """

    longitudinal = summary.longitudinal
    modules = summary.modules
    cross = summary.cross_cohort
    regulatory = summary.regulatory
    pathways = summary.pathways

    caption = f"""
Figure X. Integrated validation framework for longitudinal saliva biomarker
prioritisation.

(A) Distribution of repository-derived longitudinal protein trajectory
classes identified from paired analyses of the discovery cohort.
The analysis included
{longitudinal.total_human_proteins} human proteins,
of which
{longitudinal.robust_unique_proteins}
showed robust longitudinal responses,
{longitudinal.stringent_unique_proteins}
fulfilled stringent statistical criteria,
and
{longitudinal.persistent_unique_proteins}
displayed persistent responses across time.

(B) Module-level longitudinal responses summarising median paired
changes between baseline and post-treatment measurements.

(C) Leave-one-protein-out robustness analysis demonstrating that
module behaviour remained stable following systematic removal of
individual proteins.

(D) Cross-cohort evidence integration summarising
{cross.evidence_records}
individual evidence records,
{cross.aggregated_features}
aggregated biological features,
and
{cross.final_candidates}
final integrated candidates.

(E) Highest-ranked integrated biomarker candidates according to the
repository-derived aggregated reproducibility score.

(F) Summary of pathway enrichment and regulatory-network support,
including
{pathways.significant_pathway_terms}
significant biological pathways,
{regulatory.high_confidence_regulatory_edges}
high-confidence regulatory interactions,
and
{regulatory.matched_candidate_proteins}
candidate proteins with validated regulatory evidence.

All values shown originate directly from repository-derived analyses.
No external quantitative values or manually generated estimates were
introduced during figure generation.
"""

    return textwrap.dedent(caption).strip()



# Report generation



def build_report(
    summary: ValidationSummary,
) -> str:
    """
    Produce a human-readable report accompanying the figure.
    """

    longitudinal = summary.longitudinal
    modules = summary.modules
    cross = summary.cross_cohort
    regulatory = summary.regulatory
    pathways = summary.pathways

    lines = []

    lines.append(
        "FINAL VALIDATION FRAMEWORK REPORT"
    )
    lines.append("=" * 70)
    lines.append("")

    lines.append("LONGITUDINAL ANALYSIS")
    lines.append("-" * 70)
    lines.append(
        f"Human proteins: "
        f"{longitudinal.total_human_proteins}"
    )
    lines.append(
        f"Robust proteins: "
        f"{longitudinal.robust_unique_proteins}"
    )
    lines.append(
        f"Persistent proteins: "
        f"{longitudinal.persistent_unique_proteins}"
    )
    lines.append(
        f"Median consistency T1: "
        f"{format_percent(longitudinal.median_consistency_t1)}"
    )
    lines.append(
        f"Median consistency T6: "
        f"{format_percent(longitudinal.median_consistency_t6)}"
    )

    lines.append("")
    lines.append("MODULE ANALYSIS")
    lines.append("-" * 70)
    lines.append(
        f"Modules: {modules.total_modules}"
    )
    lines.append(
        f"Robust comparisons: "
        f"{modules.robust_module_comparisons}"
    )
    lines.append(
        f"Median coverage: "
        f"{format_percent(modules.median_module_coverage)}"
    )

    if modules.modules_robust_both_timepoints:

        lines.append("")
        lines.append(
            "Robust across both timepoints:"
        )

        for module in modules.modules_robust_both_timepoints:
            lines.append(
                f"  • {format_module_name(module)}"
            )

    lines.append("")
    lines.append("CROSS-COHORT INTEGRATION")
    lines.append("-" * 70)

    lines.append(
        f"Evidence records: "
        f"{cross.evidence_records}"
    )

    lines.append(
        f"Aggregated features: "
        f"{cross.aggregated_features}"
    )

    lines.append(
        f"Final candidates: "
        f"{cross.final_candidates}"
    )

    lines.append(
        f"High-confidence candidates: "
        f"{cross.high_confidence_candidates}"
    )

    lines.append(
        f"Median reproducibility score: "
        f"{format_float(cross.median_reproducibility_score,2)}"
    )

    lines.append("")
    lines.append("REGULATORY SUPPORT")
    lines.append("-" * 70)

    lines.append(
        f"Validated interactions: "
        f"{regulatory.all_validated_regulatory_edges}"
    )

    lines.append(
        f"High-confidence interactions: "
        f"{regulatory.high_confidence_regulatory_edges}"
    )

    lines.append(
        f"Matched candidate proteins: "
        f"{regulatory.matched_candidate_proteins}"
    )

    lines.append(
        f"Protein overlap fraction: "
        f"{format_percent(regulatory.candidate_protein_overlap_fraction)}"
    )

    lines.append("")
    lines.append("PATHWAY ENRICHMENT")
    lines.append("-" * 70)

    lines.append(
        f"Significant pathways: "
        f"{pathways.significant_pathway_terms}"
    )

    lines.append(
        f"Network nodes: "
        f"{pathways.candidate_network_nodes}"
    )

    lines.append(
        f"Network edges: "
        f"{pathways.candidate_network_edges}"
    )

    if pathways.primary_pathways:

        lines.append("")
        lines.append("Primary enriched pathways:")

        for pathway in pathways.primary_pathways:
            lines.append(f"  • {pathway}")

    return "\n".join(lines)


# Command-line entry point



def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Build the final integrated validation figure from "
            "repository-derived analysis tables."
        )
    )

    parser.add_argument(
        "--repository-root",
        type=Path,
        default=Path.cwd(),
        help=(
            "Repository root directory. Default: current working directory."
        ),
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=DEFAULT_OUTPUT_DPI,
        help=(
            f"PNG output resolution. Default: {DEFAULT_OUTPUT_DPI} DPI."
        ),
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging.",
    )

    return parser.parse_args()


def save_final_outputs(
    figure: Figure,
    caption: str,
    report: str,
    paths: RepositoryPaths,
    dpi: int,
) -> None:
    """Write all final figure and text outputs."""

    paths.output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure.savefig(
        paths.output_svg,
        format="svg",
        bbox_inches="tight",
    )

    figure.savefig(
        paths.output_pdf,
        format="pdf",
        bbox_inches="tight",
    )

    figure.savefig(
        paths.output_png,
        format="png",
        dpi=dpi,
        bbox_inches="tight",
    )

    paths.output_caption.write_text(
        caption + "\n",
        encoding="utf-8",
    )

    paths.output_report.write_text(
        report + "\n",
        encoding="utf-8",
    )

    manifest = {
        "script": SCRIPT_NAME,
        "repository_root": str(paths.repository_root),
        "output_directory": str(paths.output_directory),
        "outputs": {
            "svg": str(paths.output_svg),
            "pdf": str(paths.output_pdf),
            "png": str(paths.output_png),
            "caption": str(paths.output_caption),
            "report": str(paths.output_report),
        },
        "png_dpi": dpi,
    }

    paths.output_manifest.write_text(
        json.dumps(
            manifest,
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> int:
    """Run the complete final-figure workflow."""

    arguments = parse_arguments()
    configure_logging(arguments.verbose)

    repository_root = arguments.repository_root.expanduser().resolve()
    paths = build_repository_paths(repository_root)

    logging.info(
        "Repository root: %s",
        paths.repository_root,
    )

    logging.info("Loading repository-derived input tables.")
    tables = load_input_tables(paths)

    logging.info("Validating input tables.")
    validate_all_input_tables(tables)

    logging.info("Building integrated validation summary.")
    summary = build_validation_summary(tables)

    logging.info("Building final validation figure.")
    figure = build_final_validation_figure(
        tables=tables,
        summary=summary,
    )

    caption = build_caption(summary)
    report = build_report(summary)

    logging.info("Writing outputs to %s", paths.output_directory)
    save_final_outputs(
        figure=figure,
        caption=caption,
        report=report,
        paths=paths,
        dpi=arguments.dpi,
    )

    plt.close(figure)

    logging.info("Created: %s", paths.output_svg)
    logging.info("Created: %s", paths.output_pdf)
    logging.info("Created: %s", paths.output_png)
    logging.info("Created: %s", paths.output_caption)
    logging.info("Created: %s", paths.output_report)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        logging.exception("Final validation figure generation failed.")
        raise

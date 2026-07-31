#!/usr/bin/env python3
"""
Build a harmonised evidence registry from independent oral-fluid omics analyses.

The registry does not treat proteins, protein modules, and miRNAs as directly
equivalent biological entities. Instead, it records each result as a separate
piece of evidence in a common long-format schema.

Inputs
------
1. PXD043901 robust longitudinal protein responses
2. PXD043901 longitudinal protein-module statistics
3. GSE312585 consensus method-sensitive feature signature
4. GSE121910 exploratory salivary miRNA priorities

Outputs
-------
1. cross_cohort_evidence_registry.tsv
2. cross_cohort_evidence_registry_summary.tsv
3. cross_cohort_evidence_registry_report.txt
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


REGISTRY_COLUMNS = [
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
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a harmonised cross-cohort evidence registry."
    )

    parser.add_argument(
        "--pxd-proteins",
        type=Path,
        default=Path(
            "results/PXD043901/trajectories/"
            "robust_human_responses_preliminary.csv"
        ),
    )
    parser.add_argument(
        "--pxd-modules",
        type=Path,
        default=Path(
            "results/PXD043901/modules/"
            "module_longitudinal_statistics.csv"
        ),
    )
    parser.add_argument(
        "--gse312585-signature",
        type=Path,
        default=Path(
            "results/GSE312585/tables/"
            "feature_classification.tsv"
        ),
    )
    parser.add_argument(
        "--gse121910-priorities",
        type=Path,
        default=Path(
            "results/GSE121910/tables/"
            "exploratory_mirna_priorities.tsv"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "results/integration/cross_cohort_validation/"
            "cross_cohort_evidence_registry.tsv"
        ),
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=Path(
            "results/integration/cross_cohort_validation/"
            "cross_cohort_evidence_registry_summary.tsv"
        ),
    )
    parser.add_argument(
        "--report-output",
        type=Path,
        default=Path(
            "results/integration/cross_cohort_validation/"
            "cross_cohort_evidence_registry_report.txt"
        ),
    )

    return parser.parse_args()


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Required input file not found: {path}")


def require_columns(
    frame: pd.DataFrame,
    required: Iterable[str],
    source: Path,
) -> None:
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(
            f"{source} is missing required columns: {', '.join(missing)}"
        )


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def normalise_direction(value: object) -> str:
    if pd.isna(value):
        return "unknown"

    text = str(value).strip().lower()

    increase_terms = {
        "increase",
        "increased",
        "positive",
        "higher_in_case",
        "up",
        "upregulated",
    }
    decrease_terms = {
        "decrease",
        "decreased",
        "negative",
        "higher_in_control",
        "down",
        "downregulated",
    }

    if text in increase_terms:
        return "increase"

    if text in decrease_terms:
        return "decrease"

    return "unknown"


def direction_from_effect(effect: object) -> str:
    try:
        value = float(effect)
    except (TypeError, ValueError):
        return "unknown"

    if np.isnan(value):
        return "unknown"
    if value > 0:
        return "increase"
    if value < 0:
        return "decrease"
    return "no_change"


def classify_pxd_protein_evidence(row: pd.Series) -> str:
    consistency = row.get("direction_consistency", np.nan)
    pairs = row.get("n_observed_pairs", np.nan)
    effect = abs(row.get("median_delta_log2", np.nan))

    if (
        pd.notna(consistency)
        and pd.notna(pairs)
        and pd.notna(effect)
        and consistency >= 0.90
        and pairs >= 8
        and effect >= 0.50
    ):
        return "strong"

    if (
        pd.notna(consistency)
        and pd.notna(pairs)
        and pd.notna(effect)
        and consistency >= 0.75
        and pairs >= 8
        and effect >= 0.25
    ):
        return "moderate"

    return "supportive"


def classify_module_evidence(row: pd.Series) -> str:
    fdr = row.get("wilcoxon_fdr", np.nan)
    consistency = row.get("direction_consistency", np.nan)
    effect = abs(row.get("median_response_delta", np.nan))

    if (
        pd.notna(fdr)
        and pd.notna(consistency)
        and pd.notna(effect)
        and fdr <= 0.01
        and consistency >= 0.90
        and effect >= 0.50
    ):
        return "strong"

    if (
        pd.notna(fdr)
        and pd.notna(consistency)
        and fdr <= 0.05
        and consistency >= 0.75
    ):
        return "moderate"

    return "supportive"


def classify_gse312585_evidence(row: pd.Series) -> str:
    tier = str(row.get("method_sensitivity_tier", "")).strip().lower()
    robust = row.get("robust_contrasts", np.nan)
    consistency = row.get("mean_direction_consistency", np.nan)

    if (
        tier == "pan_contrast"
        and pd.notna(robust)
        and robust >= 6
        and pd.notna(consistency)
        and consistency >= 0.95
    ):
        return "strong"

    if (
        tier in {"pan_contrast", "high"}
        and pd.notna(robust)
        and robust >= 4
        and pd.notna(consistency)
        and consistency >= 0.90
    ):
        return "moderate"

    return "supportive"


def classify_gse121910_evidence(row: pd.Series) -> str:
    prioritised = bool(row.get("exploratory_priority", False))
    effect = abs(row.get("log2FC_case_vs_control", np.nan))
    abundance = row.get("mean_tpm", np.nan)

    if (
        prioritised
        and pd.notna(effect)
        and effect >= 2.0
        and pd.notna(abundance)
        and abundance >= 100.0
    ):
        return "strong"

    if prioritised and pd.notna(effect) and effect >= 1.0:
        return "moderate"

    return "supportive"


def build_pxd_protein_registry(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)

    required = [
        "protein_group_id",
        "comparison",
        "n_observed_pairs",
        "dominant_direction",
        "direction_consistency",
        "median_delta_log2",
        "primary_accession",
        "gene_symbol",
        "human_protein_name",
        "detected_fraction",
    ]
    require_columns(frame, required, path)

    frame = frame.copy()

    frame["effect_size"] = numeric(frame["median_delta_log2"])
    frame["direction_consistency"] = numeric(
        frame["direction_consistency"]
    )
    frame["n_observed_pairs"] = numeric(frame["n_observed_pairs"])
    frame["detected_fraction"] = numeric(frame["detected_fraction"])

    feature_name = (
        frame["gene_symbol"]
        .fillna("")
        .astype(str)
        .str.strip()
        .replace("", np.nan)
        .fillna(frame["primary_accession"].astype(str))
    )

    output = pd.DataFrame(
        {
            "feature_id": frame["primary_accession"].astype(str),
            "feature_name": feature_name,
            "feature_type": "protein",
            "modality": "proteomics",
            "dataset": "PXD043901",
            "study_role": "longitudinal_treatment_response",
            "comparison": frame["comparison"].astype(str),
            "effect_size": frame["effect_size"],
            "effect_metric": "median_delta_log2",
            "effect_direction": frame["dominant_direction"].map(
                normalise_direction
            ),
            "direction_consistency": frame["direction_consistency"],
            "robustness_value": frame["direction_consistency"],
            "robustness_metric": "paired_direction_consistency",
            "statistic": np.nan,
            "p_value": np.nan,
            "adjusted_p_value": np.nan,
            "n_observations": frame["n_observed_pairs"],
            "abundance_or_detection": frame["detected_fraction"],
            "evidence_level": frame.apply(
                classify_pxd_protein_evidence,
                axis=1,
            ),
            "is_prioritised": True,
            "source_file": str(path),
        }
    )

    return output


def build_pxd_module_registry(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)

    required = [
        "module",
        "comparison",
        "n_pairs",
        "direction_consistency",
        "median_response_delta",
        "wilcoxon_statistic",
        "wilcoxon_p",
        "wilcoxon_fdr",
    ]
    require_columns(frame, required, path)

    frame = frame.copy()

    numeric_columns = [
        "n_pairs",
        "direction_consistency",
        "median_response_delta",
        "wilcoxon_statistic",
        "wilcoxon_p",
        "wilcoxon_fdr",
    ]
    for column in numeric_columns:
        frame[column] = numeric(frame[column])

    output = pd.DataFrame(
        {
            "feature_id": "module:" + frame["module"].astype(str),
            "feature_name": frame["module"].astype(str),
            "feature_type": "protein_module",
            "modality": "proteomics",
            "dataset": "PXD043901",
            "study_role": "longitudinal_module_response",
            "comparison": frame["comparison"].astype(str),
            "effect_size": frame["median_response_delta"],
            "effect_metric": "median_module_delta",
            "effect_direction": frame["median_response_delta"].map(
                direction_from_effect
            ),
            "direction_consistency": frame["direction_consistency"],
            "robustness_value": frame["direction_consistency"],
            "robustness_metric": "paired_direction_consistency",
            "statistic": frame["wilcoxon_statistic"],
            "p_value": frame["wilcoxon_p"],
            "adjusted_p_value": frame["wilcoxon_fdr"],
            "n_observations": frame["n_pairs"],
            "abundance_or_detection": np.nan,
            "evidence_level": frame.apply(
                classify_module_evidence,
                axis=1,
            ),
            "is_prioritised": frame["wilcoxon_fdr"].le(0.05),
            "source_file": str(path),
        }
    )

    return output


def build_gse312585_registry(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")

    required = [
        "consensus_rank",
        "feature",
        "feature_class",
        "method_sensitivity_tier",
        "robust_contrasts",
        "robust_contrast_names",
        "maximum_absolute_effect",
        "median_absolute_effect",
        "mean_direction_consistency",
        "contrasts_evaluated",
    ]
    require_columns(frame, required, path)

    frame = frame.copy()

    numeric_columns = [
        "consensus_rank",
        "robust_contrasts",
        "maximum_absolute_effect",
        "median_absolute_effect",
        "mean_direction_consistency",
        "contrasts_evaluated",
    ]
    for column in numeric_columns:
        frame[column] = numeric(frame[column])

    output = pd.DataFrame(
        {
            "feature_id": frame["feature"].astype(str),
            "feature_name": frame["feature"].astype(str),
            "feature_type": frame["feature_class"].astype(str),
            "modality": "small_RNA_sequencing",
            "dataset": "GSE312585",
            "study_role": "technical_method_robustness",
            "comparison": frame["robust_contrast_names"]
            .fillna("")
            .astype(str),
            "effect_size": frame["median_absolute_effect"],
            "effect_metric": "median_absolute_method_effect",
            # This study quantifies method effects rather than disease direction.
            "effect_direction": "method_dependent",
            "direction_consistency": frame[
                "mean_direction_consistency"
            ],
            "robustness_value": (
                frame["robust_contrasts"]
                / frame["contrasts_evaluated"].replace(0, np.nan)
            ),
            "robustness_metric": "robust_contrast_fraction",
            "statistic": np.nan,
            "p_value": np.nan,
            "adjusted_p_value": np.nan,
            "n_observations": frame["contrasts_evaluated"],
            "abundance_or_detection": np.nan,
            "evidence_level": frame.apply(
                classify_gse312585_evidence,
                axis=1,
            ),
            "is_prioritised": frame[
                "method_sensitivity_tier"
            ].isin(["pan_contrast", "high"]),
            "source_file": str(path),
        }
    )

    return output


def build_gse121910_registry(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")

    required = [
        "Known miRNA",
        "n_unique_sequences",
        "log2FC_case_vs_control",
        "mean_tpm",
        "exploratory_priority",
        "direction",
    ]
    require_columns(frame, required, path)

    frame = frame.copy()

    numeric_columns = [
        "n_unique_sequences",
        "log2FC_case_vs_control",
        "mean_tpm",
    ]
    for column in numeric_columns:
        frame[column] = numeric(frame[column])

    frame["exploratory_priority"] = (
        frame["exploratory_priority"]
        .astype(str)
        .str.lower()
        .map({"true": True, "false": False})
        .fillna(False)
    )

    output = pd.DataFrame(
        {
            "feature_id": frame["Known miRNA"].astype(str),
            "feature_name": frame["Known miRNA"].astype(str),
            "feature_type": "mature_miRNA",
            "modality": "small_RNA_sequencing",
            "dataset": "GSE121910",
            "study_role": "case_control_biological_signal",
            "comparison": "case_vs_control",
            "effect_size": frame["log2FC_case_vs_control"],
            "effect_metric": "log2FC_case_vs_control",
            "effect_direction": frame["direction"].map(
                normalise_direction
            ),
            "direction_consistency": np.nan,
            "robustness_value": frame["n_unique_sequences"],
            "robustness_metric": "supporting_unique_sequences",
            "statistic": np.nan,
            "p_value": np.nan,
            "adjusted_p_value": np.nan,
            "n_observations": frame["n_unique_sequences"],
            "abundance_or_detection": frame["mean_tpm"],
            "evidence_level": frame.apply(
                classify_gse121910_evidence,
                axis=1,
            ),
            "is_prioritised": frame["exploratory_priority"],
            "source_file": str(path),
        }
    )

    return output


def add_evidence_ids(registry: pd.DataFrame) -> pd.DataFrame:
    registry = registry.reset_index(drop=True).copy()

    registry.insert(
        0,
        "evidence_id",
        [
            f"EVIDENCE_{index:05d}"
            for index in range(1, len(registry) + 1)
        ],
    )

    return registry


def build_summary(registry: pd.DataFrame) -> pd.DataFrame:
    summary = (
        registry.groupby(
            [
                "dataset",
                "study_role",
                "modality",
                "feature_type",
                "evidence_level",
            ],
            dropna=False,
        )
        .agg(
            n_evidence_records=("evidence_id", "size"),
            n_unique_features=("feature_id", "nunique"),
            n_prioritised=("is_prioritised", "sum"),
            median_absolute_effect=(
                "effect_size",
                lambda values: pd.to_numeric(
                    values,
                    errors="coerce",
                ).abs().median(),
            ),
            median_direction_consistency=(
                "direction_consistency",
                "median",
            ),
        )
        .reset_index()
        .sort_values(
            [
                "dataset",
                "study_role",
                "feature_type",
                "evidence_level",
            ]
        )
    )

    return summary


def write_report(
    registry: pd.DataFrame,
    summary: pd.DataFrame,
    output_path: Path,
) -> None:
    lines = [
        "Cross-cohort evidence registry",
        "==============================",
        "",
        f"Total evidence records: {len(registry)}",
        f"Unique feature identifiers: {registry['feature_id'].nunique()}",
        f"Datasets represented: {registry['dataset'].nunique()}",
        "",
        "Evidence records by dataset:",
    ]

    dataset_counts = registry["dataset"].value_counts().sort_index()
    for dataset, count in dataset_counts.items():
        lines.append(f"  {dataset}: {count}")

    lines.extend(
        [
            "",
            "Evidence records by modality:",
        ]
    )

    modality_counts = registry["modality"].value_counts().sort_index()
    for modality, count in modality_counts.items():
        lines.append(f"  {modality}: {count}")

    lines.extend(
        [
            "",
            "Evidence levels:",
        ]
    )

    level_counts = registry["evidence_level"].value_counts()
    for level in ["strong", "moderate", "supportive"]:
        lines.append(f"  {level}: {int(level_counts.get(level, 0))}")

    lines.extend(
        [
            "",
            "Scientific interpretation notes:",
            (
                "  - Evidence records are harmonised structurally, not treated "
                "as biologically interchangeable."
            ),
            (
                "  - GSE312585 represents technical method sensitivity and "
                "should not be interpreted as disease-direction evidence."
            ),
            (
                "  - GSE121910 is an exploratory case-control signal and does "
                "not contain inferential p-values in the source table."
            ),
            (
                "  - Cross-modality integration should subsequently occur at "
                "the target, pathway, or biological-module level."
            ),
            "",
            f"Summary rows: {len(summary)}",
            "",
        ]
    )

    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()

    input_paths = [
        args.pxd_proteins,
        args.pxd_modules,
        args.gse312585_signature,
        args.gse121910_priorities,
    ]

    for path in input_paths:
        require_file(path)

    registry_parts = [
        build_pxd_protein_registry(args.pxd_proteins),
        build_pxd_module_registry(args.pxd_modules),
        build_gse312585_registry(args.gse312585_signature),
        build_gse121910_registry(args.gse121910_priorities),
    ]

    registry = pd.concat(
        registry_parts,
        ignore_index=True,
        sort=False,
    )

    registry = add_evidence_ids(registry)
    registry = registry[REGISTRY_COLUMNS]

    registry = registry.sort_values(
        [
            "dataset",
            "modality",
            "feature_type",
            "feature_name",
            "comparison",
        ],
        na_position="last",
    ).reset_index(drop=True)

    # Rebuild stable IDs after deterministic sorting.
    registry["evidence_id"] = [
        f"EVIDENCE_{index:05d}"
        for index in range(1, len(registry) + 1)
    ]

    summary = build_summary(registry)

    for output_path in [
        args.output,
        args.summary_output,
        args.report_output,
    ]:
        output_path.parent.mkdir(parents=True, exist_ok=True)

    registry.to_csv(args.output, sep="\t", index=False)
    summary.to_csv(args.summary_output, sep="\t", index=False)
    write_report(registry, summary, args.report_output)

    print(f"Wrote registry: {args.output}")
    print(f"Wrote summary:  {args.summary_output}")
    print(f"Wrote report:   {args.report_output}")
    print(f"Evidence records: {len(registry)}")
    print(f"Unique features:  {registry['feature_id'].nunique()}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)

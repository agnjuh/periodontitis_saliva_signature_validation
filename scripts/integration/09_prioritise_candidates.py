#!/usr/bin/env python3
"""
Prioritise oral-fluid multi-omics evidence using modality-aware
reproducibility scoring.

The scoring framework preserves the distinction between:
- longitudinal protein responses,
- longitudinal protein-module responses,
- technical small-RNA method sensitivity,
- exploratory case-control miRNA signals.

Scores are therefore comparable within evidence classes, but should not be
interpreted as direct biological equivalence across modalities.

Inputs
------
cross_cohort_evidence_registry.tsv

Outputs
-------
1. scored_evidence_registry.tsv
2. prioritised_candidate_panel.tsv
3. prioritisation_by_feature_type.tsv
4. prioritisation_report.txt
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = [
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
    "adjusted_p_value",
    "n_observations",
    "abundance_or_detection",
    "evidence_level",
    "is_prioritised",
    "source_file",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Score and prioritise evidence records using "
            "modality-aware reproducibility criteria."
        )
    )

    parser.add_argument(
        "--registry",
        type=Path,
        default=Path(
            "results/integration/cross_cohort_validation/"
            "cross_cohort_evidence_registry.tsv"
        ),
    )
    parser.add_argument(
        "--scored-output",
        type=Path,
        default=Path(
            "results/integration/cross_cohort_validation/"
            "scored_evidence_registry.tsv"
        ),
    )
    parser.add_argument(
        "--panel-output",
        type=Path,
        default=Path(
            "results/integration/cross_cohort_validation/"
            "prioritised_candidate_panel.tsv"
        ),
    )
    parser.add_argument(
        "--feature-summary-output",
        type=Path,
        default=Path(
            "results/integration/cross_cohort_validation/"
            "prioritisation_by_feature_type.tsv"
        ),
    )
    parser.add_argument(
        "--report-output",
        type=Path,
        default=Path(
            "results/integration/cross_cohort_validation/"
            "prioritisation_report.txt"
        ),
    )
    parser.add_argument(
        "--top-per-feature-type",
        type=int,
        default=20,
        help="Maximum number of candidates retained per feature type.",
    )

    return parser.parse_args()


def require_columns(frame: pd.DataFrame) -> None:
    missing = sorted(set(REQUIRED_COLUMNS) - set(frame.columns))
    if missing:
        raise ValueError(
            "Evidence registry is missing required columns: "
            + ", ".join(missing)
        )


def as_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def clip_score(values: pd.Series) -> pd.Series:
    return values.clip(lower=0.0, upper=25.0)


def minmax_score(
    values: pd.Series,
    *,
    absolute: bool = False,
    upper_quantile: float = 0.95,
) -> pd.Series:
    numeric = as_numeric(values)

    if absolute:
        numeric = numeric.abs()

    observed = numeric.dropna()

    if observed.empty:
        return pd.Series(0.0, index=values.index)

    lower = observed.min()
    upper = observed.quantile(upper_quantile)

    if not np.isfinite(lower) or not np.isfinite(upper):
        return pd.Series(0.0, index=values.index)

    if upper <= lower:
        return pd.Series(
            np.where(numeric.notna(), 25.0, 0.0),
            index=values.index,
            dtype=float,
        )

    scaled = 25.0 * (numeric - lower) / (upper - lower)
    return clip_score(scaled.fillna(0.0))


def fixed_consistency_score(values: pd.Series) -> pd.Series:
    consistency = as_numeric(values).fillna(0.0)
    return clip_score(consistency * 25.0)


def evidence_level_score(values: pd.Series) -> pd.Series:
    mapping = {
        "strong": 25.0,
        "moderate": 17.5,
        "supportive": 10.0,
    }

    return (
        values.astype(str)
        .str.lower()
        .map(mapping)
        .fillna(0.0)
        .astype(float)
    )


def pvalue_quality_score(values: pd.Series) -> pd.Series:
    p_values = as_numeric(values)

    score = pd.Series(0.0, index=values.index)

    score.loc[p_values <= 0.05] = 12.5
    score.loc[p_values <= 0.01] = 18.75
    score.loc[p_values <= 0.001] = 25.0

    return score


def score_pxd_proteins(frame: pd.DataFrame) -> pd.DataFrame:
    scored = frame.copy()

    scored["effect_score"] = minmax_score(
        scored["effect_size"],
        absolute=True,
    )

    scored["consistency_score"] = fixed_consistency_score(
        scored["direction_consistency"]
    )

    scored["robustness_score"] = fixed_consistency_score(
        scored["robustness_value"]
    )

    detection = as_numeric(
        scored["abundance_or_detection"]
    ).fillna(0.0)

    observations = as_numeric(
        scored["n_observations"]
    ).fillna(0.0)

    detection_component = clip_score(detection * 12.5)
    observation_component = clip_score(
        observations / 10.0 * 12.5
    )

    scored["support_quality_score"] = (
        detection_component + observation_component
    ).clip(upper=25.0)

    return scored


def score_pxd_modules(frame: pd.DataFrame) -> pd.DataFrame:
    scored = frame.copy()

    scored["effect_score"] = minmax_score(
        scored["effect_size"],
        absolute=True,
    )

    scored["consistency_score"] = fixed_consistency_score(
        scored["direction_consistency"]
    )

    scored["robustness_score"] = fixed_consistency_score(
        scored["robustness_value"]
    )

    p_score = pvalue_quality_score(
        scored["adjusted_p_value"]
    )

    observations = as_numeric(
        scored["n_observations"]
    ).fillna(0.0)

    observation_score = clip_score(
        observations / 10.0 * 12.5
    )

    scored["support_quality_score"] = (
        0.5 * p_score + observation_score
    ).clip(upper=25.0)

    return scored


def score_gse312585(frame: pd.DataFrame) -> pd.DataFrame:
    scored = frame.copy()

    scored["effect_score"] = minmax_score(
        scored["effect_size"],
        absolute=True,
    )

    scored["consistency_score"] = fixed_consistency_score(
        scored["direction_consistency"]
    )

    scored["robustness_score"] = fixed_consistency_score(
        scored["robustness_value"]
    )

    scored["support_quality_score"] = evidence_level_score(
        scored["evidence_level"]
    )

    return scored


def score_gse121910(frame: pd.DataFrame) -> pd.DataFrame:
    scored = frame.copy()

    scored["effect_score"] = minmax_score(
        scored["effect_size"],
        absolute=True,
    )

    # No inferential direction-consistency metric exists in the source table.
    scored["consistency_score"] = 0.0

    unique_sequences = as_numeric(
        scored["robustness_value"]
    ).fillna(0.0)

    scored["robustness_score"] = clip_score(
        unique_sequences / 4.0 * 25.0
    )

    abundance = as_numeric(
        scored["abundance_or_detection"]
    )

    abundance_log = np.log10(abundance.clip(lower=0.0) + 1.0)

    scored["support_quality_score"] = minmax_score(
        abundance_log,
        absolute=False,
    )

    return scored


def score_records(registry: pd.DataFrame) -> pd.DataFrame:
    partitions = []

    masks = {
        "pxd_protein": (
            registry["dataset"].eq("PXD043901")
            & registry["feature_type"].eq("protein")
        ),
        "pxd_module": (
            registry["dataset"].eq("PXD043901")
            & registry["feature_type"].eq("protein_module")
        ),
        "gse312585": registry["dataset"].eq("GSE312585"),
        "gse121910": registry["dataset"].eq("GSE121910"),
    }

    covered = pd.Series(False, index=registry.index)

    if masks["pxd_protein"].any():
        partitions.append(
            score_pxd_proteins(
                registry.loc[masks["pxd_protein"]]
            )
        )
        covered |= masks["pxd_protein"]

    if masks["pxd_module"].any():
        partitions.append(
            score_pxd_modules(
                registry.loc[masks["pxd_module"]]
            )
        )
        covered |= masks["pxd_module"]

    if masks["gse312585"].any():
        partitions.append(
            score_gse312585(
                registry.loc[masks["gse312585"]]
            )
        )
        covered |= masks["gse312585"]

    if masks["gse121910"].any():
        partitions.append(
            score_gse121910(
                registry.loc[masks["gse121910"]]
            )
        )
        covered |= masks["gse121910"]

    uncovered = registry.loc[~covered].copy()

    if not uncovered.empty:
        uncovered["effect_score"] = minmax_score(
            uncovered["effect_size"],
            absolute=True,
        )
        uncovered["consistency_score"] = fixed_consistency_score(
            uncovered["direction_consistency"]
        )
        uncovered["robustness_score"] = minmax_score(
            uncovered["robustness_value"]
        )
        uncovered["support_quality_score"] = evidence_level_score(
            uncovered["evidence_level"]
        )
        partitions.append(uncovered)

    scored = pd.concat(
        partitions,
        axis=0,
        ignore_index=True,
        sort=False,
    )

    component_columns = [
        "effect_score",
        "consistency_score",
        "robustness_score",
        "support_quality_score",
    ]

    for column in component_columns:
        scored[column] = as_numeric(
            scored[column]
        ).fillna(0.0).clip(0.0, 25.0)

    scored["reproducibility_score"] = (
        scored[component_columns].sum(axis=1)
    ).round(3)

    scored["score_band"] = pd.cut(
        scored["reproducibility_score"],
        bins=[-np.inf, 40.0, 60.0, 80.0, np.inf],
        labels=[
            "low",
            "moderate",
            "high",
            "very_high",
        ],
        right=False,
    ).astype(str)

    scored["within_feature_type_rank"] = (
        scored.groupby("feature_type")[
            "reproducibility_score"
        ]
        .rank(
            method="min",
            ascending=False,
        )
        .astype(int)
    )

    scored["within_dataset_rank"] = (
        scored.groupby("dataset")[
            "reproducibility_score"
        ]
        .rank(
            method="min",
            ascending=False,
        )
        .astype(int)
    )

    return scored.sort_values(
        [
            "feature_type",
            "reproducibility_score",
            "feature_name",
            "comparison",
        ],
        ascending=[True, False, True, True],
    ).reset_index(drop=True)


def build_candidate_panel(
    scored: pd.DataFrame,
    top_per_feature_type: int,
) -> pd.DataFrame:
    panel = scored.loc[
        (
            scored["within_feature_type_rank"]
            <= top_per_feature_type
        )
        & (
            scored["score_band"].isin(
                ["high", "very_high"]
            )
        )
    ].copy()

    selected_columns = [
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
    ]

    return panel[selected_columns].sort_values(
        [
            "feature_type",
            "within_feature_type_rank",
            "dataset",
        ]
    )


def build_feature_type_summary(
    scored: pd.DataFrame,
) -> pd.DataFrame:
    summary = (
        scored.groupby(
            [
                "feature_type",
                "modality",
                "dataset",
            ],
            dropna=False,
        )
        .agg(
            n_records=("evidence_id", "size"),
            n_unique_features=("feature_id", "nunique"),
            mean_reproducibility_score=(
                "reproducibility_score",
                "mean",
            ),
            median_reproducibility_score=(
                "reproducibility_score",
                "median",
            ),
            maximum_reproducibility_score=(
                "reproducibility_score",
                "max",
            ),
            n_high_or_very_high=(
                "score_band",
                lambda x: x.isin(
                    ["high", "very_high"]
                ).sum(),
            ),
        )
        .reset_index()
    )

    numeric_columns = [
        "mean_reproducibility_score",
        "median_reproducibility_score",
        "maximum_reproducibility_score",
    ]

    summary[numeric_columns] = summary[
        numeric_columns
    ].round(3)

    return summary.sort_values(
        [
            "feature_type",
            "dataset",
        ]
    )


def write_report(
    scored: pd.DataFrame,
    panel: pd.DataFrame,
    output_path: Path,
) -> None:
    lines = [
        "Modality-aware candidate prioritisation",
        "=======================================",
        "",
        f"Scored evidence records: {len(scored)}",
        f"Candidate-panel records: {len(panel)}",
        "",
        "Score bands:",
    ]

    band_counts = scored["score_band"].value_counts()

    for band in [
        "very_high",
        "high",
        "moderate",
        "low",
    ]:
        lines.append(
            f"  {band}: {int(band_counts.get(band, 0))}"
        )

    lines.extend(
        [
            "",
            "Candidate-panel records by feature type:",
        ]
    )

    panel_counts = panel["feature_type"].value_counts()

    if panel_counts.empty:
        lines.append("  none")
    else:
        for feature_type, count in panel_counts.items():
            lines.append(f"  {feature_type}: {count}")

    lines.extend(
        [
            "",
            "Scoring interpretation:",
            (
                "  - Scores quantify reproducibility-oriented evidence "
                "within each evidence class."
            ),
            (
                "  - Protein, protein-module, mature-miRNA, and technical "
                "method-sensitivity scores are not direct biological "
                "equivalents."
            ),
            (
                "  - GSE312585 contributes technical robustness evidence, "
                "not disease-direction evidence."
            ),
            (
                "  - GSE121910 lacks inferential p-values and therefore "
                "receives no direction-consistency component."
            ),
            (
                "  - Cross-modality biological integration should occur "
                "later through targets, pathways, or functional modules."
            ),
            "",
        ]
    )

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()

    if not args.registry.is_file():
        raise FileNotFoundError(
            f"Evidence registry not found: {args.registry}"
        )

    registry = pd.read_csv(
        args.registry,
        sep="\t",
    )

    require_columns(registry)

    scored = score_records(registry)

    panel = build_candidate_panel(
        scored,
        top_per_feature_type=args.top_per_feature_type,
    )

    feature_summary = build_feature_type_summary(scored)

    output_paths = [
        args.scored_output,
        args.panel_output,
        args.feature_summary_output,
        args.report_output,
    ]

    for path in output_paths:
        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

    scored.to_csv(
        args.scored_output,
        sep="\t",
        index=False,
    )

    panel.to_csv(
        args.panel_output,
        sep="\t",
        index=False,
    )

    feature_summary.to_csv(
        args.feature_summary_output,
        sep="\t",
        index=False,
    )

    write_report(
        scored,
        panel,
        args.report_output,
    )

    print(f"Wrote scored registry: {args.scored_output}")
    print(f"Wrote candidate panel: {args.panel_output}")
    print(
        "Wrote feature-type summary: "
        f"{args.feature_summary_output}"
    )
    print(f"Wrote report: {args.report_output}")
    print(f"Scored records: {len(scored)}")
    print(f"Panel records:  {len(panel)}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)

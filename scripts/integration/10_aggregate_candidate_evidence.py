#!/usr/bin/env python3
"""
Aggregate record-level reproducibility scores into unique feature-level
candidate evidence.

This step prevents repeated comparisons for the same feature from occupying
multiple positions in the final candidate panel.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


GROUP_COLUMNS = [
    "feature_id",
    "feature_name",
    "feature_type",
    "modality",
    "dataset",
    "study_role",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate scored evidence into unique candidate features."
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=Path(
            "results/integration/cross_cohort_validation/"
            "scored_evidence_registry.tsv"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "results/integration/cross_cohort_validation/"
            "aggregated_candidate_evidence.tsv"
        ),
    )
    parser.add_argument(
        "--panel-output",
        type=Path,
        default=Path(
            "results/integration/cross_cohort_validation/"
            "final_candidate_panel.tsv"
        ),
    )
    parser.add_argument(
        "--report-output",
        type=Path,
        default=Path(
            "results/integration/cross_cohort_validation/"
            "aggregated_candidate_report.txt"
        ),
    )

    return parser.parse_args()


def require_columns(frame: pd.DataFrame) -> None:
    required = set(
        GROUP_COLUMNS
        + [
            "comparison",
            "effect_size",
            "effect_direction",
            "direction_consistency",
            "robustness_value",
            "reproducibility_score",
            "score_band",
            "evidence_level",
            "is_prioritised",
        ]
    )

    missing = sorted(required - set(frame.columns))

    if missing:
        raise ValueError(
            "Scored registry is missing required columns: "
            + ", ".join(missing)
        )


def join_unique(values: pd.Series) -> str:
    unique = sorted(
        {
            str(value)
            for value in values.dropna()
            if str(value).strip()
        }
    )
    return ";".join(unique)


def direction_agreement(values: pd.Series) -> float:
    directions = [
        str(value)
        for value in values.dropna()
        if str(value) in {"increase", "decrease", "no_change"}
    ]

    if not directions:
        return np.nan

    counts = pd.Series(directions).value_counts()
    return float(counts.max() / counts.sum())


def strongest_level(values: pd.Series) -> str:
    order = {
        "supportive": 1,
        "moderate": 2,
        "strong": 3,
    }

    observed = [
        str(value).lower()
        for value in values.dropna()
        if str(value).lower() in order
    ]

    if not observed:
        return "unknown"

    return max(observed, key=lambda value: order[value])


def aggregate_features(frame: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        frame.groupby(GROUP_COLUMNS, dropna=False)
        .agg(
            n_evidence_records=("evidence_id", "size"),
            n_comparisons=("comparison", "nunique"),
            comparisons=("comparison", join_unique),
            effect_directions=("effect_direction", join_unique),
            direction_agreement=(
                "effect_direction",
                direction_agreement,
            ),
            maximum_absolute_effect=(
                "effect_size",
                lambda x: pd.to_numeric(
                    x,
                    errors="coerce",
                ).abs().max(),
            ),
            median_absolute_effect=(
                "effect_size",
                lambda x: pd.to_numeric(
                    x,
                    errors="coerce",
                ).abs().median(),
            ),
            mean_direction_consistency=(
                "direction_consistency",
                "mean",
            ),
            mean_robustness_value=(
                "robustness_value",
                "mean",
            ),
            maximum_record_score=(
                "reproducibility_score",
                "max",
            ),
            mean_record_score=(
                "reproducibility_score",
                "mean",
            ),
            minimum_record_score=(
                "reproducibility_score",
                "min",
            ),
            strongest_evidence_level=(
                "evidence_level",
                strongest_level,
            ),
            any_prioritised=(
                "is_prioritised",
                "max",
            ),
        )
        .reset_index()
    )

    return grouped


def calculate_feature_scores(frame: pd.DataFrame) -> pd.DataFrame:
    scored = frame.copy()

    # Reward evidence that remains strong across all available records,
    # rather than allowing one extreme comparison to dominate.
    scored["base_feature_score"] = (
        0.50 * scored["mean_record_score"]
        + 0.30 * scored["minimum_record_score"]
        + 0.20 * scored["maximum_record_score"]
    )

    longitudinal_mask = (
        scored["dataset"].eq("PXD043901")
        & scored["n_comparisons"].ge(2)
    )

    concordant_mask = (
        longitudinal_mask
        & scored["direction_agreement"].ge(1.0)
    )

    discordant_mask = (
        longitudinal_mask
        & scored["direction_agreement"].lt(1.0)
    )

    scored["replication_adjustment"] = 0.0
    scored.loc[concordant_mask, "replication_adjustment"] = 5.0
    scored.loc[discordant_mask, "replication_adjustment"] = -5.0

    scored["aggregated_reproducibility_score"] = (
        scored["base_feature_score"]
        + scored["replication_adjustment"]
    ).clip(0.0, 100.0)

    scored["within_feature_type_percentile"] = (
        scored.groupby("feature_type")[
            "aggregated_reproducibility_score"
        ]
        .rank(
            method="average",
            ascending=True,
            pct=True,
        )
        * 100.0
    )

    scored["within_feature_type_rank"] = (
        scored.groupby("feature_type")[
            "aggregated_reproducibility_score"
        ]
        .rank(
            method="min",
            ascending=False,
        )
        .astype(int)
    )

    high_confidence = (
        scored["within_feature_type_percentile"].ge(90.0)
        & scored["aggregated_reproducibility_score"].ge(70.0)
    )

    prioritised = (
        scored["within_feature_type_percentile"].ge(75.0)
        & scored["aggregated_reproducibility_score"].ge(60.0)
    )

    scored["candidate_tier"] = np.select(
        [
            high_confidence,
            prioritised,
        ],
        [
            "high_confidence",
            "prioritised",
        ],
        default="exploratory",
    )

    numeric_columns = [
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
    ]

    scored[numeric_columns] = scored[numeric_columns].round(3)

    return scored.sort_values(
        [
            "feature_type",
            "aggregated_reproducibility_score",
            "feature_name",
        ],
        ascending=[True, False, True],
    ).reset_index(drop=True)


def build_panel(frame: pd.DataFrame) -> pd.DataFrame:
    panel = frame.loc[
        frame["candidate_tier"].isin(
            ["high_confidence", "prioritised"]
        )
    ].copy()

    return panel.sort_values(
        [
            "feature_type",
            "within_feature_type_rank",
            "feature_name",
        ]
    )


def write_report(
    aggregated: pd.DataFrame,
    panel: pd.DataFrame,
    output: Path,
) -> None:
    lines = [
        "Aggregated candidate evidence",
        "=============================",
        "",
        f"Unique candidate features evaluated: {len(aggregated)}",
        f"Features retained in final panel: {len(panel)}",
        "",
        "Candidate tiers:",
    ]

    counts = aggregated["candidate_tier"].value_counts()

    for tier in [
        "high_confidence",
        "prioritised",
        "exploratory",
    ]:
        lines.append(f"  {tier}: {int(counts.get(tier, 0))}")

    lines.extend(
        [
            "",
            "Final panel by feature type:",
        ]
    )

    panel_counts = panel["feature_type"].value_counts()

    for feature_type, count in panel_counts.items():
        lines.append(f"  {feature_type}: {count}")

    lines.extend(
        [
            "",
            "Interpretation:",
            (
                "  - Each row represents one unique feature within one "
                "dataset and evidence class."
            ),
            (
                "  - Multiple longitudinal comparisons are aggregated rather "
                "than treated as independent candidate entries."
            ),
            (
                "  - Concordant T1 and T6 protein responses receive a small "
                "replication bonus."
            ),
            (
                "  - Candidate tiers are based on both absolute score and "
                "within-feature-type percentile."
            ),
            (
                "  - Rankings remain modality-aware and should not be read "
                "as direct equivalence between proteins and miRNAs."
            ),
            "",
        ]
    )

    output.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()

    if not args.input.is_file():
        raise FileNotFoundError(
            f"Scored evidence registry not found: {args.input}"
        )

    frame = pd.read_csv(args.input, sep="\t")
    require_columns(frame)

    aggregated = aggregate_features(frame)
    aggregated = calculate_feature_scores(aggregated)
    panel = build_panel(aggregated)

    for path in [
        args.output,
        args.panel_output,
        args.report_output,
    ]:
        path.parent.mkdir(parents=True, exist_ok=True)

    aggregated.to_csv(args.output, sep="\t", index=False)
    panel.to_csv(args.panel_output, sep="\t", index=False)

    write_report(
        aggregated,
        panel,
        args.report_output,
    )

    print(f"Wrote aggregated evidence: {args.output}")
    print(f"Wrote final candidate panel: {args.panel_output}")
    print(f"Wrote report: {args.report_output}")
    print(f"Unique features evaluated: {len(aggregated)}")
    print(f"Final panel features: {len(panel)}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)

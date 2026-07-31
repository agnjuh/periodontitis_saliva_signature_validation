#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


INPUT_FILE = Path("config/dataset_registry.csv")
OUTPUT_FILE = Path("results/dataset_screening/dataset_prioritisation.csv")
SUMMARY_FILE = Path("results/dataset_screening/dataset_screening_summary.txt")


POSITIVE_VALUES = {"yes", "true", "available", "complete"}
PARTIAL_VALUES = {"partial", "limited"}
UNKNOWN_VALUES = {"", "unknown", "na", "n/a", "not reported"}


def normalise(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().lower()


def availability_score(value: object, yes_score: float = 1.0) -> float:
    value = normalise(value)

    if value in POSITIVE_VALUES:
        return yes_score
    if value in PARTIAL_VALUES:
        return yes_score * 0.5
    if value in UNKNOWN_VALUES:
        return 0.0

    return 0.0


def study_design_score(value: object) -> float:
    value = normalise(value)

    if "longitudinal" in value:
        return 3.0
    if "case-control" in value or "case control" in value:
        return 1.5
    if "cross-sectional" in value:
        return 1.0
    return 0.0


def disease_score(value: object) -> float:
    value = normalise(value)

    if "periodontitis" in value:
        return 2.0
    if "gingivitis" in value:
        return 1.0
    if "peri-implantitis" in value:
        return 0.75
    return 0.0


def sample_score(value: object) -> float:
    value = normalise(value)

    if "saliva" in value:
        return 2.0
    if "oral rinse" in value:
        return 1.5
    if "gingival crevicular fluid" in value or value == "gcf":
        return 1.0
    return 0.0


def omics_novelty_score(value: object) -> float:
    value = normalise(value)

    if "mirna" in value or "small rna" in value:
        return 2.0
    if "lncrna" in value or "transcriptomics" in value:
        return 1.5
    if "extracellular vesicle" in value or value == "ev-rna":
        return 2.0
    if "microbiome" in value or "metagenomics" in value:
        return 1.5
    if "metabolomics" in value:
        return 1.5
    if "proteomics" in value:
        return 1.0
    return 0.0


def calculate_score(row: pd.Series) -> float:
    score = 0.0

    score += study_design_score(row.get("study_design"))
    score += disease_score(row.get("disease"))
    score += sample_score(row.get("sample_type"))
    score += omics_novelty_score(row.get("omics_layer"))

    score += availability_score(
        row.get("raw_data_available"),
        yes_score=1.5,
    )
    score += availability_score(
        row.get("processed_data_available"),
        yes_score=1.0,
    )
    score += availability_score(
        row.get("clinical_metadata_available"),
        yes_score=2.0,
    )
    score += availability_score(
        row.get("patient_matching_possible"),
        yes_score=1.5,
    )
    score += availability_score(
        row.get("smoking_recorded"),
        yes_score=0.5,
    )
    score += availability_score(
        row.get("diabetes_recorded"),
        yes_score=0.5,
    )
    score += availability_score(
        row.get("antibiotic_use_recorded"),
        yes_score=0.5,
    )
    score += availability_score(
        row.get("age_sex_available"),
        yes_score=0.5,
    )

    treatment = normalise(row.get("treatment"))
    timepoints = normalise(row.get("timepoints"))

    if treatment not in UNKNOWN_VALUES:
        score += 1.0

    if timepoints not in UNKNOWN_VALUES:
        score += 1.0

    return round(score, 2)


def assign_priority(score: float) -> str:
    if score >= 15:
        return "very_high"
    if score >= 11:
        return "high"
    if score >= 7:
        return "moderate"
    return "low"


def main() -> None:
    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Input registry not found: {INPUT_FILE}"
        )

    df = pd.read_csv(INPUT_FILE)

    df["priority_score"] = df.apply(
        calculate_score,
        axis=1,
    )

    df["priority_class"] = df["priority_score"].apply(
        assign_priority
    )

    df = df.sort_values(
        by=["priority_score", "dataset_id"],
        ascending=[False, True],
    )

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    df.to_csv(
        OUTPUT_FILE,
        index=False,
    )

    summary_lines = [
        "Candidate dataset screening",
        "===========================",
        "",
        f"Datasets assessed: {len(df)}",
        "",
        "Priority distribution:",
        df["priority_class"].value_counts().to_string(),
        "",
        "Ranked candidates:",
        df[
            [
                "dataset_id",
                "omics_layer",
                "disease",
                "priority_score",
                "priority_class",
                "status",
            ]
        ].to_string(index=False),
        "",
        "Important:",
        (
            "Scores support transparent prioritisation but do not replace "
            "manual review of clinical design, treatment, sample collection, "
            "preanalytics, and publication methods."
        ),
    ]

    SUMMARY_FILE.write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    print(
        df[
            [
                "dataset_id",
                "omics_layer",
                "priority_score",
                "priority_class",
                "status",
            ]
        ].to_string(index=False)
    )

    print(f"\nSaved: {OUTPUT_FILE}")
    print(f"Saved: {SUMMARY_FILE}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Import and standardise experimentally supported miRNA-target interactions.

The script:

1. reads a miRTarBase export in XLSX, CSV, TSV, or TXT format;
2. detects common miRTarBase column-name variants;
3. standardises miRNA identifiers, target-gene symbols, species,
   experimental methods, support categories, and publication identifiers;
4. retains human miRNA-human target interactions;
5. classifies interaction evidence conservatively;
6. identifies interactions involving regulatory-source miRNAs from the
   biological-network node registry;
7. writes auditable reference tables and a summary report.

No interaction is inferred by this script. Every output edge must originate
from the supplied miRTarBase export.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd


COLUMN_ALIASES = {
    "mirna": [
        "mirna",
        "mirna name",
        "mirna_name",
        "mirna id",
        "mirna_id",
        "mirna mature id",
        "mirtarbase id",
        "mirna mature",
    ],
    "target_gene": [
        "target gene",
        "target_gene",
        "target gene symbol",
        "target_gene_symbol",
        "target",
        "target symbol",
        "gene symbol",
        "gene_symbol",
    ],
    "mirna_species": [
        "mirna species",
        "mirna_species",
        "species mirna",
        "species_mirna",
        "mirna organism",
    ],
    "target_species": [
        "target gene species",
        "target_gene_species",
        "target species",
        "target_species",
        "species target gene",
        "target organism",
    ],
    "experiments": [
        "experiments",
        "experiment",
        "experimental method",
        "experimental_method",
        "supporting experiments",
        "supporting_experiments",
        "experimental validation",
        "validation method",
    ],
    "support_type": [
        "support type",
        "support_type",
        "support category",
        "support_category",
        "evidence type",
        "evidence_type",
        "validation support",
    ],
    "publication_id": [
        "references pmid",
        "references_pmid",
        "reference pmid",
        "reference_pmid",
        "pmid",
        "pubmed id",
        "pubmed_id",
        "publication id",
        "publication_id",
        "references",
    ],
    "mti_id": [
        "mti id",
        "mti_id",
        "mirtarbase mti id",
        "mirtarbase_mti_id",
        "record id",
        "record_id",
    ],
}


STRONG_EXPERIMENT_PATTERNS = [
    r"reporter assay",
    r"luciferase",
    r"western blot",
    r"\bqpcr\b",
    r"qRT[- ]?PCR",
    r"quantitative.*PCR",
    r"immunoblot",
]

SUPPORTIVE_EXPERIMENT_PATTERNS = [
    r"\bCLIP\b",
    r"HITS[- ]?CLIP",
    r"PAR[- ]?CLIP",
    r"CLASH",
    r"microarray",
    r"\bNGS\b",
    r"RNA[- ]?seq",
    r"sequencing",
    r"proteomics",
    r"pSILAC",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Import and standardise experimentally supported "
            "miRNA-target interactions."
        )
    )

    parser.add_argument(
        "--mirtarbase",
        type=Path,
        default=Path(
            "data/reference/mirtarbase/"
            "MiRTarBase_MTI.xlsx"
        ),
        help=(
            "miRTarBase export in XLSX, XLS, CSV, TSV, or TXT format."
        ),
    )

    parser.add_argument(
        "--nodes",
        type=Path,
        default=Path(
            "results/integration/biological_network/"
            "network_node_registry.tsv"
        ),
        help="Biological-network node registry from step 11.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "results/integration/biological_network/"
            "mirtarbase"
        ),
        help="Directory for standardised miRTarBase outputs.",
    )

    parser.add_argument(
        "--sheet",
        default=None,
        help=(
            "Optional Excel sheet name. By default, the first "
            "non-empty sheet is used."
        ),
    )

    return parser.parse_args()


def normalise_column_name(value: object) -> str:
    text = str(value).strip().lower()
    text = re.sub(r"[_\-]+", " ", text)
    text = re.sub(r"\s+", " ", text)
    return text


def normalise_text(value: object) -> str:
    if pd.isna(value):
        return ""

    text = str(value).strip()

    if text.lower() in {"nan", "none", "null", "na", "n/a"}:
        return ""

    return text


def read_excel_table(
    path: Path,
    sheet: str | None,
) -> tuple[pd.DataFrame, str]:
    workbook = pd.ExcelFile(path)

    if sheet is not None:
        if sheet not in workbook.sheet_names:
            raise ValueError(
                f"Excel sheet not found: {sheet}. "
                f"Available sheets: {', '.join(workbook.sheet_names)}"
            )

        frame = pd.read_excel(
            path,
            sheet_name=sheet,
        )

        return frame, sheet

    for sheet_name in workbook.sheet_names:
        frame = pd.read_excel(
            path,
            sheet_name=sheet_name,
        )

        if not frame.empty and len(frame.columns) > 1:
            return frame, sheet_name

    raise ValueError(
        "No non-empty worksheet was found in the Excel file."
    )


def read_delimited_table(
    path: Path,
) -> pd.DataFrame:
    suffix = path.suffix.lower()

    if suffix == ".csv":
        return pd.read_csv(
            path,
            low_memory=False,
        )

    if suffix in {".tsv", ".txt"}:
        try:
            return pd.read_csv(
                path,
                sep="\t",
                low_memory=False,
            )
        except pd.errors.ParserError:
            return pd.read_csv(
                path,
                sep=None,
                engine="python",
                low_memory=False,
            )

    raise ValueError(
        f"Unsupported delimited-file extension: {suffix}"
    )


def read_mirtarbase(
    path: Path,
    sheet: str | None,
) -> tuple[pd.DataFrame, str]:
    suffix = path.suffix.lower()

    if suffix in {".xlsx", ".xls"}:
        frame, selected_sheet = read_excel_table(
            path,
            sheet,
        )

        return frame, f"Excel sheet: {selected_sheet}"

    if suffix in {".csv", ".tsv", ".txt"}:
        frame = read_delimited_table(path)
        return frame, "Delimited text file"

    raise ValueError(
        "Unsupported miRTarBase file type. Supported extensions: "
        ".xlsx, .xls, .csv, .tsv, .txt"
    )


def detect_column(
    frame: pd.DataFrame,
    logical_name: str,
    required: bool,
) -> str | None:
    aliases = {
        normalise_column_name(alias)
        for alias in COLUMN_ALIASES[logical_name]
    }

    normalised_columns = {
        normalise_column_name(column): column
        for column in frame.columns
    }

    for alias in aliases:
        if alias in normalised_columns:
            return normalised_columns[alias]

    for normalised, original in normalised_columns.items():
        for alias in aliases:
            if alias in normalised or normalised in alias:
                return original

    if required:
        raise ValueError(
            f"Could not identify required '{logical_name}' column. "
            f"Available columns: {', '.join(map(str, frame.columns))}"
        )

    return None


def detect_columns(
    frame: pd.DataFrame,
) -> dict[str, str | None]:
    required = {
        "mirna",
        "target_gene",
    }

    detected: dict[str, str | None] = {}

    for logical_name in COLUMN_ALIASES:
        detected[logical_name] = detect_column(
            frame=frame,
            logical_name=logical_name,
            required=logical_name in required,
        )

    return detected


def series_or_empty(
    frame: pd.DataFrame,
    column: str | None,
) -> pd.Series:
    if column is None:
        return pd.Series(
            [""] * len(frame),
            index=frame.index,
            dtype="object",
        )

    return frame[column].map(normalise_text)


def canonicalise_mirna(value: object) -> str:
    text = normalise_text(value)

    if not text:
        return ""

    text = text.replace("_", "-")
    text = re.sub(r"\s+", "", text)

    text = re.sub(
        r"^(HSA|hsa)-",
        "",
        text,
    )

    text = re.sub(
        r"^(MIR|Mir|mir)(?=\d)",
        "miR-",
        text,
    )

    text = re.sub(
        r"^(miRNA|mirna)-?",
        "miR-",
        text,
    )

    text = re.sub(
        r"^miR--+",
        "miR-",
        text,
    )

    match = re.match(
        r"^miR-?(\d+[a-zA-Z]?)(?:-)?([35]p)?$",
        text,
        flags=re.IGNORECASE,
    )

    if match:
        locus = match.group(1).lower()
        arm = match.group(2)

        if arm:
            return f"hsa-miR-{locus}-{arm.lower()}"

        return f"hsa-miR-{locus}"

    if text.lower().startswith("let-"):
        return "hsa-" + text.lower()

    if text.lower().startswith("mir-"):
        return "hsa-" + text.replace(
            "mir-",
            "miR-",
            1,
        )

    return text


def mirna_match_key(value: object) -> str:
    canonical = canonicalise_mirna(value)

    canonical = canonical.lower()
    canonical = re.sub(
        r"^hsa-",
        "",
        canonical,
    )
    canonical = canonical.replace("mirna", "mir")
    canonical = canonical.replace("_", "-")
    canonical = re.sub(r"-+", "-", canonical)

    return canonical


def canonicalise_gene_symbol(value: object) -> str:
    text = normalise_text(value)

    if not text:
        return ""

    text = text.split(";")[0].strip()
    text = text.split("|")[0].strip()
    text = re.sub(r"\s+", "", text)

    return text.upper()


def canonicalise_species(value: object) -> str:
    text = normalise_text(value).lower()

    if not text:
        return ""

    if (
        "homo sapiens" in text
        or text == "human"
        or text == "hsa"
        or text.startswith("hsa-")
        or "9606" in text
    ):
        return "Homo sapiens"

    if (
        "mus musculus" in text
        or text == "mouse"
        or text == "mmu"
        or "10090" in text
    ):
        return "Mus musculus"

    if (
        "rattus norvegicus" in text
        or text == "rat"
        or text == "rno"
        or "10116" in text
    ):
        return "Rattus norvegicus"

    return normalise_text(value)


def infer_mirna_species(
    mirna: object,
    species: object,
) -> str:
    standardised = canonicalise_species(species)

    if standardised:
        return standardised

    raw_mirna = normalise_text(mirna).lower()

    if raw_mirna.startswith("hsa-"):
        return "Homo sapiens"

    if raw_mirna.startswith("mmu-"):
        return "Mus musculus"

    if raw_mirna.startswith("rno-"):
        return "Rattus norvegicus"

    return ""


def normalise_publication_id(value: object) -> str:
    text = normalise_text(value)

    if not text:
        return ""

    pmids = re.findall(
        r"\b\d{5,9}\b",
        text,
    )

    if pmids:
        return ";".join(dict.fromkeys(pmids))

    return text


def split_methods(value: object) -> list[str]:
    text = normalise_text(value)

    if not text:
        return []

    parts = re.split(
        r"[;,|/]+",
        text,
    )

    methods = [
        re.sub(r"\s+", " ", part).strip()
        for part in parts
        if part.strip()
    ]

    return list(dict.fromkeys(methods))


def contains_pattern(
    text: str,
    patterns: Iterable[str],
) -> bool:
    return any(
        re.search(
            pattern,
            text,
            flags=re.IGNORECASE,
        )
        for pattern in patterns
    )


def classify_evidence_strength(
    experiments: object,
    support_type: object,
) -> str:
    experiment_text = normalise_text(experiments)
    support_text = normalise_text(support_type)

    combined = f"{experiment_text}; {support_text}"

    if re.search(
        r"\bstrong\b",
        support_text,
        flags=re.IGNORECASE,
    ):
        return "strong"

    if re.search(
        r"less strong|weak|limited",
        support_text,
        flags=re.IGNORECASE,
    ):
        return "supportive"

    if contains_pattern(
        combined,
        STRONG_EXPERIMENT_PATTERNS,
    ):
        return "strong"

    if contains_pattern(
        combined,
        SUPPORTIVE_EXPERIMENT_PATTERNS,
    ):
        return "supportive"

    if experiment_text or support_text:
        return "experimentally_supported_unspecified"

    return "unclassified"


def classify_relationship_type(
    experiments: object,
) -> str:
    text = normalise_text(experiments)

    if not text:
        return "experimentally_supported_interaction"

    if re.search(
        r"reporter|luciferase",
        text,
        flags=re.IGNORECASE,
    ):
        return "direct_regulatory_validation"

    if re.search(
        r"western blot|immunoblot|proteomics|pSILAC",
        text,
        flags=re.IGNORECASE,
    ):
        return "protein_level_validation"

    if re.search(
        r"qPCR|qRT-PCR|quantitative.*PCR|microarray|RNA-seq|NGS",
        text,
        flags=re.IGNORECASE,
    ):
        return "transcript_level_validation"

    if re.search(
        r"CLIP|CLASH",
        text,
        flags=re.IGNORECASE,
    ):
        return "physical_binding_support"

    return "experimentally_supported_interaction"


def build_standardised_table(
    raw: pd.DataFrame,
    detected: dict[str, str | None],
    source_file: Path,
) -> pd.DataFrame:
    standardised = pd.DataFrame(
        index=raw.index,
    )

    standardised["mti_id"] = series_or_empty(
        raw,
        detected["mti_id"],
    )

    standardised["mirna_raw"] = series_or_empty(
        raw,
        detected["mirna"],
    )

    standardised["mirna"] = standardised[
        "mirna_raw"
    ].map(canonicalise_mirna)

    standardised["mirna_match_key"] = standardised[
        "mirna_raw"
    ].map(mirna_match_key)

    standardised["target_gene_raw"] = series_or_empty(
        raw,
        detected["target_gene"],
    )

    standardised["target_gene"] = standardised[
        "target_gene_raw"
    ].map(canonicalise_gene_symbol)

    mirna_species_raw = series_or_empty(
        raw,
        detected["mirna_species"],
    )

    target_species_raw = series_or_empty(
        raw,
        detected["target_species"],
    )

    standardised["mirna_species"] = [
        infer_mirna_species(mirna, species)
        for mirna, species in zip(
            standardised["mirna_raw"],
            mirna_species_raw,
            strict=False,
        )
    ]

    standardised["target_species"] = (
        target_species_raw.map(canonicalise_species)
    )

    standardised["experiments"] = series_or_empty(
        raw,
        detected["experiments"],
    )

    standardised["experimental_methods"] = (
        standardised["experiments"]
        .map(split_methods)
        .map(lambda methods: ";".join(methods))
    )

    standardised["support_type_raw"] = series_or_empty(
        raw,
        detected["support_type"],
    )

    standardised["evidence_strength"] = [
        classify_evidence_strength(
            experiment,
            support,
        )
        for experiment, support in zip(
            standardised["experiments"],
            standardised["support_type_raw"],
            strict=False,
        )
    ]

    standardised["relationship_type"] = standardised[
        "experiments"
    ].map(classify_relationship_type)

    standardised["publication_id"] = series_or_empty(
        raw,
        detected["publication_id"],
    ).map(normalise_publication_id)

    standardised["evidence_database"] = "miRTarBase"
    standardised["database_version"] = (
        source_file.stem
    )
    standardised["source_file"] = str(source_file)

    standardised = standardised.loc[
        standardised["mirna"].ne("")
        & standardised["target_gene"].ne("")
    ].copy()

    standardised["mirna_species"] = standardised[
        "mirna_species"
    ].replace(
        "",
        pd.NA,
    )

    standardised["target_species"] = standardised[
        "target_species"
    ].replace(
        "",
        pd.NA,
    )

    return standardised.reset_index(drop=True)


def select_human_interactions(
    interactions: pd.DataFrame,
) -> pd.DataFrame:
    mirna_human = (
        interactions["mirna_species"]
        .fillna("")
        .eq("Homo sapiens")
        | interactions["mirna"]
        .str.startswith("hsa-", na=False)
    )

    target_species_known = interactions[
        "target_species"
    ].notna()

    target_human = (
        interactions["target_species"]
        .fillna("")
        .eq("Homo sapiens")
    )

    keep = mirna_human & (
        target_human
        | ~target_species_known
    )

    human = interactions.loc[keep].copy()

    human["mirna_species"] = (
        human["mirna_species"]
        .fillna("Homo sapiens")
    )

    human["target_species"] = (
        human["target_species"]
        .fillna("Homo sapiens_assumed")
    )

    return human.reset_index(drop=True)


def aggregate_interactions(
    interactions: pd.DataFrame,
) -> pd.DataFrame:
    if interactions.empty:
        return interactions.copy()

    grouping_columns = [
        "mirna",
        "mirna_match_key",
        "target_gene",
    ]

    aggregated = (
        interactions.groupby(
            grouping_columns,
            dropna=False,
        )
        .agg(
            mirna_species=(
                "mirna_species",
                lambda values: ";".join(
                    sorted(
                        {
                            normalise_text(value)
                            for value in values
                            if normalise_text(value)
                        }
                    )
                ),
            ),
            target_species=(
                "target_species",
                lambda values: ";".join(
                    sorted(
                        {
                            normalise_text(value)
                            for value in values
                            if normalise_text(value)
                        }
                    )
                ),
            ),
            experimental_methods=(
                "experimental_methods",
                lambda values: ";".join(
                    sorted(
                        {
                            method
                            for value in values
                            for method in split_methods(value)
                        }
                    )
                ),
            ),
            support_types=(
                "support_type_raw",
                lambda values: ";".join(
                    sorted(
                        {
                            normalise_text(value)
                            for value in values
                            if normalise_text(value)
                        }
                    )
                ),
            ),
            evidence_strengths=(
                "evidence_strength",
                lambda values: ";".join(
                    sorted(
                        {
                            normalise_text(value)
                            for value in values
                            if normalise_text(value)
                        }
                    )
                ),
            ),
            relationship_types=(
                "relationship_type",
                lambda values: ";".join(
                    sorted(
                        {
                            normalise_text(value)
                            for value in values
                            if normalise_text(value)
                        }
                    )
                ),
            ),
            publication_ids=(
                "publication_id",
                lambda values: ";".join(
                    sorted(
                        {
                            pmid
                            for value in values
                            for pmid in normalise_text(value).split(";")
                            if pmid
                        }
                    )
                ),
            ),
            n_source_records=(
                "target_gene",
                "size",
            ),
            evidence_database=(
                "evidence_database",
                "first",
            ),
            database_version=(
                "database_version",
                "first",
            ),
            source_file=(
                "source_file",
                "first",
            ),
        )
        .reset_index()
    )

    strength_rank = {
        "strong": 3,
        "experimentally_supported_unspecified": 2,
        "supportive": 1,
        "unclassified": 0,
    }

    def strongest_evidence(value: object) -> str:
        levels = [
            item
            for item in normalise_text(value).split(";")
            if item
        ]

        if not levels:
            return "unclassified"

        return max(
            levels,
            key=lambda level: strength_rank.get(
                level,
                -1,
            ),
        )

    aggregated["maximum_evidence_strength"] = (
        aggregated["evidence_strengths"]
        .map(strongest_evidence)
    )

    aggregated["n_publications"] = (
        aggregated["publication_ids"]
        .map(
            lambda value: len(
                {
                    item
                    for item in normalise_text(value).split(";")
                    if item
                }
            )
        )
    )

    aggregated["n_experimental_methods"] = (
        aggregated["experimental_methods"]
        .map(
            lambda value: len(
                {
                    item
                    for item in normalise_text(value).split(";")
                    if item
                }
            )
        )
    )

    aggregated["interaction_record_id"] = [
        f"mirtarbase:{mirna}:{target}"
        for mirna, target in zip(
            aggregated["mirna"],
            aggregated["target_gene"],
            strict=False,
        )
    ]

    strength_order = {
        "strong": 0,
        "experimentally_supported_unspecified": 1,
        "supportive": 2,
        "unclassified": 3,
    }

    aggregated["_strength_order"] = aggregated[
        "maximum_evidence_strength"
    ].map(strength_order).fillna(9)

    aggregated = (
        aggregated.sort_values(
            [
                "_strength_order",
                "n_publications",
                "n_experimental_methods",
                "mirna",
                "target_gene",
            ],
            ascending=[
                True,
                False,
                False,
                True,
                True,
            ],
        )
        .drop(columns="_strength_order")
        .reset_index(drop=True)
    )

    selected_columns = [
        "interaction_record_id",
        "mirna",
        "mirna_match_key",
        "target_gene",
        "mirna_species",
        "target_species",
        "maximum_evidence_strength",
        "evidence_strengths",
        "relationship_types",
        "experimental_methods",
        "support_types",
        "publication_ids",
        "n_publications",
        "n_experimental_methods",
        "n_source_records",
        "evidence_database",
        "database_version",
        "source_file",
    ]

    return aggregated[selected_columns]


def read_regulatory_sources(
    nodes_path: Path,
) -> pd.DataFrame:
    if not nodes_path.is_file():
        raise FileNotFoundError(
            f"Node registry not found: {nodes_path}"
        )

    nodes = pd.read_csv(
        nodes_path,
        sep="\t",
    )

    required = {
        "node_id",
        "feature_id",
        "feature_name",
        "feature_type",
        "dataset",
        "candidate_tier",
        "aggregated_reproducibility_score",
        "regulatory_source_eligible",
    }

    missing = sorted(
        required - set(nodes.columns)
    )

    if missing:
        raise ValueError(
            "Node registry is missing required columns: "
            + ", ".join(missing)
        )

    eligibility = nodes[
        "regulatory_source_eligible"
    ]

    if eligibility.dtype != bool:
        eligibility = (
            eligibility.astype(str)
            .str.lower()
            .isin({"true", "1", "yes"})
        )

    sources = nodes.loc[
        eligibility
    ].copy()

    sources["mirna_match_key"] = (
        sources["feature_name"]
        .map(mirna_match_key)
    )

    sources["canonical_mirna"] = (
        sources["feature_name"]
        .map(canonicalise_mirna)
    )

    selected = [
        "node_id",
        "feature_id",
        "feature_name",
        "canonical_mirna",
        "mirna_match_key",
        "dataset",
        "candidate_tier",
        "aggregated_reproducibility_score",
    ]

    return sources[selected].reset_index(drop=True)


def match_candidate_interactions(
    interactions: pd.DataFrame,
    sources: pd.DataFrame,
) -> pd.DataFrame:
    if interactions.empty or sources.empty:
        columns = list(interactions.columns) + [
            "mirna_node_id",
            "source_feature_id",
            "source_feature_name",
            "source_dataset",
            "source_candidate_tier",
            "source_reproducibility_score",
        ]

        return pd.DataFrame(columns=columns)

    matched = interactions.merge(
        sources[
            [
                "node_id",
                "feature_id",
                "feature_name",
                "mirna_match_key",
                "dataset",
                "candidate_tier",
                "aggregated_reproducibility_score",
            ]
        ],
        on="mirna_match_key",
        how="inner",
        validate="many_to_many",
    )

    matched = matched.rename(
        columns={
            "node_id": "mirna_node_id",
            "feature_id": "source_feature_id",
            "feature_name": "source_feature_name",
            "dataset": "source_dataset",
            "candidate_tier": "source_candidate_tier",
            "aggregated_reproducibility_score": (
                "source_reproducibility_score"
            ),
        }
    )

    selected_columns = [
        "interaction_record_id",
        "mirna_node_id",
        "source_feature_id",
        "source_feature_name",
        "mirna",
        "mirna_match_key",
        "target_gene",
        "source_dataset",
        "source_candidate_tier",
        "source_reproducibility_score",
        "maximum_evidence_strength",
        "evidence_strengths",
        "relationship_types",
        "experimental_methods",
        "support_types",
        "publication_ids",
        "n_publications",
        "n_experimental_methods",
        "n_source_records",
        "evidence_database",
        "database_version",
        "source_file",
    ]

    return (
        matched[selected_columns]
        .sort_values(
            [
                "source_reproducibility_score",
                "n_publications",
                "n_experimental_methods",
                "target_gene",
            ],
            ascending=[
                False,
                False,
                False,
                True,
            ],
        )
        .reset_index(drop=True)
    )


def build_summary(
    raw_count: int,
    standardised: pd.DataFrame,
    human: pd.DataFrame,
    aggregated: pd.DataFrame,
    sources: pd.DataFrame,
    candidate_interactions: pd.DataFrame,
) -> pd.DataFrame:
    rows = [
        {
            "metric": "raw_source_records",
            "value": raw_count,
        },
        {
            "metric": "standardised_nonempty_records",
            "value": len(standardised),
        },
        {
            "metric": "human_interaction_records",
            "value": len(human),
        },
        {
            "metric": "unique_human_interactions",
            "value": len(aggregated),
        },
        {
            "metric": "unique_human_mirnas",
            "value": aggregated["mirna"].nunique(),
        },
        {
            "metric": "unique_human_target_genes",
            "value": aggregated["target_gene"].nunique(),
        },
        {
            "metric": "regulatory_source_nodes",
            "value": len(sources),
        },
        {
            "metric": "candidate_mirna_target_interactions",
            "value": len(candidate_interactions),
        },
        {
            "metric": "candidate_target_genes",
            "value": candidate_interactions[
                "target_gene"
            ].nunique(),
        },
    ]

    evidence_counts = (
        aggregated["maximum_evidence_strength"]
        .value_counts()
    )

    for evidence_level, count in evidence_counts.items():
        rows.append(
            {
                "metric": (
                    "unique_interactions_"
                    f"{evidence_level}"
                ),
                "value": int(count),
            }
        )

    return pd.DataFrame(rows)


def write_report(
    output_path: Path,
    source_file: Path,
    source_description: str,
    detected_columns: dict[str, str | None],
    raw_count: int,
    standardised: pd.DataFrame,
    human: pd.DataFrame,
    aggregated: pd.DataFrame,
    sources: pd.DataFrame,
    candidate_interactions: pd.DataFrame,
) -> None:
    lines = [
        "miRTarBase import and standardisation",
        "====================================",
        "",
        f"Source file: {source_file}",
        f"Source format: {source_description}",
        f"Raw source records: {raw_count}",
        (
            "Standardised non-empty records: "
            f"{len(standardised)}"
        ),
        (
            "Human miRNA-target source records: "
            f"{len(human)}"
        ),
        (
            "Unique human miRNA-target interactions: "
            f"{len(aggregated)}"
        ),
        (
            "Unique human miRNAs: "
            f"{aggregated['mirna'].nunique()}"
        ),
        (
            "Unique human target genes: "
            f"{aggregated['target_gene'].nunique()}"
        ),
        "",
        "Detected source columns:",
    ]

    for logical_name, source_column in detected_columns.items():
        lines.append(
            f"  {logical_name}: "
            f"{source_column if source_column else 'not available'}"
        )

    lines.extend(
        [
            "",
            "Regulatory-source integration:",
            (
                "  Eligible miRNA source nodes: "
                f"{len(sources)}"
            ),
            (
                "  Matched candidate interactions: "
                f"{len(candidate_interactions)}"
            ),
            (
                "  Matched candidate target genes: "
                f"{candidate_interactions['target_gene'].nunique()}"
            ),
            "",
            "Eligible regulatory-source miRNAs:",
        ]
    )

    if sources.empty:
        lines.append("  none")
    else:
        for _, row in sources.iterrows():
            lines.append(
                "  "
                f"{row['feature_name']} "
                f"({row['canonical_mirna']}; "
                f"{row['candidate_tier']}; "
                f"score="
                f"{row['aggregated_reproducibility_score']})"
            )

    lines.extend(
        [
            "",
            "Evidence interpretation:",
            (
                "  - strong: strong support category or methods such as "
                "reporter assay, luciferase assay, Western blot, or qPCR."
            ),
            (
                "  - supportive: high-throughput or binding-support methods "
                "such as CLIP, sequencing, microarray, or proteomics."
            ),
            (
                "  - experimentally_supported_unspecified: an experiment "
                "or support category was provided but could not be assigned "
                "conservatively to the strong or supportive class."
            ),
            (
                "  - Evidence classes are derived annotations; raw support "
                "types and experimental methods are retained unchanged."
            ),
            (
                "  - No predicted miRNA-target relationship is generated."
            ),
            (
                "  - Target overlap with the proteomic candidate panel is "
                "performed in the next workflow step."
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

    if not args.mirtarbase.is_file():
        raise FileNotFoundError(
            "miRTarBase source file not found: "
            f"{args.mirtarbase}"
        )

    raw, source_description = read_mirtarbase(
        path=args.mirtarbase,
        sheet=args.sheet,
    )

    if raw.empty:
        raise ValueError(
            "The supplied miRTarBase table is empty."
        )

    detected = detect_columns(raw)

    standardised = build_standardised_table(
        raw=raw,
        detected=detected,
        source_file=args.mirtarbase,
    )

    if standardised.empty:
        raise ValueError(
            "No valid miRNA-target interaction remained after "
            "identifier standardisation."
        )

    human = select_human_interactions(
        standardised
    )

    if human.empty:
        raise ValueError(
            "No human miRNA-human target interaction was identified."
        )

    aggregated = aggregate_interactions(
        human
    )

    sources = read_regulatory_sources(
        args.nodes
    )

    candidate_interactions = (
        match_candidate_interactions(
            interactions=aggregated,
            sources=sources,
        )
    )

    summary = build_summary(
        raw_count=len(raw),
        standardised=standardised,
        human=human,
        aggregated=aggregated,
        sources=sources,
        candidate_interactions=candidate_interactions,
    )

    output_dir = args.output_dir
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    paths = {
        "standardised_records": (
            output_dir
            / "mirtarbase_standardised_source_records.tsv"
        ),
        "human_interactions": (
            output_dir
            / "mirtarbase_human_validated_interactions.tsv"
        ),
        "candidate_interactions": (
            output_dir
            / "candidate_mirna_validated_targets.tsv"
        ),
        "regulatory_sources": (
            output_dir
            / "regulatory_source_mirnas.tsv"
        ),
        "summary": (
            output_dir
            / "mirtarbase_import_summary.tsv"
        ),
        "report": (
            output_dir
            / "mirtarbase_import_report.txt"
        ),
    }

    standardised.to_csv(
        paths["standardised_records"],
        sep="\t",
        index=False,
    )

    aggregated.to_csv(
        paths["human_interactions"],
        sep="\t",
        index=False,
    )

    candidate_interactions.to_csv(
        paths["candidate_interactions"],
        sep="\t",
        index=False,
    )

    sources.to_csv(
        paths["regulatory_sources"],
        sep="\t",
        index=False,
    )

    summary.to_csv(
        paths["summary"],
        sep="\t",
        index=False,
    )

    write_report(
        output_path=paths["report"],
        source_file=args.mirtarbase,
        source_description=source_description,
        detected_columns=detected,
        raw_count=len(raw),
        standardised=standardised,
        human=human,
        aggregated=aggregated,
        sources=sources,
        candidate_interactions=candidate_interactions,
    )

    print(
        "Wrote standardised source records: "
        f"{paths['standardised_records']}"
    )

    print(
        "Wrote human validated interactions: "
        f"{paths['human_interactions']}"
    )

    print(
        "Wrote candidate miRNA targets: "
        f"{paths['candidate_interactions']}"
    )

    print(
        "Wrote regulatory-source miRNAs: "
        f"{paths['regulatory_sources']}"
    )

    print(
        f"Wrote import summary: "
        f"{paths['summary']}"
    )

    print(
        f"Wrote report: "
        f"{paths['report']}"
    )

    print(
        f"Raw source records: {len(raw)}"
    )

    print(
        "Unique human validated interactions: "
        f"{len(aggregated)}"
    )

    print(
        "Eligible regulatory-source miRNAs: "
        f"{len(sources)}"
    )

    print(
        "Candidate miRNA-target interactions: "
        f"{len(candidate_interactions)}"
    )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        FileNotFoundError,
        ValueError,
        KeyError,
        pd.errors.EmptyDataError,
        pd.errors.ParserError,
    ) as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        raise SystemExit(1)

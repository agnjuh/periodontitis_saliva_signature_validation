#!/usr/bin/env python3

"""
Prepare the measurable human proteomic universe for pathway enrichment.

The enrichment universe is derived from the filtered PXD043901 protein
annotation rather than from the complete human genome.

Eligible records
----------------
- organism_origin is "human" or "mixed"
- gene_symbol is present
- the gene symbol can be normalized safely
- duplicate gene symbols are collapsed

The script also verifies that all candidate-panel genes are contained in
the enrichment universe. A foreground gene set must be a subset of the
background universe for valid over-representation analysis.

Outputs
-------
- enrichment_universe_genes.txt
- enrichment_universe_mapping.tsv
- enrichment_universe_summary.tsv
- enrichment_universe_report.txt
- candidate_genes_missing_from_universe.txt
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


DEFAULT_ANNOTATION = Path(
    "data/processed/PXD043901/"
    "PXD043901_protein_annotation_extended.csv"
)

DEFAULT_CANDIDATE_GENES = Path(
    "results/integration/gene_sets/"
    "candidate_panel_genes.txt"
)

DEFAULT_OUTPUT_DIR = Path(
    "results/integration/gene_sets"
)

DEFAULT_INCLUDED_ORIGINS = (
    "human",
    "mixed",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare a measurable human proteomic background universe "
            "for pathway enrichment."
        )
    )

    parser.add_argument(
        "--annotation",
        type=Path,
        default=DEFAULT_ANNOTATION,
        help=(
            "Extended protein annotation CSV "
            f"(default: {DEFAULT_ANNOTATION})"
        ),
    )

    parser.add_argument(
        "--candidate-genes",
        type=Path,
        default=DEFAULT_CANDIDATE_GENES,
        help=(
            "Candidate-panel gene list used for universe validation "
            f"(default: {DEFAULT_CANDIDATE_GENES})"
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )

    parser.add_argument(
        "--include-origins",
        nargs="+",
        default=list(DEFAULT_INCLUDED_ORIGINS),
        help=(
            "organism_origin categories retained in the universe "
            "(default: human mixed)"
        ),
    )

    parser.add_argument(
        "--allow-missing-candidates",
        action="store_true",
        help=(
            "Do not fail when candidate-panel genes are absent from "
            "the enrichment universe. Missing genes are still reported."
        ),
    )

    return parser.parse_args()


def read_annotation(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Annotation file not found: {path}"
        )

    frame = pd.read_csv(
        path,
        dtype=str,
        keep_default_na=False,
    )

    if frame.empty:
        raise ValueError(
            f"Annotation file is empty: {path}"
        )

    return frame


def read_gene_list(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(
            f"Candidate gene file not found: {path}"
        )

    genes: list[str] = []

    for line in path.read_text(
        encoding="utf-8"
    ).splitlines():
        gene = clean_gene_symbol(line)

        if gene:
            genes.append(gene)

    return sorted(set(genes))


def clean_text(value: object) -> str:
    if value is None:
        return ""

    text = str(value).strip()

    if text.lower() in {
        "",
        "nan",
        "none",
        "na",
        "n/a",
        "null",
    }:
        return ""

    return text


def clean_gene_symbol(value: object) -> str:
    symbol = clean_text(value).upper()

    if not symbol:
        return ""

    symbol = re.sub(
        r"\s+",
        "",
        symbol,
    )

    return symbol


def split_gene_symbols(value: object) -> list[str]:
    """
    Split conservatively when an annotation cell contains multiple symbols.

    Common delimiters such as semicolon, comma, slash, and pipe are handled.
    Hyphens are not treated as delimiters because they may occur in aliases.
    """

    text = clean_text(value)

    if not text:
        return []

    parts = re.split(
        r"[;,/|]+",
        text,
    )

    genes = [
        clean_gene_symbol(part)
        for part in parts
    ]

    return sorted(
        {
            gene
            for gene in genes
            if gene
        }
    )


def write_gene_list(
    path: Path,
    genes: list[str],
) -> None:
    text = "\n".join(genes)

    if text:
        text += "\n"

    path.write_text(
        text,
        encoding="utf-8",
    )


def prepare_mapping(
    annotation: pd.DataFrame,
    included_origins: set[str],
) -> tuple[pd.DataFrame, dict[str, int]]:
    required_columns = {
        "protein_group_id",
        "primary_accession",
        "organism_origin",
        "gene_symbol",
    }

    missing = required_columns - set(
        annotation.columns
    )

    if missing:
        raise KeyError(
            "Annotation file is missing required columns: "
            + ", ".join(sorted(missing))
        )

    working = annotation.copy()

    working["organism_origin_normalized"] = (
        working["organism_origin"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    n_input_rows = len(working)

    origin_eligible = working[
        "organism_origin_normalized"
    ].isin(included_origins)

    origin_filtered = working.loc[
        origin_eligible
    ].copy()

    n_origin_eligible_rows = len(
        origin_filtered
    )

    origin_filtered["parsed_gene_symbols"] = (
        origin_filtered["gene_symbol"]
        .map(split_gene_symbols)
    )

    missing_symbol_mask = (
        origin_filtered[
            "parsed_gene_symbols"
        ].map(len).eq(0)
    )

    n_rows_without_gene_symbol = int(
        missing_symbol_mask.sum()
    )

    symbol_filtered = origin_filtered.loc[
        ~missing_symbol_mask
    ].copy()

    exploded = symbol_filtered.explode(
        "parsed_gene_symbols"
    ).rename(
        columns={
            "parsed_gene_symbols": "gene_symbol_normalized"
        }
    )

    exploded = exploded.loc[
        exploded["gene_symbol_normalized"].ne("")
    ].copy()

    optional_columns = [
        "Protein IDs",
        "Majority protein IDs",
        "Fasta headers",
        "human_protein_name",
        "detected_samples",
        "detected_fraction",
        "Unique peptides",
        "Razor + unique peptides",
        "Peptides",
        "Score",
        "Q-value",
    ]

    selected_columns = [
        "gene_symbol_normalized",
        "gene_symbol",
        "primary_accession",
        "protein_group_id",
        "organism_origin",
        "organism_origin_normalized",
    ]

    selected_columns.extend(
        column
        for column in optional_columns
        if column in exploded.columns
    )

    exploded = exploded[
        selected_columns
    ].copy()

    exploded = exploded.rename(
        columns={
            "gene_symbol_normalized": "gene_symbol",
            "gene_symbol": "original_gene_symbol",
        }
    )

    exploded = exploded.drop_duplicates(
        subset=[
            "gene_symbol",
            "primary_accession",
            "protein_group_id",
        ]
    )

    exploded = exploded.sort_values(
        [
            "gene_symbol",
            "organism_origin_normalized",
            "primary_accession",
            "protein_group_id",
        ]
    ).reset_index(drop=True)

    statistics = {
        "input_annotation_rows": n_input_rows,
        "origin_eligible_rows": n_origin_eligible_rows,
        "rows_without_gene_symbol": n_rows_without_gene_symbol,
        "mapped_protein_gene_records": len(exploded),
        "unique_universe_genes": exploded[
            "gene_symbol"
        ].nunique(),
    }

    return exploded, statistics


def build_summary(
    mapping: pd.DataFrame,
    statistics: dict[str, int],
    candidate_genes: list[str],
    missing_candidates: list[str],
    included_origins: list[str],
) -> pd.DataFrame:
    records: list[dict[str, object]] = [
        {
            "metric": "input_annotation_rows",
            "value": statistics[
                "input_annotation_rows"
            ],
        },
        {
            "metric": "origin_eligible_rows",
            "value": statistics[
                "origin_eligible_rows"
            ],
        },
        {
            "metric": "rows_without_gene_symbol",
            "value": statistics[
                "rows_without_gene_symbol"
            ],
        },
        {
            "metric": "mapped_protein_gene_records",
            "value": statistics[
                "mapped_protein_gene_records"
            ],
        },
        {
            "metric": "unique_universe_genes",
            "value": statistics[
                "unique_universe_genes"
            ],
        },
        {
            "metric": "candidate_panel_genes",
            "value": len(candidate_genes),
        },
        {
            "metric": "candidate_genes_missing_from_universe",
            "value": len(missing_candidates),
        },
        {
            "metric": "candidate_panel_coverage_fraction",
            "value": (
                (
                    len(candidate_genes)
                    - len(missing_candidates)
                )
                / len(candidate_genes)
                if candidate_genes
                else 0.0
            ),
        },
    ]

    for origin in included_origins:
        origin_genes = mapping.loc[
            mapping[
                "organism_origin_normalized"
            ].eq(origin),
            "gene_symbol",
        ].nunique()

        records.append(
            {
                "metric": (
                    f"unique_genes_with_{origin}_origin"
                ),
                "value": origin_genes,
            }
        )

    return pd.DataFrame(records)


def write_report(
    path: Path,
    *,
    annotation_path: Path,
    candidate_path: Path,
    included_origins: list[str],
    statistics: dict[str, int],
    universe_genes: list[str],
    candidate_genes: list[str],
    missing_candidates: list[str],
    mapping: pd.DataFrame,
) -> None:
    candidate_coverage = (
        (
            len(candidate_genes)
            - len(missing_candidates)
        )
        / len(candidate_genes)
        if candidate_genes
        else 0.0
    )

    lines = [
        "Proteomic enrichment-universe report",
        "====================================",
        "",
        f"Protein annotation: {annotation_path}",
        f"Candidate gene list: {candidate_path}",
        (
            "Included organism origins: "
            + ", ".join(included_origins)
        ),
        "",
        (
            "Input annotation rows: "
            f"{statistics['input_annotation_rows']}"
        ),
        (
            "Origin-eligible rows: "
            f"{statistics['origin_eligible_rows']}"
        ),
        (
            "Eligible rows without a gene symbol: "
            f"{statistics['rows_without_gene_symbol']}"
        ),
        (
            "Mapped protein-gene records: "
            f"{statistics['mapped_protein_gene_records']}"
        ),
        (
            "Unique enrichment-universe genes: "
            f"{len(universe_genes)}"
        ),
        "",
        f"Candidate-panel genes: {len(candidate_genes)}",
        (
            "Candidate genes represented in universe: "
            f"{len(candidate_genes) - len(missing_candidates)}"
        ),
        (
            "Candidate genes missing from universe: "
            f"{len(missing_candidates)}"
        ),
        (
            "Candidate-panel coverage fraction: "
            f"{candidate_coverage:.4f}"
        ),
        "",
        "Unique genes by annotation origin:",
    ]

    for origin in included_origins:
        n_genes = mapping.loc[
            mapping[
                "organism_origin_normalized"
            ].eq(origin),
            "gene_symbol",
        ].nunique()

        lines.append(
            f"  {origin}: {n_genes}"
        )

    lines.extend(
        [
            "",
            "Candidate genes missing from the universe:",
        ]
    )

    if missing_candidates:
        lines.extend(
            f"  {gene}"
            for gene in missing_candidates
        )
    else:
        lines.append("  None")

    lines.extend(
        [
            "",
            "Interpretation policy:",
            (
                "  The enrichment universe represents human genes "
                "measurable after the proteomic filtering workflow."
            ),
            (
                "  Human and mixed protein groups are retained only "
                "when a human gene symbol is available."
            ),
            (
                "  Microbiome-only protein groups are excluded from "
                "human pathway enrichment."
            ),
            (
                "  Candidate foreground sets must be subsets of this "
                "background universe."
            ),
            (
                "  This background is preferable to the complete "
                "human genome because it accounts for proteomic "
                "measurement and filtering constraints."
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
        included_origins = sorted(
            {
                clean_text(origin).lower()
                for origin in args.include_origins
                if clean_text(origin)
            }
        )

        if not included_origins:
            raise ValueError(
                "At least one organism origin must be included."
            )

        annotation = read_annotation(
            args.annotation
        )

        candidate_genes = read_gene_list(
            args.candidate_genes
        )

        mapping, statistics = prepare_mapping(
            annotation,
            set(included_origins),
        )

        universe_genes = sorted(
            mapping["gene_symbol"]
            .drop_duplicates()
            .tolist()
        )

        if not universe_genes:
            raise ValueError(
                "No genes were retained for the enrichment universe."
            )

        universe_set = set(universe_genes)

        missing_candidates = sorted(
            set(candidate_genes) - universe_set
        )

        args.output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        universe_path = (
            args.output_dir
            / "enrichment_universe_genes.txt"
        )

        mapping_path = (
            args.output_dir
            / "enrichment_universe_mapping.tsv"
        )

        summary_path = (
            args.output_dir
            / "enrichment_universe_summary.tsv"
        )

        report_path = (
            args.output_dir
            / "enrichment_universe_report.txt"
        )

        missing_path = (
            args.output_dir
            / "candidate_genes_missing_from_universe.txt"
        )

        write_gene_list(
            universe_path,
            universe_genes,
        )

        mapping.to_csv(
            mapping_path,
            sep="\t",
            index=False,
        )

        summary = build_summary(
            mapping=mapping,
            statistics=statistics,
            candidate_genes=candidate_genes,
            missing_candidates=missing_candidates,
            included_origins=included_origins,
        )

        summary.to_csv(
            summary_path,
            sep="\t",
            index=False,
        )

        write_gene_list(
            missing_path,
            missing_candidates,
        )

        write_report(
            report_path,
            annotation_path=args.annotation,
            candidate_path=args.candidate_genes,
            included_origins=included_origins,
            statistics=statistics,
            universe_genes=universe_genes,
            candidate_genes=candidate_genes,
            missing_candidates=missing_candidates,
            mapping=mapping,
        )

        print(
            f"Input annotation rows: "
            f"{statistics['input_annotation_rows']}"
        )
        print(
            f"Origin-eligible rows: "
            f"{statistics['origin_eligible_rows']}"
        )
        print(
            f"Unique enrichment-universe genes: "
            f"{len(universe_genes)}"
        )
        print(
            f"Candidate-panel genes: "
            f"{len(candidate_genes)}"
        )
        print(
            f"Candidate genes missing from universe: "
            f"{len(missing_candidates)}"
        )
        print(
            f"Wrote universe genes: {universe_path}"
        )
        print(
            f"Wrote universe mapping: {mapping_path}"
        )
        print(
            f"Wrote summary: {summary_path}"
        )
        print(
            f"Wrote report: {report_path}"
        )
        print(
            f"Wrote missing candidates: {missing_path}"
        )

        if (
            missing_candidates
            and not args.allow_missing_candidates
        ):
            print(
                "ERROR: Candidate-panel genes are missing from the "
                "enrichment universe. Inspect "
                f"{missing_path}. Re-run with "
                "--allow-missing-candidates only after confirming "
                "that the discrepancy is biologically and "
                "technically justified.",
                file=sys.stderr,
            )
            return 1

        return 0

    except Exception as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

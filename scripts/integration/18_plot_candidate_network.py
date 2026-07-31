#!/usr/bin/env python3

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch
import pandas as pd



# Input and output paths


REGULATORY_PATH = Path(
    "results/integration/regulatory_network/"
    "regulatory_edges_all.tsv"
)

ENRICHMENT_PATH = Path(
    "results/integration/pathway_enrichment/tables/"
    "significant_pathway_enrichment_results.tsv"
)

OUTPUT_DIR = Path(
    "results/integration/candidate_network/figures"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)



# Pathways retained in the focused figure


SELECTED_PATHWAYS = {
    "Antimicrobial peptides",
    "Signaling by Rho GTPases",
    "Cell-Cell communication",
    "RHO GTPases activate PAKs",
}



# Colour palette

ANTIMICROBIAL_NODE = "#9CCCF2"
ANTIMICROBIAL_EDGE = "#4F81A8"
ANTIMICROBIAL_PANEL = "#EAF5FC"

MIRNA_NODE = "#F2A59B"
MIRNA_EDGE = "#E6534B"

TARGET_NODE = "#C9B9E8"
TARGET_EDGE = "#7662AD"
TARGET_PANEL = "#F3EFFA"

PATHWAY_FILL = "#DCEFD5"
PATHWAY_BORDER = "#2E7D4A"
PATHWAY_TEXT = "#247A45"

NODE_BORDER = "#263238"
MAIN_TEXT = "#20252A"
NOTE_TEXT = "#555B60"



# Utility functions


def normalize_column_name(value: str) -> str:
    """Convert a column name into a predictable snake_case form."""

    value = str(value).strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize all DataFrame column names."""

    result = df.copy()
    result.columns = [
        normalize_column_name(column)
        for column in result.columns
    ]
    return result


def clean_string(value: object) -> str:
    """Return a stripped string, treating missing values as empty."""

    if pd.isna(value):
        return ""

    return str(value).strip()


def clean_gene(value: object) -> str:
    """Normalize a gene symbol."""

    return clean_string(value).upper()


def split_genes(value: object) -> list[str]:
    """Split semicolon-, comma-, slash-, or pipe-delimited gene symbols."""

    text = clean_string(value)

    if not text:
        return []

    genes = re.split(r"[;,/|]\s*", text)
    output: list[str] = []

    for gene in genes:
        gene = clean_gene(gene)

        if gene and gene not in output:
            output.append(gene)

    return output



# Input readers


def read_regulatory_edges() -> pd.DataFrame:
    """Read validated miRNA–candidate target relationships."""

    if not REGULATORY_PATH.exists():
        raise FileNotFoundError(
            f"Missing regulatory edge file: {REGULATORY_PATH}"
        )

    df = pd.read_csv(
        REGULATORY_PATH,
        sep="\t",
    )

    df = normalize_columns(df)

    required = {
        "source",
        "target",
        "source_type",
        "target_type",
        "candidate_gene_symbol",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            "Missing regulatory columns: "
            + ", ".join(sorted(missing))
            + "\nAvailable columns: "
            + ", ".join(df.columns)
        )

    subset = df.loc[
        df["candidate_gene_symbol"]
        .fillna("")
        .astype(str)
        .str.strip()
        .ne("")
    ].copy()

    def clean_mirna(value: object) -> str:
        text = clean_string(value)

        if ":" in text:
            text = text.split(":")[-1]

        text = text.strip()

        if text.lower().startswith("mir-"):
            text = "hsa-miR-" + text[4:]

        return text

    result = pd.DataFrame(
        {
            "mirna": subset["source"].map(clean_mirna),
            "target_gene": subset[
                "candidate_gene_symbol"
            ].map(clean_gene),
        }
    )

    result = result.loc[
        result["mirna"].ne("")
        & result["target_gene"].ne("")
    ].drop_duplicates()

    if result.empty:
        raise RuntimeError(
            "No usable miRNA-target edges were found in "
            f"{REGULATORY_PATH}"
        )

    return result


def read_pathway_memberships() -> pd.DataFrame:
    """Expand significant pathway results into gene–pathway memberships."""

    if not ENRICHMENT_PATH.exists():
        raise FileNotFoundError(
            f"Missing enrichment file: {ENRICHMENT_PATH}"
        )

    df = pd.read_csv(
        ENRICHMENT_PATH,
        sep="\t",
    )

    df = normalize_columns(df)

    required = {
        "gene_set",
        "database",
        "description",
        "gene_symbols",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            "Missing enrichment columns: "
            + ", ".join(sorted(missing))
            + "\nAvailable columns: "
            + ", ".join(df.columns)
        )

    df = df.loc[
        df["description"].isin(SELECTED_PATHWAYS)
    ].copy()

    records: list[dict[str, str]] = []

    for _, row in df.iterrows():
        for gene in split_genes(row["gene_symbols"]):
            records.append(
                {
                    "gene": gene,
                    "pathway": clean_string(
                        row["description"]
                    ),
                    "gene_set": clean_string(
                        row["gene_set"]
                    ),
                    "database": clean_string(
                        row["database"]
                    ),
                }
            )

    if not records:
        raise RuntimeError(
            "No memberships were found for the selected pathways."
        )

    return (
        pd.DataFrame(records)
        .drop_duplicates()
        .reset_index(drop=True)
    )



# Drawing helpers


def draw_protein_node(
    ax,
    x: float,
    y: float,
    label: str,
    facecolor: str,
    radius: float = 0.285,
) -> None:
    """
    Draw a larger circular protein node.

    The radius and reduced font size are selected so labels such as
    SPTBN1 and S100A8 remain inside the node.
    """

    node = Circle(
        (x, y),
        radius,
        facecolor=facecolor,
        edgecolor=NODE_BORDER,
        linewidth=1.35,
        zorder=5,
    )

    ax.add_patch(node)

    font_size = 7.8

    if len(label) >= 7:
        font_size = 7.0
    elif len(label) >= 6:
        font_size = 7.4

    ax.text(
        x,
        y,
        label,
        ha="center",
        va="center",
        fontsize=font_size,
        fontweight="semibold",
        color=MAIN_TEXT,
        zorder=6,
    )


def draw_mirna_node(
    ax,
    x: float,
    y: float,
    label: str,
) -> None:
    """Draw the central miRNA node."""

    width = 1.42
    height = 0.54

    node = FancyBboxPatch(
        (
            x - width / 2,
            y - height / 2,
        ),
        width,
        height,
        boxstyle=(
            "round,pad=0.05,"
            "rounding_size=0.07"
        ),
        facecolor=MIRNA_NODE,
        edgecolor=NODE_BORDER,
        linewidth=1.5,
        zorder=5,
    )

    ax.add_patch(node)

    ax.text(
        x,
        y,
        label,
        ha="center",
        va="center",
        fontsize=10.5,
        fontweight="bold",
        color=MAIN_TEXT,
        zorder=6,
    )


def draw_pathway_capsule(
    ax,
    x: float,
    y: float,
    label: str,
    width: float,
    height: float = 0.65,
    fontsize: float = 10.5,
) -> None:
    """
    Draw a compact rounded pathway label.

    The capsule has a defined boundary so arrows can terminate at its
    edge instead of entering the text.
    """

    pathway = FancyBboxPatch(
        (
            x - width / 2,
            y - height / 2,
        ),
        width,
        height,
        boxstyle=(
            "round,pad=0.06,"
            "rounding_size=0.20"
        ),
        facecolor=PATHWAY_FILL,
        edgecolor=PATHWAY_BORDER,
        linewidth=1.15,
        zorder=4,
    )

    ax.add_patch(pathway)

    ax.text(
        x,
        y,
        label,
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight="bold",
        color=PATHWAY_TEXT,
        linespacing=1.10,
        zorder=5,
    )


def draw_arrow(
    ax,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str,
    linewidth: float,
    mutation_scale: float = 10,
    connectionstyle: str = "arc3,rad=0.0",
    shrink_a: float = 8,
    shrink_b: float = 5,
) -> None:
    """Draw a directed edge."""

    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops={
            "arrowstyle": "-|>",
            "color": color,
            "linewidth": linewidth,
            "mutation_scale": mutation_scale,
            "shrinkA": shrink_a,
            "shrinkB": shrink_b,
            "connectionstyle": connectionstyle,
        },
        zorder=2,
    )


def add_module_panel(
    ax,
    x: float,
    y: float,
    width: float,
    height: float,
    facecolor: str,
) -> None:
    """Add a subtle rounded background panel."""

    panel = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle=(
            "round,pad=0.10,"
            "rounding_size=0.24"
        ),
        facecolor=facecolor,
        edgecolor="none",
        alpha=0.85,
        zorder=0,
    )

    ax.add_patch(panel)



# Main plot


def main() -> None:
    regulatory = read_regulatory_edges()
    memberships = read_pathway_memberships()

    antimicrobial_genes = sorted(
        memberships.loc[
            memberships["pathway"].eq(
                "Antimicrobial peptides"
            ),
            "gene",
        ].unique()
    )

    cytoskeletal_pathways = [
        "Signaling by Rho GTPases",
        "Cell-Cell communication",
        "RHO GTPases activate PAKs",
    ]

    cytoskeletal_genes = sorted(
        memberships.loc[
            memberships["pathway"].isin(
                cytoskeletal_pathways
            ),
            "gene",
        ].unique()
    )

    mirnas = sorted(
        regulatory["mirna"].unique()
    )

    if not mirnas:
        raise RuntimeError(
            "No miRNA node was detected."
        )

    mirna = mirnas[0]

    regulatory_targets = sorted(
        regulatory.loc[
            regulatory["mirna"].eq(mirna),
            "target_gene",
        ].unique()
    )

    displayed_targets = [
        gene
        for gene in regulatory_targets
        if gene in cytoskeletal_genes
    ]

    additional_targets = [
        gene
        for gene in regulatory_targets
        if gene not in displayed_targets
    ]

    fig, ax = plt.subplots(
        figsize=(15.5, 9.2)
    )

    ax.set_xlim(-6.5, 6.5)
    ax.set_ylim(-4.6, 4.6)
    ax.axis("off")

    
    # Module backgrounds
    

    add_module_panel(
        ax=ax,
        x=-6.05,
        y=-3.35,
        width=4.30,
        height=6.65,
        facecolor=ANTIMICROBIAL_PANEL,
    )

    add_module_panel(
        ax=ax,
        x=0.65,
        y=-3.35,
        width=5.40,
        height=6.65,
        facecolor=TARGET_PANEL,
    )

    
    # Main title and module headings
    

    ax.set_title(
        "Candidate miRNA–protein–pathway network",
        fontsize=18,
        fontweight="bold",
        color=MAIN_TEXT,
        pad=18,
    )

    ax.text(
        -3.90,
        3.62,
        "Innate antimicrobial defence",
        ha="center",
        fontsize=12.5,
        fontweight="bold",
        color=ANTIMICROBIAL_EDGE,
    )

    ax.text(
        3.25,
        3.62,
        "Cytoskeletal / Rho-linked module",
        ha="center",
        fontsize=12.5,
        fontweight="bold",
        color=TARGET_EDGE,
    )

    
    # Antimicrobial protein nodes
   

    antimicrobial_protein_x = -5.15

    antimicrobial_pathway_x = -2.55
    antimicrobial_pathway_y = 0.20

    antimicrobial_pathway_width = 1.60
    antimicrobial_pathway_left_edge = (
        antimicrobial_pathway_x
        - antimicrobial_pathway_width / 2
    )

    antimicrobial_y_positions = [
        2.65
        - index
        * (
            5.10
            / max(
                len(antimicrobial_genes) - 1,
                1,
            )
        )
        for index in range(
            len(antimicrobial_genes)
        )
    ]

    for gene, y in zip(
        antimicrobial_genes,
        antimicrobial_y_positions,
    ):
        draw_protein_node(
            ax=ax,
            x=antimicrobial_protein_x,
            y=y,
            label=gene,
            facecolor=ANTIMICROBIAL_NODE,
            radius=0.285,
        )

        # All arrows terminate at the left boundary of the pathway
        # capsule. They therefore cannot overlap its label.
        draw_arrow(
            ax=ax,
            start=(
                antimicrobial_protein_x
                + 0.30,
                y,
            ),
            end=(
                antimicrobial_pathway_left_edge
                - 0.04,
                antimicrobial_pathway_y,
            ),
            color=ANTIMICROBIAL_EDGE,
            linewidth=1.10,
            mutation_scale=9,
            shrink_a=3,
            shrink_b=1,
        )

    draw_pathway_capsule(
        ax=ax,
        x=antimicrobial_pathway_x,
        y=antimicrobial_pathway_y,
        label="Antimicrobial\npeptides",
        width=antimicrobial_pathway_width,
        height=0.72,
        fontsize=10.5,
    )

    
    # Central miRNA node
    

    mirna_x = -0.42
    mirna_y = 0.20

    draw_mirna_node(
        ax=ax,
        x=mirna_x,
        y=mirna_y,
        label=mirna,
    )

    
    # Regulatory target protein nodes
    

    target_x = 1.55

    target_y_positions = [
        2.60
        - index
        * (
            4.85
            / max(
                len(displayed_targets) - 1,
                1,
            )
        )
        for index in range(
            len(displayed_targets)
        )
    ]

    target_positions: dict[
        str,
        tuple[float, float],
    ] = {}

    for gene, y in zip(
        displayed_targets,
        target_y_positions,
    ):
        target_positions[gene] = (
            target_x,
            y,
        )

        draw_protein_node(
            ax=ax,
            x=target_x,
            y=y,
            label=gene,
            facecolor=TARGET_NODE,
            radius=0.285,
        )

        draw_arrow(
            ax=ax,
            start=(
                mirna_x + 0.72,
                mirna_y,
            ),
            end=(
                target_x - 0.30,
                y,
            ),
            color=MIRNA_EDGE,
            linewidth=1.55,
            mutation_scale=11,
            connectionstyle="arc3,rad=0.025",
            shrink_a=2,
            shrink_b=2,
        )

    
    # Right-hand pathway capsules
    

    pathway_x = 4.90
    pathway_width = 1.80

    pathway_positions = {
        "Signaling by Rho GTPases": (
            pathway_x,
            2.35,
        ),
        "Cell-Cell communication": (
            pathway_x,
            0.20,
        ),
        "RHO GTPases activate PAKs": (
            pathway_x,
            -1.95,
        ),
    }

    pathway_labels = {
        "Signaling by Rho GTPases":
            "Signaling by\nRho GTPases",
        "Cell-Cell communication":
            "Cell–cell\ncommunication",
        "RHO GTPases activate PAKs":
            "RHO GTPases\nactivate PAKs",
    }

    for pathway, (x, y) in pathway_positions.items():
        draw_pathway_capsule(
            ax=ax,
            x=x,
            y=y,
            label=pathway_labels[pathway],
            width=pathway_width,
            height=0.72,
            fontsize=10.5,
        )

    pathway_left_edge = (
        pathway_x - pathway_width / 2
    )

    
    # Protein-to-pathway membership edges
    

    for _, row in memberships.iterrows():
        gene = row["gene"]
        pathway = row["pathway"]

        if gene not in target_positions:
            continue

        if pathway not in pathway_positions:
            continue

        start_x, start_y = target_positions[gene]
        _, end_y = pathway_positions[pathway]

        draw_arrow(
            ax=ax,
            start=(
                start_x + 0.30,
                start_y,
            ),
            end=(
                pathway_left_edge - 0.04,
                end_y,
            ),
            color=TARGET_EDGE,
            linewidth=1.0,
            mutation_scale=9,
            connectionstyle="arc3,rad=0.015",
            shrink_a=2,
            shrink_b=1,
        )

    
    # Additional validated targets not present in selected pathways
    
    if additional_targets:
        ax.text(
            0.35,
            -3.72,
            "Additional validated targets: "
            + ", ".join(additional_targets),
            ha="center",
            va="center",
            fontsize=9,
            color=NOTE_TEXT,
        )

    
    # Legend
    

    legend_y = -4.18

    ax.scatter(
        [-5.35],
        [legend_y],
        s=145,
        facecolor=ANTIMICROBIAL_NODE,
        edgecolor=NODE_BORDER,
        linewidth=1,
        zorder=5,
    )

    ax.text(
        -5.05,
        legend_y,
        "antimicrobial protein",
        fontsize=8.5,
        va="center",
        color=MAIN_TEXT,
    )

    ax.scatter(
        [-3.20],
        [legend_y],
        s=145,
        facecolor=TARGET_NODE,
        edgecolor=NODE_BORDER,
        linewidth=1,
        zorder=5,
    )

    ax.text(
        -2.90,
        legend_y,
        "validated miRNA target",
        fontsize=8.5,
        va="center",
        color=MAIN_TEXT,
    )

    ax.plot(
        [-0.65, -0.20],
        [legend_y, legend_y],
        color=MIRNA_EDGE,
        linewidth=1.8,
    )

    ax.text(
        -0.05,
        legend_y,
        "miRNA–target support",
        fontsize=8.5,
        va="center",
        color=MAIN_TEXT,
    )

    ax.plot(
        [1.80, 2.25],
        [legend_y, legend_y],
        color=TARGET_EDGE,
        linewidth=1.5,
    )

    ax.text(
        2.40,
        legend_y,
        "pathway membership",
        fontsize=8.5,
        va="center",
        color=MAIN_TEXT,
    )

    ax.text(
        4.25,
        legend_y,
        "Enrichment does not imply activation or causality.",
        fontsize=8,
        style="italic",
        color=NOTE_TEXT,
        va="center",
    )

    plt.tight_layout()

    
    # Output files
    

    png_path = (
        OUTPUT_DIR
        / "candidate_network_coloured_v2.png"
    )

    pdf_path = (
        OUTPUT_DIR
        / "candidate_network_coloured_v2.pdf"
    )

    svg_path = (
        OUTPUT_DIR
        / "candidate_network_coloured_v2.svg"
    )

    plt.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.savefig(
        pdf_path,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.savefig(
        svg_path,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close()

    
    # Console report
    

    print(
        "Loaded regulatory targets:",
        ", ".join(regulatory_targets),
    )

    print(
        f"Regulatory edges loaded: "
        f"{len(regulatory)}"
    )

    print(
        f"Antimicrobial genes displayed: "
        f"{len(antimicrobial_genes)}"
    )

    print(
        "Regulatory pathway-linked targets displayed: "
        f"{len(displayed_targets)}"
    )

    print(
        f"Additional regulatory targets listed: "
        f"{len(additional_targets)}"
    )

    print(f"Wrote PNG: {png_path}")
    print(f"Wrote PDF: {pdf_path}")
    print(f"Wrote SVG: {svg_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env Rscript

# ==============================================================================
# Pathway over-representation analysis
# ==============================================================================
#
# This script performs pathway enrichment for standardized candidate gene sets
# using the measurable human saliva proteome as the statistical background.
#
# Gene sets
# ---------
# 1. candidate_panel
# 2. mir143_linked
# 3. nonlinked_candidate
#
# Databases
# ---------
# - Gene Ontology Biological Process
# - Reactome
# - KEGG
#
# Statistical background
# ----------------------
# The enrichment universe is the filtered, gene-annotated human proteome from
# PXD043901 rather than the complete human genome.
#
# Outputs
# -------
# results/integration/pathway_enrichment/
#   mapping/
#   tables/
#   figures/
#   pathway_enrichment_summary.tsv
#   pathway_enrichment_report.txt
#
# Important interpretation
# ------------------------
# The miR-143-3p-linked gene set contains only seven genes. Its enrichment
# results are therefore exploratory and may be unstable.
# ==============================================================================



# Utility functions


stop_with_message <- function(...) {
    stop(
        paste0(...),
        call. = FALSE
    )
}


require_package <- function(package_name) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
        stop_with_message(
            "Required R package is not installed: ",
            package_name
        )
    }
}


required_packages <- c(
    "clusterProfiler",
    "ReactomePA",
    "org.Hs.eg.db",
    "enrichplot",
    "AnnotationDbi",
    "dplyr",
    "readr",
    "ggplot2",
    "stringr",
    "tibble"
)

invisible(
    lapply(
        required_packages,
        require_package
    )
)


suppressPackageStartupMessages({
    library(clusterProfiler)
    library(ReactomePA)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(AnnotationDbi)
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(stringr)
    library(tibble)
})



# Command-line argument parsing


parse_arguments <- function(args) {
    defaults <- list(
        gene_set_dir = "results/integration/gene_sets",
        output_dir = "results/integration/pathway_enrichment",
        pvalue_cutoff = 1.0,
        qvalue_cutoff = 1.0,
        significance_fdr = 0.05,
        minimum_gene_set_size = 3,
        maximum_gene_set_size = 500,
        show_categories = 20
    )

    if (length(args) == 0) {
        return(defaults)
    }

    index <- 1

    while (index <= length(args)) {
        argument <- args[[index]]

        if (argument == "--gene-set-dir") {
            index <- index + 1
            defaults$gene_set_dir <- args[[index]]

        } else if (argument == "--output-dir") {
            index <- index + 1
            defaults$output_dir <- args[[index]]

        } else if (argument == "--pvalue-cutoff") {
            index <- index + 1
            defaults$pvalue_cutoff <- as.numeric(args[[index]])

        } else if (argument == "--qvalue-cutoff") {
            index <- index + 1
            defaults$qvalue_cutoff <- as.numeric(args[[index]])

        } else if (argument == "--significance-fdr") {
            index <- index + 1
            defaults$significance_fdr <- as.numeric(args[[index]])

        } else if (argument == "--minimum-gene-set-size") {
            index <- index + 1
            defaults$minimum_gene_set_size <- as.integer(args[[index]])

        } else if (argument == "--maximum-gene-set-size") {
            index <- index + 1
            defaults$maximum_gene_set_size <- as.integer(args[[index]])

        } else if (argument == "--show-categories") {
            index <- index + 1
            defaults$show_categories <- as.integer(args[[index]])

        } else if (argument %in% c("--help", "-h")) {
            cat(
                paste(
                    "Usage:",
                    "Rscript scripts/integration/16_run_pathway_enrichment.R",
                    "[--gene-set-dir PATH]",
                    "[--output-dir PATH]",
                    "[--pvalue-cutoff FLOAT]",
                    "[--qvalue-cutoff FLOAT]",
                    "[--significance-fdr FLOAT]",
                    "[--minimum-gene-set-size INTEGER]",
                    "[--maximum-gene-set-size INTEGER]",
                    "[--show-categories INTEGER]",
                    sep = "\n  "
                ),
                "\n"
            )

            quit(
                save = "no",
                status = 0
            )

        } else {
            stop_with_message(
                "Unknown command-line argument: ",
                argument
            )
        }

        index <- index + 1
    }

    defaults
}


config <- parse_arguments(
    commandArgs(
        trailingOnly = TRUE
    )
)



# Paths


gene_set_dir <- config$gene_set_dir
output_dir <- config$output_dir

mapping_dir <- file.path(
    output_dir,
    "mapping"
)

table_dir <- file.path(
    output_dir,
    "tables"
)

figure_dir <- file.path(
    output_dir,
    "figures"
)

dir.create(
    mapping_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    table_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


gene_set_paths <- c(
    candidate_panel = file.path(
        gene_set_dir,
        "candidate_panel_genes.txt"
    ),
    mir143_linked = file.path(
        gene_set_dir,
        "mir143_linked_genes.txt"
    ),
    nonlinked_candidate = file.path(
        gene_set_dir,
        "nonlinked_candidate_genes.txt"
    )
)

universe_path <- file.path(
    gene_set_dir,
    "enrichment_universe_genes.txt"
)



# Input validation


all_required_paths <- c(
    gene_set_paths,
    enrichment_universe = universe_path
)

missing_paths <- all_required_paths[
    !file.exists(all_required_paths)
]

if (length(missing_paths) > 0) {
    stop_with_message(
        "Missing required input files:\n",
        paste(
            paste0("  ", missing_paths),
            collapse = "\n"
        )
    )
}


read_gene_list <- function(path) {
    genes <- readLines(
        path,
        warn = FALSE,
        encoding = "UTF-8"
    )

    genes <- genes |>
        str_trim() |>
        str_to_upper()

    genes <- genes[
        !is.na(genes) &
        genes != ""
    ]

    sort(
        unique(genes)
    )
}


gene_sets_symbol <- lapply(
    gene_set_paths,
    read_gene_list
)

universe_symbol <- read_gene_list(
    universe_path
)


if (length(universe_symbol) == 0) {
    stop_with_message(
        "The enrichment universe is empty."
    )
}


for (gene_set_name in names(gene_sets_symbol)) {
    current_genes <- gene_sets_symbol[[gene_set_name]]

    if (length(current_genes) == 0) {
        stop_with_message(
            "Gene set is empty: ",
            gene_set_name
        )
    }

    missing_from_universe <- setdiff(
        current_genes,
        universe_symbol
    )

    if (length(missing_from_universe) > 0) {
        stop_with_message(
            "Gene set '",
            gene_set_name,
            "' contains genes absent from the enrichment universe:\n  ",
            paste(
                missing_from_universe,
                collapse = ", "
            )
        )
    }
}



# HGNC symbol to Entrez mapping


all_symbols <- sort(
    unique(
        c(
            universe_symbol,
            unlist(
                gene_sets_symbol,
                use.names = FALSE
            )
        )
    )
)


mapping_raw <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = all_symbols,
    keytype = "SYMBOL",
    columns = c(
        "SYMBOL",
        "ENTREZID",
        "GENENAME"
    )
)


mapping_table <- mapping_raw |>
    as_tibble() |>
    filter(
        !is.na(SYMBOL),
        SYMBOL != ""
    ) |>
    mutate(
        SYMBOL = str_to_upper(
            str_trim(SYMBOL)
        ),
        ENTREZID = as.character(
            ENTREZID
        )
    ) |>
    distinct(
        SYMBOL,
        ENTREZID,
        GENENAME
    ) |>
    arrange(
        SYMBOL,
        ENTREZID
    )


write_tsv(
    mapping_table,
    file.path(
        mapping_dir,
        "symbol_to_entrez_mapping.tsv"
    )
)


mapped_symbols <- mapping_table |>
    filter(
        !is.na(ENTREZID),
        ENTREZID != ""
    ) |>
    pull(SYMBOL) |>
    unique()


unmapped_symbols <- setdiff(
    all_symbols,
    mapped_symbols
)


write_lines <- function(values, path) {
    writeLines(
        values,
        con = path,
        useBytes = TRUE
    )
}


write_lines(
    sort(unmapped_symbols),
    file.path(
        mapping_dir,
        "unmapped_gene_symbols.txt"
    )
)


symbol_to_entrez <- function(symbols) {
    mapping_table |>
        filter(
            SYMBOL %in% symbols,
            !is.na(ENTREZID),
            ENTREZID != ""
        ) |>
        pull(ENTREZID) |>
        unique() |>
        sort()
}


universe_entrez <- symbol_to_entrez(
    universe_symbol
)


gene_sets_entrez <- lapply(
    gene_sets_symbol,
    symbol_to_entrez
)


if (length(universe_entrez) == 0) {
    stop_with_message(
        "None of the enrichment-universe genes mapped to Entrez IDs."
    )
}


for (gene_set_name in names(gene_sets_entrez)) {
    if (length(gene_sets_entrez[[gene_set_name]]) == 0) {
        stop_with_message(
            "No genes mapped to Entrez IDs for gene set: ",
            gene_set_name
        )
    }
}


mapping_summary <- tibble(
    gene_set = c(
        names(gene_sets_symbol),
        "enrichment_universe"
    ),
    input_symbols = c(
        vapply(
            gene_sets_symbol,
            length,
            integer(1)
        ),
        length(universe_symbol)
    ),
    mapped_entrez_ids = c(
        vapply(
            gene_sets_entrez,
            length,
            integer(1)
        ),
        length(universe_entrez)
    )
) |>
    mutate(
        mapping_fraction = mapped_entrez_ids / input_symbols
    )


write_tsv(
    mapping_summary,
    file.path(
        mapping_dir,
        "gene_mapping_summary.tsv"
    )
)


for (gene_set_name in names(gene_sets_entrez)) {
    write_lines(
        gene_sets_entrez[[gene_set_name]],
        file.path(
            mapping_dir,
            paste0(
                gene_set_name,
                "_entrez_ids.txt"
            )
        )
    )
}


write_lines(
    universe_entrez,
    file.path(
        mapping_dir,
        "enrichment_universe_entrez_ids.txt"
    )
)



# Enrichment helper functions


empty_enrichment_table <- function() {
    tibble(
        ID = character(),
        Description = character(),
        GeneRatio = character(),
        BgRatio = character(),
        pvalue = numeric(),
        p.adjust = numeric(),
        qvalue = numeric(),
        geneID = character(),
        Count = integer(),
        foreground_count = integer(),
        foreground_size = integer(),
        background_count = integer(),
        background_size = integer(),
        foreground_fraction = numeric(),
        background_fraction = numeric(),
        fold_enrichment = numeric(),
        gene_symbols = character(),
        significant_fdr_0_05 = logical()
    )
}


parse_ratio_numerator <- function(ratio_value) {
    as.numeric(
        sub(
            "/.*$",
            "",
            ratio_value
        )
    )
}


parse_ratio_denominator <- function(ratio_value) {
    as.numeric(
        sub(
            "^.*/",
            "",
            ratio_value
        )
    )
}


entrez_string_to_symbols <- function(entrez_string) {
    if (
        is.na(entrez_string) ||
        entrez_string == ""
    ) {
        return("")
    }

    entrez_ids <- unique(
        unlist(
            strsplit(
                entrez_string,
                "/",
                fixed = TRUE
            )
        )
    )

    mapped <- mapping_table |>
        filter(
            ENTREZID %in% entrez_ids
        ) |>
        pull(SYMBOL) |>
        unique() |>
        sort()

    paste(
        mapped,
        collapse = ";"
    )
}


standardize_enrichment_result <- function(
    enrichment_object,
    significance_fdr
) {
    if (
        is.null(enrichment_object) ||
        nrow(as.data.frame(enrichment_object)) == 0
    ) {
        return(
            empty_enrichment_table()
        )
    }

    result <- as.data.frame(
        enrichment_object
    ) |>
        as_tibble()

    required_columns <- c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "pvalue",
        "p.adjust",
        "qvalue",
        "geneID",
        "Count"
    )

    for (column_name in required_columns) {
        if (!column_name %in% colnames(result)) {
            if (
                column_name %in%
                c("pvalue", "p.adjust", "qvalue")
            ) {
                result[[column_name]] <- NA_real_
            } else if (column_name == "Count") {
                result[[column_name]] <- NA_integer_
            } else {
                result[[column_name]] <- NA_character_
            }
        }
    }

    result |>
        mutate(
            foreground_count = parse_ratio_numerator(
                GeneRatio
            ),
            foreground_size = parse_ratio_denominator(
                GeneRatio
            ),
            background_count = parse_ratio_numerator(
                BgRatio
            ),
            background_size = parse_ratio_denominator(
                BgRatio
            ),
            foreground_fraction = ifelse(
                foreground_size > 0,
                foreground_count / foreground_size,
                NA_real_
            ),
            background_fraction = ifelse(
                background_size > 0,
                background_count / background_size,
                NA_real_
            ),
            fold_enrichment = ifelse(
                background_fraction > 0,
                foreground_fraction / background_fraction,
                NA_real_
            ),
            gene_symbols = vapply(
                geneID,
                entrez_string_to_symbols,
                character(1)
            ),
            significant_fdr_0_05 = (
                !is.na(p.adjust) &
                p.adjust <= significance_fdr
            )
        ) |>
        arrange(
            p.adjust,
            pvalue,
            desc(fold_enrichment)
        )
}


safe_enrich_go <- function(
    genes,
    universe,
    config
) {
    tryCatch(
        enrichGO(
            gene = genes,
            universe = universe,
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = config$pvalue_cutoff,
            qvalueCutoff = config$qvalue_cutoff,
            minGSSize = config$minimum_gene_set_size,
            maxGSSize = config$maximum_gene_set_size,
            readable = FALSE
        ),
        error = function(error) {
            warning(
                paste0(
                    "GO enrichment failed: ",
                    conditionMessage(error)
                )
            )

            NULL
        }
    )
}


safe_enrich_reactome <- function(
    genes,
    universe,
    config
) {
    tryCatch(
        enrichPathway(
            gene = genes,
            universe = universe,
            organism = "human",
            pAdjustMethod = "BH",
            pvalueCutoff = config$pvalue_cutoff,
            qvalueCutoff = config$qvalue_cutoff,
            minGSSize = config$minimum_gene_set_size,
            maxGSSize = config$maximum_gene_set_size,
            readable = FALSE
        ),
        error = function(error) {
            warning(
                paste0(
                    "Reactome enrichment failed: ",
                    conditionMessage(error)
                )
            )

            NULL
        }
    )
}


safe_enrich_kegg <- function(
    genes,
    universe,
    config
) {
    tryCatch(
        enrichKEGG(
            gene = genes,
            universe = universe,
            organism = "hsa",
            keyType = "ncbi-geneid",
            pAdjustMethod = "BH",
            pvalueCutoff = config$pvalue_cutoff,
            qvalueCutoff = config$qvalue_cutoff,
            minGSSize = config$minimum_gene_set_size,
            maxGSSize = config$maximum_gene_set_size,
            use_internal_data = FALSE
        ),
        error = function(error) {
            warning(
                paste0(
                    "KEGG enrichment failed. This database may require ",
                    "internet access. Error: ",
                    conditionMessage(error)
                )
            )

            NULL
        }
    )
}



# Plotting


make_dotplot <- function(
    result_table,
    gene_set_name,
    database_name,
    output_path,
    show_categories
) {
    plot_data <- result_table |>
        filter(
            !is.na(p.adjust)
        ) |>
        arrange(
            p.adjust,
            desc(fold_enrichment)
        ) |>
        slice_head(
            n = show_categories
        )

    if (nrow(plot_data) == 0) {
        return(FALSE)
    }

    plot_data <- plot_data |>
        mutate(
            Description = factor(
                Description,
                levels = rev(
                    unique(Description)
                )
            ),
            minus_log10_fdr = -log10(
                pmax(
                    p.adjust,
                    .Machine$double.xmin
                )
            )
        )

    plot_object <- ggplot(
        plot_data,
        aes(
            x = fold_enrichment,
            y = Description,
            size = foreground_count,
            fill = minus_log10_fdr
        )
    ) +
        geom_point(
            shape = 21
        ) +
        scale_size_continuous(
            name = "Gene count"
        ) +
        scale_fill_viridis_c(
            name = expression(
                -log[10](FDR)
            )
        ) +
        labs(
            title = paste0(
                gene_set_name,
                ": ",
                database_name
            ),
            subtitle = paste0(
                "Background: measurable PXD043901 human proteome; n = ",
                length(universe_entrez)
            ),
            x = "Fold enrichment",
            y = NULL
        ) +
        theme_bw(
            base_size = 10
        ) +
        theme(
            plot.title = element_text(
                face = "bold"
            ),
            axis.text.y = element_text(
                size = 8
            ),
            panel.grid.minor = element_blank()
        )

    ggsave(
        filename = output_path,
        plot = plot_object,
        width = 9,
        height = max(
            5,
            0.32 * nrow(plot_data) + 2
        ),
        units = "in",
        device = cairo_pdf
    )

    TRUE
}


make_barplot <- function(
    result_table,
    gene_set_name,
    database_name,
    output_path,
    show_categories
) {
    plot_data <- result_table |>
        filter(
            !is.na(p.adjust)
        ) |>
        arrange(
            p.adjust,
            desc(fold_enrichment)
        ) |>
        slice_head(
            n = show_categories
        )

    if (nrow(plot_data) == 0) {
        return(FALSE)
    }

    plot_data <- plot_data |>
        mutate(
            Description = factor(
                Description,
                levels = rev(
                    unique(Description)
                )
            ),
            minus_log10_fdr = -log10(
                pmax(
                    p.adjust,
                    .Machine$double.xmin
                )
            )
        )

    plot_object <- ggplot(
        plot_data,
        aes(
            x = minus_log10_fdr,
            y = Description
        )
    ) +
        geom_col() +
        labs(
            title = paste0(
                gene_set_name,
                ": ",
                database_name
            ),
            subtitle = paste0(
                "Background: measurable PXD043901 human proteome; n = ",
                length(universe_entrez)
            ),
            x = expression(
                -log[10](FDR)
            ),
            y = NULL
        ) +
        theme_bw(
            base_size = 10
        ) +
        theme(
            plot.title = element_text(
                face = "bold"
            ),
            axis.text.y = element_text(
                size = 8
            ),
            panel.grid.minor = element_blank()
        )

    ggsave(
        filename = output_path,
        plot = plot_object,
        width = 9,
        height = max(
            5,
            0.32 * nrow(plot_data) + 2
        ),
        units = "in",
        device = cairo_pdf
    )

    TRUE
}



# Run enrichment


enrichment_records <- list()
summary_records <- list()
failure_records <- list()

database_functions <- list(
    GO_BP = safe_enrich_go,
    Reactome = safe_enrich_reactome,
    KEGG = safe_enrich_kegg
)


for (gene_set_name in names(gene_sets_entrez)) {
    foreground_entrez <- gene_sets_entrez[[gene_set_name]]

    message(
        "Running enrichment for ",
        gene_set_name,
        " (",
        length(foreground_entrez),
        " mapped genes)"
    )

    for (database_name in names(database_functions)) {
        message(
            "  Database: ",
            database_name
        )

        enrichment_function <- database_functions[[database_name]]

        enrichment_object <- enrichment_function(
            genes = foreground_entrez,
            universe = universe_entrez,
            config = config
        )

        result_table <- standardize_enrichment_result(
            enrichment_object,
            significance_fdr = config$significance_fdr
        )

        result_table <- result_table |>
            mutate(
                gene_set = gene_set_name,
                database = database_name,
                .before = 1
            )

        result_key <- paste(
            gene_set_name,
            database_name,
            sep = "__"
        )

        enrichment_records[[result_key]] <- result_table

        table_path <- file.path(
            table_dir,
            paste0(
                gene_set_name,
                "_",
                database_name,
                "_enrichment.tsv"
            )
        )

        write_tsv(
            result_table,
            table_path
        )

        n_tested <- nrow(result_table)

        n_significant <- sum(
            result_table$significant_fdr_0_05,
            na.rm = TRUE
        )

        minimum_fdr <- if (n_tested > 0) {
            suppressWarnings(
                min(
                    result_table$p.adjust,
                    na.rm = TRUE
                )
            )
        } else {
            NA_real_
        }

        if (
            is.infinite(minimum_fdr) ||
            is.nan(minimum_fdr)
        ) {
            minimum_fdr <- NA_real_
        }

        summary_records[[length(summary_records) + 1]] <- tibble(
            gene_set = gene_set_name,
            database = database_name,
            input_symbols = length(
                gene_sets_symbol[[gene_set_name]]
            ),
            mapped_entrez_ids = length(
                foreground_entrez
            ),
            universe_symbols = length(
                universe_symbol
            ),
            universe_entrez_ids = length(
                universe_entrez
            ),
            tested_terms = n_tested,
            significant_terms_fdr_0_05 = n_significant,
            minimum_adjusted_p_value = minimum_fdr
        )

        dotplot_path <- file.path(
            figure_dir,
            paste0(
                gene_set_name,
                "_",
                database_name,
                "_dotplot.pdf"
            )
        )

        barplot_path <- file.path(
            figure_dir,
            paste0(
                gene_set_name,
                "_",
                database_name,
                "_barplot.pdf"
            )
        )

        dotplot_created <- tryCatch(
            make_dotplot(
                result_table = result_table,
                gene_set_name = gene_set_name,
                database_name = database_name,
                output_path = dotplot_path,
                show_categories = config$show_categories
            ),
            error = function(error) {
                warning(
                    paste0(
                        "Dotplot failed for ",
                        gene_set_name,
                        " / ",
                        database_name,
                        ": ",
                        conditionMessage(error)
                    )
                )

                FALSE
            }
        )

        barplot_created <- tryCatch(
            make_barplot(
                result_table = result_table,
                gene_set_name = gene_set_name,
                database_name = database_name,
                output_path = barplot_path,
                show_categories = config$show_categories
            ),
            error = function(error) {
                warning(
                    paste0(
                        "Barplot failed for ",
                        gene_set_name,
                        " / ",
                        database_name,
                        ": ",
                        conditionMessage(error)
                    )
                )

                FALSE
            }
        )

        if (!dotplot_created || !barplot_created) {
            failure_records[[length(failure_records) + 1]] <- tibble(
                gene_set = gene_set_name,
                database = database_name,
                dotplot_created = dotplot_created,
                barplot_created = barplot_created
            )
        }
    }
}



# Combined output tables


combined_results <- bind_rows(
    enrichment_records
)


combined_path <- file.path(
    table_dir,
    "all_pathway_enrichment_results.tsv"
)

write_tsv(
    combined_results,
    combined_path
)


significant_results <- combined_results |>
    filter(
        significant_fdr_0_05
    ) |>
    arrange(
        gene_set,
        database,
        p.adjust,
        desc(fold_enrichment)
    )


write_tsv(
    significant_results,
    file.path(
        table_dir,
        "significant_pathway_enrichment_results.tsv"
    )
)


summary_table <- bind_rows(
    summary_records
) |>
    arrange(
        gene_set,
        database
    )


summary_path <- file.path(
    output_dir,
    "pathway_enrichment_summary.tsv"
)

write_tsv(
    summary_table,
    summary_path
)


if (length(failure_records) > 0) {
    figure_status <- bind_rows(
        failure_records
    )
} else {
    figure_status <- tibble(
        gene_set = character(),
        database = character(),
        dotplot_created = logical(),
        barplot_created = logical()
    )
}


write_tsv(
    figure_status,
    file.path(
        output_dir,
        "figure_generation_status.tsv"
    )
)



# Report


report_lines <- c(
    "Pathway enrichment report",
    "=========================",
    "",
    paste0(
        "Gene-set directory: ",
        gene_set_dir
    ),
    paste0(
        "Output directory: ",
        output_dir
    ),
    "",
    "Statistical design:",
    paste0(
        "  Background universe symbols: ",
        length(universe_symbol)
    ),
    paste0(
        "  Background universe Entrez IDs: ",
        length(universe_entrez)
    ),
    "  Test: over-representation analysis",
    "  Multiple-testing correction: Benjamini-Hochberg",
    paste0(
        "  Significance threshold: adjusted p-value <= ",
        config$significance_fdr
    ),
    paste0(
        "  Minimum pathway size: ",
        config$minimum_gene_set_size
    ),
    paste0(
        "  Maximum pathway size: ",
        config$maximum_gene_set_size
    ),
    "",
    "Gene mapping:"
)


for (row_index in seq_len(nrow(mapping_summary))) {
    current_row <- mapping_summary[
        row_index,
    ]

    report_lines <- c(
        report_lines,
        paste0(
            "  ",
            current_row$gene_set,
            ": ",
            current_row$mapped_entrez_ids,
            "/",
            current_row$input_symbols,
            " mapped (",
            sprintf(
                "%.4f",
                current_row$mapping_fraction
            ),
            ")"
        )
    )
}


report_lines <- c(
    report_lines,
    "",
    "Enrichment summary:"
)


for (row_index in seq_len(nrow(summary_table))) {
    current_row <- summary_table[
        row_index,
    ]

    minimum_fdr_text <- if (
        is.na(
            current_row$minimum_adjusted_p_value
        )
    ) {
        "NA"
    } else {
        format(
            current_row$minimum_adjusted_p_value,
            scientific = TRUE,
            digits = 4
        )
    }

    report_lines <- c(
        report_lines,
        paste0(
            "  ",
            current_row$gene_set,
            " / ",
            current_row$database,
            ": tested terms = ",
            current_row$tested_terms,
            "; significant terms = ",
            current_row$significant_terms_fdr_0_05,
            "; minimum FDR = ",
            minimum_fdr_text
        )
    )
}


report_lines <- c(
    report_lines,
    "",
    "Interpretation policy:",
    paste0(
        "  The primary enrichment analysis is the candidate_panel ",
        "gene set against the measurable proteomic universe."
    ),
    paste0(
        "  The mir143_linked subset contains only seven genes and ",
        "must be interpreted as exploratory."
    ),
    paste0(
        "  The nonlinked_candidate analysis is descriptive and ",
        "does not represent a formal independent control group."
    ),
    paste0(
        "  Enrichment indicates over-representation, not activation, ",
        "causation, or pathway direction."
    ),
    paste0(
        "  Supportive miRNA-target evidence must not be described as ",
        "definitive functional regulation."
    ),
    paste0(
        "  KEGG retrieval may depend on internet availability at runtime."
    ),
    "",
    "Main output tables:",
    paste0(
        "  ",
        combined_path
    ),
    paste0(
        "  ",
        file.path(
            table_dir,
            "significant_pathway_enrichment_results.tsv"
        )
    ),
    paste0(
        "  ",
        summary_path
    ),
    ""
)


report_path <- file.path(
    output_dir,
    "pathway_enrichment_report.txt"
)


writeLines(
    report_lines,
    con = report_path,
    useBytes = TRUE
)



# Console summary


cat(
    "Pathway enrichment completed.\n"
)

cat(
    "Universe symbols:",
    length(universe_symbol),
    "\n"
)

cat(
    "Universe Entrez IDs:",
    length(universe_entrez),
    "\n"
)

for (gene_set_name in names(gene_sets_symbol)) {
    cat(
        gene_set_name,
        ":",
        length(gene_sets_symbol[[gene_set_name]]),
        "symbols,",
        length(gene_sets_entrez[[gene_set_name]]),
        "Entrez IDs\n"
    )
}

cat(
    "Wrote combined results:",
    combined_path,
    "\n"
)

cat(
    "Wrote summary:",
    summary_path,
    "\n"
)

cat(
    "Wrote report:",
    report_path,
    "\n"
)

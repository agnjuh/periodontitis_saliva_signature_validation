#!/usr/bin/env Rscript

# Query experimentally validated miRNA-target interactions using multiMiR.
#
# Outputs:
#   data/reference/multimir/multimir_validated_targets.tsv
#   data/reference/multimir/multimir_query_summary.tsv
#   data/reference/multimir/multimir_query_report.txt
#
# The script does not generate predicted interactions. It requests only
# experimentally validated records from databases accessible through multiMiR.

options(
  stringsAsFactors = FALSE,
  warn = 1
)

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  position <- match(flag, args)

  if (is.na(position)) {
    return(default)
  }

  if (position == length(args)) {
    stop(
      paste("Missing value after", flag),
      call. = FALSE
    )
  }

  args[[position + 1]]
}

mirna_id <- get_arg(
  "--mirna",
  "hsa-miR-143-3p"
)

output_dir <- get_arg(
  "--output-dir",
  "data/reference/multimir"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

required_packages <- c(
  "multiMiR"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required R package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nInstall with:\n",
      "if (!requireNamespace('BiocManager', quietly = TRUE)) ",
      "install.packages('BiocManager')\n",
      "BiocManager::install('multiMiR')"
    ),
    call. = FALSE
  )
}

normalise_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

find_column <- function(
  frame,
  candidates,
  required = FALSE
) {
  existing <- names(frame)
  existing_lower <- tolower(existing)

  candidates_lower <- tolower(candidates)

  exact_index <- match(
    candidates_lower,
    existing_lower,
    nomatch = 0
  )

  exact_index <- exact_index[
    exact_index > 0
  ]

  if (length(exact_index) > 0) {
    return(existing[[exact_index[[1]]]])
  }

  for (candidate in candidates_lower) {
    partial <- grep(
      candidate,
      existing_lower,
      fixed = TRUE
    )

    if (length(partial) > 0) {
      return(existing[[partial[[1]]]])
    }
  }

  if (required) {
    stop(
      paste0(
        "Could not identify required column. Candidates: ",
        paste(candidates, collapse = ", "),
        "\nAvailable columns: ",
        paste(existing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  NULL
}

column_or_empty <- function(
  frame,
  column_name
) {
  if (is.null(column_name)) {
    return(
      rep("", nrow(frame))
    )
  }

  normalise_text(
    frame[[column_name]]
  )
}

canonicalise_gene_symbol <- function(x) {
  x <- normalise_text(x)
  x <- sub(";.*$", "", x)
  x <- sub("\\|.*$", "", x)
  x <- gsub("\\s+", "", x)
  toupper(x)
}

canonicalise_mirna <- function(x) {
  x <- normalise_text(x)
  x <- gsub("_", "-", x)
  x <- gsub("\\s+", "", x)

  missing_prefix <- (
    nzchar(x) &
    !grepl("^[a-z]{3}-", x)
  )

  x[missing_prefix] <- paste0(
    "hsa-",
    x[missing_prefix]
  )

  x
}

split_values <- function(x) {
  x <- normalise_text(x)

  values <- unlist(
    strsplit(
      x,
      split = "[;,|/]+"
    )
  )

  values <- trimws(values)
  values <- values[nzchar(values)]

  unique(values)
}

collapse_unique <- function(x) {
  values <- unique(
    unlist(
      lapply(
        x,
        split_values
      )
    )
  )

  values <- values[nzchar(values)]

  paste(
    sort(values),
    collapse = ";"
  )
}

classify_evidence_strength <- function(
  database,
  experiment,
  support_type
) {
  combined <- tolower(
    paste(
      database,
      experiment,
      support_type,
      sep = ";"
    )
  )

  if (
    grepl(
      "strong",
      combined
    )
  ) {
    return("strong")
  }

  if (
    grepl(
      paste(
        c(
          "reporter assay",
          "luciferase",
          "western blot",
          "immunoblot",
          "qrt-pcr",
          "qpcr",
          "quantitative pcr"
        ),
        collapse = "|"
      ),
      combined
    )
  ) {
    return("strong")
  }

  if (
    grepl(
      paste(
        c(
          "clip",
          "clash",
          "microarray",
          "rna-seq",
          "sequencing",
          "proteomics",
          "psilac",
          "weak",
          "less strong"
        ),
        collapse = "|"
      ),
      combined
    )
  ) {
    return("supportive")
  }

  if (
    nzchar(experiment) ||
    nzchar(support_type)
  ) {
    return(
      "experimentally_supported_unspecified"
    )
  }

  "unclassified"
}

message(
  "Querying validated interactions for: ",
  mirna_id
)

query_result <- tryCatch(
  {
    multiMiR::get_multimir(
      mirna = mirna_id,
      table = "validated",
      summary = FALSE
    )
  },
  error = function(error) {
    stop(
      paste0(
        "multiMiR query failed: ",
        conditionMessage(error),
        "\nThe query requires working internet access and ",
        "access to the multiMiR database service."
      ),
      call. = FALSE
    )
  }
)

raw_data <- NULL

if (
  isS4(query_result) &&
  "data" %in% methods::slotNames(query_result)
) {
  raw_data <- methods::slot(
    query_result,
    "data"
  )
} else if (
  is.data.frame(query_result)
) {
  raw_data <- query_result
} else if (
  is.list(query_result) &&
  !is.null(query_result[["data"]])
) {
  raw_data <- query_result[["data"]]
}

if (
  !is.null(raw_data) &&
  !is.data.frame(raw_data)
) {
  raw_data <- tryCatch(
    as.data.frame(raw_data),
    error = function(error) NULL
  )
}

if (
  is.null(raw_data) ||
  !is.data.frame(raw_data)
) {
  stop(
    paste0(
      "The multiMiR result did not contain a readable data frame. ",
      "Result class: ",
      paste(class(query_result), collapse = ", "),
      "; available S4 slots: ",
      if (isS4(query_result)) {
        paste(methods::slotNames(query_result), collapse = ", ")
      } else {
        "not an S4 object"
      }
    ),
    call. = FALSE
  )
}

if (nrow(raw_data) == 0) {
  stop(
    paste0(
      "No validated interaction was returned for ",
      mirna_id,
      "."
    ),
    call. = FALSE
  )
}

mirna_column <- find_column(
  raw_data,
  c(
    "mature_mirna_id",
    "mirna",
    "mirna_id",
    "mature_mirna"
  ),
  required = TRUE
)

target_column <- find_column(
  raw_data,
  c(
    "target_symbol",
    "target_gene",
    "target",
    "gene_symbol"
  ),
  required = TRUE
)

database_column <- find_column(
  raw_data,
  c(
    "database",
    "source_database",
    "source"
  )
)

experiment_column <- find_column(
  raw_data,
  c(
    "experiment",
    "experiments",
    "experimental_method",
    "method"
  )
)

support_column <- find_column(
  raw_data,
  c(
    "support_type",
    "support",
    "evidence_type",
    "evidence"
  )
)

pubmed_column <- find_column(
  raw_data,
  c(
    "pubmed_id",
    "pmid",
    "references_pmid",
    "publication_id"
  )
)

target_entrez_column <- find_column(
  raw_data,
  c(
    "target_entrez",
    "target_entrez_id",
    "entrez_id"
  )
)

target_ensembl_column <- find_column(
  raw_data,
  c(
    "target_ensembl",
    "target_ensembl_id",
    "ensembl_id"
  )
)

standardised <- data.frame(
  mirna = canonicalise_mirna(
    column_or_empty(
      raw_data,
      mirna_column
    )
  ),
  target_gene = canonicalise_gene_symbol(
    column_or_empty(
      raw_data,
      target_column
    )
  ),
  target_entrez_id = column_or_empty(
    raw_data,
    target_entrez_column
  ),
  target_ensembl_id = column_or_empty(
    raw_data,
    target_ensembl_column
  ),
  evidence_database = column_or_empty(
    raw_data,
    database_column
  ),
  experimental_method = column_or_empty(
    raw_data,
    experiment_column
  ),
  support_type = column_or_empty(
    raw_data,
    support_column
  ),
  publication_id = column_or_empty(
    raw_data,
    pubmed_column
  ),
  stringsAsFactors = FALSE
)

standardised <- standardised[
  nzchar(standardised$mirna) &
  nzchar(standardised$target_gene),
  ,
  drop = FALSE
]

if (nrow(standardised) == 0) {
  stop(
    "No complete miRNA-target records remained after standardisation.",
    call. = FALSE
  )
}

standardised$evidence_strength <- mapply(
  classify_evidence_strength,
  database = standardised$evidence_database,
  experiment = standardised$experimental_method,
  support_type = standardised$support_type,
  USE.NAMES = FALSE
)

standardised$query_mirna <- mirna_id
standardised$query_type <- "validated"
standardised$query_source <- "multiMiR"

group_key <- paste(
  standardised$mirna,
  standardised$target_gene,
  sep = "\t"
)

group_indices <- split(
  seq_len(nrow(standardised)),
  group_key
)

strength_rank <- c(
  strong = 4,
  experimentally_supported_unspecified = 3,
  supportive = 2,
  unclassified = 1
)

aggregate_group <- function(indices) {
  subset <- standardised[
    indices,
    ,
    drop = FALSE
  ]

  strengths <- unique(
    subset$evidence_strength
  )

  strengths <- strengths[
    nzchar(strengths)
  ]

  if (length(strengths) == 0) {
    strongest <- "unclassified"
  } else {
    ranks <- strength_rank[
      strengths
    ]

    ranks[is.na(ranks)] <- 0

    strongest <- strengths[
      which.max(ranks)
    ]
  }

  data.frame(
    interaction_record_id = paste0(
      "multimir:",
      subset$mirna[[1]],
      ":",
      subset$target_gene[[1]]
    ),
    mirna = subset$mirna[[1]],
    target_gene = subset$target_gene[[1]],
    target_entrez_id = collapse_unique(
      subset$target_entrez_id
    ),
    target_ensembl_id = collapse_unique(
      subset$target_ensembl_id
    ),
    evidence_databases = collapse_unique(
      subset$evidence_database
    ),
    maximum_evidence_strength = strongest,
    evidence_strengths = collapse_unique(
      subset$evidence_strength
    ),
    experimental_methods = collapse_unique(
      subset$experimental_method
    ),
    support_types = collapse_unique(
      subset$support_type
    ),
    publication_ids = collapse_unique(
      subset$publication_id
    ),
    n_source_records = nrow(subset),
    n_databases = length(
      split_values(
        collapse_unique(
          subset$evidence_database
        )
      )
    ),
    n_publications = length(
      split_values(
        collapse_unique(
          subset$publication_id
        )
      )
    ),
    n_experimental_methods = length(
      split_values(
        collapse_unique(
          subset$experimental_method
        )
      )
    ),
    query_mirna = mirna_id,
    query_type = "validated",
    query_source = "multiMiR",
    stringsAsFactors = FALSE
  )
}

aggregated <- do.call(
  rbind,
  lapply(
    group_indices,
    aggregate_group
  )
)

rownames(aggregated) <- NULL

order_rank <- strength_rank[
  aggregated$maximum_evidence_strength
]

order_rank[is.na(order_rank)] <- 0

aggregated <- aggregated[
  order(
    -order_rank,
    -aggregated$n_databases,
    -aggregated$n_publications,
    -aggregated$n_experimental_methods,
    aggregated$target_gene
  ),
  ,
  drop = FALSE
]

output_targets <- file.path(
  output_dir,
  "multimir_validated_targets.tsv"
)

output_raw <- file.path(
  output_dir,
  "multimir_standardised_source_records.tsv"
)

output_summary <- file.path(
  output_dir,
  "multimir_query_summary.tsv"
)

output_report <- file.path(
  output_dir,
  "multimir_query_report.txt"
)

write.table(
  standardised,
  file = output_raw,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

write.table(
  aggregated,
  file = output_targets,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

summary_table <- data.frame(
  metric = c(
    "query_mirna",
    "raw_source_records",
    "standardised_source_records",
    "unique_validated_interactions",
    "unique_target_genes",
    "represented_databases",
    "strong_interactions",
    "supportive_interactions",
    "experimentally_supported_unspecified_interactions",
    "unclassified_interactions"
  ),
  value = c(
    mirna_id,
    nrow(raw_data),
    nrow(standardised),
    nrow(aggregated),
    length(unique(aggregated$target_gene)),
    length(
      split_values(
        collapse_unique(
          aggregated$evidence_databases
        )
      )
    ),
    sum(
      aggregated$maximum_evidence_strength ==
        "strong"
    ),
    sum(
      aggregated$maximum_evidence_strength ==
        "supportive"
    ),
    sum(
      aggregated$maximum_evidence_strength ==
        "experimentally_supported_unspecified"
    ),
    sum(
      aggregated$maximum_evidence_strength ==
        "unclassified"
    )
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_table,
  file = output_summary,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

package_version <- as.character(
  utils::packageVersion("multiMiR")
)

report_lines <- c(
  "multiMiR validated-target query",
  "===============================",
  "",
  paste0(
    "Query miRNA: ",
    mirna_id
  ),
  "Query table: validated",
  paste0(
    "multiMiR package version: ",
    package_version
  ),
  paste0(
    "Raw source records: ",
    nrow(raw_data)
  ),
  paste0(
    "Standardised source records: ",
    nrow(standardised)
  ),
  paste0(
    "Unique validated interactions: ",
    nrow(aggregated)
  ),
  paste0(
    "Unique target genes: ",
    length(
      unique(
        aggregated$target_gene
      )
    )
  ),
  "",
  "Evidence policy:",
  "  - Only records returned from the multiMiR validated table are retained.",
  "  - Predicted interactions are not queried.",
  "  - Database, experiment, support-type, and publication fields are retained.",
  "  - Repeated miRNA-target records are aggregated without discarding source evidence.",
  "  - Evidence-strength labels are conservative derived annotations.",
  "",
  "Output files:",
  paste0(
    "  ",
    output_raw
  ),
  paste0(
    "  ",
    output_targets
  ),
  paste0(
    "  ",
    output_summary
  ),
  ""
)

writeLines(
  report_lines,
  con = output_report,
  useBytes = TRUE
)

message(
  "Wrote standardised source records: ",
  output_raw
)

message(
  "Wrote validated targets: ",
  output_targets
)

message(
  "Wrote query summary: ",
  output_summary
)

message(
  "Wrote query report: ",
  output_report
)

message(
  "Unique validated interactions: ",
  nrow(aggregated)
)

message(
  "Unique target genes: ",
  length(
    unique(
      aggregated$target_gene
    )
  )
)

#!/usr/bin/env Rscript
# Perform and visualize PCA on modkit pileup output -----------------------
#
# Example:
#   Rscript modkit_pca.R \
#     --bed_dir results/modkit/ \
#     --metadata sample_metadata.tsv \
#     --output_prefix modkit_pca
#

# Load libraries ----------------------------------------------------------

library(argparse)
library(tidyverse)
library(janitor)
library(glue)
library(magrittr)


# Modkit schema -----------------------------------------------------------

modkit_col_names <- c(
  "chrom", "start", "end", "name", "score", "strand", "thick_start", "thick_end", "color",
  "valid_coverage", "percent_modified", "count_modified", "count_canonical", "count_other_mod",
  "count_delete", "count_fail", "count_diff", "count_nocall"
)


# Parse arguments ---------------------------------------------------------

#' Create command-line argument parser.
#'
#' Defines CLI interface for modkit PCA pipeline including filtering, input handling,
#' and output configuration.
#'
#' @return argparse parser object.
create_parser <- function() {
  parser <- argparse::ArgumentParser(
    description = "PCA on modkit pileup with filtering and ggplot2 visualization."
  )

  parser$add_argument("--bed_files", nargs = "+", default = NULL,
    help = "Input BED or BED.GZ files"
  )

  parser$add_argument("--bed_dir", default = NULL,
    help = "Directory with modkit BED files"
  )

  parser$add_argument("--pattern", default = "\\.bed(\\.gz)?$",
    help = "Regex pattern for file discovery"
  )

  parser$add_argument("--metadata", required = TRUE,
    help = "CSV metadata file (comma-separated)"
  )

  parser$add_argument("--mod_column", default = "percent_modified",
    help = "Methylation fraction column"
  )

  parser$add_argument("--min_coverage", default = 10, type = "double",
    help = "Minimum valid coverage filter"
  )

  parser$add_argument("--min_methylation", default = 0.3, type = "double",
    help = "Minimum methylation fraction filter"
  )

  parser$add_argument("--output_dir", default = ".",
    help = "Output directory"
  )

  parser$add_argument("--output_prefix", default = "modkit_pca",
    help = "Output prefix"
  )

  parser
}


# Read metadata -----------------------------------------------------------

#' Read metadata table.
#'
#' Loads sample metadata from a comma-separated CSV file and cleans column names.
#'
#' @param metadata_file Path to metadata CSV file.
#'
#' @return Tibble containing cleaned metadata.
read_metadata <- function(metadata_file) {
  readr::read_csv(metadata_file, show_col_types = FALSE) %>%
    janitor::clean_names()
}


# Discover BED files ------------------------------------------------------

#' Discover input BED files.
#'
#' Returns a vector of file paths either from explicit input or a directory scan.
#'
#' @param bed_files Optional character vector of input BED files.
#' @param bed_dir Optional directory containing BED files.
#' @param pattern Regex pattern used to match BED files.
#'
#' @return Character vector of BED file paths.
discover_bed_files <- function(
    bed_files = NULL,
    bed_dir = NULL,
    pattern = "\\.bed(\\.gz)?$"
) {
  if (!is.null(bed_files)) return(bed_files)

  if (!is.null(bed_dir)) {
    files <- list.files(path = bed_dir, pattern = pattern, full.names = TRUE)
    return(files[!stringr::str_detect(files, "\\.gzi$|\\.tbi$|\\.csi$")])
  }

  stop("Provide bed_files or bed_dir", call. = FALSE)
}


# Read modkit file --------------------------------------------------------

#' Read a modkit pileup BED file.
#'
#' Loads modkit pileup data, applies coverage and methylation filters, and
#' converts to long format suitable for PCA.
#'
#' @param bed_file Path to BED or BED.GZ file.
#' @param mod_column Column name for methylation fraction.
#' @param min_coverage Minimum read coverage threshold.
#' @param min_methylation Minimum methylation fraction threshold.
#'
#' @return Tibble with columns: sample, locus, value.
read_modkit_bed <- function(
    bed_file,
    mod_column,
    min_coverage,
    min_methylation
) {
  sample_name <- basename(bed_file) %>%
    stringr::str_remove("_UMCUAZ......\\.pileup\\.bed(\\.gz)?$")

  bed_tbl <- readr::read_tsv(
    bed_file,
    comment = "#",
    col_names = FALSE,
    show_col_types = FALSE,
    progress = FALSE
  )

  colnames(bed_tbl) <- modkit_col_names[seq_len(ncol(bed_tbl))]
  bed_tbl %>%
    janitor::clean_names() %>%
    dplyr::filter(
      valid_coverage >= min_coverage,
      .data[[mod_column]] >= min_methylation
    ) %>%
    dplyr::mutate(
      locus = glue::glue("{chrom}_{start}_{end}"),
      sample = sample_name,
      value = .data[[mod_column]]
    ) %>%
    dplyr::select(sample, locus, value)
}


#' Read multiple modkit BED files.
#'
#' Iterates over input files and combines all filtered modkit pileup data.
#'
#' @param bed_files Character vector of BED file paths.
#' @param mod_column Methylation column name.
#' @param min_coverage Minimum coverage threshold.
#' @param min_methylation Minimum methylation threshold.
#'
#' @return Combined tibble of all samples in long format.
read_all_modkit_beds <- function(
    bed_files,
    mod_column,
    min_coverage,
    min_methylation
) {
  purrr::map_dfr(
    bed_files,
    read_modkit_bed,
    mod_column = mod_column,
    min_coverage = min_coverage,
    min_methylation = min_methylation
  )
}


# Build matrix ------------------------------------------------------------

#' Construct PCA input matrix.
#'
#' Converts long-format methylation data into a sample-by-site matrix and
#' fills missing values with zero.
#'
#' @param methylation_tbl Long-format tibble.
#'
#' @return Numeric matrix with samples as rows and loci as columns.
prepare_pca_matrix <- function(methylation_tbl) {
  wide_tbl <- methylation_tbl %>%
    tidyr::pivot_wider(names_from = locus, values_from = value, values_fill = 0)

  sample_names <- wide_tbl$sample

  mat <- wide_tbl %>%
    dplyr::select(-sample) %>%
    as.matrix()

  rownames(mat) <- sample_names
  mat
}


#' Remove zero-variance features.
#'
#' Filters loci (columns) that have no variation across samples.
#'
#' @param mat Numeric matrix (samples x loci).
#'
#' @return Filtered matrix with variable features only.
remove_zero_variance_features <- function(mat) {
  v <- apply(mat, 2, stats::var)
  mat[, v > 0, drop = FALSE]
}


# PCA ---------------------------------------------------------------------

#' Run PCA on methylation matrix.
#'
#' Performs principal component analysis using prcomp with centering and scaling.
#'
#' @param mat Numeric matrix (samples x loci).
#'
#' @return prcomp object containing PCA results.
run_pca <- function(mat) {
  stats::prcomp(mat, center = TRUE, scale. = TRUE)
}


#' Compute variance explained by PCA components.
#'
#' Uses eigenvalues derived from PCA singular values to quantify variance
#' explained per component.
#'
#' @param pca_result prcomp object.
#'
#' @return Tibble with variance explained per principal component.
compute_variance_explained <- function(pca_result) {
  eig <- pca_result$sdev^2

  tibble::tibble(
    pc = paste0("PC", seq_along(eig)),
    eigenvalue = eig,
    variance_explained = eig / sum(eig),
    cumulative_variance = cumsum(eig / sum(eig))
  )
}


# Metadata join -----------------------------------------------------------

#' Match metadata to sample names.
#'
#' Performs partial string matching between sample names and metadata table.
#'
#' @param sample_names Character vector of sample names.
#' @param metadata_tbl Metadata tibble containing partial_name column.
#'
#' @return Metadata tibble aligned to samples.
match_metadata <- function(sample_names, metadata_tbl) {
  purrr::map_dfr(sample_names, function(s) {
    m <- metadata_tbl %>%
      dplyr::filter(stringr::str_detect(s, partial_name))

    if (nrow(m) == 0) stop(paste("No metadata match:", s), call. = FALSE)

    m %>% dplyr::slice(1) %>% dplyr::mutate(sample = s)
  })
}


#' Create PCA plot input table.
#'
#' Extracts PCA scores and joins metadata annotations.
#'
#' @param pca_result prcomp object.
#' @param metadata_tbl Metadata tibble.
#'
#' @return Tibble for plotting.
create_plot_tbl <- function(pca_result, metadata_tbl) {
  scores <- tibble::as_tibble(pca_result$x, rownames = "sample")

  meta <- match_metadata(scores$sample, metadata_tbl)

  scores %>% dplyr::left_join(meta, by = "sample")
}


# Plot + save -------------------------------------------------------------

#' Create and save PCA plot.
#'
#' Generates PCA scatter plot and writes PDF to disk.
#'
#' @param plot_tbl Tibble containing PCA coordinates and metadata.
#' @param output_dir Output directory path.
#' @param output_prefix Output filename prefix.
#'
#' @return ggplot object (also saved to file).
plot_pca_and_save <- function(plot_tbl, output_dir, output_prefix, subtitle) {
  out_file <- file.path(output_dir, paste0(output_prefix, "_pca.pdf"))
  p1_p2 <- ggplot2::ggplot(
    plot_tbl,
    ggplot2::aes(x = PC1, y = PC2, color = mutation, shape = stress_condition, label=sample)
  ) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_text(vjust=-1,  nudge_x = 0.5, nudge_y = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(color = "Mutation", shape = "Stress condition (with or without)", subtitle=subtitle) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  p2_p3 <- ggplot2::ggplot(
    plot_tbl,
    ggplot2::aes(x = PC2, y = PC3, color = mutation, shape = stress_condition, label=sample)
  ) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_text(vjust=-1,  nudge_x = 0.5, nudge_y = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs( color = "Mutation", shape = "Stress condition (with or without)", subtitle=subtitle) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  p3_p4 <- ggplot2::ggplot(
    plot_tbl,
    ggplot2::aes(x = PC3, y = PC4, color = mutation, shape = stress_condition, label=sample)
  ) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_text(vjust=-1,  nudge_x = 0.5, nudge_y = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(color = "Mutation", shape = "Stress condition (with or without)", subtitle=subtitle) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    ggplot2::ggsave(out_file, gridExtra::marrangeGrob(grobs = list(p1_p2, p2_p3, p3_p4), nrow=2, ncol=2), width = 24, height = 12)

}


# Main --------------------------------------------------------------------

#' Run full modkit PCA pipeline.
#'
#' Executes file discovery, filtering, PCA computation, and visualization.
#'
#' @param args Named list of pipeline arguments (from CLI or interactive mode).
main <- function(args) {
  bed_files <- discover_bed_files(args$bed_files, args$bed_dir, args$pattern)

  meta <- read_metadata(args$metadata)

  mat_tbl <- read_all_modkit_beds(
    bed_files,
    args$mod_column,
    args$min_coverage,
    args$min_methylation
  )

  mat <- prepare_pca_matrix(mat_tbl)
  mat <- remove_zero_variance_features(mat)

  pca <- run_pca(mat)

  plot_tbl <- create_plot_tbl(pca, meta)

  subtitle = glue::glue("valid_coverage >= {args$min_coverage} and percentage modified >= {args$min_methylation}%.")
  plot_pca_and_save(plot_tbl, args$output_dir, "modkit_10_30" , subtitle)  #args$output_prefix

  variance_tbl <- compute_variance_explained(pca)

  readr::write_tsv(
    variance_tbl,
    file.path(args$output_dir, paste0(args$output_prefix, "_variance.tsv"))
  )
}


# Run ---------------------------------------------------------------------

if (interactive()) {
  args <- list(
    bed_files = NULL,
    bed_dir = "",
    metadata = "metadata.csv",
    mod_column = "percent_modified",
    min_coverage = 10,
    min_methylation = 30,
    output_dir = ".",
    output_prefix = "modkit_test"
  )
} else {
  args <- create_parser()$parse_args()
}

main(args)

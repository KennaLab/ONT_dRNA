
# 1. Load libraries ------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(GenomicFeatures)
library(ggplot2)

source("utils.R")


# 2. Input helpers -------------------------------------------------------------

#' Read a pileup BED file and append sample name
#'
#' @description
#' Reads a ModKit pileup BED file and extracts the sample name from the file
#' name (string before the first dot).
#'
#' @param file_path Character. Path to `.pileup.bed` file.
#'
#' @return Tibble with BED columns and `sample_name`.
read_pileup_bed <- function(file_path) {
  bed_dt <- fread(
    file_path,
    sep = "\t",
    header = FALSE,
    col.names = c(
      "chrom", "start", "end", "name", "score", "strand",
      "thick_start", "thick_end", "color", "valid_coverage",
      "percent_modified", "count_modified", "count_canonical",
      "count_other_mod", "count_delete", "count_fail",
      "count_diff", "count_nocall"
    )
  )

  bed_dt %>%
    mutate(
      sample_name = str_extract(basename(file_path), "^[^.]+")
    ) %>%
    as_tibble()
}


# 3. Transcript region model ---------------------------------------------------

#' Build transcript regions with lengths
#'
#' @description
#' Extracts 5'UTR, CDS, and 3'UTR exon structures grouped by transcript and
#' computes transcript lengths per region. Lengths are stored alongside the
#' ranges to avoid mismatches.
#'
#' @param gtf_path Character. Path to GTF file.
#'
#' @return Named list with elements `utr5`, `cds`, `utr3`. Each contains:
#'   - ranges: GRangesList
#'   - lengths: Named numeric vector
build_tx_regions <- function(gtf_path) {
  txdb <- makeTxDbFromGFF(gtf_path, format = "gtf")
  region_ranges <- list(
    utr5 = fiveUTRsByTranscript(txdb, use.names = TRUE),
    cds  = cdsBy(txdb, by = "tx", use.names = TRUE),
    utr3 = threeUTRsByTranscript(txdb, use.names = TRUE)
  )

  # Compute lengths per transcript within each region
  region_model <- map(region_ranges, function(grl) {
    tx_lengths <- sum(width(grl))
    tx_lengths <- setNames(tx_lengths, names(grl))
    list(
      ranges = grl,
      lengths = tx_lengths
    )
  })

  return(region_model)
}

# 4. Region mapping ------------------------------------------------------------

#' Map m6A sites to transcript region coordinates
#'
#' @description
#' Maps genomic m6A sites to transcript-relative coordinates within a specific
#' region (5'UTR, CDS, 3'UTR).
#'
#' @param bed_df Tibble. Input BED data.
#' @param region_grl GRangesList. Region exon structure.
#' @param region_lengths Named numeric vector of transcript lengths.
#' @param region_name Character. Region label.
#'
#' @return data.table with transcript-relative positions.
map_to_region <- function(
  bed_df,
  region_grl,
  region_lengths,
  region_name
) {

  # Convert BED to GRanges
  bed_gr <- GRanges(
    seqnames = bed_df$chrom,
    ranges = IRanges(start = bed_df$start + 1, end = bed_df$end),
    strand = bed_df$strand
  )
  mcols(bed_gr)$row_id <- seq_len(nrow(bed_df))

  # Flatten exon structure
  exon_gr <- unlist(region_grl)
  tx_ids <- names(exon_gr)

  # Find overlaps
  hits <- findOverlaps(bed_gr, exon_gr, ignore.strand = FALSE)
  query_idx <- queryHits(hits)
  subject_idx <- subjectHits(hits)

  # Position within exon
  exon_start <- start(exon_gr)[subject_idx]
  genomic_pos <- start(bed_gr)[query_idx]
  exon_offset <- genomic_pos - exon_start + 1

  # Cumulative exon offsets per transcript
  exon_widths <- width(exon_gr)
  exon_cumsum <- ave(exon_widths, factor(tx_ids), FUN = cumsum)
  tx_offset <- exon_cumsum[subject_idx] - exon_widths[subject_idx]
  tx_position <- tx_offset + exon_offset

  # Build result table
  result_dt <- data.table(
    sample_name = bed_df$sample_name[query_idx],
    percent_modified = bed_df$percent_modified[query_idx],
    tx_id = tx_ids[subject_idx],
    tx_pos = tx_position
  )[
    , tx_length := region_lengths[tx_id]
  ][
    !is.na(tx_length) & tx_length > 0
  ][
    , rel_pos := tx_pos / tx_length
  ][
    , region := region_name
  ][
    rel_pos >= 0 & rel_pos <= 1
  ]

  return(result_dt)
}


# 5. Metagene construction -----------------------------------------------------

#' Combine regions into a unified metagene axis
#'
#' @param dt data.table with region-relative positions.
#'
#' @return data.table with global metagene position.
combine_metagene_axis <- function(dt) {
  dt[
    , global_pos := case_when(
      region == "utr5" ~ rel_pos * (1 / 3),
      region == "cds"  ~ (1 / 3) + rel_pos * (1 / 3),
      region == "utr3" ~ (2 / 3) + rel_pos * (1 / 3)
    )
  ]

  return(dt)
}


# 6. Binning and density -------------------------------------------------------


#' Bin metagene positions
#'
#' @param dt data.table with global positions.
#' @param n_bins Integer. Number of bins.
#'
#' @return data.table with density per bin.
bin_metagene <- function(dt, n_bins = 150) {

  dt[
    , bin := floor(global_pos * n_bins)
  ][
    , bin := pmin(bin, n_bins - 1)
  ][
    , .(density = mean(percent_modified, na.rm = TRUE)),
    by = .(sample_name, bin)
  ]
}

# 7. Plotting ------------------------------------------------------------------

#' Plot HOMER-style metagene density
#'
#' @param density_dt data.table with densities.
#' @param n_bins Integer. Number of bins.
#'
#' @return ggplot object.
plot_metagene_homer <- function(density_dt, n_bins = 150) {

  density_dt %>%
    as_tibble() %>%
    mutate(position = bin / n_bins) %>%
    ggplot(aes(x = position, y = density, color = sample_name)) +
    geom_line(linewidth = 1.2) +
    geom_vline(xintercept = c(1/3, 2/3),linetype = "dashed", color = "black") +
    annotate("text", x = 0.16, y = Inf, label = "5'UTR", vjust = 2) +
    annotate("text", x = 0.5,  y = Inf, label = "CDS",  vjust = 2) +
    annotate("text", x = 0.83, y = Inf, label = "3'UTR", vjust = 2) +
    theme_minimal() +
    labs(
      title = "m6A Metagene Profile (HOMER-like)",
      x = "Transcript position",
      y = "Mean Fraction of Modified Sites",
      subtitle = paste0("Density computed per ", str(n_bins), " bins along the metagene.")
    ) +
    scale_color_viridis_d()
}

#' Plot m6A Density Along 3'UTR (Start → Middle → End)
#'
#' @description
#' Creates a bar plot showing the mean fraction of modified sites along the 3'UTR,
#' divided into start, middle, and end regions. Uses tidyverse-style pipes and
#' dplyr syntax exclusively, maintaining a tidy workflow.
#'
#' @param mapped_dt data.table. Output of `combine_metagene_axis()` with columns
#'   `region`, `rel_pos`, `percent_modified`, `sample_name`.
#'
#' @return ggplot object. Bar plot of 3'UTR density with start, middle, and end bins.
#'
#' @examples
#' utr3_plot <- plot_utr3_zoom_tidy(mapped_dt)
#' print(utr3_plot)
plot_utr3_zoom_tidy <- function(mapped_dt) {

  # Filter for 3'UTR and compute start/middle/end bins and density
  utr3_density <- mapped_dt %>%
    filter(region == "utr3") %>%
    mutate(
      utr3_bin = cut(
        rel_pos,
        breaks = c(0, 1/3, 2/3, 1),
        labels = c("Start 3'UTR", "Middle 3'UTR", "End 3'UTR"),
        include.lowest = TRUE
      )
    ) %>%
    group_by(sample_name, utr3_bin) %>%
    summarise(density = mean(percent_modified, na.rm = TRUE), .groups = "drop")

  # Generate bar plot
  utr3_density %>%
    ggplot(aes(x = utr3_bin, y = density, fill = sample_name)) +
    geom_col(position = "dodge") +
    theme_minimal() +
    labs(
      title = "m6A Density Along 3'UTR",
      x = "3'UTR Region",
      y = "Mean Fraction of Modified Sites"
    ) +
    scale_fill_viridis_d()
}


#' Plot Continuous 3'UTR Density with Finer Bins
#'
#' @description
#' Creates a line plot showing the mean fraction of modified sites along the 3'UTR.
#' The 3'UTR is divided into a user-specified number of bins (default 20) for smooth
#' visualization. Adds annotations for start, middle, and end regions of the 3'UTR.
#'
#' @param mapped_dt data.frame or tibble. Output of `combine_metagene_axis()` with
#'   columns `region`, `rel_pos`, `percent_modified`, `sample_name`.
#' @param n_bins Integer. Number of bins to divide the 3'UTR for smoother line plot. Default 20.
#'
#' @return ggplot object. Continuous line plot of 3'UTR density with annotated regions.
#'
#' @examples
#' utr3_smooth_plot <- plot_utr3_smooth(mapped_dt, n_bins = 20)
#' print(utr3_smooth_plot)
plot_utr3_smooth <- function(mapped_dt, n_bins = 20) {

  # Filter for 3'UTR and bin relative positions
  utr3_density <- mapped_dt %>%
    filter(region == "utr3") %>%
    mutate(
      bin = cut(
        rel_pos,
        breaks = seq(0, 1, length.out = n_bins + 1),
        include.lowest = TRUE,
        labels = FALSE
      )
    ) %>%
    group_by(sample_name, bin) %>%
    summarise(
      density = mean(percent_modified, na.rm = TRUE),
      bin_center = mean(rel_pos, na.rm = TRUE),
      .groups = "drop"
    )

  # Create line plot with annotated start/middle/end
  utr3_density %>%
    ggplot(aes(x = bin_center, y = density, color = sample_name)) +
    geom_line(size = 1.2) +
    geom_vline(xintercept = c(0, 0.5, 1), linetype = "dashed", color = "black") +
    annotate("text", x = 0,   y = Inf, label = "Start 3'UTR", vjust = 2, hjust = 0) +
    annotate("text", x = 0.5, y = Inf, label = "Middle 3'UTR", vjust = 2, hjust = 0.5) +
    annotate("text", x = 1,   y = Inf, label = "End 3'UTR", vjust = 2, hjust = 1) +
    theme_minimal() +
    labs(
      title = "Continuous m6A Density Along 3'UTR",
      x = "Normalized 3'UTR Position",
      y = "Mean Fraction of Modified Sites"
    ) +
    scale_color_viridis_d()
}

# 8. Pipeline ------------------------------------------------------------------

bed_dir  <- "/Users/ejong19/Documents/Werk/Kenna/ont_drna/data/modkit/pileup/"
plot_dir <- "/Users/ejong19/Documents/Werk/Kenna/ont_drna/data/plots/"
gtf_file <- "/Users/ejong19/Downloads/gencode.v49.annotation.gtf"

# 8.1 Load BED files
bed_files <- list.files(
  bed_dir,
  pattern = "\\.pileup\\.bed$",
  full.names = TRUE
)
bed_data <- map_dfr(bed_files, read_pileup_bed)
bed_data_valid_cov <- filter(bed_data, valid_coverage > 10)

# 8.2 Build transcript regions (with lengths embedded)
tx_regions <- build_tx_regions(gtf_file)

# 8.3 Map to regions
mapped_dt <- rbindlist(list(
  map_to_region(bed_data_valid_cov, tx_regions$utr5$ranges, tx_regions$utr5$lengths, "utr5"),
  map_to_region(bed_data_valid_cov, tx_regions$cds$ranges,  tx_regions$cds$lengths,  "cds"),
  map_to_region(bed_data_valid_cov, tx_regions$utr3$ranges, tx_regions$utr3$lengths, "utr3")
))

# 8.4 Build metagene axis
mapped_dt <- combine_metagene_axis(mapped_dt)

# 8.5 Bin + density
density_dt <- bin_metagene(mapped_dt, n_bins = 150)

# 8.6 Plot
metagene_plot <- plot_metagene_homer(density_dt)
utr3_plot <- plot_utr3_zoom_tidy(mapped_dt)
utr3_smooth_plot_50 <- plot_utr3_smooth(mapped_dt, n_bins = 50)

save_all_plots_to_html(
  plots = list(metagene_plot = metagene_plot, utr3_plot= utr3_plot, utr3_smooth_plot_50=utr3_smooth_plot_50),
  output_file = file.path(plot_dir, "m6a_metagene_and_utr_plots.html")
)

#!/usr/bin/env Rscript

###########################################################
# Script: visualize.R
# Author: Ellen van de Geer
# Date: 2026-03-17
# Purpose: Visualize dRNA ONT m6A data from methylartist/modkit DMR
###########################################################

# 0. Install required packages ----
cran_packages <- c("tidyverse", "ggrepel", "gt", "forcats", "Cairo")
bioc_packages <- c("DSS")

installed_cran <- cran_packages %in% rownames(installed.packages())
if (any(!installed_cran)) {
  install.packages(cran_packages[!installed_cran])
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
installed_bioc <- bioc_packages %in% rownames(installed.packages())
if (any(!installed_bioc)) {
  BiocManager::install(bioc_packages[!installed_bioc])
}

# 1. Load libraries ----
library(tidyverse)
library(ggrepel)
library(gt)

# 2. Helper functions ----

#' Pivot segment methylation data to long format
#'
#' @param data Dataframe of segment methylation
#' @return Long-format dataframe
pivot_data_to_long <- function(data) {
  data %>%
    dplyr::select(-c("seg_chrom", "seg_start", "seg_end", "seg_strand")) %>%
    tidyr::pivot_longer(
      cols = -c("seg_id", "seg_name"),
      names_to = c("sample", "suffix"),
      names_sep = ".pileup.bed_",
      values_to = "value"
    )
}

#' Prepare percentage ratio data per gene type
#'
#' @param data_long Long-format dataframe
#' @return Dataframe with percent and gene_type_factor
prepare_data_ratio <- function(data_long) {
  data_long %>%
    dplyr::filter(suffix %in% c("a_meth_calls", "a_unmeth_calls")) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(total_calls_per_sample = sum(value, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sample, seg_name, suffix) %>%
    dplyr::mutate(
      percent = (sum(value, na.rm = TRUE) / total_calls_per_sample) * 100
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(suffix), desc(percent), seg_name) %>%
    dplyr::mutate(gene_type_factor = forcats::fct_reorder(seg_name, percent))
}

#' Generate ratio plots
#'
#' @param df Dataframe from prepare_data_ratio
#' @return Named list of ggplot objects
generate_ratio_plots <- function(df) {
  base_plot <- function(df_subset, title_sub) {
    ggplot2::ggplot(
      df_subset,
      ggplot2::aes(x = gene_type_factor, y = percent + 0.5, fill = suffix)
    ) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::facet_wrap(~sample) +
      ggplot2::geom_text(
        ggplot2::aes(label = round(percent, 2)),
        position = ggplot2::position_dodge(width = 0.9),
        vjust = 0,
        size = 2.75
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      ) +
      ggplot2::labs(
        title = "Ratio of un- and methylated calls of m6A",
        subtitle = title_sub,
        x = "Gene type",
        y = "Percentage of calls"
      )
  }

  list(
    ratio_calls_m6a_plot = base_plot(df, "Per gene type and sample"),
    ratio_calls_m6a_plot_zoom20 = base_plot(df, "Per gene type and sample") +
      ggplot2::coord_cartesian(ylim = c(0, 20)),
    ratio_calls_m6a_plot_exclude_0.05p = base_plot(
      dplyr::filter(df, percent > 0.05),
      "Excluded gene types with <= 0.05% of calls"
    ),
    ratio_calls_m6a_plot_exclude_1p = base_plot(
      dplyr::filter(df, percent > 1),
      "Excluded gene types with <= 1% of calls"
    )
  )
}

#' Read BED files into named list of tibbles (memory-efficient, exclude unnecessary columns)
#'
#' @param path Directory path
#' @param strict If TRUE, abort on header mismatch
#' @return Named list of tibbles
read_bed_files <- function(path, strict = TRUE) {
  exclude_cols <- c(
    "name",
    "strand",
    "a_mod_percentages",
    "b_mod_percentages",
    "cohen_h",
    "cohen_h_low",
    "cohen_h_high"
  )
  files <- list.files(path, pattern = "\\.bed$", full.names = TRUE)

  setNames(
    lapply(files, function(f) {
      if (!file.exists(f)) {
        stop("File not found: ", f)
      }
      first_line <- readr::read_lines(f, n_max = 1)
      header <- if (stringr::str_starts(first_line, "#")) {
        strsplit(sub("^#", "", first_line), "\t")[[1]]
      } else {
        NULL
      }
      keep <- if (!is.null(header)) setdiff(header, exclude_cols) else NULL
      readr::read_tsv(f, col_names = header, comment = "#", col_select = keep)
    }),
    tools::file_path_sans_ext(basename(files))
  )
}

#' Preload GENCODE GTF once
#'
#' Loads only "gene" features and returns a GRanges object ready for annotation.
#'
#' @param gtf_file Path to GENCODE GTF file (.gtf or .gtf.gz)
#' @return GRanges object with gene_id, gene_name, gene_type metadata
preload_gencode_gtf <- function(gtf_file) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    install.packages("GenomicRanges")
  }
  if (!requireNamespace("IRanges", quietly = TRUE)) {
    install.packages("IRanges")
  }

  gtf <- readr::read_tsv(
    gtf_file,
    comment = "#",
    progress = FALSE,
    col_names = c(
      "seqname",
      "source",
      "feature",
      "start",
      "end",
      "score",
      "strand",
      "frame",
      "attribute"
    ),
    col_types = "ccccciccc"
  ) %>%
    dplyr::filter(feature == "gene") %>%
    dplyr::transmute(
      seqname,
      start = as.numeric(start),
      end = as.numeric(end),
      gene_id = stringr::str_remove_all(
        stringr::str_extract(attribute, 'gene_id "[^"]+"'),
        'gene_id "|"$'
      ),
      gene_name = stringr::str_remove_all(
        stringr::str_extract(attribute, 'gene_name "[^"]+"'),
        'gene_name "|"$'
      ),
      gene_type = stringr::str_remove_all(
        stringr::str_extract(attribute, 'gene_type "[^"]+"'),
        'gene_type "|"$'
      )
    )

  GenomicRanges::GRanges(
    seqnames = as.character(gtf$seqname),
    ranges = IRanges::IRanges(start = gtf$start, end = gtf$end),
    gene_id = gtf$gene_id,
    gene_name = gtf$gene_name,
    gene_type = gtf$gene_type
  )
}

#' Annotate genomic regions with preloaded GENCODE GTF
#'
#' Maps genomic coordinates from a data frame to gene annotations from a preloaded
#' GENCODE GRanges object. Optionally filters for significant regions based on p-value.
#' Handles multiple overlaps by concatenating gene IDs, names, and types with semicolons.
#'
#' @param df A data frame containing genomic regions with columns for chromosome, start, end, and optionally p-values.
#' @param gr_gtf A GRanges object containing preloaded GENCODE gene annotations (gene_id, gene_name, gene_type).
#' @param chrom_col Column name for chromosome (default: chrom).
#' @param start_col Column name for start position (default: start).
#' @param end_col Column name for end position (default: end).
#' @param pval_col Column name for p-values used to filter significant regions (default: map_pvalue).
#' @param alpha Numeric; significance threshold for filtering regions (default: 0.05).
#' @param filter_significant Logical; if TRUE, only regions with p-value < alpha are annotated (default: TRUE).
#'
#' @return A data frame identical to `df` with additional columns:
#'   - `gene_id`: concatenated gene IDs overlapping the region
#'   - `gene_name`: concatenated gene names overlapping the region
#'   - `gene_type`: concatenated gene types overlapping the region
#'   If no overlaps are found, the new columns contain NA.
#'
#' @export
annotate_with_gencode_preloaded <- function(
    df,
    gr_gtf,
    chrom_col = chrom,
    start_col = start,
    end_col = end,
    pval_col = map_pvalue,
    alpha = 0.05,
    filter_significant = TRUE
) {
  chrom_col <- rlang::enquo(chrom_col)
  start_col <- rlang::enquo(start_col)
  end_col <- rlang::enquo(end_col)
  pval_col <- rlang::enquo(pval_col)

  df <- df %>% dplyr::mutate(row_id = dplyr::row_number())
  df_sub <- if (filter_significant) {
    df %>% dplyr::filter(!!pval_col < alpha)
  } else {
    df
  }
  if (nrow(df_sub) == 0) {
    return(
      df %>%
        dplyr::mutate(
          gene_id = NA_character_,
          gene_name = NA_character_,
          gene_type = NA_character_
        ) %>%
        dplyr::select(-row_id)
    )
  }

  gr_df <- GenomicRanges::GRanges(
    seqnames = as.character(dplyr::pull(df_sub, !!chrom_col)),
    ranges = IRanges::IRanges(
      start = as.numeric(dplyr::pull(df_sub, !!start_col)),
      end = as.numeric(dplyr::pull(df_sub, !!end_col))
    )
  )

  hits <- GenomicRanges::findOverlaps(gr_df, gr_gtf, ignore.strand = TRUE)

  annot <- dplyr::tibble(
    row_id = df_sub$row_id[S4Vectors::queryHits(hits)],
    gene_id = gr_gtf$gene_id[S4Vectors::subjectHits(hits)],
    gene_name = gr_gtf$gene_name[S4Vectors::subjectHits(hits)],
    gene_type = gr_gtf$gene_type[S4Vectors::subjectHits(hits)]
  ) %>%
    dplyr::group_by(row_id) %>%
    dplyr::summarise(
      dplyr::across(everything(), ~ paste(unique(.x), collapse = ";")),
      .groups = "drop"
    )

  df %>%
    dplyr::left_join(annot, by = "row_id") %>%
    dplyr::select(-row_id)
}

#' Plot DMR volcano with separate effect and significance thresholds
#'
#' Generates a volcano plot for DMR (differentially methylated regions) results.
#' Classifies hits as moderate or strong based on both effect size and p-value thresholds.
#' Strong hits are highlighted and optionally labeled.
#' Effect and p-value thresholds are shown with dashed (moderate) and dotted (strong) lines,
#' and included in the plot legend for clarity.
#'
#' @param df Data frame containing DMR results.
#' @param effect_col Column name for effect size.
#' @param pval_col Column name for p-value.
#' @param chrom_col Column name for chromosome (default: chrom).
#' @param start_col Column name for start position (default: start).
#' @param end_col Column name for end position (default: end).
#' @param alpha_effect_thres Numeric, p-value threshold for moderate hits (default: 0.05).
#' @param alpha_strong_thres Numeric, p-value threshold for strong hits (default: 0.01).
#' @param effect_thresh Numeric, effect size threshold for moderate hits (default: 0.3).
#' @param strong_thresh Numeric, effect size threshold for strong hits (default: 0.5).
#' @param top_n Integer, number of top strong hits to label (default: 20).
#' @param title Character, plot title (default: "").
#' @param subtitle Character, plot subtitle (default: "").
#'
#' @return ggplot2 object representing the volcano plot.
#' @export
plot_dmr_volcano <- function(
  df,
  effect_col,
  pval_col,
  chrom_col = chrom,
  start_col = start,
  end_col = end,
  alpha_effect_thres = 0.05,
  alpha_strong_thres = 0.01,
  effect_thresh = 0.3,
  strong_thresh = 0.5,
  top_n = 20,
  title = "",
  subtitle = ""
) {
  effect_col <- rlang::enquo(effect_col)
  pval_col <- rlang::enquo(pval_col)
  chrom_col <- rlang::enquo(chrom_col)
  start_col <- rlang::enquo(start_col)
  end_col <- rlang::enquo(end_col)

  df <- df %>%
    dplyr::mutate(
      effect = !!effect_col,
      pval = pmax(!!pval_col, 1e-300),
      neglog = -log10(pval),
      region = dplyr::if_else(
        !is.na(gene_name),
        paste0(gene_name, " (", !!chrom_col, ":", !!start_col, "-", !!end_col, ")"),
        paste0(!!chrom_col, ":", !!start_col, "-", !!end_col)
      ),
      reg = dplyr::case_when(
        pval < alpha_strong_thres & effect >= strong_thresh ~ "Up (strong)",
        pval < alpha_effect_thres & effect >= effect_thresh ~ "Up (moderate)",
        pval < alpha_strong_thres & effect <= -strong_thresh ~ "Down (strong)",
        pval < alpha_effect_thres & effect <= -effect_thresh ~ "Down (moderate)",
        TRUE ~ "NS"
      )
    )

  counts <- df %>%
    dplyr::count(reg) %>%
    tidyr::complete(
      reg = c("Up (strong)", "Up (moderate)", "Down (strong)", "Down (moderate)", "NS"),
      fill = list(n = 0)
    )

  # Prepare counts as separate "dummy" data for stacking above the legend
  counts_text <- paste0(
    "Up strong: ", counts$n[counts$reg == "Up (strong)"], "\n",
    "Up moderate: ", counts$n[counts$reg == "Up (moderate)"], "\n",
    "Down strong: ", counts$n[counts$reg == "Down (strong)"], "\n",
    "Down moderate: ", counts$n[counts$reg == "Down (moderate)"], "\n",
    "Total sig: ", sum(counts$n[counts$reg != "NS"])
  )

  top_hits <- df %>%
    dplyr::filter(
      (effect >= strong_thresh & pval < alpha_strong_thres) |
      (effect <= -strong_thresh & pval < alpha_strong_thres)
    ) %>%
    dplyr::slice_min(order_by = pval, n = top_n)

  if (!requireNamespace("ggrastr", quietly = TRUE)) install.packages("ggrastr")

  gg <- ggplot2::ggplot(df, ggplot2::aes(x = effect, y = neglog)) +
    ggrastr::rasterise(
      ggplot2::geom_point(ggplot2::aes(color = reg), alpha = 0.6, size = 1.2),
      dpi = 300
    ) +
    ggplot2::geom_point(data = top_hits, color = "black", size = 2) +
    ggrepel::geom_text_repel(
      data = top_hits,
      ggplot2::aes(label = region),
      size = 3,
      max.overlaps = 20
    ) +
    # Effect size thresholds
    ggplot2::geom_vline(xintercept = c(-effect_thresh, effect_thresh), linetype = "dashed", color = "grey40") +
    ggplot2::geom_vline(xintercept = c(-strong_thresh, strong_thresh), linetype = "dotted", color = "grey20") +
    # P-value thresholds
    ggplot2::geom_hline(yintercept = -log10(alpha_effect_thres), linetype = "dashed", color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(alpha_strong_thres), linetype = "dotted", color = "grey20") +
    # Counts as separate "color" mapping (dummy) to stack with legend
    ggplot2::geom_point(aes(x = NA, y = NA, color = counts_text)) +
    ggplot2::scale_color_manual(
      values = c(
        "Up (strong)" = "#D55E00",
        "Up (moderate)" = "#F4A261",
        "Down (strong)" = "#0072B2",
        "Down (moderate)" = "#56B4E9",
        "NS" = "grey80",
        counts_text = "black"
      )
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Effect size",
      y = expression(-log[10]("p-value")),
      color = ""
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      legend.text.align = 0,
      legend.key.height = ggplot2::unit(0.8, "lines"),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(5, 120, 5, 5),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::coord_cartesian(clip = "off")

  return(gg)
}

#' Plot total coverage boxplots
#'
#' @param df_long Long dataframe with total coverage
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @return Named list of ggplot objects
plot_total_modified <- function(df_long, title = "", subtitle = "") {
  ylim_value <- stats::quantile(df_long$total_value, 0.85, na.rm = TRUE)[[1]]
  boxplot <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = sample_total, y = total_value)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Sample",
      y = "Coverage"
    )
  boxplot_zoom <- boxplot + ggplot2::coord_cartesian(ylim = c(0, ylim_value))
  list(boxplot_zoom = boxplot_zoom)
}

#' Save all plots to a single HTML page (memory-efficient)
#'
#' @param plots Nested list of ggplot objects
#' @param output_file File path for HTML output
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Rasterization DPI
#' @return NULL (writes HTML file)
save_all_plots_to_html <- function(
  plots,
  output_file,
  width = 14,
  height = 6,
  dpi = 150
) {
  if (!requireNamespace("htmltools", quietly = TRUE)) {
    install.packages("htmltools")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }

  library(htmltools)
  library(ggplot2)

  # Flatten plots + preserve names (including nested)
  flatten_plots <- function(lst, parent_name = NULL) {
    res <- list()
    nms <- names(lst)

    for (i in seq_along(lst)) {
      p <- lst[[i]]
      name <- nms[i]

      full_name <- if (!is.null(name) && nzchar(name)) {
        if (!is.null(parent_name)) paste(parent_name, name, sep = "_") else name
      } else {
        parent_name
      }

      if (is.list(p)) {
        res <- c(res, flatten_plots(p, full_name))
      } else if (!is.null(p) && inherits(p, "ggplot")) {
        res <- c(res, setNames(list(p), full_name))
      }
    }
    res
  }

  flat_plots <- flatten_plots(plots)

  # Generate names
  plot_names <- names(flat_plots)
  plot_names[is.null(plot_names) | plot_names == ""] <- paste0("plot_", seq_along(flat_plots))

  # Clean names
  clean_name <- function(x) {
    x <- tolower(x)
    x <- gsub("[^a-z0-9]+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
  }

  plot_names <- vapply(plot_names, clean_name, character(1))

  # Ensure uniqueness
  make_unique <- function(x) {
    ave(x, x, FUN = function(v) {
      if (length(v) == 1) return(v)
      paste0(v, "_", seq_along(v))
    })
  }

  plot_names <- make_unique(plot_names)

  # Create image directory
  output_dir <- dirname(output_file)
  img_dir <- file.path(output_dir, "images")
  dir.create(img_dir, showWarnings = FALSE, recursive = TRUE)

  # Build HTML body
  html_body <- lapply(seq_along(flat_plots), function(i) {
    p <- flat_plots[[i]]
    name <- plot_names[i]

    img_filename <- paste0(name, ".png")
    img_filepath <- file.path(img_dir, img_filename)

    # Save plot to persistent file
    ggplot2::ggsave(
      filename = img_filepath,
      plot = p,
      width = width,
      height = height,
      dpi = dpi,
      limitsize = FALSE
    )

    tags$div(
      tags$h3(name),
      tags$img(
        src = file.path("images", img_filename),
        style = "max-width:100%; height:auto; display:block; margin-bottom:20px;"
      )
    )
  })

  # Build HTML
  html_page <- tags$html(
    tags$head(
      tags$title("dRNA m6A Plots"),
      tags$style(HTML("
        body { font-family: Arial, sans-serif; margin: 20px; }
        h3 { margin-bottom: 5px; }
      "))
    ),
    tags$body(html_body)
  )

  htmltools::save_html(html_page, file = output_file)

  message("HTML saved to: ", normalizePath(output_file))
  message("Images saved in: ", normalizePath(img_dir))
}

# 3. Main workflow ----
main <- function() {
  segmeth_file <- "/Users/ejong19/Documents/Werk/Kenna/ont_drna/pilot1_and_pilot2.genetypes.segmeth.tsv"
  bed_dir <- "/Users/ejong19/Documents/Werk/Kenna/ont_drna/data/modkit/dmr/"
  plot_dir <- "/Users/ejong19/Documents/Werk/Kenna/ont_drna/data/plots"

  # Segment methylation plots
  segmeth_data <- readr::read_tsv(segmeth_file)
  segmeth_long <- pivot_data_to_long(segmeth_data)
  df_segmeth_percent <- prepare_data_ratio(segmeth_long)
  ratio_plots <- generate_ratio_plots(df_segmeth_percent)
  # TODO: save plot to PDF.
  rm(segmeth_data, segmeth_long, df_segmeth_percent)
  gc()

  gr_gtf <- preload_gencode_gtf(
    "/Users/ejong19/Downloads/gencode.v49.annotation.gtf"
  )
  # DMR volcano and coverage plots
  chromosomes <- paste0("chr", 1:22)
  bed_files <- read_bed_files(bed_dir)
  thres_coverage <- 3

  result_plots <- list()
  for (bed_name in names(bed_files)) {
    if (!stringr::str_ends(bed_name, "dmr_single_base_segment")) {
      next
    }
    bed_file <- bed_files[[bed_name]]

    df_chr <- dplyr::filter(bed_file, chrom %in% chromosomes)
    df_chr <- annotate_with_gencode_preloaded(df_chr, gr_gtf)

    df_chr_long <- tidyr::pivot_longer(
      df_chr,
      c(a_total, b_total),
      names_to = "sample_total",
      values_to = "total_value"
    )

    df_chr_cov <- dplyr::filter(
      df_chr,
      a_total > thres_coverage,
      b_total > thres_coverage
    )

    df_chr_cov_long <- tidyr::pivot_longer(
      df_chr_cov,
      c(a_total, b_total),
      names_to = "sample_total",
      values_to = "total_value"
    )

    volcano <- plot_dmr_volcano(
      df_chr_cov,
      effect_size,
      map_pvalue,
      top_n = 20,
      title = paste0(bed_name, ": DMR Volcano"),
      subtitle = "Selected chr1-22, coverage > 3"
    )

    volcano_balanced <- NULL
    if ("balanced_effect_size" %in% colnames(bed_file)) {
      volcano_balanced <- plot_dmr_volcano(
        df_chr_cov,
        balanced_effect_size,
        balanced_map_pvalue,
        top_n = 20,
        title = paste0(bed_name, ": Balanced DMR Volcano"),
        subtitle = "Selected chr1-22, coverage > 3"
      )
    }

    boxplots <- plot_total_modified(
      df_chr_long,
      title = paste0(bed_name, ": Coverage by sample"),
      subtitle = "Selected chr1-22"
    )
    boxplots_cov <- plot_total_modified(
      df_chr_cov_long,
      title = paste0(bed_name, ": Coverage by sample"),
      subtitle = "Selected chr1-22, coverage > 3"
    )

    result_plots[[bed_name]] <- setNames(
      list(
        volcano,
        volcano_balanced,
        boxplots,
        boxplots_cov
      ),
      c(
        paste(bed_name, "volcano", sep = "_"),
        paste(bed_name, "volcano_balanced", sep = "_"),
        paste(bed_name, "boxplots", sep = "_"),
        paste(bed_name, "boxplots_cov", sep = "_")
      )
    )
  }

  for (name in names(result_plots)) {
    save_all_plots_to_html(
      result_plots[[name]],
      file.path(plot_dir, paste0(sub("\\.[^.]*$", "",name), ".dmr_and_coverage_plots.html")),
      width = 14,
      height = 6,
      dpi = 150
    )
  }
}

# 4. Execute main ----
main()

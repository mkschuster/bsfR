#!/usr/bin/env Rscript
#
# BSF R script to create an Enhanced Volcano plot.
#
#
# Copyright 2013 - 2020 Michael K. Schuster
#
# Biomedical Sequencing Facility (BSF), part of the genomics core facility of
# the Research Center for Molecular Medicine (CeMM) of the Austrian Academy of
# Sciences and the Medical University of Vienna (MUW).
#
#
# This file is part of BSF R.
#
# BSF R is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BSF R is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with BSF R.  If not, see <http://www.gnu.org/licenses/>.

suppressPackageStartupMessages(expr = library(package = "optparse"))

argument_list <-
  optparse::parse_args(object = optparse::OptionParser(
    option_list = list(
      optparse::make_option(
        opt_str = c("--verbose", "-v"),
        action = "store_true",
        default = TRUE,
        help = "Print extra output [default]",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = c("--quiet", "-q"),
        action = "store_false",
        default = FALSE,
        dest = "verbose",
        help = "Print little output",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = c("--design-name"),
        # default = "global",
        dest = "design_name",
        help = "Design name",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--l2fc-threshold"),
        default = 1.0,
        dest = "l2fc_threshold",
        help = "Threshold for the log2(fold-change) [1.0]",
        type = "numeric"
      ),
      optparse::make_option(
        opt_str = c("--p-threshold"),
        default = 1e-05,
        dest = "p_threshold",
        help = "Threshold for the unadjusted p-value [1e-05]",
        type = "numeric"
      ),
      optparse::make_option(
        opt_str = c("--padj-threshold"),
        default = 0.1,
        dest = "padj_threshold",
        help = "Threshold for the adjusted p-value [0.1]",
        type = "numeric"
      ),
      optparse::make_option(
        opt_str = c("--padj"),
        action = "store_true",
        default = FALSE,
        dest = "plot_padj",
        help = "Plot adjusted p-values [FALSE]",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = c("--gene-path"),
        dest = "gene_path",
        help = "Gene list file path for annotation [NULL]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--genome-directory"),
        default = ".",
        dest = "genome_directory",
        help = "Genome directory path [.]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--output-directory"),
        default = ".",
        dest = "output_directory",
        help = "Output directory path [.]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--x-limits"),
        dest = "x_limits",
        help = "x-axis limits separated by a comma (lower,upper) [NULL]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--y-limits"),
        dest = "y_limits",
        help = "y-axis limits separated by a comma (lower,upper) [NULL]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--plot-dpi"),
        default = 72,
        dest = "plot_dpi",
        help = "Plot resolution in dpi [72]",
        type = "numeric"
      ),
      optparse::make_option(
        opt_str = c("--plot-width"),
        default = 7.0,
        dest = "plot_width",
        help = "Plot width in inches [7.0]",
        type = "numeric"
      ),
      optparse::make_option(
        opt_str = c("--plot-height"),
        default = 7.0,
        dest = "plot_height",
        help = "Plot height in inches [7.0]",
        type = "numeric"
      )
    )
  ))

# Check the input.

if (is.null(x = argument_list$design_name)) {
  stop("Missing --design-name option")
}

suppressPackageStartupMessages(expr = library(package = "tidyverse"))
suppressPackageStartupMessages(expr = library(package = "bsfR"))
suppressPackageStartupMessages(expr = library(package = "EnhancedVolcano"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

prefix_deseq <-
  bsfR::bsfrd_get_prefix_deseq(design_name = argument_list$design_name)

prefix_volcano <-
  bsfR::bsfrd_get_prefix_volcano(design_name = argument_list$design_name)

output_directory <-
  file.path(argument_list$output_directory, prefix_volcano)
if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

# Plot Annotation ---------------------------------------------------------


plot_labels <- character()
if (!is.null(x = argument_list$gene_path)) {
  plot_annotation_tibble <-
    bsfR::bsfrd_read_gene_set_tibble(
      genome_directory = argument_list$genome_directory,
      design_name = argument_list$design_name,
      gene_set_path = argument_list$gene_path,
      verbose = argument_list$verbose
    )
  plot_labels <- plot_annotation_tibble$gene_name
  rm(plot_annotation_tibble)
}

# Contrasts Frame ---------------------------------------------------------


contrast_tibble <-
  bsfR::bsfrd_read_contrast_tibble(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name,
    verbose = argument_list$verbose
  )
if (nrow(x = contrast_tibble) == 0L) {
  stop("No contrast remaining after selection for design name.")
}

for (contrast_index in seq_len(length.out = nrow(x = contrast_tibble))) {
  contrast_character <-
    bsfR::bsfrd_get_contrast_character(contrast_tibble = contrast_tibble, index = contrast_index)

  deseq_results_tibble <-
    bsfR::bsfrd_read_result_tibble(
      genome_directory = argument_list$genome_directory,
      design_name = argument_list$design_name,
      contrast_tibble = contrast_tibble,
      index = contrast_index,
      verbose = argument_list$verbose
    )

  # deseq_results_tibble <-
  #   dplyr::filter(.data = deseq_results_tibble,!is.na(x = log2FoldChange),!is.na(x = pvalue))

  # Enhanced Volcano plot -----------------------------------------------

  if (argument_list$plot_padj) {
    plot_paths <- file.path(output_directory,
                            paste(
                              paste(
                                prefix_volcano,
                                "contrast",
                                contrast_character,
                                "adjp",
                                sep = "_"
                              ),
                              graphics_formats,
                              sep = "."
                            ))
  } else {
    plot_paths <- file.path(output_directory,
                            paste(
                              paste(prefix_volcano,
                                    "contrast",
                                    contrast_character,
                                    sep = "_"),
                              graphics_formats,
                              sep = "."
                            ))
  }

  x <- "log2FoldChange"
  y <- if (argument_list$plot_padj) {
    "padj"
  } else {
    "pvalue"
  }

  # Replace adjusted p-values or p-values == 0.
  # The correction in EnhanceVolcano does not seem to work, as a comparison with
  # 0.0 (i.e. x == 0.0) is problematic.
  if (any(deseq_results_tibble[, y] < .Machine$double.xmin)) {
    message(
      "Adjusting some ",
      y,
      " lower than machine-specific double minimum ",
      .Machine$double.xmin
    )
    deseq_results_tibble[which(x = deseq_results_tibble[, y] < .Machine$double.xmin), y] <-
      .Machine$double.xmin
  }
  deseq_results_frame <- base::as.data.frame(deseq_results_tibble)

  ggplot_object <- EnhancedVolcano::EnhancedVolcano(
    # Without base::as.data.frame(), the log2FoldChange variable is reported to
    # be not numeric, when in fact it is.
    toptable = deseq_results_frame,
    lab = deseq_results_frame$gene_name,
    x = x,
    y = y,
    xlim = if (is.null(x = argument_list$x_limits)) {
      # Use the function default limits.
      c(min(deseq_results_frame[, x], na.rm = TRUE),
        max(deseq_results_frame[, x], na.rm = TRUE))
    } else {
      # Split the x-limits argument into a character matrix and convert into a numeric vector.
      as.numeric(x = stringr::str_split_fixed(
        string = argument_list$x_limits,
        pattern = ",",
        n = 2L
      ))
    },
    ylim = if (is.null(x = argument_list$y_limits)) {
      # Use the function default limits.
      c(0, max(-log10(deseq_results_frame[, y]), na.rm = TRUE) + 5)
    } else {
      # Split the y-limits argument into a character matrix and convert into a numeric vector.
      as.numeric(x = stringr::str_split_fixed(
        string = argument_list$y_limits,
        pattern = ",",
        n = 2L
      ))
    },
    pCutoff = if (argument_list$plot_padj) {
      argument_list$padj_threshold
    } else {
      argument_list$p_threshold
    },
    xlab = if (argument_list$plot_padj) {
      bquote(expr = ~ Log[2] ~ "fold change")
    } else {
      bquote(expr = ~ Log[2] ~ "fold change")
    },
    ylab = if (argument_list$plot_padj) {
      bquote(expr = ~ -Log[10] ~ adjusted ~ italic(P))
    } else {
      bquote(expr = ~ -Log[10] ~ italic(P))
    },
    subtitle = ggplot2::waiver(),
    caption = ggplot2::waiver(),
    legend = if (argument_list$plot_padj) {
      c("NS", "Log2 FC", "Adj. P", "Adj. P & Log2 FC")
    } else {
      c("NS", "Log2 FC", "P", "P & Log2 FC")
    },
    selectLab = plot_labels,
    drawConnectors = TRUE,
    endsConnectors = 'last'
  )

  for (plot_path in plot_paths) {
    ggplot2::ggsave(
      filename = plot_path,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height,
      units = "in",
      dpi = argument_list$plot_dpi,
      limitsize = FALSE
    )
  }
  rm(
    plot_path,
    ggplot_object,
    plot_paths,
    y,
    x,
    deseq_results_frame,
    deseq_results_tibble,
    contrast_character
  )
}
rm(
  contrast_index,
  plot_labels,
  contrast_tibble,
  output_directory,
  prefix_volcano,
  prefix_deseq,
  graphics_formats,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessionInfo())

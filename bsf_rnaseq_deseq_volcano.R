#!/usr/bin/env Rscript
#
# BSF R script to create an Enhanced Volcano plot.
#
#
# Copyright 2013 - 2019 Michael K. Schuster
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

argument_list <- parse_args(object = OptionParser(
  option_list = list(
    make_option(
      opt_str = c("--verbose", "-v"),
      action = "store_true",
      default = TRUE,
      help = "Print extra output [default]",
      type = "logical"
    ),
    make_option(
      opt_str = c("--quiet", "-q"),
      action = "store_false",
      default = FALSE,
      dest = "verbose",
      help = "Print little output",
      type = "logical"
    ),
    make_option(
      opt_str = c("--design-name"),
      # default = "global",
      dest = "design_name",
      help = "Design name",
      type = "character"
    ),
    make_option(
      opt_str = c("--l2fc-threshold"),
      default = 1.0,
      dest = "l2fc_threshold",
      help = "Threshold for the log2(fold-change) [1.0]",
      type = "numeric"
    ),
    make_option(
      opt_str = c("--p-threshold"),
      default = 1e-05,
      dest = "p_threshold",
      help = "Threshold for the unadjusted p-value [1e-05]",
      type = "numeric"
    ),
    make_option(
      opt_str = c("--padj-threshold"),
      default = 0.1,
      dest = "padj_threshold",
      help = "Threshold for the adjusted p-value [0.1]",
      type = "numeric"
    ),
    make_option(
      opt_str = c("--padj"),
      action = "store_true",
      default = FALSE,
      dest = "plot_padj",
      help = "Plot adjusted p-values [FALSE]",
      type = "logical"
    ),
    make_option(
      opt_str = c("--gene-path"),
      dest = "gene_path",
      help = "Gene list file path for annotation [NULL]",
      type = "character"
    ),
    make_option(
      opt_str = c("--genome-directory"),
      default = ".",
      dest = "genome_directory",
      help = "Genome directory path [.]",
      type = "character"
    ),
    make_option(
      opt_str = c("--output-directory"),
      default = ".",
      dest = "output_directory",
      help = "Output directory path [.]",
      type = "character"
    ),
    make_option(
      opt_str = c("--plot-dpi"),
      default = 72,
      dest = "plot_dpi",
      help = "Plot resolution in dpi [72]",
      type = "numeric"
    ),
    make_option(
      opt_str = c("--plot-width"),
      default = 7.0,
      dest = "plot_width",
      help = "Plot width in inches [7.0]",
      type = "numeric"
    ),
    make_option(
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

suppressPackageStartupMessages(expr = library(package = "bsfR"))
suppressPackageStartupMessages(expr = library(package = "EnhancedVolcano"))
suppressPackageStartupMessages(expr = library(package = "tidyverse"))

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
      gene_set_path = argument_list$gene_path
    )
  plot_labels <- plot_annotation_tibble$gene_name
  rm(plot_annotation_tibble)
}

# Contrasts Frame ---------------------------------------------------------


contrast_tibble <-
  bsfR::bsfrd_read_contrast_tibble(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name
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
      index = contrast_index
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

  ggplot_object <- EnhancedVolcano::EnhancedVolcano(
    # Without base::as.data.frame(), the log2FoldChange variable is reported to
    # be not numeric, when in fact it is.
    toptable = base::as.data.frame(deseq_results_tibble),
    lab = deseq_results_tibble$gene_name,
    x = "log2FoldChange",
    y = if (argument_list$plot_padj) {
      "padj"
    } else {
      "pvalue"
    },
    pCutoff = if (argument_list$plot_padj) {
      argument_list$padj_threshold
    } else {
      argument_list$p_threshold
    },
    xlab = if (argument_list$plot_padj) {
      bquote( ~ Log[2] ~ "fold change")
    } else {
      bquote( ~ Log[2] ~ "fold change")
    },
    ylab = if (argument_list$plot_padj) {
      bquote( ~ -Log[10] ~ adjusted ~ italic(P))
    } else {
      bquote( ~ -Log[10] ~ italic(P))
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
  rm(plot_path,
     ggplot_object,
     plot_paths,
     deseq_results_tibble,
     contrast_character)
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

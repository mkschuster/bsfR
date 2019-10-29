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

suppressPackageStartupMessages(expr = library(package = "EnhancedVolcano"))
suppressPackageStartupMessages(expr = library(package = "tidyverse"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

prefix_deseq <-
  paste("rnaseq",
        "deseq",
        argument_list$design_name,
        sep = "_")

prefix_volcano <-
  paste("rnaseq",
        "deseq",
        argument_list$design_name,
        "volcano",
        sep = "_")


output_directory <-
  file.path(argument_list$output_directory, prefix_volcano)
if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

#' Load a contrast data frame
#'
#' @param contrast_character A \code{character} scalar of the contrast
#'
#' @return A \code{tibble} of annotated DESeqResults
#' @noRd
#'
#' @examples
load_contrast_frame <- function(contrast_character) {
  # Annotated Results Frame -----------------------------------------------
  # Read the annotated data.frame with all genes.
  deseq_results_tibble <- NULL

  file_path <-
    file.path(argument_list$genome_directory,
              prefix_deseq,
              paste(
                paste(prefix_deseq,
                      "contrast",
                      contrast_character,
                      "genes",
                      sep = "_"),
                "tsv",
                sep = "."
              ))

  if (file.exists(file_path) &&
      file.info(file_path)$size > 0L) {
    message("Loading DESeqResults frame for contrast ",
            contrast_character)

    deseq_results_tibble <-
      readr::read_tsv(
        file = file_path,
        col_types = cols(
          gene_id = col_character(),
          gene_version = col_double(),
          gene_name = col_character(),
          gene_biotype = col_character(),
          gene_source = col_character(),
          location = col_character(),
          baseMean = col_double(),
          log2FoldChange = col_double(),
          lfcSE = col_double(),
          stat = col_double(),
          pvalue = col_double(),
          padj = col_double(),
          significant = col_character(),
          rank_log2_fold_change = col_double(),
          rank_base_mean = col_double(),
          rank_padj = col_double(),
          max_rank = col_double()
        )
      )

  } else {
    message("No DESeqResults frame for contrast ", contrast_character)
  }
  return(deseq_results_tibble)
}

#' Load plot labels
#'
#' @return A \code{character} vector of gene identifiers.
#' @noRd
#'
#' @examples
load_plot_labels <- function() {
  # Plot labels -----------------------------------------------------------
  if (!is.null(x = argument_list$gene_path)) {
    gene_tibble <-
      readr::read_csv(
        file = argument_list$gene_path,
        col_names = TRUE,
        col_types = cols(gene_name = col_character())
      )
    return(gene_tibble[, c("gene_name"), drop = TRUE])
  } else {
    return(character())
  }
}

plot_labels <- load_plot_labels()

# Contrasts Frame ---------------------------------------------------------


# Read a data frame of contrasts with variables "Design", "Numerator", "Denominator" and "Label".
message("Loading contrast frame")
contrast_frame <-
  read.table(
    file = file.path(
      argument_list$genome_directory,
      prefix_deseq,
      paste(
        paste(prefix_deseq, "contrasts", "summary", sep = "_"),
        "tsv",
        sep = "."
      )
    ),
    header = TRUE,
    sep = "\t",
    colClasses = c(
      "Design" = "character",
      "Numerator" = "character",
      "Denominator" = "character",
      "Label" = "character",
      "Significant" = "integer"
    ),
    stringsAsFactors = FALSE
  )
# Subset to the selected design.
contrast_frame <-
  contrast_frame[contrast_frame$Design == argument_list$design_name, , drop = FALSE]
if (nrow(x = contrast_frame) == 0L) {
  stop("No contrast remaining after selection for design name.")
}

for (i in seq_len(length.out = nrow(x = contrast_frame))) {
  # The "contrast" option of the DESeq results() function expects a list of
  # numerator and denominator.
  contrast_list <-
    list("numerator" = unlist(x = strsplit(x = contrast_frame[i, "Numerator", drop = TRUE], split = ",")),
         "denominator" = unlist(x = strsplit(x = contrast_frame[i, "Denominator", drop = TRUE], split = ",")))
  contrast_character <-
    paste(paste(contrast_list$numerator, collapse = "_"),
          "against",
          if (length(x = contrast_list$denominator) > 0L) {
            paste(contrast_list$denominator, collapse = "_")
          } else {
            "intercept"
          },
          sep = "_")
  rm(contrast_list)

  deseq_results_tibble <-
    load_contrast_frame(contrast_character = contrast_character)

  # deseq_results_tibble <-
  #   dplyr::filter(.data = deseq_results_tibble,!is.na(x = log2FoldChange),!is.na(x = pvalue))

  # Enhanced Volcano plot -----------------------------------------------

  if (argument_list$plot_padj) {
    plot_paths <- file.path(output_directory,
                            paste(
                              paste(prefix_volcano,
                                    "contrast",
                                    contrast_character,
                                    "adjp",
                                    sep = "_"),
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
      bquote(~Log[2] ~ "fold change")
    } else {
      bquote(~Log[2] ~ "fold change")
    },
    ylab = if (argument_list$plot_padj) {
      bquote(~-Log[10] ~ adjusted ~ italic(P))
    } else {
      bquote(~-Log[10] ~ italic(P))
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
rm(i)

rm(
  plot_labels,
  load_plot_labels,
  contrast_frame,
  load_contrast_frame,
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

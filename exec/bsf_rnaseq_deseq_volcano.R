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
        opt_str = c("--gene-path"),
        dest = "gene_path",
        help = "Gene set file path for custom annotation [NULL]",
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
suppressPackageStartupMessages(expr = library(package = "Nozzle.R1"))

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

# Plot Annotation Tibble --------------------------------------------------


plot_annotation_tibble <- NULL
if (!is.null(x = argument_list$gene_path)) {
  plot_annotation_tibble <-
    bsfR::bsfrd_read_gene_set_tibble(
      genome_directory = argument_list$genome_directory,
      design_name = argument_list$design_name,
      gene_set_path = argument_list$gene_path,
      verbose = argument_list$verbose
    )

  readr::write_tsv(
    x = plot_annotation_tibble,
    path = file.path(output_directory, paste0(prefix_volcano, "_gene_set.tsv")),
    col_names = TRUE
  )

  # If the tibble exists, test for NA values.
  missing_tibble <-
    dplyr::filter(.data = plot_annotation_tibble, is.na(x = .data$gene_id))
  if (nrow(x = missing_tibble) > 0L) {
    print(x = "The following genes_name values could not be resolved into gene_id values:")
    print(x = missing_tibble)
  }
  rm(missing_tibble)
}

# Contrasts Tibble --------------------------------------------------------


contrast_tibble <-
  bsfR::bsfrd_read_contrast_tibble(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name,
    verbose = argument_list$verbose
  )
if (nrow(x = contrast_tibble) == 0L) {
  stop("No contrast remaining after selection for design name.")
}

# Create a "Contrasts" report section.
nozzle_section_contrasts <-
  Nozzle.R1::newSection("Contrasts", class = SECTION.CLASS.RESULTS)
nozzle_section_contrasts <-
  addTo(parent = nozzle_section_contrasts, Nozzle.R1::newTable(table = as.data.frame(x = contrast_tibble)))

# Create a "Colcano Plots" report section.
nozzle_section_list <- list(
  "padj" = Nozzle.R1::newSection("Volcano Plots (adjusted p-value)", class = SECTION.CLASS.RESULTS),
  "pvalue" = Nozzle.R1::newSection("Volcano Plots (unadjusted p-value)", class = SECTION.CLASS.RESULTS)
)

#' Local function drawing an EnhancedVolcano object.
#'
#' @param nozzle_section_list A named \code{list} of Nozzle Report Section objects.
#' @param deseq_results_tibble A results \code{tibble}} object.
#' @param contrast_character A \code{character} scalar defining a particular contrast.
#' @param gene_labels A \code{character} vector with the gene labels named by "gene_id".
#' @param plot_index A \code{integer} index for systematic file name generation.
#' @param plot_title A \code{character} scalar with the plot title.
#'
#' @return A Nozzle Report Section.
#' @noRd
#'
#' @examples
draw_enhanced_volcano <-
  function(nozzle_section_list,
           deseq_results_tibble,
           contrast_character,
           gene_labels = character(),
           plot_index = 0L,
           plot_title = NULL) {
    for (plot_padj in c(TRUE, FALSE)) {
      plot_paths <-
        paste(
          paste(
            prefix_volcano,
            "contrast",
            contrast_character,
            if (plot_padj) {
              "padj"
            } else {
              "pvalue"
            },
            plot_index,
            sep = "_"
          ),
          graphics_formats,
          sep = "."
        )

      deseq_results_frame <-
        base::as.data.frame(deseq_results_tibble)

      x <- "log2FoldChange"
      y <- if (plot_padj) {
        "padj"
      } else {
        "pvalue"
      }

      ggplot_object <- EnhancedVolcano::EnhancedVolcano(
        # Without base::as.data.frame(), the log2FoldChange variable is reported to
        # be not numeric, when in fact it is.
        toptable = deseq_results_frame,
        lab = deseq_results_frame$gene_name,
        x = x,
        y = y,
        xlim = if (is.null(x = argument_list$x_limits)) {
          # Use the function default limits.
          c(
            min(deseq_results_frame[, x], na.rm = TRUE),
            max(deseq_results_frame[, x], na.rm = TRUE)
          )
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
          # Split the y-limits argument into a character matrix and convert into
          # a numeric vector.
          as.numeric(x = stringr::str_split_fixed(
            string = argument_list$y_limits,
            pattern = ",",
            n = 2L
          ))
        },
        pCutoff = if (plot_padj) {
          argument_list$padj_threshold
        } else {
          argument_list$p_threshold
        },
        xlab = if (plot_padj) {
          bquote(expr = ~ Log[2] ~ "fold change")
        } else {
          bquote(expr = ~ Log[2] ~ "fold change")
        },
        ylab = if (plot_padj) {
          bquote(expr = ~ -Log[10] ~ adjusted ~ italic(P))
        } else {
          bquote(expr = ~ -Log[10] ~ italic(P))
        },
        subtitle = ggplot2::waiver(),
        caption = ggplot2::waiver(),
        legend = if (plot_padj) {
          c("NS", "Log2 FC", "Adj. P", "Adj. P & Log2 FC")
        } else {
          c("NS", "Log2 FC", "P", "P & Log2 FC")
        },
        selectLab = gene_labels,
        drawConnectors = TRUE,
        endsConnectors = 'last'
      )

      for (plot_path in plot_paths) {
        ggplot2::ggsave(
          filename = file.path(output_directory, plot_path),
          plot = ggplot_object,
          width = argument_list$plot_width,
          height = argument_list$plot_height,
          units = "in",
          dpi = argument_list$plot_dpi,
          limitsize = FALSE
        )
      }

      nozzle_section_list[[if (plot_padj) {
        "padj"
      } else {
        "pvalue"
      }]] <-
        Nozzle.R1::addTo(
          parent = nozzle_section_list[[if (plot_padj) {
            "padj"
          } else {
            "pvalue"
          }]],
          Nozzle.R1::newFigure(
            file = plot_paths[2L],
            "Volcano plot for contrast ",
            Nozzle.R1::asStrong(contrast_tibble$Label[contrast_index]),
            fileHighRes = plot_paths[1L]
          )
        )

      rm(plot_path,
         ggplot_object,
         deseq_results_frame,
         x,
         y,
         plot_paths)
    }
    rm(plot_padj)

    return(nozzle_section_list)
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

  # Filter genes with NA values in either log2FoldChange, pvalue pr padj, which
  # are a consequence of Cook's distance filtering in the DESeq2::results()
  # function.
  deseq_results_tibble <-
    dplyr::filter(.data = deseq_results_tibble,!(
      is.na(x = .data$log2FoldChange) |
        is.na(x = .data$pvalue) |
        is.na(x = .data$padj)
    ))

  # Replace adjusted p-values or p-values == 0.0.
  # The correction in EnhanceVolcano does not seem to work, as testing equality
  # with double (i.e. x == 0.0) is problematic.
  #
  # FIXME: Should this block be kept to just issue the message?
  # if (any(deseq_results_tibble[, y, drop = TRUE] < .Machine$double.xmin)) {
  #   message(
  #     "Adjusting some ",
  #     y,
  #     " lower than machine-specific double minimum ",
  #     .Machine$double.xmin
  #   )
  #   deseq_results_tibble[which(x = deseq_results_tibble[, y, drop = TRUE] < .Machine$double.xmin), y] <-
  #     .Machine$double.xmin
  # }
  deseq_results_tibble <-
    dplyr::mutate(
      .data = deseq_results_tibble,
      "pvalue" = dplyr::if_else(
        condition = .data$pvalue < .Machine$double.xmin,
        true = .Machine$double.xmin,
        false = .data$pvalue
      ),
      "padj" = dplyr::if_else(
        condition = .data$padj < .Machine$double.min,
        true = .Machine$double.xmin,
        false = .data$padj
      )
    )

  # Enhanced Volcano plot -----------------------------------------------


  # Draw an empty plot with plot index 0L.
  nozzle_section_list <-
    draw_enhanced_volcano(
      nozzle_section_list = nozzle_section_list,
      deseq_results_tibble = deseq_results_tibble,
      contrast_character = contrast_character
    )

  if (!is.null(x = plot_annotation_tibble)) {
    # A plot_annotation_tibble exists to select gene labels from.
    plot_names <- unique(plot_annotation_tibble$plot_name)
    for (plot_index in seq_along(along.with = plot_names)) {
      # Filter for plot_name values.
      filtered_tibble <-
        dplyr::filter(.data = plot_annotation_tibble,
                      .data$plot_name == plot_names[plot_index])
      gene_labels <- filtered_tibble$gene_label
      names(x = gene_labels) <- filtered_tibble$gene_id
      rm(filtered_tibble)

      nozzle_section_list <-
        draw_enhanced_volcano(
          nozzle_section_list = nozzle_section_list,
          deseq_results_tibble = deseq_results_tibble,
          contrast_character = contrast_character,
          gene_labels = gene_labels,
          plot_index = plot_index,
          plot_title = plot_names[plot_index]
        )
      rm(gene_labels)
    }
    rm(plot_index, plot_names)
  }
  rm(deseq_results_tibble,
     contrast_character)
}

nozzle_report <-
  Nozzle.R1::newCustomReport("Complex Heatmap Report", version = 0)

nozzle_report <-
  Nozzle.R1::addTo(parent = nozzle_report, nozzle_section_contrasts)

nozzle_report <-
  Nozzle.R1::addTo(parent = nozzle_report, nozzle_section_list$padj)

nozzle_report <-
  Nozzle.R1::addTo(parent = nozzle_report, nozzle_section_list$pvalue)

Nozzle.R1::writeReport(report = nozzle_report,
                       filename = file.path(output_directory,
                                            paste(prefix_volcano,
                                                  "report",
                                                  sep = "_")))

rm(
  nozzle_report,
  nozzle_section_list,
  nozzle_section_contrasts,
  contrast_index,
  contrast_tibble,
  plot_annotation_tibble,
  output_directory,
  prefix_volcano,
  prefix_deseq,
  graphics_formats,
  argument_list,
  draw_enhanced_volcano
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessionInfo())

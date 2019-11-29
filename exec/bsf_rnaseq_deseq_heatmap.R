#!/usr/bin/env Rscript
#
# BSF R script to draw Complex Heatmaps from DESeq2 result tables.
#
# For each design, the bsf_rnaseq_deseq_analysis.R script saves DESeqDataSet and
# DESeqTransform objects to disk, so that running it is a prerequisite. This
# script then loads the DESeqDataSet object to get access to the column (i.e.
# sample) annotation (via colData()) and the DESeqTransform object to get access
# to the transformed counts (via assay()).
#
# The contrast summary tibble provides the number of significantly differentially
# expressed genes, per contrast.
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
      opt_str = c("--variables"),
      # default = "",
      dest = "variables",
      help = "Comma-separated list of variables [all from design]",
      type = "character"
    ),
    make_option(
      opt_str = c("--maximum-number"),
      default = 400L,
      dest = "maximum_number",
      help = "Maximum number of genes to plot or -1 for all significant [400]",
      type = "integer"
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
      opt_str = c("--plot-width"),
      default = 14.0,
      dest = "plot_width",
      help = "Plot width in inches [14.0]",
      type = "numeric"
    ),
    make_option(
      opt_str = c("--plot-height"),
      default = 36.0,
      dest = "plot_height",
      help = "Plot height in inches [36.0]",
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
suppressPackageStartupMessages(expr = library(package = "ComplexHeatmap"))
suppressPackageStartupMessages(expr = library(package = "DESeq2"))
suppressPackageStartupMessages(expr = library(package = "Nozzle.R1"))
suppressPackageStartupMessages(expr = library(package = "stringr"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

prefix_deseq <-
  bsfR::bsfrd_get_prefix_deseq(design_name = argument_list$design_name)

prefix_heatmap <-
  bsfR::bsfrd_get_prefix_heatmap(design_name = argument_list$design_name)

output_directory <-
  file.path(argument_list$output_directory, prefix_heatmap)
if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

plot_annotation_tibble <- NULL
if (!is.null(x = argument_list$gene_path)) {
  plot_annotation_tibble <-
    bsfR::bsfrd_read_gene_set_tibble(
      genome_directory = argument_list$genome_directory,
      design_name = argument_list$design_name,
      gene_set_path = argument_list$gene_path
    )

  readr::write_tsv(
    x = plot_annotation_tibble,
    path = file.path(output_directory, paste0(prefix_heatmap, "_gene_set.tsv")),
    col_names = TRUE
  )

  # If the tibble exists test for NA values.
  missing_tibble <-
    dplyr::filter(.data = plot_annotation_tibble, is.na(x = .data$gene_id))
  if (nrow(x = missing_tibble) > 0L) {
    print(x = "The following genes_name values could not be resolved into gene_id values:")
    print(x = missing_tibble)
  }
  rm(missing_tibble)
}

# DESeqDataSet ------------------------------------------------------------


# Load a pre-calculated DESeqDataSet object.
# It is required for the sample annotation.
deseq_data_set <-
  bsfR::bsfrd_read_deseq_data_set(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name
  )

# DESeqTransform ----------------------------------------------------------


# Load a previously saved "blind" or "model" DESeqTransform object
# that serves as the base for the Heatmap.
deseq_transform <-
  bsfR::bsfrd_read_deseq_transform(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name,
    model = TRUE
  )

suffix <- "model"


# Contrasts Tibble --------------------------------------------------------


# Read a contrast tibble with variables "Design", "Numerator", "Denominator" and "Label".
contrast_tibble <-
  bsfR::bsfrd_read_contrast_tibble(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name,
    summary = TRUE
  )

# Select column data variables to annotate in the heat map.
if (is.null(argument_list$variables)) {
  # Unfortunately, the DESeqDataSet desgin() accessor provides no meaningful
  # formula, if the design was not full rank. In those cases the formula is set
  # to ~1 and the Wald testing is based on a model matrix. To get the design
  # variables, load the initial design tibble and call the all.vars()
  # function on the formula object.
  design_tibble <-
    bsfR::bsfrd_read_design_tibble(
      genome_directory = argument_list$genome_directory,
      design_name = argument_list$design_name
    )
  if (nrow(x = design_tibble) == 0L) {
    stop("No design remaining after selection for design name.")
  }
  variable_names <-
    all.vars(expr = as.formula(object = design_tibble[1L, "full_formula", drop = TRUE]))
  rm(design_tibble)
} else {
  variable_names <-
    stringr::str_split(string = argument_list$variables, pattern = ",")[[1L]]
}
column_annotation_frame <-
  data.frame(SummarizedExperiment::colData(x = deseq_data_set)[, variable_names, drop = FALSE])
rm(variable_names)

# Create a "Contrasts" report section
nozzle_section_contrasts <-
  Nozzle.R1::newSection("Contrasts", class = SECTION.CLASS.RESULTS)
nozzle_section_contrasts <-
  addTo(parent = nozzle_section_contrasts, Nozzle.R1::newTable(table = as.data.frame(x = contrast_tibble)))
nozzle_section_heatmaps <-
  Nozzle.R1::newSection("Heatmaps", class = SECTION.CLASS.RESULTS)

#' Local function drawing a ComplexHeatmap object.
#'
#' @param nozzle_section A Nozzle Report Section.
#' @param deseq_results_frame A results \code{data.frame} object.
#' @param top_gene_identifiers A \code{character} vector with the top-scoring gene identifier (gene_id) values.
#' @param contrast_character A \code{character} scalar defining a particular contrast.
#' @param plot_title A \code{character} scalar with the plot title.
#' @param file_index A \code{integer} index for systematic file name generation.
#'
#' @return A Nozzle Report Section.
#' @noRd
#'
#' @examples
draw_complex_heatmap <-
  function(nozzle_section,
           deseq_results_frame,
           top_gene_identifiers,
           contrast_character,
           plot_title = NULL,
           file_index = NULL) {
    if (length(x = top_gene_identifiers) > 0L) {
      # Draw a ComplexHeatmap.
      # Select the top (gene) rows from the scaled counts matrix and calculate
      # z-scores per row to center the scale. Since base::scale() works on
      # columns, two transpositions are required.
      transformed_matrix <-
        SummarizedExperiment::assay(x = deseq_transform, i = 1L)[top_gene_identifiers, ]
      # Replace negative transformed count values with 0.
      # https://support.bioconductor.org/p/59369/
      transformed_matrix[transformed_matrix < 0] <- 1e-06
      # Add 1e-06 to the scaled counts for the log() function.

      complex_heatmap <- ComplexHeatmap::Heatmap(
        matrix = t(x = base::scale(
          x = t(x = log(x = transformed_matrix)),
          center = TRUE,
          scale = TRUE
        )),
        name = "z-score",
        row_title = "genes",
        # row_title_gp = gpar(fontsize = 7),
        column_title = "samples",
        cluster_rows = TRUE,
        show_row_dend = TRUE,
        cluster_columns = TRUE,
        show_column_dend = TRUE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 7),
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 7),
        top_annotation = ComplexHeatmap::columnAnnotation(df = column_annotation_frame)
      )
      rm(transformed_matrix)

      # Add column annotation.
      complex_heatmap <-
        complex_heatmap + ComplexHeatmap::HeatmapAnnotation(
          df = deseq_results_frame[top_gene_identifiers, c("gene_biotype", "significant"), drop = FALSE],
          which = "row",
          text = anno_text(
            x = deseq_results_frame[top_gene_identifiers, "gene_name", drop = TRUE],
            which = "row",
            gp = gpar(fontsize = 6),
            just = "left"
          )
        )

      file_path <-
        paste(if (is.null(x = file_index)) {
          # Without a file_index ...
          paste(prefix_heatmap,
                contrast_character,
                suffix,
                sep = "_")
        } else {
          # ... or with a file_index.
          paste(prefix_heatmap,
                contrast_character,
                suffix,
                file_index,
                sep = "_")
        },
        graphics_formats,
        sep = ".")

      # Draw a PDF plot (1L).
      grDevices::pdf(
        file = file.path(output_directory, file_path[1L]),
        width = argument_list$plot_width,
        height = argument_list$plot_height
      )
      if (is.null(x = plot_title)) {
        ComplexHeatmap::draw(object = complex_heatmap)
      } else {
        ComplexHeatmap::draw(object = complex_heatmap, column_title = plot_title)
      }
      base::invisible(x = grDevices::dev.off())

      # Draw a PNG plot (2L).
      grDevices::png(
        filename = file.path(output_directory, file_path[2L]),
        width = argument_list$plot_width,
        height = argument_list$plot_height,
        units = "in",
        res = 300L
      )
      if (is.null(x = plot_title)) {
        ComplexHeatmap::draw(object = complex_heatmap)
      } else {
        ComplexHeatmap::draw(object = complex_heatmap, column_title = plot_title)
      }
      base::invisible(x = grDevices::dev.off())

      nozzle_section <-
        Nozzle.R1::addTo(
          parent = nozzle_section,
          Nozzle.R1::newFigure(
            file = file_path[2L],
            "Heatmap for contrast ",
            Nozzle.R1::asStrong(contrast_tibble[contrast_index, "Label", drop = TRUE]),
            fileHighRes = file_path[1L]
          )
        )
      rm(complex_heatmap,
         file_path)
    }
    return(nozzle_section)
  }

for (contrast_index in seq_len(length.out = nrow(x = contrast_tibble))) {
  contrast_character <-
    bsfR::bsfrd_get_contrast_character(contrast_tibble = contrast_tibble, index = contrast_index)

  # Annotated Results Tibble ----------------------------------------------
  # Read the annotated results tibble with all genes for this contrast.

  deseq_results_tibble <-
    bsfR::bsfrd_read_result_tibble(
      genome_directory = argument_list$genome_directory,
      design_name = argument_list$design_name,
      contrast_tibble = contrast_tibble,
      index = contrast_index
    )
  if (is.null(x = deseq_results_tibble)) {
    rm(deseq_results_tibble, contrast_character)
    next()
  }

  # Coerce into a conventinal data.frame object and
  # reset the row names from the "gene_id" variable.
  deseq_results_frame <- as.data.frame(x = deseq_results_tibble)
  rm(deseq_results_tibble)
  row.names(x = deseq_results_frame) <-
    deseq_results_frame$gene_id

  # In case a gene_set_tibble is available, use it for filtering.
  if (is.null(x = plot_annotation_tibble)) {
    # No gene_set_tibble to select genes from.
    selected_gene_identifiers <-
      if (argument_list$maximum_number >= 0L) {
        deseq_results_frame[head(
          x = order(deseq_results_frame$max_rank, decreasing = FALSE),
          n = argument_list$maximum_number
        ), "gene_id", drop = TRUE]
      } else {
        deseq_results_frame[deseq_results_frame$significant == "yes", "gene_id", drop = TRUE]
      }
    nozzle_section_heatmaps <-
      draw_complex_heatmap(
        nozzle_section = nozzle_section_heatmaps,
        deseq_results_frame = deseq_results_frame,
        top_gene_identifiers = selected_gene_identifiers,
        contrast_character = contrast_character
      )
    rm(selected_gene_identifiers)
  } else {
    # A gene_set_tibble exists to select genes from.
    plot_names <- unique(plot_annotation_tibble$plot_name)
    for (plot_index in seq_along(along.with = plot_names)) {
      # Filter for plot_name values.
      selected_gene_identifiers <-
        dplyr::filter(.data = plot_annotation_tibble, plot_name == plot_names[plot_index])$gene_id
      nozzle_section_heatmaps <-
        draw_complex_heatmap(
          nozzle_section = nozzle_section_heatmaps,
          deseq_results_frame = deseq_results_frame,
          top_gene_identifiers = selected_gene_identifiers,
          contrast_character = contrast_character,
          plot_title = plot_names[plot_index],
          file_index = plot_index
        )
      rm(selected_gene_identifiers)
    }
    rm(plot_index, plot_names)
  }
  rm(contrast_character,
     deseq_results_frame)
}
rm(contrast_index)

nozzle_report <-
  Nozzle.R1::newCustomReport("Complex Heatmap Report", version = 0)

nozzle_report <-
  Nozzle.R1::addTo(parent = nozzle_report, nozzle_section_contrasts)

nozzle_report <-
  Nozzle.R1::addTo(parent = nozzle_report, nozzle_section_heatmaps)

Nozzle.R1::writeReport(report = nozzle_report,
                       filename = file.path(
                         output_directory,
                         paste(prefix_heatmap,
                               "report",
                               suffix,
                               sep = "_")
                       ))

rm(
  nozzle_report,
  nozzle_section_heatmaps,
  nozzle_section_contrasts,
  column_annotation_frame,
  output_directory,
  contrast_tibble,
  suffix,
  deseq_transform,
  deseq_data_set,
  prefix_heatmap,
  prefix_deseq,
  graphics_formats,
  argument_list,
  plot_annotation_tibble,
  draw_complex_heatmap
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessionInfo())
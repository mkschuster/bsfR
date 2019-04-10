#!/usr/bin/env Rscript
#
# BSF R script to draw a Heatmap from DESeq2 result tables.
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
      help = "Comma-separated list of variables",
      type = "character"
    ),
    make_option(
      opt_str = c("--maximum-number"),
      default = 400L,
      dest = "maximum_number",
      help = "Maximum number of genes to plot [400]",
      type = "integer"
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

suppressPackageStartupMessages(expr = library(package = "ComplexHeatmap"))
suppressPackageStartupMessages(expr = library(package = "DESeq2"))
suppressPackageStartupMessages(expr = library(package = "ggplot2"))
suppressPackageStartupMessages(expr = library(package = "stringr"))

# Save plots in the following formats.

# FIXME: The graphics_formats are not used for non ggplot2 plots.
graphics_formats <- c("pdf", "png")

prefix <-
  paste("rnaseq",
        "deseq",
        argument_list$design_name,
        sep = "_")

output_directory <-
  file.path(argument_list$output_directory, prefix)
if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

# DESeqDataSet ------------------------------------------------------------


# Load a pre-calculated DESeqDataSet object.
deseq_data_set <- NULL

file_path <-
  file.path(argument_list$genome_directory,
            prefix,
            paste0(prefix, "_deseq_data_set.Rdata"))
if (file.exists(file_path) &&
    file.info(file_path)$size > 0L) {
  message("Loading a DESeqDataSet object")
  load(file = file_path)
} else {
  stop(paste0(
    "Require a pre-calculated DESeqDataSet object in file: ",
    file_path
  ))
}
rm(file_path)


# DESeqTransform ----------------------------------------------------------
# Load a previously saved "blind" or "model" DESeqTransform object
# that serves as the base for the Heatmap.
deseq_transform <- NULL
suffix <- "model"

file_path <-
  file.path(argument_list$genome_directory,
            prefix,
            paste(
              paste(prefix, "deseq", "transform", suffix, sep = "_"),
              "Rdata",
              sep = "."
            ))
if (file.exists(file_path) &&
    file.info(file_path)$size > 0L) {
  message(paste("Loading a", suffix, "DESeqTransform object"))
  load(file = file_path)
} else {
  stop(paste0(
    "Require a pre-calculated DESeqTransform object in file: ",
    file_path
  ))
}
rm(file_path)


# Contrasts Frame ---------------------------------------------------------
# Read a data frame of contrasts with variables "Design", "Numerator", "Denominator" and "Label".
message("Loading contrast frame")
contrast_frame <-
  read.table(
    file = file.path(
      argument_list$genome_directory,
      prefix,
      paste(prefix, "contrasts.tsv", sep = "_")
    ),
    header = TRUE,
    sep = "\t",
    colClasses = c(
      "Design" = "character",
      "Numerator" = "character",
      "Denominator" = "character",
      "Label" = "character"
    ),
    stringsAsFactors = FALSE
  )
# Subset to the selected design.
contrast_frame <-
  contrast_frame[contrast_frame$Design == argument_list$design_name, ]
if (nrow(contrast_frame) == 0L) {
  stop("No design remaining after selection for design name.")
}
contrast_frame$Significant <- 0L

# Select column data variables to annotate in the heat map.
column_annotation_frame <-
  data.frame(colData(x = deseq_data_set)[, stringr::str_split(string = argument_list$variables, pattern = ",")[[1L]], drop = FALSE])

for (i in seq_len(length.out = nrow(x = contrast_frame))) {
  # The "contrast" option of the DESeq results() function expects a list of numerator and denominator.
  contrast_list <-
    list("numerator" = unlist(x = strsplit(x = contrast_frame[i, "Numerator"], split = ",")),
         "denominator" = unlist(x = strsplit(x = contrast_frame[i, "Denominator"], split = ",")))
  contrast_character <-
    paste(paste(contrast_list$numerator, collapse = "_"),
          "against",
          if (length(x = contrast_list$denominator) > 0L) {
            paste(contrast_list$denominator, collapse = "_")
          } else {
            "intercept"
          },
          sep = "_")
  
  # Annotated Results Frame -----------------------------------------------
  # Read the annotated data.frame with all genes.
  
  file_path <-
    file.path(argument_list$genome_directory,
              prefix,
              paste(
                paste(prefix,
                      "contrast",
                      contrast_character,
                      "genes",
                      sep = "_"),
                "tsv",
                sep = "."
              ))
  
  if (file.exists(file_path) &&
      file.info(file_path)$size > 0L) {
    message(paste0("Loading DESeqResults frame for ", contrast_character))
    
    deseq_results_frame <-
      read.table(
        file = file_path,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
      )
    # Reset the row names from the gene_id variable.
    row.names(x = deseq_results_frame) <-
      deseq_results_frame$gene_id
    
    # FIXME: Make this configurable so that either top-ranked genes or just the significant genes could be used.
    # Now try a heat map with the top 400 genes
    
    top_gene_identifiers <-
      deseq_results_frame[head(
        x = order(deseq_results_frame$max_rank, decreasing = FALSE),
        n = argument_list$maximum_number
      ), "gene_id"]
    
    # Draw a ComplexHeatmap.
    complex_heatmap <- ComplexHeatmap::Heatmap(
      # Select the top (gene) rows from the scaled counts matrix and
      # calculate z-scores per row to center the scale.
      # Since base::scale() works on columns, two transpositions are required.
      # The following publication uses log2(TPM + 1) before scaling.
      # http://www.nature.com.ez.srv.meduniwien.ac.at/articles/nature24489/figures/1
      matrix = t(x = base::scale(
        x = t(x = log(x = SummarizedExperiment::assay(
          x = deseq_transform, i = 1L
        )[top_gene_identifiers, ])),
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
    # Add column annotation.
    
    grDevices::pdf(
      file = file.path(output_directory,
                       paste(
                         paste(
                           prefix,
                           "contrast",
                           contrast_character,
                           "heatmap",
                           suffix,
                           sep = "_"
                         ),
                         "pdf",
                         sep = "."
                       )),
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ComplexHeatmap::draw(object = complex_heatmap)
    base::invisible(x = grDevices::dev.off())
    rm(complex_heatmap, top_gene_identifiers)
  } else {
    message(paste0("Missing DESeqResults data.frame for ",
                   contrast_character))
  }
  rm(contrast_character, contrast_list, deseq_results_frame, file_path)
}

rm(
  i,
  column_annotation_frame,
  output_directory,
  contrast_frame,
  suffix,
  deseq_transform,
  deseq_data_set,
  prefix,
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

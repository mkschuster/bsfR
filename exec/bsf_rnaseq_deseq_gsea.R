#!/usr/bin/env Rscript
#
# BSF R script to export normalised counts for the BROAD Gene Set Enrichment
# Analysis GSEA programme.
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
      )
    )
  ))

# Check the input.

if (is.null(x = argument_list$design_name)) {
  stop("Missing --design-name option")
}

suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
suppressPackageStartupMessages(expr = library(package = "tidyverse"))
suppressPackageStartupMessages(expr = library(package = "BiocVersion"))
suppressPackageStartupMessages(expr = library(package = "bsfR"))
suppressPackageStartupMessages(expr = library(package = "DESeq2"))

prefix_gsea <-
  paste("rnaseq", "deseq", "gsea", argument_list$design_name, sep = "_")

output_directory <-
  file.path(argument_list$output_directory, prefix_gsea)
if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

# Load a pre-calculated annotation tibble.
annotation_tibble <- bsfR::bsfrd_read_annotation_tibble(
  genome_directory = argument_list$genome_directory,
  design_name = argument_list$design_name,
  feature_types = "gene",
  gtf_file_path = argument_list$gtf_reference,
  genome = argument_list$genome_version,
  verbose = argument_list$verbose
)

# Load a pre-calculated DESeqDataSet object.
deseq_data_set <-
  bsfR::bsfrd_read_deseq_data_set(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name,
    verbose = argument_list$verbose
  )

# Export the DESeq2 normalised counts in Gene Cluster Text file format (*.gct).

file_path <- file.path(output_directory,
                       paste(
                         paste(prefix_gsea,
                               "counts",
                               "normalised",
                               sep = "_"),
                         "gct",
                         sep = "."
                       ))

if (argument_list$verbose) {
  message("Exporting GCT header ...")
}

# Write the header consisting of a first fixed line ("#1.2") and a second line
# with the number of features (rows) and the number of samples (columns).

readr::write_lines(x = c("#1.2", paste(dim(x = deseq_data_set), collapse = "\t")), file = file_path)

if (argument_list$verbose) {
  message("Exporting GCT normalised counts table ...")
}

# Extract the "gene_name" (gene symbol) and "gene_id" (Ensembl gene identifier)
# variables form the annotation tibble and left join to the normalised count
# matrix.

count_tibble <-
  dplyr::right_join(
    x = dplyr::select(.data = annotation_tibble, .data$gene_name, .data$gene_id),
    y = tibble::as_tibble(
      x = DESeq2::counts(object = deseq_data_set, normalized = TRUE),
      rownames = "gene_id"
    ),
    by = "gene_id"
  )

# Rename the "gene_name" and "gene_id" variables to "NAME" and "Description",
# respectively, to match the GCT file format.

count_tibble <-
  dplyr::rename(
    .data = count_tibble,
    "NAME" = .data$gene_name,
    "Description" = .data$gene_id
  )

# Append the tibble to the header, but make sure the header line gets also exported.

readr::write_tsv(
  x = count_tibble,
  file = file_path,
  append = TRUE,
  col_names = TRUE
)
rm(count_tibble,
   file_path,
   deseq_data_set,
   annotation_tibble)

# Read a contrast tibble with variables "Design", "Numerator", "Denominator" and "Label".
contrast_tibble <-
  bsfR::bsfrd_read_contrast_tibble(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name,
    summary = TRUE,
    verbose = argument_list$verbose
  )

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
      index = contrast_index,
      verbose = argument_list$verbose
    )

  if (is.null(x = deseq_results_tibble)) {
    warning("Could not read result tibble: ", contrast_character)
    rm(deseq_results_tibble, contrast_character)
    next()
  }

  if (argument_list$verbose) {
    message("Exporting ranked list: ", contrast_character)
  }

  readr::write_tsv(
    x = dplyr::select(.data = deseq_results_tibble, .data$gene_name, .data$log2FoldChange),
    file = file.path(output_directory, paste(
      paste(prefix_gsea, contrast_character, sep = "_"), "rnk", sep = "."
    )),
    col_names = FALSE
  )

  rm(deseq_results_tibble, contrast_character)
}

rm(contrast_index,
   contrast_tibble,
   output_directory,
   prefix_gsea,
   argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessioninfo::session_info())

#!/usr/bin/env Rscript
#
# BSF R script to extract results of a DESeq2 analysis.
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
        help = "design name",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--lfc-threshold"),
        default = 0.0,
        dest = "lfc_threshold",
        help = "Log-fold change threshold [0.0]",
        type = "numeric"
      ),
      optparse::make_option(
        opt_str = c("--padj-threshold"),
        default = 0.1,
        dest = "padj_threshold",
        help = "Adjusted p-value threshold [0.1]",
        type = "numeric"
      ),
      optparse::make_option(
        opt_str = c("--threads"),
        default = 1L,
        dest = "threads",
        help = "Number of parallel processing threads",
        type = "integer"
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
suppressPackageStartupMessages(expr = library(package = "BiocParallel"))
suppressPackageStartupMessages(expr = library(package = "DESeq2"))
suppressPackageStartupMessages(expr = library(package = "EnhancedVolcano"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

message("Processing design '", argument_list$design_name, "'")

# Set the number of parallel threads in the MulticoreParam instance.

BiocParallel::register(BPPARAM = MulticoreParam(workers = argument_list$threads))

# The working directory is the analyis genome directory.
# Create a new sub-directory for results if it does not exist.

prefix <-
  bsfR::bsfrd_get_prefix_deseq(design_name = argument_list$design_name)

output_directory <-
  file.path(argument_list$output_directory, prefix)
if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

# Plot Annotation ---------------------------------------------------------


# Load a pre-calculated annotation tibble.
annotation_tibble <- bsfR::bsfrd_read_annotation_tibble(
  genome_directory = argument_list$genome_directory,
  design_name = argument_list$design_name,
  verbose = argument_list$verbose
)

# DESeqDataSet ------------------------------------------------------------


# Load a pre-calculated DESeqDataSet object.
deseq_data_set <-
  bsfR::bsfrd_read_deseq_data_set(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name,
    verbose = argument_list$verbose
  )


# Contrasts Tibble --------------------------------------------------------


# Read a tibble of contrasts with variables "Design", "Numerator", "Denominator"
# and "Label".
contrast_tibble <-
  bsfR::bsfrd_read_contrast_tibble(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name,
    summary = FALSE,
    verbose = argument_list$verbose
  )
if (nrow(x = contrast_tibble) == 0L) {
  stop("No contrast remaining after selection for design name.")
}

# Add new variables for the summary.
contrast_tibble <-
  tibble::add_column(
    .data = contrast_tibble,
    Significant = 0L,
    SignificantUp = 0L,
    SignificantDown = 0L
  )

for (contrast_index in seq_len(length.out = nrow(x = contrast_tibble))) {
  contrast_list <-
    bsfR::bsfrd_get_contrast_list(contrast_tibble = contrast_tibble, index = contrast_index)

  contrast_character <-
    bsfR::bsfrd_get_contrast_character(contrast_tibble = contrast_tibble, index = contrast_index)

  # Results ---------------------------------------------------------------


  # Check for the significant genes table and if it exists already, read it to
  # get the number of significant genes for the contrasts summary tibble.
  file_path <-
    file.path(output_directory,
              paste(
                paste(prefix,
                      "contrast",
                      contrast_character,
                      "significant",
                      sep = "_"),
                "tsv",
                sep = "."
              ))

  if (file.exists(file_path) &&
      file.info(file_path)$size > 0L) {
    # Read a DESeqResults Tibble ------------------------------------------

    message("Skipping DESeqResults for ", contrast_character)

    deseq_significant_tibble <-
      bsfR::bsfrd_read_result_tibble(
        genome_directory = argument_list$genome_directory,
        design_name = argument_list$design_name,
        contrast_tibble = contrast_tibble,
        index = contrast_index,
        verbose = argument_list$verbose
      )

    # Record the number of significant genes.
    contrast_tibble[contrast_index, "Significant"] <-
      nrow(x = deseq_significant_tibble)

    contrast_tibble[contrast_index, "SignificantUp"] <-
      nrow(x = dplyr::filter(.data = deseq_significant_tibble,
                             .data$log2FoldChange > 0.0))

    contrast_tibble[contrast_index, "SignificantDown"] <-
      nrow(x = dplyr::filter(.data = deseq_significant_tibble,
                             .data$log2FoldChange < 0.0))

    rm(deseq_significant_tibble)
  } else {
    # Create a DESeqResults Tibble ----------------------------------------
    message("Creating DESeqResults for ", contrast_character)

    deseq_results_default <-
      DESeq2::results(
        object = deseq_data_set,
        contrast = contrast_list,
        lfcThreshold = argument_list$lfc_threshold,
        alpha = argument_list$padj_threshold,
        format = "DataFrame",
        # If option "tidy" is TRUE, a classical data.frame is returned.
        tidy = FALSE,
        parallel = TRUE
      )

    # Shrink log2-fold change values.
    #
    # NOTE: The "ashr" method shrinks log2-fold changes on the basis of the
    # DESeqResults object and can thus deal with model matrices not being full
    # rank. If a "res" option is provided for type "ashr", then both options,
    # "coef" and "contrast" are ignored. Type "ashr" also ignores option
    # "lfcThreshold". The "format" option "DataFrame" implies that a (modified)
    # DESeqResults object is returned where variable "stat" is missing, but
    # which extends "DataFrame".
    #
    # Method "normal" should no longer be used, while "apeglm" can apparently
    # deal with contrasts.

    message("Shrinking log2-fold changes for ", contrast_character)
    deseq_results_shrunk <- DESeq2::lfcShrink(
      dds = deseq_data_set,
      res = deseq_results_default,
      type = "ashr",
      format = "DataFrame",
      parallel = TRUE
    )

    # DESeqResults Tibble -------------------------------------------------


    # Coerce the DataFrame into a tibble and add lots of useful variables.
    deseq_results_tibble <-
      tibble::as_tibble(x = as.data.frame(x = deseq_results_shrunk))
    deseq_results_tibble <- dplyr::mutate(
      .data = deseq_results_tibble,
      "gene_id" = row.names(x = deseq_results_shrunk),
      # Calculate a factor indicating significance.
      "significant" = factor(
        x = if_else(
          condition = .data$padj <= argument_list$padj_threshold,
          true = "yes",
          false = "no",
          missing = "no"
        ),
        levels = c("no", "yes")
      ),
      # Calculate ranks for ...
      # (1) the effect size (log2FoldChange), ...
      "rank_log2_fold_change" = base::rank(
        x = -abs(x = .data$log2FoldChange),
        ties.method = c("min")
      ),
      # (2) the absolute level (baseMean) and ...
      "rank_base_mean" = base::rank(
        x = -.data$baseMean,
        ties.method = c("min")
      ),
      # (3) the statistical significance (padj).
      "rank_padj" = base::rank(x = .data$padj, ties.method = c("min")),
      # Calculate the maximum of the three ranks.
      "max_rank" = base::pmax(
        .data$rank_log2_fold_change,
        .data$rank_base_mean,
        .data$rank_padj
      )
    )
    # Left join with the reference transcriptome annotation tibble.
    deseq_results_tibble <- dplyr::left_join(x = annotation_tibble,
                                             y = deseq_results_tibble,
                                             by = "gene_id")

    readr::write_tsv(x = deseq_results_tibble,
                     path = file.path(output_directory,
                                      paste(
                                        paste(prefix,
                                              "contrast",
                                              contrast_character,
                                              "genes",
                                              sep = "_"),
                                        "tsv",
                                        sep = "."
                                      )))

    # MA Plot -------------------------------------------------------------


    # Create a MA plot.
    message("Creating a MA plot for ", contrast_character)
    plot_paths <-
      file.path(output_directory,
                paste(
                  paste(prefix,
                        "contrast",
                        contrast_character,
                        "ma",
                        sep = "_"),
                  graphics_formats,
                  sep = "."
                ))
    names(x = plot_paths) <- names(x = graphics_formats)

    grDevices::pdf(
      file = plot_paths["pdf"],
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    DESeq2::plotMA(object = deseq_results_shrunk)
    base::invisible(x = grDevices::dev.off())

    grDevices::png(
      filename = plot_paths["png"],
      width = argument_list$plot_width,
      height = argument_list$plot_height,
      units = "in",
      res = 300L
    )
    DESeq2::plotMA(object = deseq_results_shrunk)
    base::invisible(x = grDevices::dev.off())
    rm(plot_paths)

    # Enhanced Volcano Plot -----------------------------------------------


    plot_paths <- file.path(output_directory,
                            paste(
                              paste(prefix,
                                    "contrast",
                                    contrast_character,
                                    "volcano",
                                    sep = "_"),
                              graphics_formats,
                              sep = "."
                            ))
    message("Creating an EnhancedVolcano plot for ", contrast_character)

    # EnhancedVolcano needs coercing the tibble into a data.frame.
    deseq_results_frame <- as.data.frame(x = deseq_results_tibble)

    ggplot_object <- EnhancedVolcano::EnhancedVolcano(
      toptable = deseq_results_frame,
      lab = deseq_results_frame$gene_name,
      x = "log2FoldChange",
      y = "pvalue"
    )
    rm(deseq_results_frame)

    for (plot_path in plot_paths) {
      ggplot2::ggsave(
        filename = plot_path,
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height,
        limitsize = FALSE
      )
    }
    rm(plot_path,
       ggplot_object,
       plot_paths)

    # Significant DESeqResults Tibble -------------------------------------

    # Filter for significan genes last, as they are the criterion for running
    # DESeq2::results() above.
    deseq_significant_tibble <-
      dplyr::filter(.data = deseq_results_tibble, .data$padj <= argument_list$padj_threshold)

    readr::write_tsv(x = deseq_significant_tibble, path = file.path(output_directory,
                                                                    paste(
                                                                      paste(
                                                                        prefix,
                                                                        "contrast",
                                                                        contrast_character,
                                                                        "significant",
                                                                        sep = "_"
                                                                      ),
                                                                      "tsv",
                                                                      sep = "."
                                                                    )))

    # Record the number of significant genes.
    contrast_tibble[contrast_index, "Significant"] <-
      nrow(x = deseq_significant_tibble)

    contrast_tibble[contrast_index, "SignificantUp"] <-
      nrow(x = dplyr::filter(.data = deseq_significant_tibble,
                             .data$log2FoldChange > 0.0))

    contrast_tibble[contrast_index, "SignificantDown"] <-
      nrow(x = dplyr::filter(.data = deseq_significant_tibble,
                             .data$log2FoldChange < 0.0))

    rm(
      deseq_significant_tibble,
      deseq_results_tibble,
      deseq_results_shrunk,
      deseq_results_default
    )
  }
  rm(file_path,
     contrast_character,
     contrast_list)
}
rm(contrast_index)

# Write Summary Tibble ----------------------------------------------------


readr::write_tsv(x = contrast_tibble,
                 path = file.path(
                   output_directory,
                   paste(prefix, "contrasts", "summary.tsv", sep = "_")
                 ))

rm(
  contrast_tibble,
  deseq_data_set,
  annotation_tibble,
  output_directory,
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

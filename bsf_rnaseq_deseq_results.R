#!/usr/bin/env Rscript
#
# BSF R script to extract results of a DESeq2 analysis.
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
      help = "design name",
      type = "character"
    ),
    make_option(
      opt_str = c("--lfc-threshold"),
      default = 0.0,
      dest = "lfc_threshold",
      help = "Log-fold change threshold [0.0]",
      type = "numeric"
    ),
    make_option(
      opt_str = c("--padj-threshold"),
      default = 0.1,
      dest = "padj_threshold",
      help = "Adjusted p-value threshold [0.1]",
      type = "numeric"
    ),
    make_option(
      opt_str = c("--threads"),
      default = 1L,
      dest = "threads",
      help = "Number of parallel processing threads",
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

suppressPackageStartupMessages(expr = library(package = "tidyverse"))
suppressPackageStartupMessages(expr = library(package = "bsfR"))
suppressPackageStartupMessages(expr = library(package = "BiocParallel"))
suppressPackageStartupMessages(expr = library(package = "DESeq2"))
suppressPackageStartupMessages(expr = library(package = "EnhancedVolcano"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

message(paste0("Processing design '", argument_list$design_name, "'"))

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


# Load a pre-calculated annotation frame.
annotation_tibble <- bsfR::bsfrd_read_annotation_tibble(
  genome_directory = argument_list$genome_directory,
  design_name = argument_list$design_name
)

# DESeqDataSet ------------------------------------------------------------


# Load a pre-calculated DESeqDataSet object.
deseq_data_set <-
  bsfR::bsfrd_read_deseq_data_set(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name
  )


# Contrasts Frame ---------------------------------------------------------


# Read a data frame of contrasts with variables "Design", "Numerator", "Denominator" and "Label".
contrast_tibble <-
  bsfR::bsfrd_read_contrast_tibble(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name,
    summary = FALSE
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

  # Check for the significant genes table and if it exist already,
  # read it to get the number of significant genes for the summary data frame.
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
    message(paste0("Skipping DESeqResults for ", contrast_character))

    deseq_merge_significant <-
      read.table(
        file = file_path,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
      )

    # Record the number of significant genes.
    contrast_tibble[contrast_index, "Significant"] <-
      nrow(x = deseq_merge_significant)

    contrast_tibble[contrast_index, "SignificantUp"] <-
      nrow(x = deseq_merge_significant[deseq_merge_significant$log2FoldChange > 0.0, , drop = FALSE])

    contrast_tibble[contrast_index, "SignificantDown"] <-
      nrow(x = deseq_merge_significant[deseq_merge_significant$log2FoldChange < 0.0, , drop = FALSE])

    rm(deseq_merge_significant)
  } else {
    message(paste0("Creating DESeqResults for ", contrast_character))

    deseq_results_default <-
      DESeq2::results(
        object = deseq_data_set,
        contrast = contrast_list,
        lfcThreshold = argument_list$lfc_threshold,
        alpha = argument_list$padj_threshold,
        format = "DataFrame",
        tidy = FALSE,
        # If tidy is TRUE, a classical data.frame is returned.
        parallel = TRUE
      )

    # Run lfcShrink() if possible.
    deseq_results_shrunk <- NULL
    # NOTE: The "ashr" method shrinks log2-fold changes on the basis of the
    # results table and can thus deal with model matrices not being full rank.
    #
    # if (any(attr(x = deseq_data_set, which = "full_rank"))) {
    # The original DESEqDataSet object is full rank, so that log2-fold changes
    # can be shrunk. The any() function returns FALSE for NULL values.
    message(paste0("Shrinking log2-fold changes for ", contrast_character))
    deseq_results_shrunk <- DESeq2::lfcShrink(
      dds = deseq_data_set,
      # For "ashr", if "res" is provided, then "coef" and "contrast" are ignored.
      contrast = contrast_list,
      res = deseq_results_default,
      type = "ashr",
      lfcThreshold = argument_list$lfc_threshold,
      format = "DataFrame",
      parallel = TRUE
    )
    # }

    # MA Plot ---------------------------------------------------------------


    # Create a MA plot.
    message(paste0("Creating a MA plot for ", contrast_character))
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
    if (is.null(x = deseq_results_shrunk)) {
      DESeq2::plotMA(object = deseq_results_default)
    } else {
      DESeq2::plotMA(object = deseq_results_shrunk)
    }
    base::invisible(x = grDevices::dev.off())

    grDevices::png(
      filename = plot_paths["png"],
      width = argument_list$plot_width,
      height = argument_list$plot_height,
      units = "in",
      res = 300L
    )
    if (is.null(x = deseq_results_shrunk)) {
      DESeq2::plotMA(object = deseq_results_default)
    } else {
      DESeq2::plotMA(object = deseq_results_shrunk)
    }
    base::invisible(x = grDevices::dev.off())
    rm(plot_paths)

    # Enhanced Volcano plot -----------------------------------------------


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
    # To annotate gene symbols rather than Ensembl gene identifiers, the
    # DESeqResults DataFrame needs subsetting into a data.frame and merging with
    # the annotation tibble.
    message(paste0("Creating an EnhancedVolcano plot for ", contrast_character))
    deseq_results_frame <-
      if (is.null(x = deseq_results_shrunk))
        data.frame(
          gene_id = rownames(x = deseq_results_default),
          deseq_results_default[, c("baseMean",
                                    "log2FoldChange",
                                    "lfcSE",
                                    "stat",
                                    "pvalue",
                                    "padj"), drop = FALSE]
        )
    else
      data.frame(
        gene_id = rownames(x = deseq_results_shrunk),
        deseq_results_shrunk[, c("baseMean",
                                 "log2FoldChange",
                                 "lfcSE",
                                 "pvalue", # no "stat" column
                                 "padj"), drop = FALSE]
      )

    deseq_results_merged <-
      base::merge(x = annotation_tibble, y = deseq_results_frame, by = "gene_id")

    ggplot_object <- EnhancedVolcano::EnhancedVolcano(
      toptable = deseq_results_merged,
      lab = deseq_results_merged$gene_name,
      x = "log2FoldChange",
      y = "pvalue"
    )
    rm(deseq_results_merged, deseq_results_frame)

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
       ggplot_object, plot_paths)

    # DESeqResults DataFrame ----------------------------------------------


    # Adjust the DESeqResults DataFrame for merging with the annotation DataFrame,
    # by setting the rownames() as "gene_id" variable.
    deseq_results_frame <-
      S4Vectors::DataFrame(
        gene_id = rownames(x = deseq_results_default),
        deseq_results_default[, c("baseMean",
                                  "log2FoldChange",
                                  "lfcSE",
                                  "stat",
                                  "pvalue",
                                  "padj"), drop = FALSE],
        significant = factor(x = "no", levels = c("no", "yes"))
      )

    # Assign the "significant" factor on the basis of the adjusted p-value threshold.
    deseq_results_frame[!is.na(x = deseq_results_frame$padj) &
                          deseq_results_frame$padj <= argument_list$padj_threshold, "significant"] <-
      "yes"

    # Calculate ranks for ...
    # (1) the effect size (log2FoldChange), ...
    if (is.null(x = deseq_results_shrunk)) {
      deseq_results_frame$rank_log2_fold_change <-
        base::rank(
          x = -abs(x = deseq_results_frame$log2FoldChange),
          ties.method = c("min")
        )
    } else {
      deseq_results_frame$rank_log2_fold_change <-
        base::rank(
          x = -abs(x = deseq_results_shrunk$log2FoldChange),
          ties.method = c("min")
        )
    }

    # (2) the absolute level (baseMean) and ...
    deseq_results_frame$rank_base_mean <-
      base::rank(x = -deseq_results_frame$baseMean,
                 ties.method = c("min"))

    # (3) the statistical significance (padj).
    deseq_results_frame$rank_padj <-
      base::rank(x = deseq_results_frame$padj, ties.method = c("min"))

    # Calculate the maximum of the three ranks.
    deseq_results_frame$max_rank <-
      base::pmax(
        deseq_results_frame$rank_log2_fold_change,
        deseq_results_frame$rank_base_mean,
        deseq_results_frame$rank_padj
      )

    deseq_merge_complete <-
      merge(x = annotation_tibble, y = deseq_results_frame, by = "gene_id")

    write.table(
      x = deseq_merge_complete,
      file = file.path(output_directory,
                       paste(
                         paste(prefix,
                               "contrast",
                               contrast_character,
                               "genes",
                               sep = "_"),
                         "tsv",
                         sep = "."
                       )),
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE
    )

    # Significant DESeqResults DataFrame ------------------------------------


    deseq_merge_significant <-
      subset(x = deseq_merge_complete, padj <= argument_list$padj_threshold)

    # Record the number of significant genes.
    contrast_tibble[contrast_index, "Significant"] <-
      nrow(x = deseq_merge_significant)

    contrast_tibble[contrast_index, "SignificantUp"] <-
      nrow(x = deseq_merge_significant[deseq_merge_significant$log2FoldChange > 0.0, , drop = FALSE])

    contrast_tibble[contrast_index, "SignificantDown"] <-
      nrow(x = deseq_merge_significant[deseq_merge_significant$log2FoldChange < 0.0, , drop = FALSE])

    write.table(
      x = deseq_merge_significant,
      file = file.path(output_directory,
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
                       )),
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE
    )

    rm(
      deseq_merge_significant,
      deseq_merge_complete,
      deseq_results_frame,
      deseq_results_shrunk,
      deseq_results_default,
      contrast_character,
      contrast_list
    )
  }
  rm(file_path)
}
rm(contrast_index)

# Write summary frame -----------------------------------------------------


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

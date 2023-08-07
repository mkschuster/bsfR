#!/usr/bin/env Rscript
#
# Copyright 2013 - 2022 Michael K. Schuster
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

# Description -------------------------------------------------------------


# BSF R script to extract results of a DESeq2 analysis.

# Option Parsing ----------------------------------------------------------


suppressPackageStartupMessages(expr = library(package = "optparse"))

argument_list <-
  optparse::parse_args(object = optparse::OptionParser(
    option_list = list(
      optparse::make_option(
        opt_str = "--verbose",
        action = "store_true",
        default = TRUE,
        help = "Print extra output [default]",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = "--quiet",
        action = "store_false",
        default = FALSE,
        dest = "verbose",
        help = "Print little output",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = "--design-name",
        # default = "global",
        dest = "design_name",
        help = "design name",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--threads",
        default = 1L,
        dest = "threads",
        help = "Number of parallel processing threads",
        type = "integer"
      ),
      optparse::make_option(
        opt_str = "--genome-directory",
        default = ".",
        dest = "genome_directory",
        help = "Genome directory path [.]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--output-directory",
        default = ".",
        dest = "output_directory",
        help = "Output directory path [.]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--plot-width",
        default = 7.0,
        dest = "plot_width",
        help = "Plot width in inches [7.0]",
        type = "double"
      ),
      optparse::make_option(
        opt_str = "--plot-height",
        default = 7.0,
        dest = "plot_height",
        help = "Plot height in inches [7.0]",
        type = "double"
      )
    )
  ))

if (is.null(x = argument_list$design_name)) {
  stop("Missing --design-name option")
}

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "dplyr"))
suppressPackageStartupMessages(expr = library(package = "ggplot2"))
suppressPackageStartupMessages(expr = library(package = "readr"))
suppressPackageStartupMessages(expr = library(package = "tibble"))
# Bioconductor
suppressPackageStartupMessages(expr = library(package = "BiocVersion"))
suppressPackageStartupMessages(expr = library(package = "BiocParallel"))
suppressPackageStartupMessages(expr = library(package = "DESeq2"))
suppressPackageStartupMessages(expr = library(package = "EnhancedVolcano"))
suppressPackageStartupMessages(expr = library(package = "IHW"))
# BSF
suppressPackageStartupMessages(expr = library(package = "bsfR"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

message("Processing design '", argument_list$design_name, "'")

# Set the number of parallel threads in the MulticoreParam instance.

BiocParallel::register(BPPARAM = BiocParallel::MulticoreParam(workers = argument_list$threads))

# The working directory is the analysis genome directory.
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

# Design List -------------------------------------------------------------


design_list <- bsfR::bsfrd_read_design_list(
  genome_directory = argument_list$genome_directory,
  design_name = argument_list$design_name,
  verbose = argument_list$verbose
)

# Feature Annotation ------------------------------------------------------


# Load a pre-calculated feature annotation S4Vectors::DataFrame.
feature_dframe <- bsfR::bsfrd_initialise_feature_dframe(
  genome_directory = argument_list$genome_directory,
  design_name = argument_list$design_name,
  feature_types = "gene",
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
if (base::nrow(x = contrast_tibble) == 0L) {
  stop("No contrast remaining after selection for design name.")
}

# Add new variables to the summary.
contrast_tibble <-
  tibble::add_column(
    .data = contrast_tibble,
    "Significant" = 0L,
    "SignificantUp" = 0L,
    "SignificantDown" = 0L,
    "Effective" = 0L,
    "EffectiveUp" = 0L,
    "EffectiveDown" = 0L
  )

for (contrast_index in seq_len(length.out = base::nrow(x = contrast_tibble))) {
  contrast_list <-
    bsfR::bsfrd_get_contrast_list(contrast_tibble = contrast_tibble, index = contrast_index)

  contrast_character <-
    bsfR::bsfrd_get_contrast_character(contrast_tibble = contrast_tibble, index = contrast_index)

  # DESeqResults ----------------------------------------------------------


  deseq_results <- NULL
  file_path <-
    file.path(output_directory,
              paste(
                paste(prefix,
                      "contrast",
                      contrast_character,
                      "results",
                      sep = "_"),
                "rds",
                sep = "."
              ))

  if (file.exists(file_path) &&
      file.info(file_path)$size > 0L) {
    message("Reading DESeqResults for ", contrast_character)

    deseq_results <-
      readr::read_rds(file = file_path)
  } else {
    # Create DESeqResults -------------------------------------------------


    # Use independent hypothesis weighting from Bioconductor package IHW. The
    # IHW::hw.DESeqResults function sets the padj variable in the DESeqResults
    # object, while the IHW::ihwResult object is still available as meta data.
    #
    # Always use the default value 0.0 for lfcThreshold to test for the
    # log2-fold changes being 0.0 with the default altHypothesis greaterAbs.
    #
    # If option "tidy" is TRUE, a classical data.frame is returned.

    message("Creating DESeqResults for ", contrast_character)
    deseq_results_raw <-
      DESeq2::results(
        object = deseq_data_set,
        contrast = contrast_list,
        alpha = design_list$padj_threshold,
        filterFun = IHW::ihw,
        parallel = TRUE
      )

    # Shrink log2-fold change values --------------------------------------


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
    deseq_results <- DESeq2::lfcShrink(
      dds = deseq_data_set,
      res = deseq_results_raw,
      type = "ashr",
      parallel = TRUE
    )

    # Serialise the shrunken DESeqResults object to disk.
    # It still contains the IHW::ihwResult object as meta data.

    readr::write_rds(x = deseq_results,
                     file = file_path,
                     compress = "gz")

    rm(deseq_results_raw)
  }
  rm(file_path)

  # Record the number of significant features in the summary tibble regardless.

  contrast_tibble$Significant[contrast_index] <-
    sum(deseq_results$padj <= design_list$padj_threshold, na.rm = TRUE)

  contrast_tibble$SignificantUp[contrast_index] <-
    sum(
      deseq_results$padj <= design_list$padj_threshold &
        deseq_results$log2FoldChange > 0.0,
      na.rm = TRUE
    )

  contrast_tibble$SignificantDown[contrast_index] <-
    sum(
      deseq_results$padj <= design_list$padj_threshold &
        deseq_results$log2FoldChange < 0.0,
      na.rm = TRUE
    )

  # Record the number of effective features in the summary tibble regardless.

  contrast_tibble$Effective[contrast_index] <-
    sum(base::abs(x = deseq_results$log2FoldChange) >= design_list$l2fc_threshold,
        na.rm = TRUE)

  contrast_tibble$EffectiveUp[contrast_index] <-
    sum(
      base::abs(x = deseq_results$log2FoldChange) >= design_list$l2fc_threshold &
        deseq_results$log2FoldChange > 0.0,
      na.rm = TRUE
    )

  contrast_tibble$EffectiveDown[contrast_index] <-
    sum(
      base::abs(x = deseq_results$log2FoldChange) >= design_list$l2fc_threshold &
        deseq_results$log2FoldChange < 0.0,
      na.rm = TRUE
    )


  # IHW weights plots -----------------------------------------------------


  plot_paths <-
    file.path(output_directory,
              paste(
                paste(
                  prefix,
                  "contrast",
                  contrast_character,
                  "ihw",
                  "weights",
                  sep = "_"
                ),
                graphics_formats,
                sep = "."
              ))
  base::names(x = plot_paths) <- base::names(x = graphics_formats)

  if (all(file.exists(plot_paths)) &&
      all(file.info(plot_paths)$size > 0L)) {
    message("Skipping IHW weights plots for ", contrast_character)
  } else {
    message("Creating IHW weights plots for ", contrast_character)

    ggplot_object <-
      IHW::plot(x = S4Vectors::metadata(x = deseq_results)$ihwResult,
                what = "weights")

    ggplot_object <-
      ggplot_object + ggplot2::ggtitle(label = "Weights", subtitle = "Independent Hypothesis Weighting")

    for (plot_path in plot_paths) {
      ggplot2::ggsave(
        filename = plot_path,
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height,
        limitsize = FALSE
      )
    }
    rm(plot_path, ggplot_object)
  }
  rm(plot_paths)

  # IHW decision boundaries plots -----------------------------------------


  plot_paths <-
    file.path(output_directory,
              paste(
                paste(
                  prefix,
                  "contrast",
                  contrast_character,
                  "ihw",
                  "decision",
                  "boundaries",
                  sep = "_"
                ),
                graphics_formats,
                sep = "."
              ))
  base::names(x = plot_paths) <- base::names(x = graphics_formats)

  if (all(file.exists(plot_paths)) &&
      all(file.info(plot_paths)$size > 0L)) {
    message("Skipping IHW decision boundaries plots for ",
            contrast_character)
  } else {
    message("Creating IHW decision boundaries plots for ",
            contrast_character)

    ggplot_object <-
      IHW::plot(x = S4Vectors::metadata(x = deseq_results)$ihwResult,
                what = "decisionboundary")

    ggplot_object <-
      ggplot_object + ggplot2::ggtitle(label = "Decision Boundaries", subtitle = "Independent Hypothesis Weighting")

    for (plot_path in plot_paths) {
      ggplot2::ggsave(
        filename = plot_path,
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height,
        limitsize = FALSE
      )
    }
    rm(plot_path, ggplot_object)
  }
  rm(plot_paths)

  # IHW raw versus adjusted p-values plots --------------------------------


  plot_paths <-
    file.path(output_directory,
              paste(
                paste(
                  prefix,
                  "contrast",
                  contrast_character,
                  "ihw",
                  "pvalues",
                  sep = "_"
                ),
                graphics_formats,
                sep = "."
              ))
  base::names(x = plot_paths) <- base::names(x = graphics_formats)

  if (all(file.exists(plot_paths)) &&
      all(file.info(plot_paths)$size > 0L)) {
    message("Skipping IHW p-values plots for ",
            contrast_character)
  } else {
    message("Creating IHW p-values plots for ", contrast_character)

    ggplot_object <-
      ggplot2::ggplot(data = IHW::as.data.frame(x = S4Vectors::metadata(x = deseq_results)$ihwResult))

    ggplot_object <- ggplot_object +
      ggplot2::geom_point(
        mapping = ggplot2::aes(
          x = .data$pvalue,
          y = .data$adj_pvalue,
          colour = .data$group
        ),
        size = 0.25
      )

    ggplot_object <-
      ggplot_object + ggplot2::scale_colour_hue(c = 150, l = 70, drop = FALSE)

    ggplot_object <-
      ggplot_object + ggplot2::labs(
        x = "p-value",
        y = "adjusted p-value",
        title = "Raw versus Adjusted p-Values",
        subtitle = "Independent Hypothesis Weighting"
      )

    for (plot_path in plot_paths) {
      ggplot2::ggsave(
        filename = plot_path,
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height,
        limitsize = FALSE
      )
    }
    rm(plot_path, ggplot_object)
  }
  rm(plot_paths)

  # MA Plots --------------------------------------------------------------


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
  base::names(x = plot_paths) <- base::names(x = graphics_formats)

  if (all(file.exists(plot_paths)) &&
      all(file.info(plot_paths)$size > 0L)) {
    message("Skipping MA plots for ",
            contrast_character)
  } else {
    message("Creating MA plots for ", contrast_character)

    grDevices::pdf(
      file = plot_paths["pdf"],
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    DESeq2::plotMA(object = deseq_results)
    base::invisible(x = grDevices::dev.off())

    grDevices::png(
      filename = plot_paths["png"],
      width = argument_list$plot_width,
      height = argument_list$plot_height,
      units = "in",
      res = 300L
    )
    DESeq2::plotMA(object = deseq_results)
    base::invisible(x = grDevices::dev.off())
  }
  rm(plot_paths)

  # Annotated DESeqResults Tibble -----------------------------------------


  # Combine columns of the feature annotation DataFrame with the DESeqResults
  # DataFrame and coerce both into a Tibble.

  deseq_results_tibble <-
    tibble::as_tibble(x = S4Vectors::as.data.frame(x = S4Vectors::combineCols(x = feature_dframe, deseq_results)))

  # Add lots of useful variables.

  deseq_results_tibble <- dplyr::mutate(
    .data = deseq_results_tibble,
    # Calculate a logical indicating significance.
    "significant" = dplyr::if_else(
      condition = .data$padj <= .env$design_list$padj_threshold,
      true = TRUE,
      false = FALSE,
      missing = FALSE
    ),
    # Calculate a logical indicating effectiveness.
    "effective" = dplyr::if_else(
      condition = base::abs(x = .data$log2FoldChange) >= .env$design_list$l2fc_threshold,
      true = TRUE,
      false = FALSE,
      missing = FALSE
    ),
    # Calculate ranks for ...
    # (1) the effect size (log2FoldChange), ...
    "rank_log2_fold_change" = base::rank(
      x = -base::abs(x = .data$log2FoldChange),
      ties.method = c("min")
    ),
    # (2) the absolute level (baseMean) and ...
    "rank_base_mean" = base::rank(x = -.data$baseMean,
                                  ties.method = c("min")),
    # (3) the statistical significance (padj).
    "rank_padj" = base::rank(x = .data$padj, ties.method = c("min")),
    # Calculate the maximum of the three ranks.
    "max_rank" = base::pmax(
      .data$rank_log2_fold_change,
      .data$rank_base_mean,
      .data$rank_padj
    )
  )

  # The annotated results Tibble may contain a variable set of feature
  # annotation variables. Consequently, the readr::read_tsv() column types are
  # no longer predictable, so that the annotated results tibble is also
  # serialised to disk.
  readr::write_rds(x = deseq_results_tibble,
                   file = file.path(output_directory,
                                    paste(
                                      paste(prefix,
                                            "contrast",
                                            contrast_character,
                                            "differential",
                                            sep = "_"),
                                      "rds",
                                      sep = "."
                                    )))

  readr::write_tsv(x = deseq_results_tibble,
                   file = file.path(output_directory,
                                    paste(
                                      paste(prefix,
                                            "contrast",
                                            contrast_character,
                                            "differential",
                                            sep = "_"),
                                      "tsv",
                                      sep = "."
                                    )))

  # Significant DESeqResults Tibble -------------------------------------


  # Filter for significant features, which also removes observations with NA
  # values.
  readr::write_tsv(
    x = dplyr::filter(.data = deseq_results_tibble,
                      .data$padj <= .env$design_list$padj_threshold),
    file = file.path(output_directory,
                     paste(
                       paste(prefix,
                             "contrast",
                             contrast_character,
                             "significant",
                             sep = "_"),
                       "tsv",
                       sep = "."
                     ))
  )


  # Effective DESeqResults Tibble ---------------------------------------


  # Filter for effective features, which also removes observations with NA
  # values.
  readr::write_tsv(
    x = dplyr::filter(
      .data = deseq_results_tibble,
      base::abs(x = .data$log2FoldChange) >= .env$design_list$l2fc_threshold
    ),
    file = file.path(output_directory,
                     paste(
                       paste(prefix,
                             "contrast",
                             contrast_character,
                             "effective",
                             sep = "_"),
                       "tsv",
                       sep = "."
                     ))
  )

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
  base::names(x = plot_paths) <- base::names(x = graphics_formats)

  if (all(file.exists(plot_paths)) &&
      all(file.info(plot_paths)$size > 0L)) {
    message("Skipping EnhancedVolcano plots for ", contrast_character)
  } else {
    message("Creating EnhancedVolcano plots for ", contrast_character)

    # EnhancedVolcano needs the annotated tibble from above, but filtered for NA
    # values in the padj variable and coerced into a data.frame.
    deseq_results_frame <-
      base::as.data.frame(x = dplyr::filter(.data = deseq_results_tibble,!rlang::are_na(x = .data$padj)))

    ggplot_object <- EnhancedVolcano::EnhancedVolcano(
      toptable = deseq_results_frame,
      lab = deseq_results_frame$gene_name,
      x = "log2FoldChange",
      y = "padj",
      xlab = bquote(expr = ~ Log[2] ~ "fold change"),
      ylab = bquote(expr = ~ -Log[10] ~ adjusted ~ italic(P)),
      subtitle = ggplot2::waiver(),
      caption = ggplot2::waiver(),
      pCutoff = design_list$padj_threshold,
      FCcutoff = design_list$l2fc_threshold,
      legendLabels = c(
        "NS",
        expression(log[2] ~ FC),
        expression("adj." ~ italic(P)),
        expression("adj." ~ italic(P) ~ "and" ~ log[2] ~ FC)
      )
    )

    for (plot_path in plot_paths) {
      ggplot2::ggsave(
        filename = plot_path,
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height,
        limitsize = FALSE
      )
    }
    rm(plot_path, ggplot_object, deseq_results_frame)
  }
  rm(plot_paths)

  rm(deseq_results_tibble,
     deseq_results,
     contrast_character,
     contrast_list)
}
rm(contrast_index)

# Write Summary Tibble ----------------------------------------------------


readr::write_tsv(x = contrast_tibble,
                 file = file.path(
                   output_directory,
                   paste(prefix, "contrasts", "summary.tsv", sep = "_")
                 ))

rm(
  contrast_tibble,
  deseq_data_set,
  feature_dframe,
  design_list,
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

print(x = sessioninfo::session_info())

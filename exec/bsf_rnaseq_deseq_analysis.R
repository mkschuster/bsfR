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


# BSF R script to run a DESeq2 analysis.
#
# Reads are counted on the basis of a reference (Ensembl) transcriptome supplied
# as a GTF file and imported via rtracklayer::import() into exon
# GenomicRanges::GRanges objects. The exon GenomicRanges::GRanges are
# subsequently converted into a GenomicRanges::GRangesList object by (Ensembl)
# gene identifiers. A SummarizedExperiment::RangedSummarizedExperiment object is
# created by the GenomicAlignments::summarizeOverlaps() function. Reads in
# secondary alignments and reads failing vendor quality filtering are thereby
# dismissed.
#
# Sample Annotation DataFrame Description ---------------------------------
#
#
# See the bsfrd_read_sample_dframe() function for variable names.
#
# Design Annotation Tibble Description ------------------------------------
#
#
# See the bsfrd_read_design_tibble() function for variable names.
#
# Contrast Annotation Tibble Description ----------------------------------
#
#
# See the bsfrd_read_contrast_tibble() function for variable names.
#
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
        help = "Design name",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--gtf-reference",
        default = NULL,
        dest = "gtf_reference",
        help = "GTF file specifying a reference transcriptome",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--genome-version",
        default = NULL,
        dest = "genome_version",
        help = "Genome version",
        type = "character"
      ),
      optparse::make_option(
        # This option is required for Likelihood Ratio Testing (LRT)
        opt_str = "--padj-threshold",
        default = 0.1,
        dest = "padj_threshold",
        help = "Adjusted p-value threshold [0.1]",
        type = "numeric"
      ),
      optparse::make_option(
        # This option is required for PCA plots
        opt_str = "--pca-dimensions",
        default = 4L,
        dest = "pca_dimensions",
        help = "Principal components to plot [4]",
        type = "integer"
      ),
      optparse::make_option(
        # This option is required for PCA plots
        opt_str = "--pca-top-number",
        default = 500L,
        dest = "pca_top_number",
        help = "Number of most variable genes for PCA [500]",
        type = "integer"
      ),
      optparse::make_option(
        opt_str = "--threads",
        default = 1L,
        dest = "threads",
        help = "Number of parallel processing threads [1]",
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
      # optparse::make_option(
      #   opt_str = "--plot-factor",
      #   default = 0.5,
      #   dest = "plot_factor",
      #   help = "Plot width increase per 24 samples [0.5]",
      #   type = "numeric"
      # ),
      optparse::make_option(
        opt_str = "--plot-width",
        default = 7.0,
        dest = "plot_width",
        help = "Plot width in inches [7.0]",
        type = "numeric"
      ),
      optparse::make_option(
        opt_str = "--plot-height",
        default = 7.0,
        dest = "plot_height",
        help = "Plot height in inches [7.0]",
        type = "numeric"
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
suppressPackageStartupMessages(expr = library(package = "purrr"))
suppressPackageStartupMessages(expr = library(package = "readr"))
suppressPackageStartupMessages(expr = library(package = "stringr"))
suppressPackageStartupMessages(expr = library(package = "tibble"))
suppressPackageStartupMessages(expr = library(package = "tidyr"))
# CRAN
suppressPackageStartupMessages(expr = library(package = "caret"))
suppressPackageStartupMessages(expr = library(package = "pheatmap"))
# Bioconductor
suppressPackageStartupMessages(expr = library(package = "BiocVersion"))
suppressPackageStartupMessages(expr = library(package = "DESeq2"))
suppressPackageStartupMessages(expr = library(package = "BiocParallel"))
suppressPackageStartupMessages(expr = library(package = "GenomicAlignments"))
suppressPackageStartupMessages(expr = library(package = "RColorBrewer"))
suppressPackageStartupMessages(expr = library(package = "Rsamtools"))
suppressPackageStartupMessages(expr = library(package = "genefilter"))
suppressPackageStartupMessages(expr = library(package = "rtracklayer"))
# BSF
suppressPackageStartupMessages(expr = library(package = "bsfR"))

# Global Variables --------------------------------------------------------


# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

# Define a design-specific prefix.

prefix <-
  bsfR::bsfrd_get_prefix_deseq(design_name = argument_list$design_name)

# Define a design-specific output directory.

output_directory <-
  file.path(argument_list$output_directory, prefix)
if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

# Initialise a DESeqDataSet object ----------------------------------------


#' Fix a Model Matrix.
#'
#' Attempt to fix a model matrix by removing empty columns or by removing linear
#' combinations.
#'
#' @param model_matrix_local A model \code{matrix}.
#'
#' @return A model \code{matrix}.
#'
#' @examples
#' \dontrun{
#' model_matrix
#'   <- fix_model_matrix(model_matrix_local = model_matrix)
#' }
#' @noRd
fix_model_matrix <- function(model_matrix_local) {
  # Check, whether the model matrix is full rank.
  # This is based on the DESeq2::checkFullRank() function.
  full_rank <- TRUE

  if (qr(x = model_matrix_local)$rank < base::ncol(x = model_matrix_local)) {
    message("The model matrix is not full rank.")
    full_rank <- FALSE

    model_all_zero <-
      apply(
        X = model_matrix_local,
        MARGIN = 2,
        FUN = function(model_matrix_column) {
          return(all(model_matrix_column == 0))
        }
      )

    if (any(model_all_zero)) {
      message(
        "Levels or combinations of levels without any samples have resulted in\n",
        "column(s) of zeros in the model matrix:\n  ",
        paste(base::colnames(x = model_matrix_local)[model_all_zero], collapse = "\n  "),
        "\n",
        "Attempting to fix the model matrix by removing empty columns."
      )

      model_matrix_local <-
        model_matrix_local[, -which(x = model_all_zero)]
    } else {
      linear_combinations_list <-
        caret::findLinearCombos(x = model_matrix_local)

      message_character <-
        purrr::map_chr(
          .x = linear_combinations_list$linearCombos,
          .f = function(x) {
            return(paste0(
              "  Linear combinations:\n    ",
              paste(base::colnames(x = model_matrix_local)[x], collapse = "\n    "),
              "\n"
            ))
          }
        )

      message(
        "One or more variables or interaction terms in the design formula are\n",
        "linear combinations of the others.\n",
        message_character,
        "\n",
        "Attempting to fix the model by removing linear combinations:\n  ",
        paste(base::colnames(x = model_matrix_local)[linear_combinations_list$remove], collapse = "\n  ")
      )

      model_matrix_local <-
        model_matrix_local[,-linear_combinations_list$remove]
    }
    rm(model_all_zero)
  }

  return(list("model_matrix" = model_matrix_local, "full_rank" = full_rank))
}

#' Check a Model Matrix.
#'
#' Check a model matrix for being full rank and attempt to fix it by removing
#' empty columns or by removing linear combinations.
#'
#' @param model_matrix A model \code{matrix}.
#'
#' @return A named \code{list}.
#' \describe{
#' \item{model_matrix}{A model \code{matrix}.}
#' \item{formula_full_rank}{A \code{logical} indicating that the original
#' formula was full rank.}
#' }
#'
#' @seealso fix_model_matrix
#' @examples
#' \dontrun{
#' result_list <-
#'   check_model_matrix(model_matrix = model_matrix)
#' }
#' @noRd
check_model_matrix <- function(model_matrix) {
  formula_full_rank <- NULL

  matrix_full_rank <- FALSE

  while (!matrix_full_rank) {
    result_list <- fix_model_matrix(model_matrix_local = model_matrix)

    model_matrix <- result_list$model_matrix

    matrix_full_rank <- result_list$full_rank

    if (is.null(x = formula_full_rank)) {
      # Capture the intial state of the model matrix, which represents the formula.
      formula_full_rank <- result_list$full_rank
    }
  }

  # Return the model matrix and indicate whether the original formula was full
  # rank and the current model matrix is.
  return(
    list(
      "model_matrix" = model_matrix,
      "formula_full_rank" = formula_full_rank,
      "matrix_full_rank" = matrix_full_rank
    )
  )
}

#' Initialise a DESeqDataSet Object.
#'
#' Initialise or load a \code{DESeq2::DESeqDataSet} object.
#'
#' @param design_list A named \code{list} of design information.
#' @references argument_list
#' @references output_directory
#' @references prefix
#'
#' @return A \code{DESeq2::DESeqDataSet} object.
#'
#' @examples
#' \dontrun{
#' deseq_data_set <-
#'   initialise_deseq_data_set(design_list = design_list)
#' }
#' @noRd
initialise_deseq_data_set <- function(design_list) {
  deseq_data_set <- NULL

  file_path <-
    file.path(output_directory,
              paste0(prefix, "_deseq_data_set.rds"))
  if (file.exists(file_path) &&
      file.info(file_path)$size > 0L) {
    message("Loading a DESeqDataSet object")
    deseq_data_set <- base::readRDS(file = file_path)
  } else {
    ranged_summarized_experiment <-
      bsfR::bsfrd_initialise_rse(
        genome_directory = argument_list$genome_directory,
        design_list = design_list,
        gtf_path = argument_list$gtf_reference,
        genome_version = argument_list$genome_version,
        verbose = argument_list$verbose)

    # Create a model matrix based on the model formula and column (sample
    # annotation) data and check whether it is full rank.
    model_matrix <- stats::model.matrix.default(
      object = stats::as.formula(object = design_list$full_formula),
      data = SummarizedExperiment::colData(x = ranged_summarized_experiment)
    )

    result_list <-
      check_model_matrix(model_matrix = model_matrix)

    if (argument_list$verbose) {
      message("Writing initial model matrix")
      utils::write.table(
        x = base::as.data.frame(x = model_matrix),
        file = file.path(
          output_directory,
          paste0(prefix, "_model_matrix_full_initial.tsv")
        ),
        sep = "\t",
        row.names = TRUE,
        col.names = TRUE
      )

      message("Writing modified model matrix")
      utils::write.table(
        x = base::as.data.frame(x = result_list$model_matrix),
        file = file.path(
          output_directory,
          paste0(prefix, "_model_matrix_full_modified.tsv")
        ),
        sep = "\t",
        row.names = TRUE,
        col.names = TRUE
      )
    }
    rm(model_matrix)

    if (result_list$formula_full_rank) {
      # The design formula *is* full rank, so set it as "design" option directly.
      message("Creating a DESeqDataSet object with a model formula")
      deseq_data_set <-
        DESeq2::DESeqDataSet(se = ranged_summarized_experiment,
                             design = stats::as.formula(object = design_list$full_formula))

      # Start DESeq2 Wald testing.
      #
      # Set betaPrior = FALSE for consistent result names for designs,
      # regardless of interaction terms. DESeq2 seems to set betaPrior = FALSE
      # upon interaction terms, automatically.
      #
      # See: https://support.bioconductor.org/p/84366/
      #
      # betaPrior also has to be FALSE in case a user-supplied full model matrix is specified.
      message("Started DESeq Wald testing with a model formula")
      deseq_data_set <-
        DESeq2::DESeq(
          object = deseq_data_set,
          test = "Wald",
          betaPrior = FALSE,
          parallel = TRUE
        )
      attr(x = deseq_data_set, which = "full_rank") <- TRUE
    } else {
      # The original design formula was not full rank. Unfortunately, to
      # initialise the DESeqDataSet, a model matrix can apparently not be used
      # directly. Thus, use the simplest possible design (i.e. ~ 1) for
      # initialisation and perform Wald testing with the full model matrix.
      message("Creating a DESeqDataSet object with design formula ~ 1")
      deseq_data_set <-
        DESeq2::DESeqDataSet(se = ranged_summarized_experiment,
                             design = ~ 1)
      attr(x = deseq_data_set, which = "full_rank") <- FALSE

      message("Started DESeq Wald testing with a model matrix")
      deseq_data_set <-
        DESeq2::DESeq(
          object = deseq_data_set,
          test = "Wald",
          betaPrior = FALSE,
          full = result_list$model_matrix,
          parallel = TRUE
        )
    }

    base::saveRDS(object = deseq_data_set, file = file_path)

    rm(result_list, ranged_summarized_experiment)
  }
  rm(file_path)

  return(deseq_data_set)
}

# Initialise a DESeqTransform object --------------------------------------


#' Initialise a DESeqTransform Object.
#'
#' Initialise or load a \code{DESeq2::DESeqTransform} object.
#'
#' @param blind A \code{logical} scalar to create the
#'   \code{DESeq2::DESeqTransform} object blindly or based on the model.
#' @references output_directory
#' @references prefix
#'
#' @return A \code{DESeq2::DESeqTransform} object.
#'
#' @examples
#' \dontrun{
#' deseq_transform <-
#'   initialise_deseq_transform(deseq_data_set = deseq_data_set, blind = FALSE)
#' }
#' @noRd
initialise_deseq_transform <-
  function(deseq_data_set, blind = FALSE) {
    deseq_transform <- NULL

    suffix <- if (blind)
      "blind"
    else
      "model"

    file_path <-
      file.path(output_directory,
                paste(
                  paste(prefix, "deseq", "transform", suffix, sep = "_"),
                  "rds",
                  sep = "."
                ))
    if (file.exists(file_path) &&
        file.info(file_path)$size > 0L) {
      message("Loading a ", suffix, " DESeqTransform object")
      deseq_transform <- base::readRDS(file = file_path)
    } else {
      message("Creating a ", suffix, " DESeqTransform object")
      # Run variance stabilizing transformation (VST) to get homoscedastic data
      # for PCA plots.
      deseq_transform <-
        DESeq2::varianceStabilizingTransformation(object = deseq_data_set, blind = blind)

      base::saveRDS(object = deseq_transform, file = file_path)
    }
    rm(file_path, suffix)

    return(deseq_transform)
  }

# Plot FPKM Values --------------------------------------------------------


#' Create an FPKM Density Plot.
#'
#' Create a density plot of log10(FPKM) values and save PDF and PNG documents.
#'
#' @param object A \code{matrix} object returned by \code{DESeq2::fpkm}.
#' @references argument_list
#' @references graphics_formats
#' @references output_directory
#' @references prefix
#'
#' @return \code{NULL}
#'
#' @examples
#' \dontrun{
#' plot_fpkm_values(object = DESeq2::fpkm(object = dds))
#' }
#' @noRd
plot_fpkm_values <- function(object) {
  plot_paths <- file.path(output_directory,
                          paste(
                            paste(prefix,
                                  "fpkm",
                                  "density",
                                  sep = "_"),
                            graphics_formats,
                            sep = "."
                          ))

  if (all(file.exists(plot_paths) &&
          file.info(plot_paths)$size > 0L)) {
    message("Skipping a FPKM density plot")
  } else {
    message("Create a FPKM density plot")

    # Pivot the tibble to get just a "name" and a "value" variable.
    ggplot_object <-
      ggplot2::ggplot(data = tidyr::pivot_longer(
        data = tibble::as_tibble(x = object),
        cols = tidyselect::everything()
      ))

    ggplot_object <-
      ggplot_object +
      ggplot2::geom_density(
        mapping = ggplot2::aes(
          x = log10(x = .data$value),
          y = ggplot2::after_stat(x = .data$density),
          colour = .data$name
        ),
        alpha = I(1 / 3),
        na.rm = TRUE
      )

    ggplot_object <-
      ggplot_object +
      ggplot2::labs(
        x = "log10(FPKM)",
        y = "Density",
        colour = "Sample",
        title = "FPKM Density",
        subtitle = paste("Design", argument_list$design_name)
      )

    # Increase the plot width per 24 samples.
    # The number of samples is the number of columns of the matrix.
    # Rather than argument_list$plot_factor, a fixed number of 0.25 is used here.
    plot_width <-
      argument_list$plot_width + (ceiling(x = base::ncol(x = object) / 24L) - 1L) * argument_list$plot_width * 0.25

    for (plot_path in plot_paths) {
      ggplot2::ggsave(
        filename = plot_path,
        plot = ggplot_object,
        width = plot_width,
        height = argument_list$plot_height,
        limitsize = FALSE
      )
    }
    rm(plot_path, ggplot_object)
  }
  rm(plot_paths)
}

# Plot Cook's Distances ---------------------------------------------------


#' Create a Cook's Distances Box Plot.
#'
#' Create a box plot of Cook's distances and save PDF and PNG documents.
#'
#' @param object A \code{DESeq2::DESeqDataSet} object.
#' @references argument_list
#' @references graphics_formats
#' @references output_directory
#' @references prefix
#'
#' @return \code{NULL}
#'
#' @examples
#' \dontrun{
#' }
#' @noRd
plot_cooks_distances <- function(object) {
  plot_paths <-
    file.path(output_directory, paste(
      paste(prefix, "cooks", "distances", sep = "_"),
      graphics_formats,
      sep = "."
    ))

  if (all(file.exists(plot_paths) &&
          file.info(plot_paths)$size > 0L)) {
    message("Skipping a Cook's distances box plot")
  } else {
    message("Creating a Cook's distances box plot")

    ggplot_object <-
      ggplot2::ggplot(data = tidyr::pivot_longer(
        data = tibble::as_tibble(x = SummarizedExperiment::assays(x = object)$cooks),
        cols = tidyselect::everything()
      ))

    ggplot_object <-
      ggplot_object +
      ggplot2::geom_boxplot(mapping = ggplot2::aes(x = .data$name,
                                                   y = .data$value),
                            na.rm = TRUE)

    ggplot_object <-
      ggplot_object +
      ggplot2::labs(
        x = "Sample",
        y = "Cook's Distance",
        title = "Cook's Distance per Sample",
        subtitle = paste("Design", argument_list$design_name)
      )

    ggplot_object <-
      ggplot_object +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          size = ggplot2::rel(x = 0.8),
          hjust = 0.0,
          vjust = 0.5,
          angle = 90.0
        )
      )

    # Increase the plot width per 24 samples.
    # The number of samples is the number of rows of the
    # SummarizedExperiment::colData() S4Vectors::DataFrame. Rather than
    # argument_list$plot_factor, a fixed number of 0.33 is used here.
    plot_width <-
      argument_list$plot_width + (ceiling(x = S4Vectors::nrow(x = SummarizedExperiment::colData(x = object)) / 24L) - 1L) * argument_list$plot_width * 0.33

    for (plot_path in plot_paths) {
      ggplot2::ggsave(
        filename = plot_path,
        plot = ggplot_object,
        width = plot_width,
        height = argument_list$plot_height,
        limitsize = FALSE
      )
    }
    rm(plot_path, plot_width, ggplot_object)
  }
  rm(plot_paths)
}

# Plot RIN Scores ---------------------------------------------------------


#' Create a RIN Score Density Plot.
#'
#' Create a density plot of the "RIN" variable in the sample annotation
#' DataFrame. Add vertical lines at 1.0, 4.0 and 7.0 in red, yellow and green,
#' respectively, to indicate quality tiers. Save PDF and PNG documents.
#'
#' @param object A \code{DESeq2::DESeqDataSet} object.
#' @references argument_list
#' @references graphics_formats
#' @references output_directory
#' @references prefix
#'
#' @return \code{NULL}
#'
#' @examples
#' \dontrun{
#' plot_rin_scores(object = dds)
#' }
#' @noRd
plot_rin_scores <- function(object) {
  plot_paths <- file.path(output_directory,
                          paste(
                            paste(prefix,
                                  "rin",
                                  "density",
                                  sep = "_"),
                            graphics_formats,
                            sep = "."
                          ))
  if (all(file.exists(plot_paths) &&
          file.info(plot_paths)$size > 0L)) {
    message("Skipping a RIN score density plot")
  } else {
    if ("RIN" %in% S4Vectors::colnames(x = SummarizedExperiment::colData(x = object))) {
      message("Creating a RIN score density plot")
      ggplot_object <-
        ggplot2::ggplot(data = S4Vectors::as.data.frame(x = SummarizedExperiment::colData(x = object)))

      ggplot_object <-
        ggplot_object +
        ggplot2::xlim(RIN = c(0.0, 10.0))

      ggplot_object <-
        ggplot_object +
        ggplot2::geom_vline(xintercept = 1.0,
                            colour = "red",
                            linetype = 2L)

      ggplot_object <-
        ggplot_object +
        ggplot2::geom_vline(xintercept = 4.0,
                            colour = "yellow",
                            linetype = 2L)

      ggplot_object <-
        ggplot_object +
        ggplot2::geom_vline(xintercept = 7.0,
                            colour = "green",
                            linetype = 2L)

      ggplot_object <-
        ggplot_object +
        ggplot2::geom_density(mapping = ggplot2::aes(
          x = .data$RIN,
          y = ggplot2::after_stat(x = .data$density)
        ))

      ggplot_object <-
        ggplot_object +
        ggplot2::labs(
          x = "RIN score",
          y = "density",
          title = "RNA Integrity Number (RIN) Density Plot",
          subtitle = paste("Design", argument_list$design_name)
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
    } else {
      message("Skipping a RIN score density plot. No RIN variable in column annotation.")
    }
  }
  rm(plot_paths)
}

# Plot Multi-Dimensional Scaling (MDS) ------------------------------------


#' Create a Multi-Dimensional Scaling (MDS) Plot.
#'
#' Create an MDS plot based on Classical (Metric) Multi-Dimensional Scaling of
#' Euclidean distances of transformed counts for each gene. Save PDF and PNG
#' documents.
#'
#' @param object A \code{DESeq2::DESeqTransform} object.
#' @param plot_list A \code{list} of \code{list} objects configuring plots and
#'   their \code{ggplot2} aesthetic mappings.
#' @param blind A \code{logical} scalar to indicate a blind or model-based
#'   \code{DESeq2::DESeqTransform} object.
#' @references argument_list
#' @references graphics_formats
#' @references output_directory
#' @references prefix
#'
#' @return \code{NULL}
#'
#' @examples
#' \dontrun{
#' plot_mds(object = dds, plot_list = plot_list, blind = FALSE)
#' }
#' @noRd
plot_mds <- function(object,
                     plot_list = list(),
                     blind = FALSE) {
  suffix <- if (blind)
    "blind"
  else
    "model"

  message("Creating ", suffix, " MDS plots:")

  dist_object <-
    stats::dist(x = t(x = SummarizedExperiment::assay(x = object, i = 1L)))

  dist_matrix <- as.matrix(x = dist_object)

  # Convert the Multi-Dimensional Scaling matrix into a S4Vectors::DataFrame and
  # bind its columns to the sample annotation S4Vectors::DataFrame.
  mds_frame <-
    base::cbind(
      base::data.frame(stats::cmdscale(d = dist_matrix)),
      S4Vectors::as.data.frame(x = SummarizedExperiment::colData(x = object))
    )

  purrr::walk(
    .x = plot_list,
    .f = function(geom_list) {
      aes_character <-
        bsfR::bsfrd_geometrics_list_to_character(geom_list = geom_list)

      plot_paths <- file.path(output_directory,
                              paste(
                                paste(prefix,
                                      "mds",
                                      aes_character,
                                      suffix,
                                      sep = "_"),
                                graphics_formats,
                                sep = "."
                              ))

      if (all(file.exists(plot_paths) &&
              file.info(plot_paths)$size > 0L)) {
        message("  Skipping MDS plot: ", aes_character)
      } else {
        message("  Creating MDS plot: ", aes_character)

        ggplot_object <-
          ggplot2::ggplot(data = mds_frame)

        # geom_line
        if (!is.null(x = geom_list$geom_line)) {
          mapping_list <-
            ggplot2::aes(x = .data$X1, y = .data$X2)
          if (!is.null(x = geom_list$geom_line$colour)) {
            mapping_list <-
              utils::modifyList(x = mapping_list, val = ggplot2::aes(colour = .data[[geom_list$geom_line$colour]]))
          }
          if (!is.null(x = geom_list$geom_line$group)) {
            mapping_list <-
              utils::modifyList(x = mapping_list, val = ggplot2::aes(group = .data[[geom_list$geom_line$group]]))
          }
          if (!is.null(x = geom_list$geom_line$linetype)) {
            mapping_list <-
              utils::modifyList(x = mapping_list,
                                val = ggplot2::aes(linetype = .data[[geom_list$geom_line$linetype]]))
          }
          ggplot_object <-
            ggplot_object +
            ggplot2::geom_line(mapping = mapping_list,
                               alpha = I(1 / 3))
          rm(mapping_list)
        }

        # geom_point
        if (!is.null(x = geom_list$geom_point)) {
          mapping_list <-
            ggplot2::aes(x = .data$X1, y = .data$X2)
          if (!is.null(x = geom_list$geom_point$colour)) {
            mapping_list <-
              utils::modifyList(x = mapping_list, val = ggplot2::aes(colour = .data[[geom_list$geom_point$colour]]))
          }
          if (!is.null(x = geom_list$geom_point$shape)) {
            mapping_list <-
              utils::modifyList(x = mapping_list, val = ggplot2::aes(shape = .data[[geom_list$geom_point$shape]]))
          }
          ggplot_object <-
            ggplot_object +
            ggplot2::geom_point(mapping = mapping_list,
                                size = 2.0,
                                alpha = I(1 / 3))
          if (!is.null(x = geom_list$geom_point$shape)) {
            # For more than six shapes (scale_shape()), a manual scale
            # (scale_shape_manual()) needs setting up.
            # https://ggplot2.tidyverse.org/reference/scale_shape.html
            ggplot_object <-
              ggplot_object +
              ggplot2::scale_shape_manual(values = seq_len(length.out = nlevels(x = mds_frame[, geom_list$geom_point$shape])))
          }
          rm(mapping_list)
        }

        # geom_text
        if (!is.null(x = geom_list$geom_text)) {
          mapping_list <-
            ggplot2::aes(x = .data$X1, y = .data$X2)
          if (!is.null(x = geom_list$geom_text$label)) {
            mapping_list <-
              utils::modifyList(x = mapping_list, val = ggplot2::aes(label = .data[[geom_list$geom_text$label]]))
          }
          if (!is.null(x = geom_list$geom_text$colour)) {
            mapping_list <-
              utils::modifyList(x = mapping_list, val = ggplot2::aes(colour = .data[[geom_list$geom_text$colour]]))
          }
          ggplot_object <-
            ggplot_object +
            ggplot2::geom_text(mapping = mapping_list,
                               size = 2.0,
                               alpha = I(1 / 3))
          rm(mapping_list)
        }

        # geom_path
        if (!is.null(x = geom_list$geom_path)) {
          mapping_list <-
            ggplot2::aes(x = .data$X1, y = .data$X2)
          if (!is.null(x = geom_list$geom_path$colour)) {
            mapping_list <-
              utils::modifyList(x = mapping_list, val = ggplot2::aes(colour = .data[[geom_list$geom_path$colour]]))
          }
          if (!is.null(x = geom_list$geom_path$group)) {
            mapping_list <-
              utils::modifyList(x = mapping_list, val = ggplot2::aes(group = .data[[geom_list$geom_path$group]]))
          }
          if (!is.null(x = geom_list$geom_path$linetype)) {
            mapping_list <-
              utils::modifyList(x = mapping_list,
                                val = ggplot2::aes(linetype = .data[[geom_list$geom_path$linetype]]))
          }
          ggplot_object <-
            ggplot_object +
            ggplot2::geom_path(
              mapping = mapping_list,
              arrow = if (is.null(x = geom_list$geom_path$arrow))
                NULL
              else
                grid::arrow(
                  length = grid::unit(x = 0.08, units = "inches"),
                  type = "closed"
                )
            )
          rm(mapping_list)
        }

        ggplot_object <-
          ggplot_object +
          ggplot2::theme_bw() +
          ggplot2::coord_fixed()

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
      rm(plot_paths, aes_character)
    }
  )

  rm(mds_frame,
     dist_matrix,
     dist_object,
     suffix)
}


# Plot Heatmap ------------------------------------------------------------


#' Create a Heat Map Plot.
#'
#' Create a heat map plot based on hierarchical clustering of Euclidean
#' distances of transformed counts for each gene. Save PDF and PNG documents.
#'
#' @param object A \code{DESeq2::DESeqTransform} object.
#' @param plot_list A \code{list} of \code{list} objects configuring plots and
#'   their \code{ggplot2} aesthetic mappings.
#' @param blind A \code{logical} scalar to indicate a blind or model-based
#'   \code{DESeq2::DESeqTransform} object.
#' @references argument_list
#' @references graphics_formats
#' @references output_directory
#' @references prefix
#'
#' @return
#'
#' @examples
#' \dontrun{
#' plot_heatmap(object = dds, plot_list = plot_list, blind = FALSE)
#' }
#' @noRd
plot_heatmap <- function(object,
                         plot_list = list(),
                         blind = FALSE) {
  suffix <- if (blind)
    "blind"
  else
    "model"

  message("Creating ", suffix, " Heatmap plots:")

  # Transpose the counts table, since dist() works with columns and
  # assign the sample names as column and row names to the resulting matrix.
  dist_object <-
    dist(x = t(x = SummarizedExperiment::assay(x = object, i = 1L)))

  dist_matrix <- as.matrix(x = dist_object)

  base::colnames(x = dist_matrix) <- object$sample
  base::rownames(x = dist_matrix) <- object$sample

  purrr::walk(
    .x = plot_list,
    .f = function(geom_list) {
      aes_character <-
        bsfR::bsfrd_geometrics_list_to_character(geom_list = geom_list)
      message("  Creating heatmap plot: ", aes_character)

      aes_factors <-
        unique(x = unlist(x = geom_list, use.names = TRUE))

      plotting_frame <-
        S4Vectors::as.data.frame(x = SummarizedExperiment::colData(x = object))[, aes_factors, drop = FALSE]

      pheatmap_object <- NULL

      if (FALSE) {
        # Add the aes_factors together to create a new grouping factor
        group_factor <- if (length(x = aes_factors) > 1) {
          factor(x = apply(
            X = plotting_frame,
            MARGIN = 1,
            FUN = paste,
            collapse = " : "
          ))
        } else {
          SummarizedExperiment::colData(x = object)[[aes_factors]]
        }

        # Assign the grouping factor to the distance matrix row names.
        base::colnames(x = dist_matrix) <- NULL
        base::rownames(x = dist_matrix) <-
          paste(object$sample, group_factor, sep = " - ")

        pheatmap_object <-
          pheatmap::pheatmap(
            mat = dist_matrix,
            color = grDevices::colorRampPalette(colors = base::rev(x = RColorBrewer::brewer.pal(
              n = 9, name = "Blues"
            )))(255),
            clustering_distance_rows = dist_object,
            clustering_distance_cols = dist_object,
            fontsize_row = 6,
            silent = TRUE
          )

        rm(group_factor)
      } else {
        # Draw a heatmap with covariate column annotation.
        pheatmap_object <-
          pheatmap::pheatmap(
            mat = dist_matrix,
            color = grDevices::colorRampPalette(colors = base::rev(x = RColorBrewer::brewer.pal(
              n = 9, name = "Blues"
            )))(255),
            clustering_distance_rows = dist_object,
            clustering_distance_cols = dist_object,
            annotation_col = plotting_frame,
            show_rownames = TRUE,
            show_colnames = FALSE,
            fontsize_row = 6,
            silent = TRUE
          )
      }

      plot_paths <- file.path(output_directory,
                              paste(
                                paste(prefix,
                                      "heatmap",
                                      aes_character,
                                      suffix,
                                      sep = "_"),
                                graphics_formats,
                                sep = "."
                              ))
      base::names(x = plot_paths) <- names(x = graphics_formats)

      # PDF output
      grDevices::pdf(
        file = plot_paths["pdf"],
        width = argument_list$plot_width,
        height = argument_list$plot_height
      )
      grid::grid.draw(pheatmap_object$gtable)
      base::invisible(x = grDevices::dev.off())

      # PNG output
      grDevices::png(
        filename = plot_paths["png"],
        width = argument_list$plot_width,
        height = argument_list$plot_height,
        units = "in",
        res = 300L
      )
      grid::grid.draw(pheatmap_object$gtable)
      base::invisible(x = grDevices::dev.off())

      rm(plot_paths,
         pheatmap_object,
         plotting_frame,
         aes_factors,
         aes_character)
    }
  )

  rm(dist_matrix,
     dist_object,
     suffix)
}

# Plot Principal Component Analysis (PCA) ---------------------------------


#' Create a Principal Component Analysis (PCA) Plot.
#'
#' Create a principal component analysis plot based on the top-most variant rows
#' (i.e. genes) calculated via \code{genefilter::rowVars()}. Save PDF and PNG
#' documents.
#'
#' @param object A \code{DESeq2::DESeqTransform} object.
#' @param plot_list A \code{list} of \code{list} objects configuring plots and
#'   their \code{ggplot2} aesthetic mappings.
#' @param blind A \code{logical} scalar to indicate a blind or model-based
#'   \code{DESeq2::DESeqTransform} object.
#' @references argument_list
#' @references graphics_formats
#' @references output_directory
#' @references prefix
#'
#' @return
#'
#' @examples
#' \dontrun{
#' }
#' @noRd
plot_pca <- function(object,
                     plot_list = list(),
                     blind = FALSE) {
  suffix <- if (blind)
    "blind"
  else
    "model"

  message("Creating ", suffix, " PCA plots:")

  # Calculate the variance for each row (i.e., gene).
  row_variance <-
    genefilter::rowVars(x = SummarizedExperiment::assay(x = object, i = 1L))

  # Order by decreasing variance and select the top number of rows (i.e., genes).
  selected_rows <-
    order(row_variance, decreasing = TRUE)[seq_len(length.out = min(argument_list$pca_top_number, length(x = row_variance)))]

  # Perform a PCA on the (count) matrix returned by
  # SummarizedExperiment::assay() for the selected genes.
  pca_object <-
    stats::prcomp(x = t(x = SummarizedExperiment::assay(x = object, i = 1L)[selected_rows, ]))

  rm(selected_rows)

  # Plot the variance for a maximum of 100 components.
  plot_paths <- file.path(output_directory,
                          paste(
                            paste(prefix,
                                  "pca",
                                  "variance",
                                  suffix,
                                  sep = "_"),
                            graphics_formats,
                            sep = "."
                          ))

  plotting_tibble <-
    tibble::tibble(
      "component" = seq_along(along.with = pca_object$sdev),
      "variance" = pca_object$sdev ^ 2 / sum(pca_object$sdev ^ 2)
    )

  ggplot_object <-
    ggplot2::ggplot(data = plotting_tibble[seq_len(length.out = min(100L, length(x = pca_object$sdev))), , drop = FALSE])

  ggplot_object <-
    ggplot_object +
    ggplot2::geom_point(mapping = ggplot2::aes(x = .data$component, y = .data$variance))

  ggplot_object <-
    ggplot_object +
    ggplot2::labs(
      x = "Principal Component",
      y = "Variance",
      title = "Variance by Principal Component",
      subtitle = paste("Design", argument_list$design_name)
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
  rm(plot_path, ggplot_object, plotting_tibble, plot_paths)

  # The pca_object$x matrix has as many columns and rows as there are samples.
  pca_dimensions <-
    min(argument_list$pca_dimensions, base::ncol(x = pca_object$x))

  # Create combinations of all possible principal component pairs.
  pca_pair_matrix <-
    combn(x = seq_len(length.out = pca_dimensions), m = 2L)

  # Calculate the contribution to the total variance for each principal
  # component. Establish a label list with principal components and their
  # respective percentage of the total variance.
  label_list <-
    as.list(x = sprintf(
      fmt = "PC%i (%.3f%%)",
      seq_along(along.with = pca_object$sdev),
      100 * (pca_object$sdev ^ 2 / sum(pca_object$sdev ^ 2))
    ))

  base::names(x = label_list) <-
    paste0("PC", seq_along(along.with = pca_object$sdev))

  label_function <- function(value) {
    return(label_list[value])
  }

  purrr::walk(
    .x = plot_list,
    .f = function(geom_list) {
      aes_character <-
        bsfR::bsfrd_geometrics_list_to_character(geom_list = geom_list)
      message("  Creating PCA plot: ", aes_character)

      plot_paths <- file.path(output_directory,
                              paste(
                                paste(prefix,
                                      "pca",
                                      aes_character,
                                      suffix,
                                      sep = "_"),
                                graphics_formats,
                                sep = "."
                              ))

      # Assemble the data for the plot from the rotated data matrix.
      pca_frame <- base::as.data.frame(x = pca_object$x)

      plotting_frame <-
        base::data.frame(
          "component_1" = factor(levels = paste0(
            "PC", seq_len(length.out = pca_dimensions)
          )),
          "component_2" = factor(levels = paste0(
            "PC", seq_len(length.out = pca_dimensions)
          )),
          "x" = numeric(),
          "y" = numeric(),
          # Also initialise all variables of the column data, but do not
          # include any data (i.e., 0L rows).
          S4Vectors::as.data.frame(x = SummarizedExperiment::colData(x = object)[0L,])
        )

      for (column_number in seq_len(length.out = base::ncol(x = pca_pair_matrix))) {
        pca_label_1 <-
          paste0("PC", pca_pair_matrix[1L, column_number])

        pca_label_2 <-
          paste0("PC", pca_pair_matrix[2L, column_number])

        plotting_frame <- base::rbind(
          plotting_frame,
          base::data.frame(
            "component_1" = pca_label_1,
            "component_2" = pca_label_2,
            "x" = pca_frame[, pca_label_1],
            "y" = pca_frame[, pca_label_2],
            S4Vectors::as.data.frame(x = SummarizedExperiment::colData(x = object))
          )
        )

        rm(pca_label_1, pca_label_2)
      }
      rm(column_number)

      ggplot_object <- ggplot2::ggplot(data = plotting_frame)

      # geom_line
      if (!is.null(x = geom_list$geom_line)) {
        mapping_list <-
          ggplot2::aes(x = .data$x, y = .data$y)
        if (!is.null(x = geom_list$geom_line$colour)) {
          mapping_list <-
            utils::modifyList(x = mapping_list, val = ggplot2::aes(colour = .data[[geom_list$geom_line$colour]]))
        }
        if (!is.null(x = geom_list$geom_line$group)) {
          mapping_list <-
            utils::modifyList(x = mapping_list, val = ggplot2::aes(group = .data[[geom_list$geom_line$group]]))
        }
        ggplot_object <-
          ggplot_object +
          ggplot2::geom_line(mapping = mapping_list,
                             alpha = I(1 / 3))
        rm(mapping_list)
      }

      # geom_point
      if (!is.null(x = geom_list$geom_point)) {
        mapping_list <-
          ggplot2::aes(x = .data$x, y = .data$y)
        if (!is.null(x = geom_list$geom_point$colour)) {
          mapping_list <-
            utils::modifyList(x = mapping_list, val = ggplot2::aes(colour = .data[[geom_list$geom_point$colour]]))
        }
        if (!is.null(x = geom_list$geom_point$shape)) {
          mapping_list <-
            utils::modifyList(x = mapping_list, val = ggplot2::aes(shape = .data[[geom_list$geom_point$shape]]))
        }
        ggplot_object <-
          ggplot_object +
          ggplot2::geom_point(mapping = mapping_list,
                              size = 2.0,
                              alpha = I(1 / 3))
        if (!is.null(x = geom_list$geom_point$shape)) {
          # For more than six shapes (scale_shape()), a manual scale
          # (scale_shape_manual()) needs setting up.
          # https://ggplot2.tidyverse.org/reference/scale_shape.html
          ggplot_object <-
            ggplot_object +
            ggplot2::scale_shape_manual(values = seq_len(length.out = nlevels(x = plotting_frame[, geom_list$geom_point$shape])))
        }
        rm(mapping_list)
      }

      # geom_text
      if (!is.null(x = geom_list$geom_text)) {
        mapping_list <-
          ggplot2::aes(x = .data$x, y = .data$y)
        if (!is.null(x = geom_list$geom_text$label)) {
          mapping_list <-
            utils::modifyList(x = mapping_list, val = ggplot2::aes(label = .data[[geom_list$geom_text$label]]))
        }
        if (!is.null(x = geom_list$geom_text$colour)) {
          mapping_list <-
            utils::modifyList(x = mapping_list, val = ggplot2::aes(colour = .data[[geom_list$geom_text$colour]]))
        }
        ggplot_object <-
          ggplot_object +
          ggplot2::geom_text(mapping = mapping_list,
                             size = 2.0,
                             alpha = I(1 / 3))
        rm(mapping_list)
      }

      # geom_path
      if (!is.null(x = geom_list$geom_path)) {
        mapping_list <-
          ggplot2::aes(x = .data$x, y = .data$y)
        if (!is.null(x = geom_list$geom_path$colour)) {
          mapping_list <-
            utils::modifyList(x = mapping_list, val = ggplot2::aes(colour = .data[[geom_list$geom_path$colour]]))
        }
        if (!is.null(x = geom_list$geom_path$group)) {
          mapping_list <-
            utils::modifyList(x = mapping_list, val = ggplot2::aes(group = .data[[geom_list$geom_path$group]]))
        }
        if (!is.null(x = geom_list$geom_path$linetype)) {
          mapping_list <-
            utils::modifyList(x = mapping_list, val = ggplot2::aes(linetype = .data[[geom_list$geom_path$linetype]]))
        }
        ggplot_object <-
          ggplot_object +
          ggplot2::geom_path(
            mapping = mapping_list,
            arrow = grid::arrow(
              length = grid::unit(x = 0.08, units = "inches"),
              type = "closed"
            )
          )
        rm(mapping_list)
      }

      ggplot_object <-
        ggplot_object +
        ggplot2::facet_grid(
          rows = ggplot2::vars(.data$component_1),
          cols = ggplot2::vars(.data$component_2),
          labeller = ggplot2::labeller(component_1 = label_function, component_2 = label_function)
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

      if (argument_list$verbose) {
        # Write the PCA plot data frame.
        utils::write.table(
          x = plotting_frame,
          file = file.path(output_directory,
                           paste(
                             paste(prefix,
                                   "pca",
                                   aes_character,
                                   suffix,
                                   sep = "_"),
                             "tsv",
                             sep = "."
                           )),
          sep = "\t",
          row.names = FALSE,
          col.names = TRUE
        )
      }
      rm(pca_frame,
         plotting_frame,
         aes_character)
    }
  )

  rm(
    label_function,
    label_list,
    pca_pair_matrix,
    pca_object,
    row_variance,
    pca_dimensions,
    suffix
  )
}

# Start of main script ----------------------------------------------------


message("Processing design '", argument_list$design_name, "'")

# Set the number of parallel threads in the MulticoreParam instance.
BiocParallel::register(BPPARAM = BiocParallel::MulticoreParam(workers = argument_list$threads))

global_design_list <- bsfR::bsfrd_read_design_list(
  genome_directory = argument_list$genome_directory,
  design_name = argument_list$design_name,
  verbose = argument_list$verbose
)

annotation_tibble <- bsfR::bsfrd_read_annotation_tibble(
  genome_directory = argument_list$genome_directory,
  design_name = argument_list$design_name,
  feature_types = "gene",
  gtf_file_path = argument_list$gtf_reference,
  genome = argument_list$genome_version,
  verbose = argument_list$verbose
)

deseq_data_set <-
  initialise_deseq_data_set(design_list = global_design_list)

# TODO: Write the matrix of gene versus model coefficients as a data.frame to
# disk.

# Cooks Distances Plot ----------------------------------------------------

plot_cooks_distances(object = deseq_data_set)

# RIN Score Plot ----------------------------------------------------------


# If RIN scores are annotated in the sample frame, plot their distribution.
plot_rin_scores(object = deseq_data_set)

# Likelihood Ratio Test (LRT) ---------------------------------------------


#' Run an Likelihood Ratio Test (LRT) on a Reduced Formula.
#'
#' @param reduced_formula_character A \code{character} scalar of a reduced
#'   formula with a "reduced_name" attribute.
#' @references global_design_list
#' @references output_directory
#' @references prefix
#' @references global_design_list
#'
#' @return A \code{tbl_df} of LRT summary information.
#' @noRd
#'
#' @examples
lrt_reduced_formula_test <-
  function(reduced_formula_character) {
    lrt_tibble <- NULL

    # Skip NA or empty character vectors.
    if (is.na(x = reduced_formula_character) ||
        !base::nzchar(x = reduced_formula_character)) {
      # Return NULL instead of a tibble, which can still be processed by
      # purrr::map_dfr().
      return(lrt_tibble)
    }

    file_path_all <-
      file.path(output_directory,
                paste(paste(
                  prefix,
                  "lrt",
                  attr(x = reduced_formula_character, which = "reduced_name"),
                  sep = "_"
                ),
                "tsv",
                sep = "."))

    file_path_significant <-
      file.path(output_directory,
                paste(
                  paste(
                    prefix,
                    "lrt",
                    attr(x = reduced_formula_character, which = "reduced_name"),
                    "significant",
                    sep = "_"
                  ),
                  "tsv",
                  sep = "."
                ))

    if (file.exists(file_path_significant) &&
        file.info(file_path_significant)$size > 0L) {
      message(
        "Skipping reduced formula: ",
        attr(x = reduced_formula_character, which = "reduced_name")
      )

      # Read the existing table to count the number of significant genes after LRT.
      deseq_results_lrt_frame <-
        utils::read.table(file = file_path_significant,
                          header = TRUE,
                          sep = "\t")

      lrt_tibble <- tibble::tibble(
        "design" = global_design_list$design,
        "full_formula" = global_design_list$full_formula,
        "reduced_name" = attr(x = reduced_formula_character, which = "reduced_name"),
        "reduced_formula" = reduced_formula_character,
        "significant" = base::nrow(x = deseq_results_lrt_frame)
      )

      rm(deseq_results_lrt_frame)
    } else {
      message(
        "Processing reduced formula: ",
        attr(x = reduced_formula_character, which = "reduced_name")
      )

      # DESeq LRT requires either two model formulas or two model matrices.
      # Create a reduced model matrix and check whether it is full rank.
      formula_full <-
        stats::as.formula(object = global_design_list$full_formula)

      result_list_full <-
        check_model_matrix(
          model_matrix = stats::model.matrix.default(
            object = formula_full,
            data = SummarizedExperiment::colData(x = deseq_data_set)
          )
        )

      formula_reduced <-
        stats::as.formula(object = reduced_formula_character)

      model_matrix_reduced <-
        stats::model.matrix.default(object = formula_reduced,
                                    data = SummarizedExperiment::colData(x = deseq_data_set))

      result_list_reduced <-
        check_model_matrix(model_matrix = model_matrix_reduced)

      if (argument_list$verbose) {
        message("Writing initial reduced model matrix")

        utils::write.table(
          x = base::as.data.frame(x = model_matrix_reduced),
          file = file.path(
            output_directory,
            paste0(
              prefix,
              "_model_matrix_",
              attr(x = reduced_formula_character, which = "reduced_name"),
              "_initial.tsv"
            )
          ),
          sep = "\t",
          row.names = TRUE,
          col.names = TRUE
        )

        message("Writing modified reduced model matrix")

        utils::write.table(
          x = base::as.data.frame(x = result_list_reduced$model_matrix),
          file = file.path(
            output_directory,
            paste0(
              prefix,
              "_model_matrix_",
              attr(x = reduced_formula_character, which = "reduced_name"),
              "_modified.tsv"
            )
          ),
          sep = "\t",
          row.names = TRUE,
          col.names = TRUE
        )
      }

      full_rank <-
        result_list_full$formula_full_rank &
        result_list_reduced$formula_full_rank

      deseq_data_set_lrt <-
        DESeq2::DESeq(
          object = deseq_data_set,
          test = "LRT",
          full = if (full_rank)
            formula_full
          else
            result_list_full$model_matrix,
          reduced = if (full_rank)
            formula_reduced
          else
            result_list_reduced$model_matrix
        )

      rm(
        full_rank,
        result_list_reduced,
        result_list_full,
        formula_reduced,
        formula_full
      )
      # print(x = paste("DESeqDataSet LRT result names for", attr(x = reduced_formula_character, which = "reduced_name")))
      # print(x = DESeq2::resultsNames(object = deseq_data_set_lrt))

      # Convert the DESeqResults object that extends the S4Vectors::DataFrame
      # with additional meta data annotation into a tibble.
      deseq_results_lrt_tibble <-
        tibble::as_tibble(x = S4Vectors::as.data.frame(
          x = DESeq2::results(
            object = deseq_data_set_lrt,
            format = "DataFrame",
            # If tidy is TRUE, a classical data.frame is returned.
            tidy = FALSE,
            parallel = TRUE
          )
        ),
        rownames = "gene_id")

      # Set a "significant" variable to "yes" or "no".
      deseq_results_lrt_tibble <-
        dplyr::mutate(.data = deseq_results_lrt_tibble,
                      "significant" = factor(
                        x = dplyr::if_else(
                          condition = .data$padj <= .env$argument_list$padj_threshold,
                          true = "yes",
                          false = "no",
                          missing = "no"
                        ),
                        levels = c("no", "yes")
                      ))

      # Join with the annotation tibble.
      deseq_results_lrt_tibble <-
        dplyr::right_join(x = annotation_tibble,
                          y = deseq_results_lrt_tibble,
                          by = "gene_id")

      # Write all genes.
      readr::write_tsv(x = deseq_results_lrt_tibble,
                       file = file_path_all)

      # Filter only significant genes.
      deseq_results_lrt_tibble <-
        dplyr::filter(.data = deseq_results_lrt_tibble,
                      .data$padj <= .env$argument_list$padj_threshold)

      # Write only significant genes.
      readr::write_tsv(x = deseq_results_lrt_tibble,
                       file = file_path_significant)

      lrt_tibble <- tibble::tibble(
        "design" = global_design_list$design,
        "full_formula" = global_design_list$full_formula,
        "reduced_name" = attr(x = reduced_formula_character, which = "reduced_name"),
        "reduced_formula" = reduced_formula_character,
        "significant" = base::nrow(deseq_results_lrt_tibble)
      )

      rm(deseq_results_lrt_tibble,
         deseq_data_set_lrt)
    }
    rm(file_path_all, file_path_significant)

    return(lrt_tibble)
  }

#' Convert a Reduced Formula Character Vector into a Scalar.
#'
#' @param reduced_formula_character A two component \code{character} vector
#'   after splitting the reduced formula name [1L] and the reduced formula
#'   [2L].
#'
#' @return A \code{character} scalar of the reduced formula with a
#'   "reduced_name" attribute.
#' @noRd
#'
#' @examples
lrt_reduced_formula_scalar <-
  function(reduced_formula_character) {
    character_scalar <- reduced_formula_character[2L]

    attr(x = character_scalar, which = "reduced_name") <-
      reduced_formula_character[1L]

    return(character_scalar)
  }

# The "reduced_formulas" variable of the design list encodes one or more named
# reduced formulas for LRT (e.g., "name_1:~genotype + gender;name_2:~1"). Split
# reduced formulas on ";" characters, then select the formula [2L] and assign
# the name [1L] as "reduced_name" attribute. Finally run LRT tests on the
# reduced formulas and combine the resulting LRT tibbles into a LRT summary tibble.

lrt_summary_tibble <-
  purrr::map_dfr(
    .x = purrr::map(
      .x = stringr::str_split(
        string = stringr::str_split(
          string = global_design_list$reduced_formulas[1L],
          pattern = stringr::fixed(pattern = ";")
        )[[1L]],
        pattern = stringr::fixed(pattern = ":")
      ),
      .f = lrt_reduced_formula_scalar
    ),
    .f = lrt_reduced_formula_test
  )

# Write the LRT summary tibble to disk.
utils::write.table(
  x = lrt_summary_tibble,
  file = file.path(output_directory,
                   paste(
                     paste(prefix,
                           "lrt",
                           "summary",
                           sep = "_"),
                     "tsv",
                     sep = "."
                   )),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

rm(lrt_summary_tibble)


# Plot Aesthetics ---------------------------------------------------------


# The plot_aes variable of the design data frame supplies a pipe-separated list
# of plot definitions. Each plot definition contains a semi-colon-separated list
# of geometric objects and their associated aesthetics for each plot, which is a
# comma-separated list of aesthetics=variable mappings.
#
# geom_point:colour=test_1,shape=test_2;geom_line:colour=test_3,group=test_4|
# geom_point:colour=test_a,shape=test_b;geom_line:colour=test_c,group=test_d
#
# Convert into a list of list objects with variables and aesthetics as names.
plot_list <-
  bsfR::bsfrd_plots_character_to_list(plots_character = global_design_list$plot_aes)

# DESeqDataSet Results ----------------------------------------------------


print(x = "DESeqDataSet result names:")
print(x = DESeq2::resultsNames(object = deseq_data_set))

# Export RAW counts -------------------------------------------------------


# Export the raw counts from the DESeqDataSet object.
readr::write_tsv(
  x = dplyr::right_join(
    x = annotation_tibble,
    y = tibble::as_tibble(
      x = DESeq2::counts(object = deseq_data_set, normalized = FALSE),
      rownames = "gene_id"
    ),
    by = "gene_id"
  ),
  file = file.path(output_directory,
                   paste(
                     paste(prefix,
                           "counts",
                           "raw",
                           sep = "_"),
                     "tsv",
                     sep = "."
                   ))
)

# Export normalised counts ------------------------------------------------


# Export the normalised counts from the DESeq2::DESeqDataSet object.
readr::write_tsv(
  x = dplyr::right_join(
    x = annotation_tibble,
    y = tibble::as_tibble(
      x = DESeq2::counts(object = deseq_data_set, normalized = TRUE),
      rownames = "gene_id"
    ),
    by = "gene_id"
  ),
  file = file.path(output_directory,
                   paste(
                     paste(prefix,
                           "counts",
                           "normalised",
                           sep = "_"),
                     "tsv",
                     sep = "."
                   ))
)

# Export FPKM values ------------------------------------------------------


# Retrieve and plot FPKM values.
fpkm_matrix <- DESeq2::fpkm(object = deseq_data_set)
plot_fpkm_values(object = fpkm_matrix)

# Export FPKM values from the DESeq2::DESeqDataSet object.
readr::write_tsv(
  x = dplyr::right_join(
    x = annotation_tibble,
    y = tibble::as_tibble(x = fpkm_matrix, rownames = "gene_id"),
    by = "gene_id"
  ),
  file = file.path(output_directory,
                   paste(
                     paste(prefix,
                           "fpkms",
                           sep = "_"),
                     "tsv",
                     sep = "."
                   ))
)

rm(fpkm_matrix)

# DESeqTransform ----------------------------------------------------------
# MDS Plot ----------------------------------------------------------------
# PCA plot ----------------------------------------------------------------
# Heatmap Plot ------------------------------------------------------------
# Export VST counts -------------------------------------------------------


for (blind in c(FALSE, TRUE)) {
  suffix <- if (blind)
    "blind"
  else
    "model"

  deseq_transform <-
    initialise_deseq_transform(deseq_data_set = deseq_data_set, blind = blind)

  plot_mds(object = deseq_transform,
           plot_list = plot_list,
           blind = blind)

  plot_pca(object = deseq_transform,
           plot_list = plot_list,
           blind = blind)

  plot_heatmap(object = deseq_transform,
               plot_list = plot_list,
               blind = blind)

  # Export the vst counts from the DESeq2::DESeqTransform object
  readr::write_tsv(
    x = dplyr::right_join(
      x = annotation_tibble,
      y = tibble::as_tibble(
        x = SummarizedExperiment::assay(x = deseq_transform, i = 1L),
        rownames = "gene_id"
      ),
      by = "gene_id"
    ),
    file = file.path(output_directory,
                     paste(
                       paste(prefix,
                             "counts",
                             "vst",
                             suffix,
                             sep = "_"),
                       "tsv",
                       sep = "."
                     ))
  )

  rm(deseq_transform, suffix)
}
rm(blind, plot_list)

rm(
  annotation_tibble,
  deseq_data_set,
  global_design_list,
  output_directory,
  prefix,
  graphics_formats,
  argument_list,
  plot_pca,
  plot_heatmap,
  plot_mds,
  plot_rin_scores,
  plot_cooks_distances,
  plot_fpkm_values,
  lrt_reduced_formula_test,
  lrt_reduced_formula_scalar,
  initialise_deseq_transform,
  initialise_deseq_data_set,
  check_model_matrix,
  fix_model_matrix
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessioninfo::session_info())

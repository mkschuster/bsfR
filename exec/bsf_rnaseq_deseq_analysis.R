#!/usr/bin/env Rscript
#
# BSF R script to run a DESeq2 analysis.
#
# Reads are counted on the basis of a reference (Ensembl) transcriptome supplied
# as a GTF file and imported via rtracklayer::import() into exon GRanges
# objects. The exon GRanges are subsequently converted into a GRangesList object
# by (Ensembl) gene identifiers. A RangedSummarizedExpriment object is created
# by the GenomicAlignments::summarizeOverlaps() function. Reads in secondary
# alignments and reads failing vendor quality filtering are thereby dismissed.
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
#
# Sample Annotation DataFrame Description ---------------------------------
#
#
# The rnaseq_deseq_PREFIX_samples.tsv sample annotation DataFrame supports the
# following variables.
#
#   "bam_path":
#      A "character" vector of BAM file paths.
#
#   "bai_path":
#      A "character" vector of BAI file paths.
#
#   "sample":
#      The sample name.
#
#   "run":
#     Optional. If present, indicates that technical replicates should be
#     collapsed according to information in the "sample" variable. The "run"
#     variable provides the original sample name before collapsing technical
#     replicates.
#
#   "designs":
#      A "character" vector of comma-separated values of designs,
#      a particular sample should be part of.
#
#   "library_type":
#      A "factor" vector with levels "unstranded", "first" and "second" to
#      indicate the strandedness of the RNA-seq protocol and whether the first
#      or second strand gets sequenced. Illumina TruSeq standed mRNA sequences
#      the second strand so that reads need inverting before counting
#      strand-specifically.
#
#   "sequencing_type":
#      A "factor" vector with levels "SE" and "PE" indicating single-end or
#      paired-end sequencing, repsectivley and thus counting as read pairs or
#      not.
#
#   "total_counts":
#      An "integer" vector with total counts per sample. Calculated
#      automatically based on the colSums() of the counts() function.
#
#   "RIN":
#      A "numeric" vector providing the RNA integrity number (RIN) score per
#      sample. If available, the RIN score distribution will be plotted.
#
#

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
      opt_str = c("--gtf-reference"),
      default = NULL,
      dest = "gtf_reference",
      help = "GTF file specifying a reference transcriptome",
      type = "character"
    ),
    make_option(
      opt_str = c("--genome-version"),
      default = NULL,
      dest = "genome_version",
      help = "Genome version",
      type = "character"
    ),
    make_option(
      # This option is required for Likelihood Ratio Testing (LRT)
      opt_str = c("--padj-threshold"),
      default = 0.1,
      dest = "padj_threshold",
      help = "Adjusted p-value threshold [0.1]",
      type = "numeric"
    ),
    make_option(
      # This option is required for PCA plots
      opt_str = c("--pca-dimensions"),
      default = 4L,
      dest = "pca_dimensions",
      help = "Principal components to plot [4]",
      type = "integer"
    ),
    make_option(
      # This option is required for PCA plots
      opt_str = c("--pca-top-number"),
      default = 500L,
      dest = "pca_top_number",
      help = "Number of most variable genes for PCA [500]",
      type = "integer"
    ),
    make_option(
      opt_str = c("--threads"),
      default = 1L,
      dest = "threads",
      help = "Number of parallel processing threads [1]",
      type = "integer"
    ),
    # make_option(
    #   opt_str = c("--plot-factor"),
    #   default = 0.5,
    #   dest = "plot_factor",
    #   help = "Plot width increase per 24 samples [0.5]",
    #   type = "numeric"
    # ),
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

if (is.null(x = argument_list$design_name)) {
  stop("Missing --design-name option")
}

suppressPackageStartupMessages(expr = library(package = "DESeq2"))
suppressPackageStartupMessages(expr = library(package = "BiocParallel"))
suppressPackageStartupMessages(expr = library(package = "GenomicAlignments"))
suppressPackageStartupMessages(expr = library(package = "RColorBrewer"))
suppressPackageStartupMessages(expr = library(package = "Rsamtools"))
suppressPackageStartupMessages(expr = library(package = "caret"))
suppressPackageStartupMessages(expr = library(package = "plyr"))
suppressPackageStartupMessages(expr = library(package = "genefilter"))
suppressPackageStartupMessages(expr = library(package = "tidyverse"))
suppressPackageStartupMessages(expr = library(package = "grid"))
suppressPackageStartupMessages(expr = library(package = "pheatmap"))
suppressPackageStartupMessages(expr = library(package = "rtracklayer"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

prefix <-
  paste("rnaseq",
        "deseq",
        argument_list$design_name,
        sep = "_")

output_directory <- prefix

# Global variable for the Design list, assigned by the initialise_design_list() function.
# FIXME: Remove the global variable and pass the list into R functions().
global_design_list <- NULL

# Initialise Gene Annotation DataFrame object -----------------------------


#' Initialise or load a gene annotation DataFrame.
#'
#' @references argument_list
#' @references output_directory
#' @references prefix
#' @return A \code{DataFrame} object.
#'
#' @examples
initialise_annotation_frame <- function() {
  # Load pre-existing gene annotation data frame or create it from the reference GTF file.
  data_frame <- NULL

  file_path <-
    file.path(output_directory,
              paste(prefix, "annotation.tsv", sep = "_"))
  if (file.exists(file_path) &&
      file.info(file_path)$size > 0L) {
    message("Loading an annotation frame")
    data_frame <-
      read.table(
        file = file_path,
        header = TRUE,
        sep = "\t",
        comment.char = "",
        stringsAsFactors = FALSE
      )
  } else {
    # Extracting a list of gene names from the GrangesList object
    # inside the DESeqDataSet object seemingly takes forever.
    # Therefore, import the GTF file once more, but this time only the "gene" features.
    # gene_name_list <- lapply(X = rowRanges(x = deseq_data_set), FUN = function(x) { S4Vectors::mcols(x = x)[1L, "gene_name"] })
    message("Reading reference GTF gene features")
    gene_ranges <-
      rtracklayer::import(
        con = argument_list$gtf_reference,
        format = "gtf",
        genome = argument_list$genome_version,
        feature.type = "gene"
      )
    message("Creating annotation frame")
    data_frame <-
      S4Vectors::mcols(x = gene_ranges)[, c("gene_id",
                                            "gene_version",
                                            "gene_name",
                                            "gene_biotype",
                                            "gene_source")]
    # Add the location as an Ensembl-like location, lacking the coordinate system name and version.
    data_frame$location <-
      as(object = gene_ranges, Class = "character")
    rm(gene_ranges)
    write.table(
      x = data_frame,
      file = file_path,
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE
    )
  }
  rm(file_path)
  return(data_frame)
}

# Initialise Sample Annotation Data Frame object --------------------------


#' Initialise Sample Annotation
#'
#' @param factor_levels A \code{character} vector with a packed string to assign
#'   factor levels.
#' @references \code{argument_list$design_name}
#' @references \code{output_directory}
#' @references \code{prefix}
#' @return A \code{DataFrame} with sample annotation.
#'
#' @examples initialise_sample_frame(factor_levels="factor_1:level_1,level_2;factor_2:level_A,level_B")
initialise_sample_frame <- function(factor_levels) {
  # Read the BSF Python sample TSV file as a data.frame and convert into a DataFrame.
  # Import strings as factors and cast to character vectors where required.
  message("Loading sample DataFrame")
  data_frame <-
    as(
      object = read.table(
        file = file.path(output_directory, paste(prefix, 'samples.tsv', sep = '_')),
        header = TRUE,
        sep = "\t",
        comment.char = "",
        stringsAsFactors = TRUE
      ),
      "DataFrame"
    )
  rownames(x = data_frame) <- data_frame$sample

  # Select only those samples, which have the design name annotated in the designs variable.
  index_logical <-
    unlist(x = lapply(
      X = stringr::str_split(
        string = as.character(data_frame$designs),
        pattern = stringr::fixed(pattern = ",")
      ),
      FUN = function(character_1) {
        # character_1 is a character vector resulting from the split on ",".
        return(any(argument_list$design_name %in% character_1))
      }
    ))
  data_frame <- data_frame[index_logical,]
  rm(index_logical)

  if (nrow(x = data_frame) == 0L) {
    stop("No sample remaining after selection for design name.")
  }

  # The sequencing_type and library_type variables are required to set options
  # for the GenomicAlignments::summarizeOverlaps() read counting function.

  if (!any("sequencing_type" %in% names(x = data_frame))) {
    stop("A sequencing_type variable is missing from the sample annotation frame.")
  }

  if (!any("library_type" %in% names(x = data_frame))) {
    stop("A library_type variable is missing from the sample annotation frame.")
  }

  # Re-level the library_type and sequencing_type variables.
  data_frame$library_type <-
    factor(x = data_frame$library_type,
           levels = c("unstranded", "first", "second"))
  data_frame$sequencing_type <-
    factor(x = data_frame$sequencing_type,
           levels = c("SE", "PE"))

  # The 'factor_levels' variable of the design data frame specifies the order of factor levels.
  # Turn the factor_levels variable into a list of character vectors, where the factor names are set
  # as attributes of the list components.
  # factor_levels='factor_1:level_1,level_2;factor_2:level_A,level_B'
  factor_list <- lapply(
    # The "factor_levels" variable is a character vector, always with just a single component.
    # Split by ';' then by ':' character.
    X = stringr::str_split(
      string = stringr::str_split(
        string = factor_levels[1L],
        pattern = stringr::fixed(pattern = ";")
      )[[1L]],
      pattern = stringr::fixed(pattern = ":")
    ),
    FUN = function(character_1) {
      # Split the second component of character_1, the factor levels, on ",".
      character_2 <-
        unlist(x = stringr::str_split(
          string = character_1[2L],
          pattern = stringr::fixed(pattern = ",")
        ))
      # Set the first component of character_1, the factor name, as attribute.
      attr(x = character_2, which = "factor") <-
        character_1[1L]
      return(character_2)
    }
  )

  # Apply the factor levels to each factor.
  design_variables <- names(x = data_frame)
  for (i in seq_along(along.with = factor_list)) {
    factor_name <- attr(x = factor_list[[i]], which = "factor")
    if (factor_name != "") {
      if (factor_name %in% design_variables) {
        data_frame[, factor_name] <-
          factor(x = as.character(x = data_frame[, factor_name]),
                 levels = factor_list[[i]])
        # Check for NA values in case a factor level was missing.
        if (any(is.na(x = data_frame[, factor_name]))) {
          stop("Missing values after assigning factor levels for factor name ",
               factor_name)
        }
      } else {
        stop("Factor name ",
             factor_name,
             " does not resemble a variable of the design frame.")
      }
    }
    rm(factor_name)
  }
  rm(i, design_variables, factor_list)

  # Drop any unused levels from the sample data frame before retrning it.
  return(droplevels(x = data_frame))
}


# Initialise a Design list object -----------------------------------------

#' Initialise a named \code{list} for the selected design.
#'
#' @references \code{argument_list$design_name}
#' @references \code{output_directory}
#' @references \code{prefix}
#' @return A named \code{list} for the selected design.
#'
#' @examples
initialise_design_list <- function() {
  # Read the BSF Python design TSV file as a data.frame and convert into a DataFrame.
  message("Loading a design DataFrame")
  data_frame <-
    as(
      object = read.table(
        file = file.path(output_directory, paste(prefix, 'designs.tsv', sep = '_')),
        header = TRUE,
        sep = "\t",
        colClasses = c(
          "design" = "character",
          "full_formula" = "character",
          "reduced_formulas" = "character",
          "factor_levels" = "character",
          "plot_aes" = "character"
        ),
        fill = TRUE,
        comment.char = "",
        stringsAsFactors = FALSE
      ),
      Class = "DataFrame"
    )
  data_frame <-
    data_frame[data_frame$design == argument_list$design_name,]

  if (nrow(data_frame) == 0L) {
    stop("No design remaining after selection for design name.")
  }

  return(as.list(x = data_frame))
}


# Initialise a RangedSummarizedExperiment object --------------------------


#' Initialise or load a \code{RangedSummarizedExperiment} object.
#'
#' @param design_list A named \code{list} of design information.
#' @references argument_list
#' @references output_directory
#' @references prefix
#' @return A \code{RangedSummarizedExperiment} object.
#'
#' @examples
initialise_ranged_summarized_experiment <- function(design_list) {
  # Load a pre-existing RangedSummarizedExperiment object or create it by counting BAM files.
  ranged_summarized_experiment <- NULL

  file_path <-
    file.path(output_directory,
              paste0(prefix, "_ranged_summarized_experiment.Rdata"))
  if (file.exists(file_path) &&
      file.info(file_path)$size > 0L) {
    message("Loading a RangedSummarizedExperiment object")
    load(file = file_path)
  } else {
    sample_frame <-
      initialise_sample_frame(factor_levels = design_list$factor_levels)

    message("Reading reference GTF exon features")
    # The DESeq2 and RNA-seq vignettes suggest using TcDB objects, but for the moment,
    # we need extra annotation provided by Ensembl GTF files.
    exon_ranges <-
      rtracklayer::import(
        con = argument_list$gtf_reference,
        format = "gtf",
        genome = argument_list$genome_version,
        feature.type = "exon"
      )
    # Convert (i.e. split) the GRanges object into a GRangesList object
    # by gene identifiers.
    gene_ranges_list <-
      GenomicRanges::split(x = exon_ranges,
                           f = S4Vectors::mcols(x = exon_ranges)$gene_id)

    # Process per library_type and sequencing_type and merge the RangedSummarizedExperiment objects.
    ranged_summarized_experiment <- NULL

    for (library_type in levels(x = sample_frame$library_type)) {
      for (sequencing_type in levels(x = sample_frame$sequencing_type)) {
        message(
          "Processing library_type: ",
          library_type,
          " sequencing_type: ",
          sequencing_type
        )
        sub_sample_frame <-
          sample_frame[(sample_frame$library_type == library_type) &
                         (sample_frame$sequencing_type == sequencing_type), ]

        if (nrow(x = sub_sample_frame) == 0L) {
          rm(sub_sample_frame)
          next()
        }

        # Create a BamFileList object and set the samples as names.
        message("Creating a BamFileList object")
        bam_file_list <- Rsamtools::BamFileList(
          file = as.character(x = sub_sample_frame$bam_path),
          index = as.character(x = sub_sample_frame$bai_path),
          yieldSize = 2000000L,
          asMates = (sequencing_type == "PE")
        )
        # If a "run" variable is defined, technical replicates need collapsing
        # and the "sample" variable has duplicate values. Hence, use "run"
        # instead of "sample" for naming.
        if ("run" %in% names(x = sub_sample_frame)) {
          names(x = bam_file_list) <-
            as.character(x = sub_sample_frame$run)
        } else {
          names(x = bam_file_list) <-
            as.character(x = sub_sample_frame$sample)
        }

        message("Creating a RangedSummarizedExperiment object")
        sub_ranged_summarized_experiment <-
          GenomicAlignments::summarizeOverlaps(
            features = gene_ranges_list,
            reads = bam_file_list,
            mode = "Union",
            ignore.strand = (library_type == "unstranded"),
            # Exclude reads that represent secondary alignments or fail the vendor quality filter.
            param = Rsamtools::ScanBamParam(
              flag = Rsamtools::scanBamFlag(
                isSecondaryAlignment = FALSE,
                isNotPassingQualityControls = FALSE
              )
            ),
            # Invert the strand for protocols that sequence the second strand.
            preprocess.reads = if (library_type == "second")
              GenomicAlignments::invertStrand
          )
        SummarizedExperiment::colData(x = sub_ranged_summarized_experiment) <-
          sub_sample_frame
        # Combine RangedSummarizedExperiment objects with the same ranges,
        # but different samples via cbind().
        if (is.null(x = ranged_summarized_experiment)) {
          ranged_summarized_experiment <- sub_ranged_summarized_experiment
        } else {
          ranged_summarized_experiment <-
            cbind(ranged_summarized_experiment,
                  sub_ranged_summarized_experiment)
        }
        rm(sub_ranged_summarized_experiment,
           sub_sample_frame,
           bam_file_list)
      }
      rm(sequencing_type)
    }
    rm(library_type)

    # Collapse technical replicates if variable "run" is defined.
    sample_frame <-
      SummarizedExperiment::colData(x = ranged_summarized_experiment)
    if ("run" %in% names(x = sample_frame)) {
      message("Collapsing technical replicates.")
      # To avoid mismatching column and row names between
      # assay matrices and the column data annotation, variables
      # "sample" and "run" should be used. So set the original samples as runs
      # and rename the collapsed_sample variable into the sample variable.
      ranged_summarized_experiment <- DESeq2::collapseReplicates(
        object = ranged_summarized_experiment,
        groupby = sample_frame$sample,
        run = sample_frame$run,
        renameCols = TRUE
      )
    }
    rm(sample_frame)

    # Calculate colSums() of SummarizedExperiment::assays()$counts and add as
    # total_count into the colData data frame.
    sample_frame <-
      SummarizedExperiment::colData(x = ranged_summarized_experiment)
    sample_frame$total_counts <-
      base::colSums(
        x = SummarizedExperiment::assays(x = ranged_summarized_experiment)$counts,
        na.rm = TRUE
      )
    SummarizedExperiment::colData(x = ranged_summarized_experiment) <-
      sample_frame

    rm(gene_ranges_list,
       exon_ranges,
       sample_frame)
    save(ranged_summarized_experiment, file = file_path)
  }
  rm(file_path)

  return(ranged_summarized_experiment)
}

# Initialise a DESeqDataSet object ----------------------------------------


#' Fix a model matrix
#'
#' Attempt to fix a model matrix by removing empty columns or
#' by removing linear combinations.
#'
#' @param model_matrix_local A model \code{matrix}.
#'
#' @return A model \code{matrix}.
#'
#' @examples
fix_model_matrix <- function(model_matrix_local) {
  # Check, whether the model matrix is full rank.
  # This is based on the DESeq2::checkFullRank() function.
  full_rank <- TRUE
  if (qr(x = model_matrix_local)$rank < ncol(x = model_matrix_local)) {
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
        paste(colnames(x = model_matrix_local)[model_all_zero], collapse = "\n  "),
        "\n",
        "Attempting to fix the model matrix by removing empty columns."
      )
      model_matrix_local <-
        model_matrix_local[,-which(x = model_all_zero)]
    } else {
      linear_combinations_list <-
        caret::findLinearCombos(x = model_matrix_local)
      message_character <-
        unlist(x = lapply(
          X = linear_combinations_list$linearCombos,
          FUN = function(x) {
            return(paste0(
              "  Linear combinations:\n    ",
              paste(colnames(x = model_matrix_local)[x], collapse = "\n    "),
              "\n"
            ))
          }
        ))
      message(
        "One or more variables or interaction terms in the design formula are\n",
        "linear combinations of the others.\n",
        message_character,
        "\n",
        "Attempting to fix the model by removing linear combinations:\n  ",
        paste(colnames(x = model_matrix_local)[linear_combinations_list$remove], collapse = "\n  ")
      )
      model_matrix_local <-
        model_matrix_local[, -linear_combinations_list$remove]
    }
    rm(model_all_zero)
  }
  return(list("model_matrix" = model_matrix_local, "full_rank" = full_rank))
}

#' Check a model matrix for being full rank.
#'
#' @param model_matrix A model \code{matrix}.
#'
#' @return A named \code{list} of "model_matrix" and "formula_full_rank", a
#'   \code{locical} to indicate whether the original formula was already full
#'   rank.
#'
#' @examples
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

  # Return the model matrix and indicate whether the original formula was full rank and
  # the current model matrix is.
  return(
    list(
      "model_matrix" = model_matrix,
      "formula_full_rank" = formula_full_rank,
      "matrix_full_rank" = matrix_full_rank
    )
  )
}

#' Initialise or load a DESeqDataSet object.
#'
#' @param design_list A named \code{list} of design information.
#'
#' @return A \code{DESeqDataSet} object.
#'
#' @examples
initialise_deseq_data_set <- function(design_list) {
  deseq_data_set <- NULL

  file_path <-
    file.path(output_directory,
              paste0(prefix, "_deseq_data_set.Rdata"))
  if (file.exists(file_path) &&
      file.info(file_path)$size > 0L) {
    message("Loading a DESeqDataSet object")
    load(file = file_path)
  } else {
    ranged_summarized_experiment <-
      initialise_ranged_summarized_experiment(design_list = design_list)

    # Create a model matrix based on the model formula and column (sample annotation) data
    # and check whether it is full rank.
    model_matrix <- stats::model.matrix.default(
      object = as.formula(object = design_list$full_formula),
      data = SummarizedExperiment::colData(x = ranged_summarized_experiment)
    )
    result_list <-
      check_model_matrix(model_matrix = model_matrix)
    if (argument_list$verbose) {
      message("Writing initial model matrix")
      write.table(
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
      write.table(
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
                             design = as.formula(object = design_list$full_formula))
      # Start DESeq2 Wald testing.
      # Set betaPrior = FALSE for consistent result names for designs, regardless of interaction terms.
      # DESeq2 seems to set betaPrior = FALSE upon interaction terms, automatically.
      # See: https://support.bioconductor.org/p/84366/
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
      # The orignal design formula was not full rank.
      # Unfortunately, to initialise the DESeqDataSet,
      # a model matrix can apparently not be used directly.
      # Thus, use the simplest possible design (i.e. ~ 1) for initialisation and
      # peform Wald testing with the full model matrix.
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

    save(deseq_data_set, file = file_path)
    rm(result_list, ranged_summarized_experiment)
  }
  rm(file_path)

  return(deseq_data_set)
}

# Initialise a DESeqTransform object --------------------------------------


#' Initialise or load a DESeqTransform object.
#'
#' @param blind A \code{logical} scalar to create the DESeqTransform object
#'   blindly or based on the model.
#' @return A \code{DESeqTransform} object.
#'
#' @examples
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
                  "Rdata",
                  sep = "."
                ))
    if (file.exists(file_path) &&
        file.info(file_path)$size > 0L) {
      message("Loading a ", suffix, " DESeqTransform object")
      load(file = file_path)
    } else {
      message("Creating a ", suffix, " DESeqTransform object")
      # Run variance stabilizing transformation (VST) to get homoskedastic data for PCA plots.
      deseq_transform <-
        DESeq2::varianceStabilizingTransformation(object = deseq_data_set, blind = blind)
      save(deseq_transform, file = file_path)
    }
    rm(file_path, suffix)

    return(deseq_transform)
  }

# Convert aes_list into character -----------------------------------------


#' Convert the aes_list into a simple character string for file and plot naming.
#'
#' @param aes_list A \code{list} of aestethics.
#'
#' @return A \code{character} scalar.
#'
#' @examples
aes_list_to_character <- function(aes_list) {
  aes_character <-
    unlist(x = lapply(
      X = aes_list,
      FUN = function(list_1) {
        character_1 <- unlist(x = list_1)
        # Use paste() twice rather than c() to retain the order of name and value.
        return(paste(paste(
          names(x = character_1), character_1, sep = "_"
        ), collapse = "_"))
      }
    ))
  return(paste(c(
    sub(
      pattern = 'geom_',
      replacement = '',
      x = names(x = aes_character)
    ),
    aes_character
  ), collapse = "__"))
}


# Plot FPKM Values --------------------------------------------------------


#' Plot FPKM values.
#'
#' @param object A \code{matrix} object returned by \code{DESeq2::fpkm}.
#'
#' @return \code{NULL}
#'
#' @examples
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
    # Pivot the data frame to get just a "name" and a "value" variable.
    ggplot_object <-
      ggplot2::ggplot(data = tidyr::pivot_longer(
        data = tibble::as_tibble(x = object),
        cols = tidyselect::everything()
      ))

    ggplot_object <-
      ggplot_object + ggplot2::geom_density(
        mapping = ggplot2::aes(
          x = log10(.data$value),
          y = ..density..,
          colour = .data$name
        ),
        alpha = I(1 / 3)
      )

    ggplot_object <-
      ggplot_object + ggplot2::labs(
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
      argument_list$plot_width + (ceiling(x = ncol(x = object) / 24L) - 1L) * argument_list$plot_width * 0.25

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


#' Plot Cook's distances as box plot.
#'
#' @param object A \code{DESeqDataSet} object.
#' @return \code{NULL}
#'
#' @examples
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
      ggplot_object + ggplot2::geom_boxplot(mapping = ggplot2::aes(x = .data$name,
                                                                   y = .data$value))
    ggplot_object <-
      ggplot_object + ggplot2::labs(
        x = "Sample",
        y = "Cook's Distance",
        title = "Cook's Distance per Sample",
        subtitle = paste("Design", argument_list$design_name)
      )
    ggplot_object <-
      ggplot_object + ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 0,
        size = ggplot2::rel(x = 0.8)
      ))

    # Increase the plot width per 24 samples.
    # The number of samples is the number of rows of the colData() DataFrame.
    # Rather than argument_list$plot_factor, a fixed number of 0.33 is used here.
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


#' Plot RIN scores.
#'
#' @param object A \code{DESeqDataSet} object.
#'
#' @return \code{NULL}
#'
#' @examples
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
    if ("RIN" %in% colnames(x = SummarizedExperiment::colData(x = object))) {
      message("Creating a RIN score density plot")
      ggplot_object <-
        ggplot2::ggplot(data = BiocGenerics::as.data.frame(x = SummarizedExperiment::colData(x = object)))

      ggplot_object <-
        ggplot_object + ggplot2::xlim(RIN = c(0.0, 10.0))

      ggplot_object <-
        ggplot_object + ggplot2::geom_vline(xintercept = 1.0,
                                            colour = "red",
                                            linetype = 2L)

      ggplot_object <-
        ggplot_object + ggplot2::geom_vline(xintercept = 4.0,
                                            colour = "yellow",
                                            linetype = 2L)

      ggplot_object <-
        ggplot_object + ggplot2::geom_vline(xintercept = 7.0,
                                            colour = "green",
                                            linetype = 2L)

      ggplot_object <-
        ggplot_object + ggplot2::geom_density(mapping = ggplot2::aes(x = .data$RIN, y = ..density..))

      ggplot_object <-
        ggplot_object + ggplot2::labs(
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


#' Plot a Multi-Dimensional Scaling (MDS) analysis
#'
#' The MDS plot is based on Classical (Metric) Multi-Dimensional Scaling of
#'"Euclidean" distances of transformed counts for each gene.
#'
#' @param object DESeqTransform object
#' @param plot_list List of lists configuring plots and their ggplot2 aesthetic mappings
#' @param blind bool to indicate a blind or model-based DESeqTransform object
#'
#' @return \code{NULL}
#'
#' @examples
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
  # Convert the Mulitdimensional Scaling matrix into a DataFrame and
  # bind its columns to the sample annotation DataFrame.
  mds_frame <-
    base::cbind(
      base::data.frame(stats::cmdscale(d = dist_matrix)),
      BiocGenerics::as.data.frame(x = SummarizedExperiment::colData(x = object))
    )

  dummy_list <-
    lapply(
      X = plot_list,
      FUN = function(aes_list) {
        aes_character <-
          unique(x = unlist(x = aes_list, use.names = TRUE))

        plot_paths <- file.path(output_directory,
                                paste(
                                  paste(
                                    prefix,
                                    "mds",
                                    aes_list_to_character(aes_list = aes_list),
                                    suffix,
                                    sep = "_"
                                  ),
                                  graphics_formats,
                                  sep = "."
                                ))

        if (all(file.exists(plot_paths) &&
                file.info(plot_paths)$size > 0L)) {
          message(paste(c("  Skipping MDS plot:", aes_character), collapse = " "))
        } else {
          message(paste(c("  Creating MDS plot:", aes_character), collapse = " "))

          ggplot_object <-
            ggplot2::ggplot(data = mds_frame)

          # geom_line
          if (!is.null(x = aes_list$geom_line)) {
            ggplot_object <-
              ggplot_object +
              ggplot2::geom_line(
                mapping = ggplot2::aes_(
                  x = quote(expr = X1),
                  y = quote(expr = X2),
                  colour = if (is.null(x = aes_list$geom_line$colour))
                    NULL
                  else
                    as.name(x = aes_list$geom_line$colour),
                  group = if (is.null(x = aes_list$geom_line$group))
                    NULL
                  else
                    as.name(x = aes_list$geom_line$group),
                  linetype = if (is.null(x = aes_list$geom_path$linetype))
                    NULL
                  else
                    as.name(x = aes_list$geom_path$linetype)
                ),
                alpha = I(1 / 3)
              )
          }

          # geom_point
          if (!is.null(x = aes_list$geom_point)) {
            ggplot_object <- ggplot_object +
              ggplot2::geom_point(
                mapping = ggplot2::aes_(
                  x = quote(expr = X1),
                  y = quote(expr = X2),
                  colour = if (is.null(x = aes_list$geom_point$colour))
                    NULL
                  else
                    as.name(x = aes_list$geom_point$colour),
                  shape = if (is.null(x = aes_list$geom_point$shape))
                    NULL
                  else
                    as.name(x = aes_list$geom_point$shape)
                ),
                size = 2.0,
                alpha = I(1 / 3)
              )
            if (!is.null(x = aes_list$geom_point$shape)) {
              # For more than six shapes (scale_shape()), a manual scale
              # (scale_shape_manual()) needs setting up.
              # https://ggplot2.tidyverse.org/reference/scale_shape.html
              ggplot_object <-
                ggplot_object + ggplot2::scale_shape_manual(values = seq_len(length.out = nlevels(x = mds_frame[, aes_list$geom_point$shape])))
            }
          }

          # geom_text
          if (!is.null(x = aes_list$geom_text)) {
            ggplot_object <- ggplot_object +
              ggplot2::geom_text(
                mapping = ggplot2::aes_(
                  x = quote(expr = X1),
                  y = quote(expr = X2),
                  label = if (is.null(x = aes_list$geom_text$label))
                    "x"
                  else
                    as.name(x = aes_list$geom_text$label),
                  colour = if (is.null(x = aes_list$geom_text$colour))
                    NULL
                  else
                    as.name(x = aes_list$geom_text$colour)
                ),
                size = 2.0,
                alpha = I(1 / 3)
              )
          }

          # geom_path
          if (!is.null(x = aes_list$geom_path)) {
            ggplot_object <-
              ggplot_object + ggplot2::geom_path(
                mapping = ggplot2::aes_(
                  x = quote(expr = X1),
                  y = quote(expr = X2),
                  colour = if (is.null(x = aes_list$geom_path$colour))
                    NULL
                  else
                    as.name(x = aes_list$geom_path$colour),
                  group = if (is.null(x = aes_list$geom_path$group))
                    NULL
                  else
                    as.name(x = aes_list$geom_path$group),
                  linetype = if (is.null(x = aes_list$geom_path$linetype))
                    NULL
                  else
                    as.name(x = aes_list$geom_path$linetype)
                ),
                arrow = if (is.null(x = aes_list$geom_path$arrow))
                  NULL
                else
                  arrow(
                    length = unit(x = 0.08, units = "inches"),
                    type = "closed"
                  )
              )
          }

          ggplot_object <- ggplot_object +
            ggplot2::theme_bw() +
            coord_fixed()

          # ggplot_object <- ggplot_object + ggplot2::xlim(min(mds_frame$X1, mds_frame$X2), max(mds_frame$X1, mds_frame$X2))
          # ggplot_object <- ggplot_object + ggplot2::ylim(min(mds_frame$X1, mds_frame$X2), max(mds_frame$X1, mds_frame$X2))

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

  rm(dummy_list,
     mds_frame,
     dist_matrix,
     dist_object,
     suffix)
}


# Plot Heatmap ------------------------------------------------------------


#' Plot a heatmap
#'
#' The heatmap plot is based on hierarchical clustering of "Euclidean" distances
#' of transformed counts for each gene.
#'
#' @param object A \code{DESeqTransform} object.
#' @param plot_list A \code{list} of \code{list} objects configuring plots and
#'   their \code{ggplot2} aesthetic mappings.
#' @param blind A \code{logical} scalar to indicate a blind or model-based
#'   \code{DESeqTransform} object.
#'
#' @return
#'
#' @examples
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

  dummy_list <-
    lapply(
      X = plot_list,
      FUN = function(aes_list) {
        aes_character <-
          unique(x = unlist(x = aes_list, use.names = TRUE))
        message(paste(c(
          "  Creating heat map plot:", aes_character
        ), collapse = " "))

        plotting_frame <-
          BiocGenerics::as.data.frame(x = SummarizedExperiment::colData(x = object))[, aes_character, drop = FALSE]

        pheatmap_object <- NULL

        if (FALSE) {
          # Add the aes_character factors together to create a new grouping factor
          group_factor <- if (length(x = aes_character) > 1) {
            factor(x = apply(
              X = plotting_frame,
              MARGIN = 1,
              FUN = paste,
              collapse = " : "
            ))
          } else {
            SummarizedExperiment::colData(x = object)[[aes_character]]
          }

          # Assign the grouping factor to the distance matrix row names.
          base::colnames(x = dist_matrix) <- NULL
          base::rownames(x = dist_matrix) <-
            paste(object$sample, group_factor, sep = " - ")

          pheatmap_object <-
            pheatmap::pheatmap(
              mat = dist_matrix,
              color = grDevices::colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(
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
              color = grDevices::colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(
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
                                  paste(
                                    prefix,
                                    "heatmap",
                                    paste(aes_character, collapse = "_"),
                                    suffix,
                                    sep = "_"
                                  ),
                                  graphics_formats,
                                  sep = "."
                                ))
        names(x = plot_paths) <- names(x = graphics_formats)

        # PDF output
        grDevices::pdf(
          file = plot_paths["pdf"],
          width = argument_list$plot_width,
          height = argument_list$plot_height
        )
        # grid::grid.newpage()
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
        # grid::grid.newpage()
        grid::grid.draw(pheatmap_object$gtable)
        base::invisible(x = grDevices::dev.off())

        rm(plot_paths,
           pheatmap_object,
           plotting_frame,
           aes_character)
      }
    )

  rm(dummy_list,
     dist_matrix,
     dist_object,
     suffix)
}

# Plot Principal Component Analysis (PCA) ---------------------------------


#' Plot a Principal Component Analysis (PCA)
#'
#' @param object A \code{DESeqTransform} object.
#' @param plot_list A \code{list} of \code{list} objects configuring plots and
#'   their \code{ggplot2} aesthetic mappings.
#' @param blind A \code{logical} scalar to indicate a blind or model-based
#'   \code{DESeqTransform} object.
#'
#' @return
#'
#' @examples
plot_pca <- function(object,
                     plot_list = list(),
                     blind = FALSE) {
  suffix <- if (blind)
    "blind"
  else
    "model"

  message("Creating ", suffix, " PCA plots:")

  # Calculate the variance for each gene.
  row_variance <-
    genefilter::rowVars(x = SummarizedExperiment::assay(x = object, i = 1L))
  # Select the top number of genes by variance.
  selected_rows <-
    order(row_variance, decreasing = TRUE)[seq_len(length.out = min(argument_list$pca_top_number, length(x = row_variance)))]

  # Perform a PCA on the (count) matrix returned by SummarizedExperiment::assay() for the selected genes.
  pca_object <-
    stats::prcomp(x = t(x = SummarizedExperiment::assay(x = object, i = 1L)[selected_rows,]))
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
    ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(x = .data$component, y = .data$variance))
  ggplot_object <-
    ggplot_object + ggplot2::labs(
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
    min(argument_list$pca_dimensions, ncol(x = pca_object$x))
  # Create combinations of all possible principal component pairs.
  pca_pair_matrix <-
    combn(x = seq_len(length.out = pca_dimensions), m = 2L)

  # Calculate the contribution to the total variance for each principal component.
  # Establish a label list with principal components and their respective percentage of the total variance.
  label_list <-
    as.list(sprintf(
      fmt = "PC%i (%.3f%%)",
      seq_along(along.with = pca_object$sdev),
      100 * (pca_object$sdev ^ 2 / sum(pca_object$sdev ^ 2))
    ))
  names(x = label_list) <-
    paste0("PC", seq_along(along.with = pca_object$sdev))

  label_function <- function(value) {
    return(label_list[value])
  }

  dummy_list <-
    lapply(
      X = plot_list,
      FUN = function(aes_list) {
        aes_character <-
          unique(x = unlist(x = aes_list, use.names = TRUE))
        message(paste(c("  Creating PCA plot:", aes_character), collapse = " "))

        plot_paths <- file.path(output_directory,
                                paste(
                                  paste(
                                    prefix,
                                    "pca",
                                    aes_list_to_character(aes_list = aes_list),
                                    suffix,
                                    sep = "_"
                                  ),
                                  graphics_formats,
                                  sep = "."
                                ))

        # Assemble the data for the plot from the rotated data matrix.
        pca_frame <- base::as.data.frame(x = pca_object$x)
        plotting_frame <-
          base::data.frame(
            component_1 = factor(levels = paste0(
              "PC", seq_len(length.out = pca_dimensions)
            )),
            component_2 = factor(levels = paste0(
              "PC", seq_len(length.out = pca_dimensions)
            )),
            x = numeric(),
            y = numeric(),
            # Also initalise all variables of the column data, but do not include data (i.e. 0L rows).
            BiocGenerics::as.data.frame(x = SummarizedExperiment::colData(x = object)[0L, ])
          )

        for (column_number in seq_len(length.out = ncol(x = pca_pair_matrix))) {
          pca_label_1 <-
            paste0("PC", pca_pair_matrix[1L, column_number])
          pca_label_2 <-
            paste0("PC", pca_pair_matrix[2L, column_number])
          plotting_frame <- base::rbind(
            plotting_frame,
            base::data.frame(
              component_1 = pca_label_1,
              component_2 = pca_label_2,
              x = pca_frame[, pca_label_1],
              y = pca_frame[, pca_label_2],
              BiocGenerics::as.data.frame(x = SummarizedExperiment::colData(x = object))
            )
          )
          rm(pca_label_1, pca_label_2)
        }
        rm(column_number)

        ggplot_object <- ggplot2::ggplot(data = plotting_frame)

        # geom_line
        if (!is.null(x = aes_list$geom_line)) {
          ggplot_object <-
            ggplot_object +
            ggplot2::geom_line(
              mapping = ggplot2::aes_(
                x = quote(expr = x),
                y = quote(expr = y),
                colour = if (is.null(x = aes_list$geom_line$colour))
                  NULL
                else
                  as.name(x = aes_list$geom_line$colour),
                group = if (is.null(x = aes_list$geom_line$group))
                  NULL
                else
                  as.name(x = aes_list$geom_line$group)
              ),
              alpha = I(1 / 3)
            )
        }

        # geom_point
        if (!is.null(x = aes_list$geom_point)) {
          ggplot_object <- ggplot_object +
            ggplot2::geom_point(
              mapping = ggplot2::aes_(
                x = quote(expr = x),
                y = quote(expr = y),
                colour = if (is.null(x = aes_list$geom_point$colour))
                  NULL
                else
                  as.name(x = aes_list$geom_point$colour),
                shape = if (is.null(x = aes_list$geom_point$shape))
                  NULL
                else
                  as.name(x = aes_list$geom_point$shape)
              ),
              size = 2.0,
              alpha = I(1 / 3)
            )
          if (!is.null(x = aes_list$geom_point$shape)) {
            # For more than six shapes (scale_shape()), a manual scale
            # (scale_shape_manual()) needs setting up.
            # https://ggplot2.tidyverse.org/reference/scale_shape.html
            ggplot_object <-
              ggplot_object + ggplot2::scale_shape_manual(values = seq_len(length.out = nlevels(x = plotting_frame[, aes_list$geom_point$shape])))
          }
        }

        # geom_text
        if (!is.null(x = aes_list$geom_text)) {
          ggplot_object <- ggplot_object +
            ggplot2::geom_text(
              mapping = ggplot2::aes_(
                x = quote(expr = x),
                y = quote(expr = y),
                label = if (is.null(x = aes_list$geom_text$label))
                  "x"
                else
                  as.name(x = aes_list$geom_text$label),
                colour = if (is.null(x = aes_list$geom_text$colour))
                  NULL
                else
                  as.name(x = aes_list$geom_text$colour)
              ),
              size = 2.0,
              alpha = I(1 / 3)
            )
        }

        # geom_path
        if (!is.null(x = aes_list$geom_path)) {
          ggplot_object <-
            ggplot_object + ggplot2::geom_path(
              mapping = ggplot2::aes_(
                x = quote(expr = x),
                y = quote(expr = y),
                colour = if (is.null(x = aes_list$geom_path$colour))
                  NULL
                else
                  as.name(x = aes_list$geom_path$colour),
                group = if (is.null(x = aes_list$geom_path$group))
                  NULL
                else
                  as.name(x = aes_list$geom_path$group),
                linetype = if (is.null(x = aes_list$geom_path$linetype))
                  NULL
                else
                  as.name(x = aes_list$geom_path$linetype)
              ),
              arrow = arrow(
                length = unit(x = 0.08, units = "inches"),
                type = "closed"
              )
            )
        }

        ggplot_object <-
          ggplot_object + facet_grid(
            rows = ggplot2::vars(component_1),
            cols = ggplot2::vars(component_2),
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
          write.table(
            x = plotting_frame,
            file = file.path(output_directory,
                             paste(
                               paste(
                                 prefix,
                                 "pca",
                                 aes_list_to_character(aes_list = aes_list),
                                 suffix,
                                 sep = "_"
                               ),
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
    dummy_list,
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
BiocParallel::register(BPPARAM = MulticoreParam(workers = argument_list$threads))

# The working directory is the analyis genome directory.
# Create a new sub-directory for results if it does not exist.
if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

global_design_list <- initialise_design_list()
annotation_frame <- initialise_annotation_frame()
deseq_data_set <-
  initialise_deseq_data_set(design_list = global_design_list)

# TODO: Write the matrix of gene versus model coefficients as a data.frame to disk.

# Cooks Distances Plot ----------------------------------------------------

plot_cooks_distances(object = deseq_data_set)

# RIN Score Plot ----------------------------------------------------------


# If RIN scores are annotated in the sample frame, plot their distribution.
plot_rin_scores(object = deseq_data_set)

# Likelihood Ratio Test (LRT) ---------------------------------------------


# The "reduced_formulas" variable of the "design" data frame encodes reduced
# formulas for LRT. Example: "name_1:~genotype + gender;name_2:~1"
reduced_formula_list <-
  lapply(
    X = stringr::str_split(
      string = stringr::str_split(
        string = global_design_list$reduced_formulas[1L],
        pattern = stringr::fixed(pattern = ";")
      )[[1L]],
      pattern = stringr::fixed(pattern = ":")
    ),
    FUN = function(character_1) {
      reduced_formula_character <- character_1[2L]
      attr(x = reduced_formula_character, which = "reduced_name") <-
        character_1[1L]
      return(reduced_formula_character)
    }
  )

# The plyr::ldply() function interates over a list and returns a data frame.
# Each list element should yield a data frame.
reduced_formula_frame <- plyr::ldply(
  .data = reduced_formula_list,
  .fun = function(reduced_formula_character) {
    summary_frame <- NULL
    # Skip NA or empty character vectors.
    if (is.na(x = reduced_formula_character) ||
        !base::nzchar(x = reduced_formula_character)) {
      # Return NULL instead of a data.frame, which can still be processed by rbind().
      return(summary_frame)
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
      deseq_merge_significant <-
        read.table(file = file_path_significant,
                   header = TRUE,
                   sep = "\t")
      summary_frame <- data.frame(
        "design" = global_design_list$design,
        "full_formula" = global_design_list$full_formula,
        "reduced_name" = attr(x = reduced_formula_character, which = "reduced_name"),
        "reduced_formula" = reduced_formula_character,
        "significant" = nrow(deseq_merge_significant),
        stringsAsFactors = FALSE
      )
      rm(deseq_merge_significant)
    } else {
      message(
        "Processing reduced formula: ",
        attr(x = reduced_formula_character, which = "reduced_name")
      )
      # DESeq LRT testing requires either two model formulas or two model matrices.
      # Create a reduced model matrix and check whether it is full rank.
      formula_full <-
        as.formula(object = global_design_list$full_formula)
      result_list_full <-
        check_model_matrix(
          model_matrix = stats::model.matrix.default(
            object = formula_full,
            data = SummarizedExperiment::colData(x = deseq_data_set)
          )
        )

      formula_reduced <-
        as.formula(object = reduced_formula_character)
      model_matrix_reduced <-
        stats::model.matrix.default(object = formula_reduced,
                                    data = SummarizedExperiment::colData(x = deseq_data_set))
      result_list_reduced <-
        check_model_matrix(model_matrix = model_matrix_reduced)
      if (argument_list$verbose) {
        message("Writing initial reduced model matrix")
        write.table(
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
        write.table(
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
        DESeq(
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
      # print(x = resultsNames(object = deseq_data_set_lrt))
      deseq_results_lrt <-
        DESeq2::results(
          object = deseq_data_set_lrt,
          format = "DataFrame",
          tidy = FALSE,
          # If tidy is TRUE, a classical data.frame is returned.
          parallel = TRUE
        )
      # Re-adjust the DESeqResults DataFrame for merging with the annotation DataFrame.
      deseq_results_lrt_frame <-
        DataFrame(
          gene_id = rownames(x = deseq_results_lrt),
          deseq_results_lrt[, c("baseMean",
                                "log2FoldChange",
                                "lfcSE",
                                "stat",
                                "pvalue",
                                "padj")],
          significant = factor(x = "no", levels = c("no", "yes"))
        )
      deseq_results_lrt_frame[!is.na(x = deseq_results_lrt_frame$padj) &
                                deseq_results_lrt_frame$padj <= argument_list$padj_threshold, "significant"] <-
        "yes"

      # Write all genes.

      deseq_merge_complete <-
        merge(x = annotation_frame, y = deseq_results_lrt_frame, by = "gene_id")

      write.table(
        x = deseq_merge_complete,
        file = file_path_all,
        sep = "\t",
        col.names = TRUE,
        row.names = FALSE
      )

      # Write only significant genes.
      deseq_merge_significant <-
        subset(x = deseq_merge_complete, padj <= argument_list$padj_threshold)

      write.table(
        x = deseq_merge_significant,
        file = file_path_significant,
        sep = "\t",
        col.names = TRUE,
        row.names = FALSE
      )
      summary_frame <- data.frame(
        "design" = global_design_list$design,
        "full_formula" = global_design_list$full_formula,
        "reduced_name" = attr(x = reduced_formula_character, which = "reduced_name"),
        "reduced_formula" = reduced_formula_character,
        "significant" = nrow(deseq_merge_significant),
        stringsAsFactors = FALSE
      )
      rm(
        deseq_merge_significant,
        deseq_merge_complete,
        deseq_results_lrt_frame,
        deseq_results_lrt,
        deseq_data_set_lrt
      )
    }
    rm(file_path_all, file_path_significant)
    return(summary_frame)
  }
)

# Write the reduced formula (LRT) summary frame to disk.
write.table(
  x = reduced_formula_frame,
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
rm(reduced_formula_frame, reduced_formula_list)


# Plot Aesthetics ---------------------------------------------------------


# The plot_aes variable of the design data frame supplies a semi-colon-separated
# list of geometric objects and their associated aestethics for each plot, which
# is a comma-separared list of aestethics=variable mappings.
# plot_aes='colour=group,shape=gender;colour=group,shape=extraction'
# geom_point:colour=test_1,shape=test_2;geom_line:colour=test_3,group=test_4|
# geom_point:colour=test_a,shape=test_b;geom_line:colour=test_c,group=test_d
# Convert into a list of list objects with variables and aesthetics as names.
plot_list <-
  lapply(
    X = stringr::str_split(
      string = global_design_list$plot_aes[1L],
      pattern = stringr::fixed(pattern = "|")
    )[[1L]],
    FUN = function(plot_character) {
      single_geom_list <- lapply(
        X = stringr::str_split(
          string = stringr::str_split(
            string = plot_character[1L],
            pattern = stringr::fixed(pattern = ";")
          )[[1L]],
          pattern = stringr::fixed(pattern = ":")
        ),
        FUN = function(geom_aes_character) {
          # Split on "," characters and assign names (geometric names) to the list components (aesthetic list).
          geom_aes_list <-
            lapply(
              X = stringr::str_split(
                string = geom_aes_character[2L],
                pattern = stringr::fixed(pattern = ",")
              ),
              FUN = function(aes_character) {
                # Split on "=" characters and assign names (aestetic names) to the list components (variable names).
                temporary_list <-
                  stringr::str_split(string = aes_character,
                                     pattern = stringr::fixed(pattern = "="))
                aes_list <-
                  lapply(
                    X = temporary_list,
                    FUN = function(temporary_character) {
                      return(temporary_character[2L])
                    }
                  )
                names(x = aes_list) <-
                  lapply(
                    X = temporary_list,
                    FUN = function(temporary_character) {
                      return(temporary_character[1L])
                    }
                  )
                rm(temporary_list)
                return(aes_list)
              }
            )
          names(x = geom_aes_list) <- geom_aes_character[1L]
          return(geom_aes_list)
        }
      )
      # Flatten the single geometric list and assign geometric object names to the list components.
      geom_list <-
        lapply(
          X = single_geom_list,
          FUN = function(temporary_list) {
            return(temporary_list[[1L]])
          }
        )
      names(geom_list) <-
        lapply(
          X = single_geom_list,
          FUN = function(temporary_list) {
            return(names(x = temporary_list[1L]))
          }
        )
      rm(single_geom_list)
      return(geom_list)
    }
  )


# DESeqDataSet Results ----------------------------------------------------


print(x = "DESeqDataSet result names:")
print(x = DESeq2::resultsNames(object = deseq_data_set))

# Export RAW counts -------------------------------------------------------


# Export the raw counts from the DESeqDataSet object.
counts_frame <-
  as(object = SummarizedExperiment::assays(x = deseq_data_set)$counts,
     Class = "DataFrame")
counts_frame$gene_id <- row.names(x = counts_frame)

write.table(
  x = merge(x = annotation_frame, y = counts_frame, by = "gene_id"),
  file = file.path(output_directory,
                   paste(prefix,
                         "counts",
                         "raw.tsv",
                         sep = "_")),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)
rm(counts_frame)

# Export FPKM values ------------------------------------------------------


# Retrieve and plot FPKM values.
fpkm_matrix <- DESeq2::fpkm(object = deseq_data_set)
plot_fpkm_values(object = fpkm_matrix)

# Export FPKM values from the DESeqDataSet object
fpkm_frame <-
  as(object = fpkm_matrix,
     Class = "DataFrame")
fpkm_frame$gene_id <- row.names(x = fpkm_frame)

write.table(
  x = merge(x = annotation_frame, y = fpkm_frame, by = "gene_id"),
  file = file.path(output_directory,
                   paste(prefix,
                         "fpkms.tsv",
                         sep = "_")),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)
rm(fpkm_frame, fpkm_matrix)

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

  # Export the vst counts from the DESeqTransform object
  counts_frame <-
    as(object = SummarizedExperiment::assay(x = deseq_transform, i = 1L),
       Class = "DataFrame")
  counts_frame$gene_id <- row.names(x = counts_frame)

  write.table(
    x = merge(x = annotation_frame, y = counts_frame, by = "gene_id"),
    file = file.path(output_directory,
                     paste(
                       paste(prefix,
                             "counts",
                             "vst",
                             suffix,
                             sep = "_"), "tsv", sep = "."
                     )),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
  )

  rm(counts_frame, deseq_transform, suffix)
}
rm(blind, plot_list)

# Save an R image for project-specific post-processing later.
# save.image()

rm(
  annotation_frame,
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
  aes_list_to_character,
  initialise_deseq_transform,
  initialise_deseq_data_set,
  check_model_matrix,
  fix_model_matrix,
  initialise_ranged_summarized_experiment,
  initialise_design_list,
  initialise_sample_frame,
  initialise_annotation_frame
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessionInfo())

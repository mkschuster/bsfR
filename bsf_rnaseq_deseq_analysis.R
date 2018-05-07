#! /usr/bin/env Rscript
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --distribution=block
#SBATCH --mem=131072
#SBATCH --time=2-00:00:00
#SBATCH --partition=mediumq
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --error=.bsf_rnaseq_deseq_analysis_%j.err
#SBATCH --output=.bsf_rnaseq_deseq_analysis_%j.out
#
# BSF R script to run a DESeq2 analysis.
#
# Copyright 2013 - 2017 Michael K. Schuster
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
# The rnaseq_deseq_PREFIX_samples.tsv sample annotation DataFrame supports
# the following variables.
#
#   "bam_path":
#      A "character" vector of BAM file paths.
#
#   "bai_path":
#      A "character" vector of BAI file paths.
#
#   "designs":
#      A "character" vector of comma-separated values of designs,
#      a particular sample should be part of.
#
#   "library_type":
#      A "factor" vector with levels "unstranded", "first" and "second" to
#      indicate the strandedness of the RNA-seq protocol and whether the
#      first or second strand get sequenced. Illumina TruSeq standed mRNA
#      sequences the second strand so that reads need inverting before
#      counting strand-specifically.
#
#   "sequencing_type":
#      A factor with levels "SE" and "PE" to indicate single-end or
#      paired-end sequencing and thus counting as read pairs or not.
#
#   "total_counts":
#      Total counts per sample.
#      Calculated automatically based on the colSums() of the counts() function.
#
#   "RIN":
#      A numeric vector providing the RNA integrity number (RIN) score
#      per sample. If available, the RIN score distribution will be plotted.
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

suppressPackageStartupMessages(expr = library(package = "DESeq2"))  # for DESeq2::DESeqDataSet() and DESeq2::DESeq()
suppressPackageStartupMessages(expr = library(package = "BiocParallel"))  # for BiocParallel::register()
suppressPackageStartupMessages(expr = library(package = "GenomicAlignments"))  # for GenomicAlignments::summarizeOverlaps()
suppressPackageStartupMessages(expr = library(package = "RColorBrewer"))  # for RColorBrewer::brewer.pal()
suppressPackageStartupMessages(expr = library(package = "Rsamtools"))  # for Rsamtools::BamFileList()
suppressPackageStartupMessages(expr = library(package = "caret"))  # For caret::findLinearCombos()
suppressPackageStartupMessages(expr = library(package = "genefilter"))  # for genefilter::rowVars()
suppressPackageStartupMessages(expr = library(package = "ggplot2"))  # for ggplot2::ggplot()
suppressPackageStartupMessages(expr = library(package = "grid"))  # for grid::grid.newpage() and grid::grid.draw()
suppressPackageStartupMessages(expr = library(package = "pheatmap"))  # for pheatmap::pheatmap()
# suppressPackageStartupMessages(expr = library(package = "reshape2"))  # for reshape2::melt()
suppressPackageStartupMessages(expr = library(package = "rtracklayer"))  # for rtracklayer::import() GTF import
suppressPackageStartupMessages(expr = library(package = "stringi"))  # For stringi::stri_split_fixed()

# Save plots in the following formats.

graphics_formats <- c("pdf", "png")

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
#' @return DataFrame
#' @export
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
    # gene_name_list <- lapply(X = rowRanges(x = deseq_data_set), FUN = function(x) { mcols(x = x)[1L, "gene_name"] })
    message("Reading reference GTF gene features")
    gene_ranges <-
      import(
        con = argument_list$gtf_reference,
        format = "gtf",
        genome = argument_list$genome_version,
        feature.type = "gene"
      )
    message("Creating annotation frame")
    data_frame <-
      mcols(x = gene_ranges)[, c("gene_id",
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
#' @param factor_levels character vector with a packed string to assign factor levels
#' @references argument_list
#' @references output_directory
#' @references prefix
#' @return Sample DataFrame
#' @export
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
      X = stri_split_fixed(str = as.character(data_frame$designs),
                           pattern = ","),
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
  # for the summarizeOverlaps() read counting function.
  
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
    X = stri_split_fixed(
      str = stri_split_fixed(str = factor_levels[1L],
                             pattern = ";")[[1L]],
      pattern = ":"
    ),
    FUN = function(character_1) {
      # Split the second component of character_1, the factor levels, on ",".
      character_2 <-
        unlist(x = stri_split_fixed(str = character_1[2L],
                                    pattern = ","))
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
    if (factor_name %in% design_variables) {
      data_frame[, factor_name] <-
        factor(x = as.character(x = data_frame[, factor_name]),
               levels = factor_list[[i]])
      # Check for NA values in case a factor level was missing.
      if (any(is.na(x = data_frame[, factor_name]))) {
        stop(
          paste0(
            "Missing values after assigning factor levels for factor name ",
            factor_name
          )
        )
      }
    } else {
      stop(
        paste0(
          "Factor name ",
          factor_name,
          " does not resemble a variable of the design frame."
        )
      )
    }
    rm(factor_name)
  }
  rm(i, design_variables, factor_list)
  
  # Drop any unused levels from the sample data frame before retrning it.
  return(droplevels(x = data_frame))
}


# Initialise a Design list object -----------------------------------------

#' Initialise a Design list
#'
#' @references argument_list
#' @references output_directory
#' @references prefix
#' @return Named list of the selected design
#' @export
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


#' Initialise or load a RangedSummarizedExperiment object.
#'
#' @param design_list Named list of design information
#' @references argument_list
#' @references output_directory
#' @references prefix
#' @return RangedSummarizedExperiment object
#' @export
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
      import(
        con = argument_list$gtf_reference,
        format = "gtf",
        genome = argument_list$genome_version,
        feature.type = "exon"
      )
    # Convert (i.e. split) the GRanges object into a GRangesList object
    # by gene identifiers.
    gene_ranges_list <-
      split(x = exon_ranges, f = mcols(x = exon_ranges)$gene_id)
    
    # Process per library_type and sequencing_type and merge the RangedSummarizedExperiment objects.
    ranged_summarized_experiment <- NULL
    
    for (library_type in levels(x = sample_frame$library_type)) {
      for (sequencing_type in levels(x = sample_frame$sequencing_type)) {
        message(
          paste0(
            "Processing library_type: ",
            library_type,
            " sequencing_type: ",
            sequencing_type
          )
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
        bam_file_list <- BamFileList(
          file = as.character(x = sub_sample_frame$bam_path),
          index = as.character(x = sub_sample_frame$bai_path),
          yieldSize = 2000000L,
          asMates = (sequencing_type == "PE")
        )
        names(x = bam_file_list) <-
          as.character(x = sub_sample_frame$sample)
        
        message("Creating a RangedSummarizedExperiment object")
        sub_ranged_summarized_experiment <-
          summarizeOverlaps(
            features = gene_ranges_list,
            reads = bam_file_list,
            mode = "Union",
            ignore.strand = (library_type == "unstranded"),
            # Exclude reads that represent secondary alignments or fail the vendor quality filter.
            param = ScanBamParam(
              flag = scanBamFlag(
                isSecondaryAlignment = FALSE,
                isNotPassingQualityControls = FALSE
              )
            ),
            # Invert the strand for protocols that sequence the second strand.
            preprocess.reads = if (library_type == "second")
              invertStrand
          )
        colData(x = sub_ranged_summarized_experiment) <-
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
    
    # Calculate colSums() of assay() and add as total_count into the colData data frame.
    sample_frame <- colData(x = ranged_summarized_experiment)
    sample_frame$total_counts <- colSums(x = assay(x = ranged_summarized_experiment), na.rm = TRUE)
    colData(x = ranged_summarized_experiment) <- sample_frame
    
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
#' @param model_matrix_local Model matrix
#'
#' @return
#' @export
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
        "Levels or combinations of levels without any samples have resulted in column(s) of zeros in the model matrix."
      )
      message(paste(colnames(x = model_matrix_local)[model_all_zero], collapse = ", "))
      message("Attempting to fix the model matrix by removing empty columns.")
      model_matrix_local <-
        model_matrix_local[,-which(x = model_all_zero)]
    } else {
      message(
        "One or more variables or interaction terms in the design formula are linear combinations of the others."
      )
      message("Attempting to fix the model by removing linear combinations.")
      linear_combinations_list <-
        findLinearCombos(x = model_matrix_local)
      # print(x = linear_combinations_list)
      # print(x = colnames(x = model_matrix_local)[linear_combinations_list$remove])
      model_matrix_local <-
        model_matrix_local[, -linear_combinations_list$remove]
    }
    rm(model_all_zero)
  }
  return(list("model_matrix" = model_matrix_local, "full_rank" = full_rank))
}

#' Check a model matrix for being full rank.
#'
#' @param model_matrix Model matrix
#'
#' @return Named list of "model_matrix" and "formula_full_rank", a boolean to indicate
#' whether the original formula was already full rank.
#' @export
#'
#' @examples
check_model_matrix <- function(model_matrix) {
  # Write the unmodified model matrix to disk if argument "--verbose" was set.
  if (FALSE) {
    # FIXME: Writing out original model matrices no longer works,
    # because this function is now used on full and reduced matrices.
    if (argument_list$verbose) {
      message("Writing model matrix")
      model_frame <- as.data.frame(x = model_matrix)
      model_path <-
        file.path(output_directory,
                  paste0(prefix, "_model_matrix_initial.tsv"))
      write.table(
        x = model_frame,
        file = model_path,
        sep = "\t",
        row.names = TRUE,
        col.names = TRUE
      )
      rm(model_path, model_frame)
    }
  }
  
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
  
  if (FALSE) {
    # FIXME: Writing out modified model matrices no longer works,
    # because this function is now used on full and reduced matrices.
    if (argument_list$verbose) {
      message("Writing modified model matrix")
      model_frame <- as.data.frame(x = model_matrix)
      model_path <-
        file.path(output_directory,
                  paste0(prefix, "_model_matrix_modified.tsv"))
      write.table(
        x = model_frame,
        file = model_path,
        sep = "\t",
        row.names = TRUE,
        col.names = TRUE
      )
      rm(model_path, model_frame)
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
#' @param design_list Named list of design information
#'
#' @return DESeqDataSet object
#' @export
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
    result_list <-
      check_model_matrix(model_matrix = stats::model.matrix.default(
        object = as.formula(object = design_list$full_formula),
        data = colData(x = ranged_summarized_experiment)
      ))
    
    if (result_list$formula_full_rank) {
      # The design formula *is* full rank, so set it as "design" option directly.
      message("Creating a DESeqDataSet object with a model formula")
      deseq_data_set <-
        DESeqDataSet(se = ranged_summarized_experiment,
                     design = as.formula(object = design_list$full_formula))
      # Start DESeq2 Wald testing.
      # Set betaPrior = FALSE for consistent result names for designs, regardless of interaction terms.
      # DESeq2 seems to set betaPrior = FALSE upon interaction terms, automatically.
      # See: https://support.bioconductor.org/p/84366/
      # betaPrior also has to be FALSE in case a user-supplied full model matrix is specified.
      message("Started DESeq Wald testing with a model formula")
      deseq_data_set <-
        DESeq(
          object = deseq_data_set,
          test = "Wald",
          betaPrior = FALSE,
          parallel = TRUE
        )
    } else {
      # The orignal design formula was not full rank.
      # Unfortunately, to initialise the DESeqDataSet,
      # a model matrix can apparently not be used directly.
      # Thus, use the simplest possible design (i.e. ~ 1) for initialisation and
      # peform Wald testing with the full model matrix.
      message("Creating a DESeqDataSet object with design formula ~ 1")
      deseq_data_set <-
        DESeqDataSet(se = ranged_summarized_experiment,
                     design = ~ 1)
      message("Started DESeq Wald testing with a model matrix")
      deseq_data_set <-
        DESeq(
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
#' @param blind bool to create the DESeqTransform object blindly or based on
#'              the model
#' @return DESeqTransform object
#' @export
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
      message(paste("Loading a", suffix, "DESeqTransform object"))
      load(file = file_path)
    } else {
      message(paste("Creating a", suffix, "DESeqTransform object"))
      # Run variance stabilizing transformation (VST) to get homoskedastic data for PCA plots.
      deseq_transform <-
        varianceStabilizingTransformation(object = deseq_data_set, blind = blind)
      save(deseq_transform, file = file_path)
    }
    rm(file_path, suffix)
    
    return(deseq_transform)
  }

# Convert aes_list into character -----------------------------------------


#' Convert the aes_list into a simple character string for file and plot naming.
#'
#' @param aes_list
#'
#' @return
#' @export
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


# Plot RIN scores ---------------------------------------------------------


#' Plot RIN scores.
#'
#' @param object DESeqDataSet object
#'
#' @return
#' @export
#'
#' @examples
plot_rin_scores <- function(object) {
  plot_paths <- file.path(output_directory,
                          paste(
                            paste(prefix,
                                  "rin_density",
                                  sep = "_"),
                            graphics_formats,
                            sep = "."
                          ))
  if (all(file.exists(plot_paths) &&
          file.info(plot_paths)$size > 0L)) {
    message("Skipping a RIN score density plot")
  } else {
    if ("RIN" %in% colnames(x = colData(x = object))) {
      message("Creating a RIN score density plot")
      ggplot_object <-
        ggplot(data = as.data.frame(x = colData(x = object)))
      ggplot_object <-
        ggplot_object + ggtitle(label = "RNA Integry Number (RIN) Density Plot")
      ggplot_object <- ggplot_object + xlim(RIN = c(0.0, 10.0))
      ggplot_object <-
        ggplot_object + geom_vline(xintercept = 1.0,
                                   colour = "red",
                                   linetype = 2L)
      ggplot_object <-
        ggplot_object + geom_vline(xintercept = 4.0,
                                   colour = "yellow",
                                   linetype = 2L)
      ggplot_object <-
        ggplot_object + geom_vline(xintercept = 7.0,
                                   colour = "green",
                                   linetype = 2L)
      ggplot_object <-
        ggplot_object + geom_density(mapping = aes(x = RIN, y = ..density..))
      for (plot_path in plot_paths) {
        ggsave(
          filename = plot_path,
          width = argument_list$plot_width,
          height = argument_list$plot_height,
          plot = ggplot_object
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
#' @return
#' @export
#'
#' @examples
plot_mds <- function(object,
                     plot_list = list(),
                     blind = FALSE) {
  suffix <- if (blind)
    "blind"
  else
    "model"
  
  message(paste("Creating", suffix, "MDS plots"))
  
  dist_object <- dist(x = t(x = assay(x = object)))
  dist_matrix <- as.matrix(x = dist_object)
  mds_frame <-
    cbind(data.frame(cmdscale(d = dist_matrix)), as.data.frame(colData(x = object)))
  
  dummy_list <-
    lapply(
      X = plot_list,
      FUN = function(aes_list) {
        # plot_mds(object = deseq_transform, aes_list = aes_list)
        ggplot_object <-
          ggplot(data = mds_frame)
        
        # geom_line
        if (!is.null(x = aes_list$geom_line)) {
          ggplot_object <-
            ggplot_object +
            geom_line(
              mapping = aes_(
                x = quote(expr = X1),
                y = quote(expr = X2),
                colour = if (is.null(x = aes_list$geom_line$colour))
                  "black"
                else
                  as.name(x = aes_list$geom_line$colour),
                group = if (is.null(x = aes_list$geom_line$group))
                  1L
                else
                  as.name(x = aes_list$geom_line$group),
                linetype = if (is.null(x = aes_list$geom_path$linetype))
                  "solid"
                else
                  as.name(x = aes_list$geom_path$linetype)
              ),
              alpha = I(1 / 3)
            )
        }
        
        # geom_point
        if (!is.null(x = aes_list$geom_point)) {
          ggplot_object <- ggplot_object +
            geom_point(
              mapping = aes_(
                x = quote(expr = X1),
                y = quote(expr = X2),
                colour = if (is.null(x = aes_list$geom_point$colour))
                  "black"
                else
                  as.name(x = aes_list$geom_point$colour),
                shape = if (is.null(x = aes_list$geom_point$shape))
                  15L
                else
                  as.name(x = aes_list$geom_point$shape)
              ),
              size = 2.0,
              alpha = I(1 / 3)
            )
          if (is.null(x = aes_list$geom_point$shape)) {
            # In case the shape is not mapped, use values without scaling.
            ggplot_object <- ggplot_object + scale_shape_identity()
          }
        }
        
        # geom_text
        if (!is.null(x = aes_list$geom_text)) {
          ggplot_object <- ggplot_object +
            geom_text(
              mapping = aes_(
                x = quote(expr = X1),
                y = quote(expr = X2),
                label = if (is.null(x = aes_list$geom_text$label))
                  "x"
                else
                  as.name(x = aes_list$geom_text$label),
                colour = if (is.null(x = aes_list$geom_text$colour))
                  "black"
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
            ggplot_object + geom_path(
              mapping = aes_(
                x = quote(expr = X1),
                y = quote(expr = X2),
                colour = if (is.null(x = aes_list$geom_path$colour))
                  "black"
                else
                  as.name(x = aes_list$geom_path$colour),
                group = if (is.null(x = aes_list$geom_path$group))
                  1L
                else
                  as.name(x = aes_list$geom_path$group),
                linetype = if (is.null(x = aes_list$geom_path$linetype))
                  "solid"
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
          theme_bw() +
          coord_fixed()
        
        # ggplot_object <- ggplot_object + xlim(min(mds_frame$X1, mds_frame$X2), max(mds_frame$X1, mds_frame$X2))
        # ggplot_object <- ggplot_object + ylim(min(mds_frame$X1, mds_frame$X2), max(mds_frame$X1, mds_frame$X2))
        
        for (graphics_format in graphics_formats) {
          ggsave(filename = file.path(
            output_directory,
            paste(
              paste(
                prefix,
                "mds",
                aes_list_to_character(aes_list = aes_list),
                suffix,
                sep = "_"
              ),
              graphics_format,
              sep = "."
            )
          ),
          plot = ggplot_object)
        }
        rm(graphics_format,
           ggplot_object)
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
#' The heatmap plot is based on hierarchical clustering of
#' "Euclidean" distances of transformed counts for each gene.
#'
#' @param object DESeqTransform object
#' @param aes_list ggplot2 aesthethics list
#' @param blind bool to indicate a blind or model-based DESeqTransform object
#'
#' @return
#' @export
#'
#' @examples
plot_heatmap <- function(object,
                         aes_list = list(),
                         blind = FALSE) {
  suffix <- if (blind)
    "blind"
  else
    "model"
  
  message(paste("Creating a", suffix, "Heatmap plot"))
  
  aes_character <-
    unique(x = unlist(x = aes_list, use.names = TRUE))
  message(paste("  Heat map plot:", aes_character))
  if (!all(aes_character %in% names(x = colData(x = object)))) {
    stop("the argument 'aes_character' should specify columns of colData(dds)")
  }
  
  plotting_frame <-
    as.data.frame(colData(x = object)[, aes_character, drop = FALSE])
  
  # Add the aes_character factors together to create a new grouping factor
  group_factor <- if (length(x = aes_character) > 1) {
    factor(x = apply(
      X = plotting_frame,
      MARGIN = 1,
      FUN = paste,
      collapse = " : "
    ))
  } else {
    colData(x = object)[[aes_character]]
  }
  
  dist_object <- dist(x = t(x = assay(x = object)))
  dist_matrix <- as.matrix(x = dist_object)
  colnames(x = dist_matrix) <- NULL
  rownames(x = dist_matrix) <-
    paste(object$sample, group_factor, sep = "-")
  # TODO: Rather use the ComplexHeatmap package?
  pheatmap_object <-
    pheatmap(
      mat = dist_matrix,
      clustering_distance_rows = dist_object,
      clustering_distance_cols = dist_object,
      color = colorRampPalette(colors = rev(x = brewer.pal(
        n = 9, name = "Blues"
      )))(255),
      fontsize_row = 6
    )
  
  # PDF output
  pdf(
    file = file.path(output_directory,
                     paste(
                       paste(
                         prefix,
                         "heatmap",
                         paste(aes_character, collapse = "_"),
                         suffix,
                         sep = "_"
                       ), "pdf", sep = "."
                     )),
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  # grid::grid.newpage()
  grid::grid.draw(pheatmap_object$gtable)
  device_number <- dev.off()
  
  # PNG output
  png(
    file = file.path(output_directory,
                     paste(
                       paste(
                         prefix,
                         "heatmap",
                         paste(aes_character, collapse = "_"),
                         suffix,
                         sep = "_"
                       ), "png", sep = "."
                     )),
    width = argument_list$plot_width,
    height = argument_list$plot_height,
    units = "in",
    res = 300L
  )
  # grid::grid.newpage()
  grid::grid.draw(pheatmap_object$gtable)
  device_number <- dev.off()
  
  rm(
    device_number,
    dist_matrix,
    dist_object,
    group_factor,
    plotting_frame,
    aes_character,
    suffix
  )
}

# Plot Principal Component Analysis (PCA) ---------------------------------


#' Plot a Principal Component Analysis (PCA)
#'
#' @param object DESeqTransform object
#' @param plot_list List of lists configuring plots and their ggplot2 aesthetic mappings
#' @param blind bool to indicate a blind or model-based DESeqTransform object
#'
#' @return
#' @export
#'
#' @examples
plot_pca <- function(object,
                     plot_list = list(),
                     blind = FALSE) {
  suffix <- if (blind)
    "blind"
  else
    "model"
  
  message(paste("Creating", suffix, "PCA plots"))
  
  # Calculate the variance for each gene.
  row_variance <- rowVars(assay(x = object))
  # Select the top number of genes by variance.
  selected_rows <-
    order(row_variance, decreasing = TRUE)[seq_len(length.out = min(argument_list$pca_top_number, length(x = row_variance)))]
  
  # Perform a PCA on the data in assay(x) for the selected genes
  pca_object <-
    prcomp(x = t(x = assay(x = object)[selected_rows,]))
  rm(selected_rows)
  
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
        message(paste("  PCA plot:", aes_character, collapse = " "))
        if (!all(aes_character %in% names(x = colData(x = object)))) {
          stop("the argument 'aes_character' should specify columns of colData(dds)")
        }
        # Assemble the data for the plot.
        pca_frame <- as.data.frame(pca_object$x)
        plotting_frame <-
          data.frame(
            component1 = factor(levels = paste0(
              "PC", seq_len(length.out = pca_dimensions)
            )),
            component2 = factor(levels = paste0(
              "PC", seq_len(length.out = pca_dimensions)
            )),
            x = numeric(),
            y = numeric(),
            # Also initalise all variables of the column data, but do not include data (i.e. 0L rows).
            colData(x = object)[0L, ]
          )
        
        for (column_number in seq_len(length.out = ncol(x = pca_pair_matrix))) {
          pca_label_1 <-
            paste0("PC", pca_pair_matrix[1L, column_number])
          pca_label_2 <-
            paste0("PC", pca_pair_matrix[2L, column_number])
          plotting_frame <- rbind.data.frame(
            plotting_frame,
            data.frame(
              component_1 = pca_label_1,
              component_2 = pca_label_2,
              x = pca_frame[, pca_label_1],
              y = pca_frame[, pca_label_2],
              colData(x = object)
            )
          )
          rm(pca_label_1, pca_label_2)
        }
        rm(column_number)
        
        ggplot_object <- ggplot(data = plotting_frame)
        
        # geom_line
        if (!is.null(x = aes_list$geom_line)) {
          ggplot_object <-
            ggplot_object +
            geom_line(
              mapping = aes_(
                x = quote(expr = x),
                y = quote(expr = y),
                colour = if (is.null(x = aes_list$geom_line$colour))
                  "black"
                else
                  as.name(x = aes_list$geom_line$colour),
                group = if (is.null(x = aes_list$geom_line$group))
                  1L
                else
                  as.name(x = aes_list$geom_line$group)
              ),
              alpha = I(1 / 3)
            )
        }
        
        # geom_point
        if (!is.null(x = aes_list$geom_point)) {
          ggplot_object <- ggplot_object +
            geom_point(
              mapping = aes_(
                x = quote(expr = x),
                y = quote(expr = y),
                colour = if (is.null(x = aes_list$geom_point$colour))
                  "black"
                else
                  as.name(x = aes_list$geom_point$colour),
                shape = if (is.null(x = aes_list$geom_point$shape))
                  15L
                else
                  as.name(x = aes_list$geom_point$shape)
              ),
              size = 2.0,
              alpha = I(1 / 3)
            )
          if (is.null(x = aes_list$geom_point$shape)) {
            # In case the shape is not mapped, use values without scaling.
            ggplot_object <- ggplot_object + scale_shape_identity()
          }
        }
        
        # geom_text
        if (!is.null(x = aes_list$geom_text)) {
          ggplot_object <- ggplot_object +
            geom_text(
              mapping = aes_(
                x = quote(expr = x),
                y = quote(expr = y),
                label = if (is.null(x = aes_list$geom_text$label))
                  "x"
                else
                  as.name(x = aes_list$geom_text$label),
                colour = if (is.null(x = aes_list$geom_text$colour))
                  "black"
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
            ggplot_object + geom_path(
              mapping = aes_(
                x = quote(expr = x),
                y = quote(expr = y),
                colour = if (is.null(x = aes_list$geom_path$colour))
                  "black"
                else
                  as.name(x = aes_list$geom_path$colour),
                group = if (is.null(x = aes_list$geom_path$group))
                  "black"
                else
                  as.name(x = aes_list$geom_path$group),
                linetype = if (is.null(x = aes_list$geom_path$linetype))
                  "solid"
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
            facets = component_1 ~ component_2,
            labeller = labeller(component_1 = label_function, component_2 = label_function)
          )
        for (graphics_format in graphics_formats) {
          ggsave(filename = file.path(
            output_directory,
            paste(
              paste(
                prefix,
                "pca",
                paste(aes_character, collapse = "_"),
                suffix,
                sep = "_"
              ),
              graphics_format,
              sep = "."
            )
          ))
        }
        rm(
          graphics_format,
          ggplot_object,
          pca_frame,
          plotting_frame,
          aes_character
        )
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


message(paste0("Processing design '", argument_list$design_name, "'"))

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

# RIN Score Plot ----------------------------------------------------------


# If RIN scores are annotated in the sample frame, plot their distribution.
plot_rin_scores(object = deseq_data_set)


# Likelihood Ratio Test (LRT) ---------------------------------------------


# The "reduced_formulas" variable of the "design" data frame encodes reduced formulas for LRT.
# Example: "name_1:~genotype + gender;name_2:~1"
reduced_formula_list <-
  lapply(
    X = stri_split_fixed(
      str = stri_split_fixed(str = global_design_list$reduced_formulas[1L], pattern = ";")[[1L]],
      pattern = ":"
    ),
    FUN = function(character_1) {
      reduced_formula_character <- character_1[2L]
      attr(x = reduced_formula_character, which = "reduced_name") <-
        character_1[1L]
      return(reduced_formula_character)
    }
  )

temporary_list <- lapply(
  X = reduced_formula_list,
  FUN = function(reduced_formula_character) {
    # Skip empty character vectors.
    if (!nzchar(x = reduced_formula_character)) {
      return()
    }
    file_path <-
      file.path(output_directory,
                paste(paste(
                  prefix,
                  "lrt",
                  attr(x = reduced_formula_character, which = "reduced_name"),
                  sep = "_"
                ),
                "tsv",
                sep = "."))
    
    if (file.exists(file_path) &&
        file.info(file_path)$size > 0L) {
      message(paste0(
        "Skipping reduced formula: ",
        attr(x = reduced_formula_character, which = "reduced_name")
      ))
    } else {
      message(paste0(
        "Processing reduced formula: ",
        attr(x = reduced_formula_character, which = "reduced_name")
      ))
      # DESeq LRT testing requires either two model formulas or two model matrices.
      # Create a reduced model matrix and check whether it is full rank.
      formula_full <-
        as.formula(object = global_design_list$full_formula)
      result_list_full <-
        check_model_matrix(model_matrix = stats::model.matrix.default(object = formula_full,
                                                                      data = colData(x = deseq_data_set)))
      formula_reduced <-
        as.formula(object = reduced_formula_character)
      result_list_reduced <-
        check_model_matrix(model_matrix = stats::model.matrix.default(object = formula_reduced,
                                                                      data = colData(x = deseq_data_set)))
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
        results(
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
      
      deseq_merge <-
        merge(x = annotation_frame, y = deseq_results_lrt_frame, by = "gene_id")
      
      write.table(
        x = deseq_merge,
        file = file_path,
        sep = "\t",
        col.names = TRUE,
        row.names = FALSE
      )
      # Write only significant.
      deseq_merge_significant <-
        subset(x = deseq_merge, padj <= argument_list$padj_threshold)
      
      file_path <-
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
      
      write.table(
        x = deseq_merge_significant,
        file = file_path,
        sep = "\t",
        col.names = TRUE,
        row.names = FALSE
      )
      rm(
        deseq_merge_significant,
        deseq_merge,
        deseq_results_lrt_frame,
        deseq_results_lrt,
        deseq_data_set_lrt
      )
    }
    rm(file_path)
  }
)
rm(temporary_list, reduced_formula_list)


# Plot Aesthetics ---------------------------------------------------------


# The plot_aes variable of the design data frame supplies a semi-colon-separated list of
# geometric objects and their associated aestethics for each plot,
# which is a comma-separared list of aestethics=variable mappings.
# plot_aes='colour=group,shape=gender;colour=group,shape=extraction'
# geom_point:colour=test_1,shape=test_2;geom_line:colour=test_3,group=test_4|
# geom_point:colour=test_a,shape=test_b;geom_line:colour=test_c,group=test_d
# Convert into a list of list objects with variables and aesthetics as names.
plot_list <-
  lapply(
    X = stri_split_fixed(str = global_design_list$plot_aes[1L], pattern = "|")[[1L]],
    FUN = function(plot_character) {
      single_geom_list <- lapply(
        X = stri_split_fixed(
          str = stri_split_fixed(str = plot_character[1L], pattern = ";")[[1L]],
          pattern = ":"
        ),
        FUN = function(geom_aes_character) {
          # Split on "," characters and assign names (geometric names) to the list components (aesthetic list).
          geom_aes_list <-
            lapply(
              X = stri_split_fixed(str = geom_aes_character[2L], pattern = ","),
              FUN = function(aes_character) {
                # Split on "=" characters and assign names (aestetic names) to the list components (variable names).
                temporary_list <-
                  stri_split_fixed(str = aes_character, pattern = "=")
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
print(x = resultsNames(object = deseq_data_set))

# Export RAW counts -------------------------------------------------------


# Export the raw counts from the DESeqDataSet object.
counts_frame <-
  as(object = assays(x = deseq_data_set)$counts,
     Class = "DataFrame")
counts_frame$gene_id <- row.names(x = counts_frame)
deseq_merge <-
  merge(x = annotation_frame, y = counts_frame, by = "gene_id")
file_path <-
  file.path(output_directory,
            paste(prefix,
                  "counts",
                  "raw.tsv",
                  sep = "_"))
write.table(
  x = deseq_merge,
  file = file_path,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)
rm(file_path, counts_frame)

# Export FPKM values ------------------------------------------------------


# Export FPKM values from the DESeqDataSet object
counts_frame <-
  as(object = fpkm(object = deseq_data_set), Class = "DataFrame")
counts_frame$gene_id <- row.names(x = counts_frame)
deseq_merge <-
  merge(x = annotation_frame, y = counts_frame, by = "gene_id")
file_path <-
  file.path(output_directory,
            paste(prefix,
                  "fpkms.tsv",
                  sep = "_"))
write.table(
  x = deseq_merge,
  file = file_path,
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)
rm(file_path, counts_frame)

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
  
  dummy_list <-
    lapply(
      X = plot_list,
      FUN = function(aes_list) {
        plot_heatmap(object = deseq_transform,
                     aes_list = aes_list,
                     blind = blind)
      }
    )
  rm(dummy_list)
  
  # Export the vst counts from the DESeqTransform object
  counts_frame <-
    as(object = assay(x = deseq_transform, i = 1),
       Class = "DataFrame")
  counts_frame$gene_id <- row.names(x = counts_frame)
  deseq_merge <-
    merge(x = annotation_frame, y = counts_frame, by = "gene_id")
  file_path <-
    file.path(output_directory,
              paste(paste(prefix,
                          "counts",
                          "vst",
                          suffix,
                          sep = "_"), "tsv", sep = "."))
  write.table(
    x = deseq_merge,
    file = file_path,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
  )
  rm(file_path, counts_frame)
  
  rm(suffix)
}
rm(blind, plot_list)

# Save an R image for project-specific post-processing later.
# save.image()

rm(
  deseq_merge,
  annotation_frame,
  deseq_transform,
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

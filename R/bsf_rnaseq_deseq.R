#
# Library module of bsfR RNA-seq DESeq2 functions.
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

#' Get a DESeq2 Prefix.
#'
#' Get a design-specific \code{DESeq2} analysis prefix.
#'
#' @param design_name A \code{character} scalar with the design name.
#'
#' @return A \code{character} scalar with the \code{DESeq2} analysis prefix.
#' @export
#'
#' @examples
#' \dontrun{
#' design_name <- "global"
#' prefix_deseq <- bsfrd_get_prefix_deseq(
#'   design_name = design_name)
#' }
bsfrd_get_prefix_deseq <- function(design_name) {
  return(paste("rnaseq",
               "deseq",
               design_name,
               sep = "_"))
}

#' Get a DESeq2 Enrichr Prefix.
#'
#' Get a design-specific \code{DESeq2} analysis Enrichr prefix.
#'
#' @param design_name A \code{character} scalar with the design name.
#'
#' @return A \code{character} scalar with the \code{DESeq2} Enrichr prefix.
#' @export
#'
#' @examples
#' \dontrun{
#' design_name <- "global"
#' prefix_enrichr <- bsfrd_get_prefix_enrichr(
#'   design_name = design_name)
#' }
bsfrd_get_prefix_enrichr <- function(design_name) {
  return(paste("rnaseq",
               "deseq",
               design_name,
               "enrichr",
               sep = "_"))
}

#' Get a DESeq2 Gene Ontology Prefix.
#'
#' Get a design-specific \code{DESeq2} analysis Gene Ontology prefix.
#'
#' @param design_name A \code{character} scalar with the design name.
#'
#' @return A \code{character} scalar with the \code{DESeq2} Gene Ontology prefix.
#' @export
#'
#' @examples
#' \dontrun{
#' design_name <- "global"
#' prefix_go <- bsfrd_get_prefix_go(
#'   design_name = design_name)
#' }
bsfrd_get_prefix_go <- function(design_name) {
  return(paste("rnaseq",
               "deseq",
               design_name,
               "go",
               sep = "_"))
}

#' Get a DESeq2 Heat Map Prefix.
#'
#' Get a design-specific \code{DESeq2} analysis heat map prefix.
#'
#' @param design_name A \code{character} scalar with the design name.
#'
#' @return A \code{character} scalar with the \code{DESeq2} heat map prefix.
#' @export
#'
#' @examples
#' \dontrun{
#' design_name <- "global"
#' prefix_deseq_heatmap <- bsfrd_get_prefix_heatmap(
#'   design_name = design_name)
#' }
bsfrd_get_prefix_heatmap <- function(design_name) {
  return(paste("rnaseq",
               "deseq",
               design_name,
               "heatmap",
               sep = "_"))
}

#' Get a DESeq2 Volcano Prefix.
#'
#' Get a design-specific \code{DESeq2} analysis volcano prefix.
#'
#' @param design_name A \code{character} scalar with the design name.
#'
#' @return A \code{character} scalar with the \code{DESeq2} volcano prefix.
#' @export
#'
#' @examples
#' \dontrun{
#' design_name <- "global"
#' prefix_deseq_volcano <- bsfrd_get_prefix_volcano(
#'   design_name = design_name)
#' }
bsfrd_get_prefix_volcano <- function(design_name) {
  return(paste("rnaseq",
               "deseq",
               design_name,
               "volcano",
               sep = "_"))
}

#' Read a DESeq2 Contrast Tibble.
#'
#' Read a \code{DESeq2} analysis contrasts \code{tibble} from a tab-separated
#' value file and automatically sub-set to a particular design.
#'
#' Contrast \code{tibble}:
#' \describe{
#' \item{Design}{A \code{character} with the design name.}
#' \item{Numerator}{A \code{character} with the numerator as of \code{DESeq2::resultNames()}.}
#' \item{Denominator}{A \code{character} with the denominator as of \code{DESeq2::resultNames()}.}
#' \item{Label}{A \code{character} with a human-readable label.}
#' \item{Exclude}{A \code{logical} to exclude the design from reporting.}
#' }
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param summary A \code{logical} scalar to load a contrast summary rather than
#'   a contrast \code{tibble}.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{tibble} with contrast information.
#' @export
#' @importFrom rlang .data .env
#'
#' @examples
#' \dontrun{
#' contrast_tibble <- bsfrd_read_contrast_tibble(
#'   genome_directory = genome_directory,
#'   design_name = design_name,
#'   summary = FALSE,
#'   verbose = FALSE)
#' }
bsfrd_read_contrast_tibble <-
  function(genome_directory,
           design_name,
           summary = FALSE,
           verbose = FALSE) {
    prefix_deseq <- bsfrd_get_prefix_deseq(design_name = design_name)

    file_path_components <- c(prefix_deseq, "contrasts")

    col_types <- readr::cols(
      Design = readr::col_character(),
      Numerator = readr::col_character(),
      Denominator = readr::col_character(),
      Label = readr::col_character(),
      Exclude = readr::col_logical()
    )

    if (summary) {
      file_path_components <- c(file_path_components, "summary")

      col_types$cols <-
        c(
          col_types$cols,
          readr::cols(
            Significant = readr::col_integer(),
            SignificantUp = readr::col_integer(),
            SignificantDown = readr::col_integer()
          )$cols
        )
    }

    if (verbose) {
      message("Loading a contrast tibble ...")
    }

    contrast_tibble <- readr::read_tsv(
      file = file.path(
        genome_directory,
        prefix_deseq,
        paste(paste(file_path_components, collapse = "_"),
              "tsv",
              sep = ".")
      ),
      col_names = TRUE,
      col_types = col_types
    )
    rm(col_types, file_path_components, prefix_deseq)

    # Subset to the selected design.

    return(dplyr::filter(.data = contrast_tibble, .data$Design == .env$design_name))
  }

#' Get a DESeq2 Contrast List.
#'
#' Get a named \code{list} describing a particular contrast of a \code{DESeq2}
#' analysis.
#'
#' @param contrast_tibble A \code{tibble} with "Numerator" and "Denominator"
#'   variables.
#' @param index An \code{integer} scalar pointing at a particular \code{tibble}
#'   row.
#'
#' @return A named \code{list} of "numerator" and "denominator" \code{character}
#'   vectors. \code{NA} values in the contrast \code{tibble} are replaced by
#'   empty \code{character} vectors.
#' @export
#'
#' @examples
#' \dontrun{
#' contrast_list <- bsfrd_get_contrast_list(
#'   contrast_tibble = contrast_tibble,
#'   index = 1L)
#' }
bsfrd_get_contrast_list <- function(contrast_tibble, index) {
  numerator_character <-
    stringr::str_split(string = contrast_tibble$Numerator[index],
                       pattern = stringr::fixed(pattern = ","))[[1L]]

  denomintor_character <-
    stringr::str_split(string = contrast_tibble$Denominator[index],
                       pattern = stringr::fixed(pattern = ","))[[1L]]

  character_list <-
    list("numerator" = if (length(x = numerator_character == 1L) &&
                           is.na(x = numerator_character)) {
      character()
    } else {
      numerator_character
    },
    "denominator" = if (length(x = denomintor_character) == 1L &&
                        is.na(x = denomintor_character)) {
      character()
    } else {
      denomintor_character
    })
  rm(denomintor_character, numerator_character)

  return(character_list)
}

#' Get a DESeq2 Contrast Character Scalar.
#'
#' Get a \code{character} scalar describing a particular contrast of a
#' \code{DESeq2} analysis.
#'
#' @param contrast_tibble A \code{tibble} with Numerator and Denominator
#'   variables.
#' @param index An \code{integer} scalar pointing at a particular \code{tibble}
#'   row.
#'
#' @return A \code{character} scalar describing a particular contrast.
#' @export
#'
#' @examples
#' \dontrun{
#' contrast_character <- bsfrd_get_contrast_character(
#'   contrast_tibble = contrast_tibble,
#'   index = 1L)
#' }
bsfrd_get_contrast_character <- function(contrast_tibble, index) {
  contrast_list <-
    bsfrd_get_contrast_list(contrast_tibble = contrast_tibble, index = index)
  return(paste(
    paste(contrast_list$numerator, collapse = "_"),
    "against",
    if (length(x = contrast_list$denominator) == 0L ||
        is.na(contrast_list$denominator)) {
      "intercept"
    } else {
      paste(contrast_list$denominator, collapse = "_")
    },
    sep = "_"
  ))
}

#' Read a DESeq2 Design Tibble.
#'
#' Read a \code{DESeq2} analysis design \code{tibble} from a tab-separated value
#' file and automatically sub-set to a particular design.
#'
#' Design tibble:
#' \describe{
#' \item{design}{A \code{character} with the design name.}
#' \item{exclude}{A \code{logical} to exclude the design from reporting.}
#' \item{full_formula}{A \code{character} with the full model formula.}
#' \item{reduced_formulas}{A \code{character} with a comma-separated list of reduced model formulas.}
#' \item{factor_levels}{A \code{character} with semicolon-separated factors and their levels.
#' e.g. factor_levels="factor_1:level_1,level_2;factor_2:level_A,level_B"}
#' \item{plot_aes}{A \code{character} with \code{ggplot2} aesthetics.}
#' }
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{tibble} with design information.
#' @export
#' @importFrom rlang .data .env
#'
#' @examples
#' \dontrun{
#' design_tibble <- bsfrd_read_design_tibble(
#'   genome_directory = genome_directory,
#'   design_name = design_name,
#'   verbose = FALSE)
#' }
bsfrd_read_design_tibble <-
  function(genome_directory, design_name, verbose = FALSE) {
    prefix_deseq <- bsfrd_get_prefix_deseq(design_name = design_name)

    if (verbose) {
      message("Loading a design tibble ...")
    }

    design_tibble <- readr::read_tsv(
      file = file.path(
        genome_directory,
        prefix_deseq,
        paste(paste(prefix_deseq, "designs", sep = "_"), "tsv", sep = ".")
      ),
      col_names = TRUE,
      col_types = readr::cols(
        design = readr::col_character(),
        exclude = readr::col_logical(),
        full_formula = readr::col_character(),
        reduced_formulas = readr::col_character(),
        factor_levels = readr::col_character(),
        plot_aes = readr::col_character()
      )
    )
    rm(prefix_deseq)

    return(dplyr::filter(.data = design_tibble, .data$design == .env$design_name))
  }

#' Read a DESeq2 Design List.
#'
#' Read a \code{DESeq2} analysis design \code{tibble} from a tab-separated value file,
#' automatically sub-set to a particular design and return as a named \code{list}.
#'
#' @param genome_directory A \code{character} scalar with the genome directory.
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A named \code{list} with design information.
#' \describe{
#' \item{design}{A \code{character} with the design name.}
#' \item{exclude}{A \code{logical} to exclude the design from reporting.}
#' \item{full_formula}{A \code{character} with the full model formula.}
#' \item{reduced_formulas}{A \code{character} with a comma-separated list of reduced model formulas.}
#' \item{factor_levels}{A \code{character} with semicolon-separated factors and their levels.
#' e.g. factor_levels="factor_1:level_1,level_2;factor_2:level_A,level_B"}
#' \item{plot_aes}{A \code{character} with \code{ggplot2} aesthetics.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' design_list <- bsfrd_read_design_list(
#'   genome_directory = genome_directory,
#'   design_name = design_name,
#'   verbose = FALSE)
#' }
bsfrd_read_design_list <-
  function(genome_directory, design_name, verbose = FALSE) {
    return(as.list(
      x = bsfR::bsfrd_read_design_tibble(
        genome_directory = genome_directory,
        design_name = design_name,
        verbose = verbose
      )
    ))
  }

#' Process a Factor Specification.
#'
#' Private function to process a factor specification by splitting factor levels
#' and setting a "factor_name" attribute.
#'
#' @param factor_specification A \code{character} vector of exactly two
#'   components, the factor name [1L] and comma-separated factor levels [2L].
#'
#' @return A \code{character} vector with factor levels and an attribute
#'   "factor" specifying the factor name.
#' @noRd
#'
#' @examples
#' \dontrun{
#' factor_levels <- .bsfrd_process_factor_specification(
#'   factor_specification = "factor_name:level_1,level_2")
#' }
.bsfrd_process_factor_specification <-
  function(factor_specification) {
    # Split the second component of factor_specification, the factor levels, on ",".
    factor_levels <-
      stringr::str_split(string = factor_specification[2L],
                         pattern = stringr::fixed(pattern = ","))[[1L]]

    # Set the first component of factor_specification, the factor name, as attribute.
    attr(x = factor_levels, which = "factor_name") <-
      factor_specification[1L]

    return(factor_levels)
  }

#' Read a DESeq2 Sample Annotation DataFrame.
#'
#' Read a \code{DESeq2} analysis sample annotation \code{S4Vectors::DataFrame}
#' from a tab-separated value file, select only those samples that match the
#' \code{design_name} and re-level \code{factor} vectors according to the
#' \code{factor_level} specification or the design \code{list}, automatically
#' loaded via \code{bsfR::bsfrd_read_design_list()}.
#'
#' Sample DataFrame:
#' \describe{
#' \item{bam_path}{A \code{character} vector of BAM file paths.}
#' \item{bai_path}{A \code{character} vector of BAI file paths.}
#' \item{sample}{A \code{character} vector of sample names.}
#' \item{run}{A \code{character} vector of original sample names. Optional. If
#' present, indicates that technical replicates should be collapsed according to
#' information in the "sample" variable. The "run"variable provides the original
#' sample name before collapsing technical replicates.}
#' \item{designs}{A \code{character} vector of comma-separated values of
#' designs, a particular sample should be part of.}
#' \item{library_type}{A \code{factor} vector with levels "unstranded", "first"
#' and "second" to indicate the strand-orientation of the RNA-seq protocol and
#' whether the first or second strand gets sequenced. Illumina TruSeq stranded
#' mRNA sequences the second strand so that reads need inverting before counting
#' strand-specifically.}
#' \item{sequencing_type}{A \code{factor} vector with levels "SE" and "PE"
#' indicating single-end or paired-end sequencing, respectively and thus
#' counting as read pairs or not.}
#' \item{total_counts}{An \code{integer} vector with total counts per sample.
#' Calculated automatically based on the colSums() of the counts() function.}
#' \item{RIN}{A \code{numeric} vector providing the RNA integrity number (RIN)
#' score per sample. If available, the RIN score distribution will be plotted.}
#' }
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param factor_levels A \code{character} vector with a packed string to assign
#'   factor levels. If \code{NULL} will be read via \code{bsfR::bsfrd_read_design_list()}.
#' @param verbose A \code{logical} scalar to emit messages.
#' @return A \code{S4Vectors::DataFrame} with sample annotation.
#' @export
#'
#' @examples
#' \dontrun{
#' sample_frame <- bsfrd_read_sample_frame(
#'   genome_directory = genome_directory,
#'   design_name = design_name,
#'   factor_levels="factor_1:level_1,level_2;factor_2:level_A,level_B",
#'   verbose = TRUE
#' )
#' }
bsfrd_read_sample_frame <-
  function(genome_directory,
           design_name,
           factor_levels = NULL,
           verbose = FALSE) {
    prefix_deseq <-
      bsfrd_get_prefix_deseq(design_name = design_name)

    # If factor_levels is not defined, it can be read from the design_list.
    if (is.null(x = factor_levels)) {
      design_list <-
        bsfR::bsfrd_read_design_list(
          genome_directory = genome_directory,
          design_name = design_name,
          verbose = verbose
        )
      factor_levels <- design_list$factor_levels
      rm(design_list)
    }

    # Read the BSF Python sample TSV file as a data.frame and convert into a
    # S4Vectors::DataFrame. Import strings as factors and cast to character
    # vectors where required.
    if (verbose) {
      message("Loading a sample S4Vectors::DataFrame ...")
    }

    mcols_frame <-
      methods::as(
        object = utils::read.table(
          file = file.path(
            genome_directory,
            prefix_deseq,
            paste(paste(prefix_deseq, "samples", sep = "_"), "tsv", sep = ".")
          ),
          header = TRUE,
          sep = "\t",
          comment.char = "",
          stringsAsFactors = TRUE
        ),
        "DataFrame"
      )
    rownames(x = mcols_frame) <- mcols_frame$sample

    # Select only those samples, which have the design name annotated in the
    # designs variable split into a character vector.
    index_logical <-
      purrr::map_lgl(
        .x = stringr::str_split(
          string = as.character(mcols_frame$designs),
          pattern = stringr::fixed(pattern = ",")
        ),
        .f = ~ design_name %in% .
      )
    mcols_frame <- mcols_frame[index_logical, , drop = FALSE]
    rm(index_logical)

    if (nrow(x = mcols_frame) == 0L) {
      stop("No sample remaining after selection for design name.")
    }

    # The sequencing_type and library_type variables are required to set options
    # for the GenomicAlignments::summarizeOverlaps() read counting function.

    if (!"sequencing_type" %in% names(x = mcols_frame)) {
      stop("A sequencing_type variable is missing from the sample annotation frame.")
    }

    if (!"library_type" %in% names(x = mcols_frame)) {
      stop("A library_type variable is missing from the sample annotation frame.")
    }

    # Re-level the library_type and sequencing_type variables.
    mcols_frame$library_type <-
      factor(x = mcols_frame$library_type,
             levels = c("unstranded", "first", "second"))

    mcols_frame$sequencing_type <-
      factor(x = mcols_frame$sequencing_type,
             levels = c("SE", "PE"))

    # The "factor_levels" variable of the design data frame specifies the order
    # of factor levels.
    #
    # factor_levels="factor_1:level_1,level_2;factor_2:level_A,level_B"
    #
    # Turn the factor_levels character scalar into a list of character vectors
    # by splitting factor specifications by ";" then factor names and levels by
    # ":" characters. Set the factor names  as "factor_name" attribute of the
    # character vectors of factor levels list components.
    factor_list <- purrr::map(
      .x = stringr::str_split(
        string = stringr::str_split(
          string = factor_levels[1L],
          pattern = stringr::fixed(pattern = ";")
        )[[1L]],
        pattern = stringr::fixed(pattern = ":")
      ),
      .f = .bsfrd_process_factor_specification
    )

    # Apply the factor levels to each factor.
    design_variables <- names(x = mcols_frame)
    for (i in seq_along(along.with = factor_list)) {
      factor_name <- attr(x = factor_list[[i]], which = "factor_name")
      if (!is.na(x = factor_name) && factor_name != "") {
        if (factor_name %in% design_variables) {
          mcols_frame[, factor_name] <-
            factor(x = as.character(x = mcols_frame[, factor_name]),
                   levels = factor_list[[i]])
          # Check for NA values in case a factor level was missing.
          if (any(is.na(x = mcols_frame[, factor_name]))) {
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

    # Drop any unused levels from the sample data frame before returning it.
    return(droplevels(x = mcols_frame))
  }

#' Read a RangedSummarizedExperiment Object.
#'
#' Read a pre-calculated \code{SummarizedExperiment::RangedSummarizedExperiment}
#' object or initialise it from a sample annotation sheet loaded via
#' \code{bsfR::bsfrd_read_sample_frame()}, the \code{gtf_path} and
#' \code{genome_version}.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param gtf_path A \code{character} scalar with the reference transcriptome
#'   GTF file path. Only required, if the
#'   \code{SummarizedExperiment::RangedSummarizedExperiment} object needs
#'   initialising.
#' @param genome_version A \code{character} scalar with the genome version of
#'   the reference GTF file. Only required, if the
#'   \code{SummarizedExperiment::RangedSummarizedExperiment} object needs
#'   initialising.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{SummarizedExperiment::RangedSummarizedExperiment} object or
#'   \code{NULL}.
#' @export
#'
#' @examples
#' \dontrun{
#' bsfrd_read_summarized_experiment <- bsfrd_read_summarized_experiment(
#'   genome_directory = genome_directory,
#'   design_name = "global",
#'   gtf_path = "Homo_sapiens.e100.gtf.gz",
#'   genome_version = "hg38",
#'   verbose = FALSE)
#' }
bsfrd_read_summarized_experiment <-
  function(genome_directory,
           design_name,
           gtf_path = NULL,
           genome_version = NULL,
           verbose = FALSE) {
    ranged_summarized_experiment <- NULL

    prefix_deseq <-
      bsfrd_get_prefix_deseq(design_name = design_name)

    file_path <-
      file.path(
        genome_directory,
        prefix_deseq,
        paste0(prefix_deseq, "_ranged_summarized_experiment.rds")
      )

    if (file.exists(file_path) &&
        file.info(file_path)$size > 0L) {
      if (verbose) {
        message("Loading a RangedSummarizedExperiment object ...")
      }

      ranged_summarized_experiment <-
        base::readRDS(file = file_path)
    } else {
      # Get a S4Vectors::DataFrame with sample annotation.

      sample_frame <-
        bsfR::bsfrd_read_sample_frame(
          genome_directory = genome_directory,
          design_name = design_name,
          factor_levels = NULL,
          verbose = verbose
        )

      if (verbose) {
        message("Reading reference GTF exon features ...")
      }

      # The DESeq2 and RNA-seq vignettes suggest using TxDB objects, but for the
      # moment, we need extra annotation provided by Ensembl GTF files.

      exon_granges <-
        rtracklayer::import(
          con = gtf_path,
          format = "gtf",
          genome = genome_version,
          feature.type = "exon"
        )

      # Convert (i.e. split) the GenomicRanges::GRanges object into a
      # GenomicRanges::GRangesList object by gene identifiers.

      gene_granges_list <-
        GenomicRanges::split(x = exon_granges,
                             f = S4Vectors::mcols(x = exon_granges)$gene_id)

      # Process per library_type and sequencing_type and merge the
      # SummarizedExperiment::RangedSummarizedExperiment objects.

      for (library_type in levels(x = sample_frame$library_type)) {
        for (sequencing_type in levels(x = sample_frame$sequencing_type)) {
          if (verbose) {
            message(
              "Processing library_type: ",
              library_type,
              " sequencing_type: ",
              sequencing_type
            )
          }

          sub_sample_frame <-
            sample_frame[(sample_frame$library_type == library_type) &
                           (sample_frame$sequencing_type == sequencing_type),]

          if (nrow(x = sub_sample_frame) == 0L) {
            rm(sub_sample_frame)
            next()
          }

          # Create a BamFileList object and set the samples as names.

          if (verbose) {
            message("Creating a BamFileList object ...")
          }

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

          if (verbose) {
            message("Creating a RangedSummarizedExperiment object ...")
          }

          sub_ranged_summarized_experiment <-
            GenomicAlignments::summarizeOverlaps(
              features = gene_granges_list,
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

          # Combine SummarizedExperiment::RangedSummarizedExperiment objects
          # with the same GenomicRanges::GRanges, but different samples via
          # SummarizedExperiment::cbind().

          if (is.null(x = ranged_summarized_experiment)) {
            ranged_summarized_experiment <- sub_ranged_summarized_experiment
          } else {
            ranged_summarized_experiment <-
              SummarizedExperiment::cbind(ranged_summarized_experiment,
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
        if (verbose) {
          message("Collapsing technical replicates ...")
        }

        # To avoid mismatching column and row names between assay matrices and
        # the column data annotation, variables "sample" and "run" should be
        # used. So set the original samples as runs and rename the
        # collapsed_sample variable into the sample variable.

        ranged_summarized_experiment <- DESeq2::collapseReplicates(
          object = ranged_summarized_experiment,
          groupby = sample_frame$sample,
          run = sample_frame$run,
          renameCols = TRUE
        )
      }
      rm(sample_frame)

      # Calculate colSums() of SummarizedExperiment::assays()$counts and add as
      # total_count into the SummarizedExperiment::colData()
      # S4Vectors::DataFrame.

      sample_frame <-
        SummarizedExperiment::colData(x = ranged_summarized_experiment)

      sample_frame$total_counts <-
        base::colSums(
          x = SummarizedExperiment::assays(x = ranged_summarized_experiment)$counts,
          na.rm = TRUE
        )

      SummarizedExperiment::colData(x = ranged_summarized_experiment) <-
        sample_frame

      rm(gene_granges_list,
         exon_granges,
         sample_frame)

      base::saveRDS(object = ranged_summarized_experiment, file = file_path)
    }
    rm(file_path, prefix_deseq)

    return(ranged_summarized_experiment)
  }

#' Read a DESeqDataSet Object.
#'
#' Read a pre-calculated \code{DESeq2::DESeqDataSet} object.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{DESeq2::DESeqDataSet} object or \code{NULL}.
#' @export
#'
#' @examples
#' \dontrun{
#' deseq_data_set <- bsfrd_read_deseq_data_set(
#'   genome_directory = genome_directory,
#'   design_name = design_name,
#'   verbose = FALSE)
#' }
bsfrd_read_deseq_data_set <-
  function(genome_directory, design_name, verbose = FALSE) {
    deseq_data_set <- NULL
    prefix_deseq <-
      bsfrd_get_prefix_deseq(design_name = design_name)

    file_path <-
      file.path(genome_directory,
                prefix_deseq,
                paste0(prefix_deseq, "_deseq_data_set.rds"))
    if (file.exists(file_path) &&
        file.info(file_path)$size > 0L) {
      if (verbose) {
        message("Loading a DESeqDataSet object ...")
      }
      deseq_data_set <- base::readRDS(file = file_path)
    } else {
      warning("Require a pre-calculated DESeqDataSet object in file: ",
              file_path)
    }
    rm(file_path, prefix_deseq)

    return(deseq_data_set)
  }

#' Read a DESeqTransform Object.
#'
#' Read a previously saved "blind" or "model" \code{DESeq2::DESeqTransform} object.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param model A \code{logical} scalar to retrieve a model aware (\code{TRUE})
#'   or a blind (\code{FALSE}) \code{DESeq2::DESeqTransform} object.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{DESeq2::DESeqTransform} object or \code{NULL}.
#' @export
#'
#' @examples
#' \dontrun{
#' deseq_transform <- bsfrd_read_deseq_transform(
#'   genome_directory = genome_directory,
#'   design_name = design_name,
#'   model = TRUE,
#'   verbose = FALSE)
#' }
bsfrd_read_deseq_transform <-
  function(genome_directory,
           design_name,
           model = TRUE,
           verbose = FALSE) {
    deseq_transform <- NULL
    suffix <- if (model) {
      "model"
    } else {
      "blind"
    }
    prefix_deseq <-
      bsfrd_get_prefix_deseq(design_name = design_name)
    file_path <-
      file.path(genome_directory,
                prefix_deseq,
                paste(
                  paste(prefix_deseq, "deseq", "transform", suffix, sep = "_"),
                  "rds",
                  sep = "."
                ))
    if (file.exists(file_path) &&
        file.info(file_path)$size > 0L) {
      if (verbose) {
        message("Loading a ", suffix, " DESeqTransform object ...")
      }
      deseq_transform <- base::readRDS(file = file_path)
    } else {
      warning("Require a pre-calculated DESeqTransform object in file: ",
              file_path)
    }
    rm(file_path, prefix_deseq, suffix)

    return(deseq_transform)
  }

#' Read a DESeqResults Object.
#'
#' Read a previously saved \code{DESeq2::DESeqResults} object.
#'
#' Either \code{contrast_tibble} and \code{index} or just a (valid)
#' \code{contrast_character} option are required. The \code{contrast_character}
#' takes precedence.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param contrast_tibble A \code{tibble} with "Numerator" and "Denominator"
#'   variables.
#' @param index An \code{integer} scalar pointing at a particular \code{tibble}
#'   row.
#' @param contrast_character A \code{character} scalar specifying the contrast.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{DESeq2::DESeqResults} object for a particular contrast or
#'   \code{NULL}.
#' @export
#'
#' @examples
#' \dontrun{
#' deseq_results <- bsfrd_read_deseq_results(
#'   genome_directory = genome_directory,
#'   design_name = design_name,
#'   contrast_tibble = contrast_tibble,
#'   index = index,
#'   verbose = FALSE)
#'
#' deseq_results <- bsfrd_read_deseq_results(
#'   genome_directory = genome_directory,
#'   design_name = design_name,
#'   contrast_character = contrast_character,
#'   verbose = FALSE)
#' }
bsfrd_read_deseq_results <-
  function(genome_directory,
           design_name,
           contrast_tibble = NULL,
           index = NULL,
           contrast_character = NULL,
           verbose = FALSE) {
    deseq_results <- NULL
    prefix_deseq <-
      bsfrd_get_prefix_deseq(design_name = design_name)
    if (is.null(x = contrast_character)) {
      if (is.null(x = contrast_tibble) || is.null(x = index)) {
        warning(
          "Either a contrast_tibble and index or a (valid) contrast_character option are required."
        )
        return(NULL)
      }
      contrast_character <-
        bsfrd_get_contrast_character(contrast_tibble = contrast_tibble, index = index)
    }
    file_path <-
      file.path(genome_directory,
                prefix_deseq,
                paste(
                  paste(
                    prefix_deseq,
                    "contrast",
                    contrast_character,
                    "results",
                    sep = "_"
                  ),
                  "rds",
                  sep = "."
                ))

    if (file.exists(file_path) &&
        file.info(file_path)$size > 0L) {
      if (verbose) {
        message("Loading a DESeq2::DESeqResults object for contrast: ",
                contrast_character)
      }

      deseq_results <- base::readRDS(file = file_path)
    } else {
      warning("Missing DESeq2::DESeqResults object for contrast: ",
              contrast_character)
    }
    rm(file_path, prefix_deseq)

    return(deseq_results)
  }

#' Read a DESeqResults Tibble.
#'
#' Read a previously saved \code{DESeq2::DESeqResults} \code{tibble}.
#'
#' Either \code{contrast_tibble} and \code{index} or just a (valid)
#' \code{contrast_character} option are required. The \code{contrast_character}
#' takes precedence.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param contrast_tibble A \code{tibble} with "Numerator" and "Denominator"
#'   variables.
#' @param index An \code{integer} scalar pointing at a particular \code{tibble}
#'   row.
#' @param contrast_character A \code{character} scalar specifying the contrast.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{tibble} of \code{DESeq2::DESeqResults} for a particular contrast.
#' @export
#'
#' @examples
#' \dontrun{
#' result_tibble <- bsfrd_read_result_tibble(
#'   genome_directory = genome_directory,
#'   design_name = design_name,
#'   contrast_tibble = contrast_tibble,
#'   index = index,
#'   verbose = FALSE)
#'
#' result_tibble <- bsfrd_read_result_tibble(
#'   genome_directory = genome_directory,
#'   design_name = design_name,
#'   contrast_character = contrast_character,
#'   verbose = FALSE)
#' }
bsfrd_read_result_tibble <-
  function(genome_directory,
           design_name,
           contrast_tibble = NULL,
           index = NULL,
           contrast_character = NULL,
           verbose = FALSE) {
    deseq_results_tibble <- NULL
    prefix_deseq <-
      bsfrd_get_prefix_deseq(design_name = design_name)
    if (is.null(x = contrast_character)) {
      if (is.null(x = contrast_tibble) || is.null(x = index)) {
        warning(
          "Either a contrast_tibble and index or a (valid) contrast_character option are required."
        )
        return(NULL)
      }
      contrast_character <-
        bsfrd_get_contrast_character(contrast_tibble = contrast_tibble, index = index)
    }
    file_path <-
      file.path(genome_directory,
                prefix_deseq,
                paste(
                  paste(prefix_deseq,
                        "contrast",
                        contrast_character,
                        "genes",
                        sep = "_"),
                  "tsv",
                  sep = "."
                ))

    if (file.exists(file_path) &&
        file.info(file_path)$size > 0L) {
      if (verbose) {
        message("Loading a DESeqResults tibble for contrast: ",
                contrast_character)
      }

      deseq_results_tibble <-
        readr::read_tsv(
          file = file_path,
          col_types = readr::cols(
            gene_id = readr::col_character(),
            gene_version = readr::col_double(),
            gene_name = readr::col_character(),
            gene_biotype = readr::col_character(),
            gene_source = readr::col_character(),
            location = readr::col_character(),
            baseMean = readr::col_double(),
            log2FoldChange = readr::col_double(),
            lfcSE = readr::col_double(),
            # stat = readr::col_double(), # No longer available after lfcShrink().
            pvalue = readr::col_double(),
            padj = readr::col_double(),
            significant = readr::col_character(),
            rank_log2_fold_change = readr::col_double(),
            rank_base_mean = readr::col_double(),
            rank_padj = readr::col_double(),
            max_rank = readr::col_double()
          )
        )
    } else {
      warning("Missing DESeqResults tibble for contrast: ",
              contrast_character)
    }
    rm(file_path, prefix_deseq)

    return(deseq_results_tibble)
  }

#' Read a Feature Annotation Tibble.
#'
#' Read a previously saved feature annotation \code{tibble} from a tab-separated
#' value file or import it from a GTF file.
#'
#' Features are extracted from a transcriptome reference GTF file. For the
#' moment, Ensembl-specific files with "gene", "transcript" and "exon" features
#' are supported.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param feature_types A \code{character} vector of GTF feature types to be
#'   imported. Defaults to "genes".
#' @param gtf_file_path A \code{character} scalar with the reference GTF file
#'   path.
#' @param genome A \code{character} scalar with the genome version or a
#'   \code{GenomeInfoDb::Seqinfo} object.
#' @param verbose A \code{logical} scalar to emit messages.
#' @return A \code{tibble} with feature annotation.
#' @export
#'
#' @examples
#' \dontrun{
#' annotation_tibble <- bsfR::bsfrd_read_annotation_tibble(
#'   genome_directory = genome_directory,
#'   design_name = design_name,
#'   feature_types = "gene",
#'   gtf_file_path = gtf_file_path,
#'   genome = "hg38",
#'   verbose = TRUE
#' )
#' }
bsfrd_read_annotation_tibble <-
  function(genome_directory,
           design_name,
           feature_types = "gene",
           gtf_file_path = NULL,
           genome = NA,
           verbose = FALSE) {
    stopifnot(all(feature_types %in% c("gene", "transcript", "exon")))

    prefix_deseq <-
      bsfrd_get_prefix_deseq(design_name = design_name)

    # Load pre-existing gene annotation tibble or create it from the reference
    # GTF file.
    annotation_tibble <- NULL

    file_path <-
      file.path(genome_directory,
                prefix_deseq,
                paste(
                  paste(prefix_deseq, "annotation", paste(sort(feature_types), collapse = "_"), sep = "_"),
                  "tsv",
                  sep = "."
                ))
    if (file.exists(file_path) &&
        file.info(file_path)$size > 0L) {
      if (verbose) {
        message("Loading an annotation tibble ...")
      }

      col_types <- readr::cols(
        "gene_id" = readr::col_character(),
        "gene_version" = readr::col_character(),
        "gene_name" = readr::col_character(),
        "gene_biotype" = readr::col_character(),
        "gene_source" = readr::col_character()
      )

      if (any(c("transcript", "exon") %in% feature_types)) {
        col_types$cols <-
          c(
            col_types$cols,
            readr::cols(
              "transcript_id" = readr::col_character(),
              "transcript_version" = readr::col_character(),
              "transcript_name" = readr::col_character(),
              "transcript_biotype" = readr::col_character(),
              "transcript_source" = readr::col_character()
            )$cols
          )
      }

      if ("exon" %in% feature_types) {
        col_types$cols <-
          c(
            col_types$cols,
            readr::cols(
              "exon_id" = readr::col_character(),
              "exon_version" = readr::col_character(),
              "exon_number" = readr::col_character()
            )$cols
          )
      }

      annotation_tibble <-
        readr::read_tsv(file = file_path,
                        col_names = TRUE,
                        col_types = col_types)

      rm(col_types)
    } else {
      if (verbose) {
        message("Reading reference GTF features ...")
      }

      if (is.null(x = gtf_file_path)) {
        stop("Missing gtf_file_path option")
      }

      granges_object <-
        rtracklayer::import(
          con = gtf_file_path,
          format = "gtf",
          genome = genome,
          feature.type = feature_types
        )

      if (verbose) {
        message("Creating an annotation tibble ...")
      }
      variable_names <- c("gene_id",
                          "gene_version",
                          "gene_name",
                          "gene_biotype",
                          "gene_source")

      if (any(c("transcript", "exon") %in% feature_types)) {
        variable_names <- c(
          variable_names,
          "transcript_id",
          "transcript_version",
          "transcript_name",
          "transcript_biotype",
          "transcript_source"
        )
      }

      if ("exon" %in% feature_types) {
        variable_names <-
          c(variable_names,
            "exon_id",
            "exon_version",
            "exon_number")
      }

      mcols_frame <-
        S4Vectors::mcols(x = granges_object)[, variable_names]
      rm(variable_names)

      # Add the location as an Ensembl-like location, lacking the coordinate
      # system name and version.
      mcols_frame$location <-
        methods::as(object = granges_object, Class = "character")
      rm(granges_object)

      annotation_tibble <-
        tibble::as_tibble(x = as.data.frame(x = mcols_frame))
      rm(mcols_frame)

      readr::write_tsv(x = annotation_tibble, file = file_path)
    }
    rm(file_path)

    return(annotation_tibble)
  }

#' Read a Gene Set Tibble.
#'
#' Read a gene set \code{tibble} for gene annotation or selection.
#'
#' The \code{tibble} should have the following variables:
#' \describe{
#' \item{gene_id}{The Ensembl gene identifier from the annotation \code{tibble}.}
#' \item{gene_name}{The offical gene symbol.}
#' \item{gene_label}{The gene label to be plotted instead of the official symbol.}
#' \item{plot_name}{The sub-plot (i.e. heat map) to apply this label to.}
#' }
#'
#' Missing "gene_id" values are filled in on the basis of "gene_name" values and
#' the annotation \code{tibble}, which in turn is based on the reference GTF
#' file. Missing "gene_label" values are then filled in on the basis of the
#' "gene_name" variable.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param gene_set_path A \code{character} scalar with the gene set file path.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{tibble} with gene set information.
#' @export
#'
#' @examples
#' \dontrun{
#' gene_set_tibble <- bsfrd_read_gene_set_tibble(
#'   genome_directory = genome_directory,
#'   design_name = design_name,
#'   gene_set_path = gene_set_path,
#'   verbose = FALSE)
#' }
bsfrd_read_gene_set_tibble <-
  function(genome_directory,
           design_name,
           gene_set_path,
           verbose = FALSE) {
    if (verbose) {
      message("Loading a gene set tibble ...")
    }

    gene_set_tibble <-
      readr::read_csv(
        file = gene_set_path,
        col_names = TRUE,
        col_types = readr::cols(
          gene_id = readr::col_character(),
          gene_name = readr::col_character(),
          gene_label = readr::col_character(),
          plot_name = readr::col_character()
        )
      )

    # Find all those observations in "gene_id" that are NA or empty.
    missing_ids <-
      is.na(x = gene_set_tibble$gene_id) |
      gene_set_tibble$gene_id == ""
    missing_indices <- which(x = missing_ids)
    if (length(x = missing_indices) > 0L) {
      # Read the central transcriptome annotation tibble.
      annotation_tibble <-
        bsfrd_read_annotation_tibble(genome_directory = genome_directory, design_name = design_name)
      # Associate empty "gene_id" values with corresponding "gene_name" values.
      missing_names <-
        gene_set_tibble$gene_name[missing_indices]
      # Reset the missing "gene_id" values, by matching missing names in the annotation_tibble.
      gene_set_tibble$gene_id[missing_indices] <-
        annotation_tibble$gene_id[match(x = missing_names, table = annotation_tibble$gene_name)]
      rm(missing_names, annotation_tibble)
    }

    # Find all those observations in "gene_label" that are NA or empty and
    # replace with "gene_name" values.
    missing_ids <-
      is.na(x = gene_set_tibble$gene_label) |
      gene_set_tibble$gene_label == ""
    gene_set_tibble$gene_label[missing_ids] <-
      gene_set_tibble$gene_name[missing_ids]

    rm(missing_indices, missing_ids)

    return(gene_set_tibble)
  }

#' Convert Aesthetic and Variable Specifications.
#'
#' Split a \code{character} vector of multiple aesthetic and variable pairs
#' separated by "=" characters. Convert into a named \code{list} with the
#' variables as components and the aesthetics as names.
#'
#' @param aes_var_character A \code{character} vector of multiple aesthetic and
#'   variable pairs separated by "=" characters.
#'
#' @return A named \code{list} of variables as components and aesthetics as
#'   names.
#' @noRd
#'
#' @examples
#' \dontrun{
#' aes_var_list <-
#'   .bsfrd_convert_aesthetic_variable(
#'      aes_var_character = c("colour=variable_A", "shape=variable_B"))
#' }
.bsfrd_convert_aesthetic_variable <- function(aes_var_character) {
  # Split the character vector on "=" characters.
  character_list <-
    stringr::str_split(string = aes_var_character,
                       pattern = stringr::fixed(pattern = "="))

  # Assign the variables [2L] as list components.
  aes_var_list <- purrr::map(.x = character_list, .f = ~ .[2L])

  # Assign the aesthetics [1L] as names.
  names(x = aes_var_list) <-
    purrr::map_chr(.x = character_list, .f = ~ .[1L])

  rm(character_list)

  return(aes_var_list)
}

#' Convert a Geometric and Aesthetic Specification.
#'
#' Split the second component of a two-component \code{character} vector with
#' geometric ([1L]) and multiple aesthetic and variable definitions ([2L]) into
#' a named list of aesthetic components named by the geometric.
#'
#' @param geom_aes_character A two-component \code{character} vector of
#'   geometric ([1L]) and multiple aesthetics and variable definitions
#'   ([2L]).
#'
#' @return A named \code{list} of \code{list} objects of aesthetic and variable
#'   information named by the geometric. The \code{list} contains a single
#'   component to allow for naming of the component.
#' @noRd
#'
#' @examples
#' \dontrun{
#' geom_aes_list <-
#'   .bsfrd_convert_geometric_aesthetic(
#'     geom_aes_character = c("geom_point", "colour=variable_A,shape=variable_B"))
#' }
.bsfrd_convert_geometric_aesthetic <- function(geom_aes_character) {
  # Split geometric definitions on "," characters and assign names (geometric
  # names) to the list components (aesthetic list).

  # Assign the aesthetics definitions [2L] as the list component.
  # Since only component [2L] gets split, the list contains just one component.
  geom_aes_list <- purrr::map(
    .x = stringr::str_split(
      string = geom_aes_character[2L],
      pattern = stringr::fixed(pattern = ",")
    ),
    .f = .bsfrd_convert_aesthetic_variable
  )

  # Assign the geometric [1L] as name to the only list component.
  names(x = geom_aes_list) <- geom_aes_character[1L]

  return(geom_aes_list)
}

#' Convert a Plot Specification.
#'
#' Split a \code{character} vector of a single plot specification into a
#' two-component \code{character} vector of geometric ([1L]) and multiple
#' aesthetic and variable definitions ([2L]).
#'
#' @param plot_character
#'
#' @return A named \code{list} of aesthetic specification \code{list} objects
#'   named by geometrics.
#' @noRd
#'
#' @examples
#' \dontrun{
#' plot_list <-
#'   .bsfrd_convert_plot(
#'     plot_character =
#'     "geom_point:colour=variable_A,shape=variable_B;geom_text:colour=variable_C,label=variable_D")
#' }
.bsfrd_convert_plot <- function(plot_character) {
  # Split the character vector with multiple plot specifications on ";"
  # characters into single plot specifications. Split each plot specification
  # into two-component character vectors with geometric [1L] and aesthetic [2L]
  # information. The resulting list contains a single component named by the geometric.
  single_geom_aes_list <- purrr::map(
    .x = stringr::str_split(
      string = stringr::str_split(
        string = plot_character[1L],
        pattern = stringr::fixed(pattern = ";")
      )[[1L]],
      pattern = stringr::fixed(pattern = ":")
    ),
    .f = .bsfrd_convert_geometric_aesthetic
  )

  # Flatten the list of single-component geometric lists into a single,
  # plot-specific list. Select the single list-component and assign geometric
  # names to the new list components.

  # Assign the first (and only) list component [[1L]] as new list components.
  plot_list <- purrr::map(.x = single_geom_aes_list, .f = ~ .[[1L]])

  # Assign the geometric names of the first (and only) list component [1L] as
  # new list component names.
  names(x = plot_list) <-
    purrr::map_chr(.x = single_geom_aes_list, .f = ~ names(x = .[1L]))

  rm(single_geom_aes_list)

  return(plot_list)
}


#' Convert Multiple Plot Specifications.
#'
#' Convert a \code{character} scalar encoding mutliple plot specificatios into a
#' \code{list}.
#'
#' @param plots_character A \code{character} scalar encoding multiple plot
#'   specifications.
#'
#' @return A \code{list} of plot specification \code{list} objects.
#' @export
#' @seealso bsfrd_plots_list_to_character
#'
#' @examples
#' \dontrun{
#' plot_list <-
#'   bsfR::bsfrd_plots_character_to_list(
#'     plots_character =
#'     "geom_point:colour=variable_A,shape=variable_B;geom_text:colour=variable_C,label=variable_D")
#' }
bsfrd_plots_character_to_list <- function(plots_character) {
  # Split multiple plot definitions separated by "|" characters.
  plot_list <- purrr::map(
    .x = stringr::str_split(
      string = plots_character[1L],
      pattern = stringr::fixed(pattern = "|")
    )[[1L]],
    .f = .bsfrd_convert_plot
  )

  return(plot_list)
}


#' Convert an Aesthetics List into a Character Scalar.
#'
#' Convert a \code{list} of variables as components and aesthetics as names into
#' a "_"_separated \code{character} scalar.
#'
#' @param aes_list A \code{list} of variables as components and aesthetics as
#'   names.
#'
#' @return A "_"_separated \code{character} scalar.
#' @noRd
#'
#' @examples
#' \dontrun{
#'  aes_character <-
#'    bsfR::.bsfrd_convert_aesthetics_list(
#'      aes_list = list("colour" = "variable_A", "shape" = "variable_B")
#'    )
#' }
.bsfrd_convert_aesthetics_list <- function(aes_list) {
  # This could also be achieved via unlist().
  aes_character <- purrr::map_chr(.x = aes_list, .f = ~ .)

  # Combine the aesthetics and the variables, before collapsing into a single
  # "_"-separated character scalar.
  return(paste(paste(
    names(x = aes_character), aes_character, sep = "_"
  ), collapse = "_"))
}

#' Convert a List of Geometrics and Aesthetics into a Character Scalar.
#'
#' Convert a \code{list} of \code{list} objects with aesthetics and variable
#' information named by geometrics into a "__"-separated \code{character}
#' scalar. Remove the "geom_" prefix to get a cleaner character representation
#' for file naming.
#'
#' @param geom_list A \code{list} of aesthetics and variables \code{list}
#'   objects as components and geometrics as names.
#'
#' @return A "__"-separated \code{character} scalar.
#' @noRd
#'
#' @examples
#' \dontrun{
#'  geom_character <-
#'    .bsfrd_convert_geometrics_list(
#'      geom_list = list(
#'        "geom_point" = list("colour" = "variable_A", "shape" = "variable_B"),
#'        "geom_text" = list("colour" = "variable_C", "shape" = "variable_D")
#'      )
#'    )
#' }
.bsfrd_convert_geometrics_list <- function(geom_list) {
  geom_character <-
    purrr::map_chr(.x = geom_list, .f = .bsfrd_convert_aesthetics_list)

  # Remove the "geom_" prefix from the geometrics, then combine them with the
  # aesthetics, before collapsing into a single "__"-separated character scalar.

  return(paste(paste(
    sub(
      pattern = 'geom_',
      replacement = '',
      x = names(x = geom_character)
    ),
    geom_character,
    sep = "_"
  ), collapse = "__"))
}

#' Convert a Plot Specification List into a Character Scalar.
#'
#' Convert a \code{list} of plot specifications into a \code{character} scalar.
#'
#' @param plots_list A \code{list} of plot specifications.
#'
#' @return A \code{character} scalar.
#' @export
#' @seealso bsfrd_plots_character_to_list
#'
#' @examples
#' \dontrun{
#'  plot_character <-
#'    bsfR::bsfrd_plots_list_to_character(
#'      plots_list = list(
#'        list(
#'          "geom_point" = list("colour" = "variable_A", "shape" = "variable_B"),
#'          "geom_text" = list("colour" = "variable_C", "shape" = "variable_D")
#'        ),
#'        list(
#'          "geom_point" = list("colour" = "variable_1", "shape" = "variable_2"),
#'          "geom_text" = list("colour" = "variable_3", "shape" = "variable_4")
#'        )
#'      )
#'    )
#' }
bsfrd_plots_list_to_character <- function(plots_list) {
  return(purrr::map_chr(.x = plots_list, .f = .bsfrd_convert_geometrics_list))
}

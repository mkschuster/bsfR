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

#' Get a design-specific DESeq2 analysis prefix.
#'
#' @param design_name A \code{character} scalar with the design name.
#'
#' @return A \code{character} scalar with the DESeq2 analysis prefix.
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

#' Get a design-specific Enrichr analysis prefix.
#'
#' @param design_name A \code{character} scalar with the design name.
#'
#' @return A \code{character} scalar with the DESeq2 Enrichr prefix.
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

#' Get a design-specific GO analysis prefix.
#'
#' @param design_name A \code{character} scalar with the design name.
#'
#' @return A \code{character} scalar with the DESeq2 GO prefix.
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

#' Get a design-specific DESeq2 heatmap prefix.
#'
#' @param design_name A \code{character} scalar with the design name.
#'
#' @return A \code{character} scalar with the DESeq2 heatmap prefix.
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

#' Get a design-specific DESeq2 volcano prefix.
#'
#' @param design_name A \code{character} scalar with the design name.
#'
#' @return A \code{character} scalar with the DESeq2 volcano prefix.
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

#' Read a DESeq2 analysis contrasts tibble and automatically sub-set to a
#' particular design.
#'
#' Contrast tibble:
#' \describe{
#' \item{Design}{A \code{character} with the design name}
#' \item{Numerator}{A \code{character} with the numerator as of \code{DESeq2::resultNames()}}
#' \item{Denominator}{A \code{character} with the denominator as of \code{DESeq2::resultNames()}}
#' \item{Label}{A \code{character} with a human-readable label}
#' \item{Exclude}{A \code{logical} to exclude the design from reporting}
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

#' Get a named \code{list} describing a particular contrast of a DESeq2
#' analysis.
#'
#' @param contrast_tibble A \code{tibble} with Numerator and Denominator
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

#' Get a \code{character} scalar describing a particular contrast of a DESeq2 analysis.
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

#' Read a DESeq2 analysis design tibble and automatically sub-set to a
#' particular design.
#'
#' Design tibble:
#' \describe{
#' \item{design}{A \code{character} with the design name}
#' \item{exclude}{A \code{logical} to exclude the design from reporting}
#' \item{full_formula}{A \code{character} with the full model formula}
#' \item{reduced_formulas}{A \code{character} with a comma-separated list of reduced model formulas}
#' \item{factor_levels}{A \code{character} with semicolon-separated factors and their levels}
#' \item{plot_aes}{A \code{character} with \code{ggplot2} aesthetics}
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

#' Read a DESeq2 analysis design tibble, automatically sub-set to a
#' particular design and return as a list.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A named \code{list} with design information.
#' \describe{
#' \item{design}{A \code{character} with the design name}
#' \item{exclude}{A \code{logical} to exclude the design from reporting}
#' \item{full_formula}{A \code{character} with the full model formula}
#' \item{reduced_formulas}{A \code{character} with a comma-separated list of reduced model formulas}
#' \item{factor_levels}{A \code{character} with semicolon-separated factors and their levels}
#'   factor_levels="factor_1:level_1,level_2;factor_2:level_A,level_B")
#' \item{plot_aes}{A \code{character} with \code{ggplot2} aesthetics}
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

#' Private function to match a design name in the vector of sample-specific
#' design names.
#'
#' @param design_names A \code{character} vector of sample-specific design names.
#' @param design_name A \code{character} scalar with a design name
#'
#' @return A \code{logical} scalar if the design name scalar is in the design
#'   names vector.
#' @noRd
#'
#' @examples
#' \dontrun{
#' design_logical <- .bsfrd_match_design_name(
#'   design_names=c("global", "test"),
#'   design_name="global")
#' }
.bsfrd_match_design_name <- function(design_names, design_name) {
  return(design_name %in% design_names)
}

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

#' Read a Sample Annotation DataFrame.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param factor_levels A \code{character} vector with a packed string to assign
#'   factor levels. If \code{NULL} will be read via \code{bsfrd_read_design_list()}.
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
      message("Loading sample DataFrame")
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
    # designs variable.
    index_logical <-
      unlist(
        x = lapply(
          X = stringr::str_split(
            string = as.character(mcols_frame$designs),
            pattern = stringr::fixed(pattern = ",")
          ),
          FUN = .bsfrd_match_design_name,
          design_name = design_name
        )
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
    # of factor levels. Turn the factor_levels variable into a list of character
    # vectors, where the factor names are set as attributes of the character
    # vectors of factor levels list components.
    #
    # factor_levels="factor_1:level_1,level_2;factor_2:level_A,level_B"
    factor_list <- lapply(
      # The "factor_levels" variable is a character vector, always with just a
      # single component.
      #
      # Split by the ";" then by the ":" character.
      X = stringr::str_split(
        string = stringr::str_split(
          string = factor_levels[1L],
          pattern = stringr::fixed(pattern = ";")
        )[[1L]],
        pattern = stringr::fixed(pattern = ":")
      ),
      FUN = .bsfrd_process_factor_specification
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

    # Drop any unused levels from the sample data frame before retrning it.
    return(droplevels(x = mcols_frame))
  }

#' Read a pre-calculated RangedSummarizedExperiment object.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{RangedSummarizedExperiment} object or \code{NULL}.
#' @export
#'
#' @examples
#' \dontrun{
#' bsfrd_read_summarized_experiment <- bsfrd_read_summarized_experiment(
#'   genome_directory = genome_directory,
#'   design_name = design_name,
#'   verbose = FALSE)
#' }
bsfrd_read_summarized_experiment <-
  function(genome_directory, design_name, verbose = FALSE) {
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
      warning("Require a pre-calculated RangedSummarizedExperiment object in file: ",
              file_path)
    }
    rm(file_path, prefix_deseq)

    return(ranged_summarized_experiment)
  }

#' Read a pre-calculated DESeqDataSet object.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{DESeqDataSet} object or \code{NULL}.
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

#' Read a previously saved "blind" or "model" DESeqTransform object.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param model A \code{logical} scalar to retrieve a model aware (\code{TRUE})
#'   or a blind (\code{FALSE}) \code{DESeqTransform} object.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{DESeqTransform} object or \code{NULL}.
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

#' Read a previously saved DESeqResults object.
#'
#' Either contrast_tibble and index or just a (valid) contrast_character option
#' are required. The contrast_character takes precedence.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param contrast_tibble A \code{tibble} with Numerator and Denominator
#'   variables.
#' @param index An \code{integer} scalar pointing at a particular \code{tibble}
#'   row.
#' @param contrast_character A \code{character} scalar specifying the contrast.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{DESeqResults} object for a particular contrast or
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
        message("Loading a DESeqResults object for contrast: ",
                contrast_character)
      }

      deseq_results <- base::readRDS(file = file_path)
    } else {
      warning("Missing DESeqResults object for contrast: ",
              contrast_character)
    }
    rm(file_path, prefix_deseq)

    return(deseq_results)
  }

#' Read a previously saved DESeqResults tibble.
#'
#' Either contrast_tibble and index or just a (valid) contrast_character option
#' are required. The contrast_character takes precedence.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param contrast_tibble A \code{tibble} with Numerator and Denominator
#'   variables.
#' @param index An \code{integer} scalar pointing at a particular \code{tibble}
#'   row.
#' @param contrast_character A \code{character} scalar specifying the contrast.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{tibble} of DESeqResults for a particular contrast.
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

#' Read a feature annotation tibble or import it from a GTF file.
#'
#' Features are extracted from a transcriptome reference GTF file. For the
#' moment, Ensembl-specific files with "gene", "transcript" and "exon" features
#' are supported.
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

      readr::write_tsv(x = annotation_tibble, path = file_path)
    }
    rm(file_path)

    return(annotation_tibble)
  }

#' Read a gene set tibble for gene annotation or selection.
#'
#' The tibble should have the following variables:
#'   gene_id:    The Ensembl gene identifier from the annotation tibble.
#'   gene_name:  The offical gene symbol.
#'   gene_label: The gene label to be plotted instead of the official symbol.
#'   plot_name:  The sub-plot (i.e. heatmap) to apply this label to.
#'
#' Missing 'gene_id' values are filled in on the basis of 'gene_name' values and
#' the annotation tibble, which in turn is based on the reference GTF file.
#' Missing 'gene_label' values are then filled in on the basis of the
#' 'gene_name' variable.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param design_name A \code{character} scalar with the design name.
#' @param gene_set_path A \code{character} scalar with the gene set file path.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{tibble} object of gene set information.
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

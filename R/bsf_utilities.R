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


# A BSF R library module of general utility functions.

# Functions ---------------------------------------------------------------


#' Summarise reports in subdirectories.
#'
#' List all directories under a directory_path and filter only those, which base
#' name matches the directory_pattern. Parse the sample name from the base name
#' and use the report_function to concatenate the directory path and the
#' constant report name before reading and parsing the report and return a list
#' of results to be bound into a tibble by rows.
#'
#' This function applies to reports with constant names in sample-specific
#' sub-directories, while the bsfu_summarise_report_files() function below
#' applies to sample-specific files.
#'
#' @param directory_path A \code{character} scalar with the the base directory
#'   path.
#' @param directory_pattern A \code{character} scalar with a regular expression
#'   to filter directories and extract sample names.
#' @param report_function A \code{function} reading and parsing report files
#'   before returning a \code{list} of results.
#' @param variable_name A \code{character} scalar with the variable name, which
#'   is parsed form the file_pattern and added to the results \code{tbl_df}.
#'   Defaults to "sample_name".
#' @param ... A arguments passed to the \code{purrr::map_dfr} function.
#'
#' @return A \code{tbl_df} of results depending on the result \code{list}
#'   objects returned by the report_function.
#' @export
#'
#' @examples
#' \dontrun{
#' result_tibble <- bsfu_summarise_report_directories(
#'   directory_path = ".",
#'   directory_pattern = "^rnaseq_tophat_(.*)$",
#'   report_function = parse_report,
#'   variable_name = "sample_name"
#' )
#' }
bsfu_summarise_report_directories <-
  function(directory_path,
           directory_pattern,
           report_function,
           variable_name = "sample_name",
           ...) {
    directory_paths <-
      base::list.dirs(path = directory_path,
                      full.names = TRUE,
                      recursive = FALSE)

    directory_names <- base::basename(path = directory_paths)

    # Find those directory paths matching the directory pattern.
    directory_indices <-
      base::grep(pattern = directory_pattern, x = directory_names)

    result_tibble <-
      dplyr::bind_rows(purrr::map(.x = directory_paths[directory_indices],
                                  .f = report_function,
                                  ... = ...))

    # Prepend the "variable_name" variable and append the "directory_path"
    # variable.

    result_tibble <- dplyr::bind_cols(
      "{variable_name}" := base::gsub(
        pattern = directory_pattern,
        replacement = "\\1",
        x = directory_names[directory_indices]
      ),
      result_tibble,
      "directory_path" = directory_paths[directory_indices]
    )

    rm(directory_indices, directory_names, directory_paths)

    return(result_tibble)
  }

#' Summarise reports in file paths.
#'
#' List all file paths under a directory_path and filter only those, which base
#' name matches the file_pattern. Parse the sample name from the base name
#' and use the report_function to concatenate the directory path an the constant
#' report name before reading and parsing the report and return a list of
#' results to be bound into a tibble by rows.
#'
#' This function applies to sample-specific reports, while the
#' bsfu_summarise_report_directories() function above applies to reports with
#' constant names in sample-specific sub-directories.
#'
#' @param directory_path A \code{character} scalar with the the base directory
#'   path.
#' @param file_pattern A \code{character} scalar with a regular expression to
#'   filter files and extract sample names.
#' @param report_function A \code{function} reading and parsing report files
#'   before returning a \code{list} of results.
#' @param variable_name A \code{character} scalar with the variable name, which
#'   is parsed form the file_pattern and added to the results \code{tbl_df}.
#'   Defaults to "sample_name".
#' @param ... A arguments passed to the \code{purrr::map_dfr} function.
#'
#' @return A \code{tbl_df} of results depending on the result \code{list}
#'   objects returned by the report_function.
#' @export
#'
#' @examples
#' \dontrun{
#' result_tibble <- bsfu_summarise_report_files(
#'   directory_path = ".",
#'   file_pattern = "^star_align_(.*)_Log\\.final\\.out$",
#'   report_function = parse_report,
#'   variable_name = "read_group"
#' )
#' }
bsfu_summarise_report_files <-
  function(directory_path,
           file_pattern,
           report_function,
           variable_name = "sample_name",
           ...) {
    file_paths <-
      base::list.files(
        path = directory_path,
        pattern = file_pattern,
        full.names = TRUE,
        recursive = TRUE
      )

    result_tibble <- dplyr::bind_rows(purrr::map(.x = file_paths,
                                                 .f = report_function))

    # Prepend the "variable_name" variable and append the "file_path" variable.

    result_tibble <- dplyr::bind_cols(
      "{variable_name}" := base::gsub(
        pattern = file_pattern,
        replacement = "\\1",
        x = base::basename(path = file_paths)
      ),
      result_tibble,
      "file_path" = file_paths
    )

    rm(file_paths)

    return(result_tibble)
  }

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
#' and use the report_function to concatenate the directory path an the constant
#' report name before reading and parsing the report and return a list of
#' results to be bound into a tibble by rows.
#'
#' @param directory_path A \code{character} scalar with the the base directory
#'   path.
#' @param directory_pattern A \code{character} scalar witha regular expression
#'   to filter directories and extract sample names.
#' @param report_function A \code{function} reading and parsing report files
#'   before returning a \code{list} of results.
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
#'   report_function = parse_report
#' )
#' }
bsfu_summarise_report_directories <-
  function(directory_path,
           directory_pattern,
           report_function,
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
      purrr::map_dfr(.x = directory_paths[directory_indices],
                     .f = report_function,
                     ... = ...)

    # Prepend the "sample_name" variable and append the "directory_path" variable.
    result_tibble <- dplyr::bind_cols(
      "sample_name" = base::gsub(
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


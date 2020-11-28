#
# Library module of bsfR RNA-seq Tuxedo functions.
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

#' Get a comparison-specific Tuxedo analysis prefix.
#'
#' @param comparison_name A \code{character} scalar with the comparison name.
#'
#' @return A \code{character} scalar with the Tuxedo analysis prefix.
#' @export
#'
#' @examples
#' comparison_name <- "global"
#' prefix_tuxedo <- bsfrt_get_prefix_tuxedo(
#'   comparison_name = comparison_name)
bsfrt_get_prefix_tuxedo <- function(comparison_name) {
  return(paste("rnaseq", "tuxedo", comparison_name, sep = "_"))
}

#' Get a comparison-specific Cuffdiff analysis prefix.
#'
#' @param comparison_name A \code{character} scalar with the comparison name.
#'
#' @return A \code{character} scalar with the Cuffdiff analysis prefix.
#' @export
#'
#' @examples
#' comparison_name <- "global"
#' prefix_cuffdiff <- bsfrt_get_prefix_cuffdiff(
#'   comparison_name = comparison_name)
bsfrt_get_prefix_cuffdiff <- function(comparison_name) {
  return(paste("rnaseq", "cuffdiff", comparison_name, sep = "_"))
}

#' Get a comparison-specific Process Cuffdiff analysis prefix.
#'
#' @param comparison_name A \code{character} scalar with the comparison name.
#'
#' @return A \code{character} scalar with the Process Cuffdiff analysis prefix.
#' @export
#'
#' @examples
#' comparison_name <- "global"
#' prefix_process_cuffdiff <- bsfrt_get_prefix_process_cuffdiff(
#'   comparison_name = comparison_name)
bsfrt_get_prefix_process_cuffdiff <- function(comparison_name) {
  return(paste("rnaseq", "process", "cuffdiff", comparison_name, sep = "_"))
}

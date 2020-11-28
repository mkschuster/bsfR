#
# Library module of bsfR RNA-seq Cufflinks functions.
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

#' Get a comparison-specific Cuffdiff analysis prefix.
#'
#' @param comparison_name A \code{character} scalar with the comparison name.
#'
#' @return A \code{character} scalar with the Cuffdiff analysis prefix.
#' @export
#'
#' @examples
#' comparison_name <- "global"
#' prefix_cuffdiff <- bsfrc_get_prefix_cuffdiff(
#'   comparison_name = comparison_name)
bsfrc_get_prefix_cuffdiff <- function(comparison_name) {
  return(paste("rnaseq",
               "cuffdiff",
               comparison_name,
               sep = "_"))
}

#' Get a comparison-specific process Cuffdiff analysis prefix.
#'
#' @param comparison_name A \code{character} scalar with the comparison name.
#'
#' @return A \code{character} scalar with the process Cuffdiff analysis prefix.
#' @export
#'
#' @examples
#' comparison_name <- "global"
#' prefix_process_cuffdiff <- bsfrc_get_prefix_process_cuffdiff(
#'   comparison_name = comparison_name)
bsfrc_get_prefix_process_cuffdiff <- function(comparison_name) {
  return(paste("rnaseq",
               "process",
               "cuffdiff",
               comparison_name,
               sep = "_"))
}

#' Read a comparison-specific gene annotation tibble.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param comparison_name A \code{character} scalar with the comparison name.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{tibble} with gene annotation information.
#' @export
#'
#' @examples
#' comparison_name <- "global"
#' genome_directory <- "."
#' gene_annotation_tibble <- bsfrc_read_gene_annotation_tibble(
#'   genome_directory = genome_directory,
#'   comparison_name = comparison_name,
#'   verbose = FALSE)
bsfrc_read_gene_annotation_tibble <-
  function(genome_directory,
           comparison_name,
           verbose = FALSE) {
    gene_annotation_tibble <- NULL

    prefix_process_cuffdiff <-
      bsfrc_get_prefix_process_cuffdiff(comparison_name = comparison_name)

    file_path <- file.path(genome_directory,
                           prefix_process_cuffdiff,
                           paste(
                             paste(prefix_process_cuffdiff, "genes", "annotation", sep = "_"),
                             "tsv",
                             sep = "."
                           ))

    if (file.exists(file_path) &&
        file.info(file_path)$size > 0L) {
      if (verbose) {
        message("Loading a gene annotation tibble ...")
      }

      gene_annotation_tibble <-
        readr::read_tsv(
          file = file_path,
          col_types = readr::cols(
            gene_id = readr::col_character(),
            gene_short_name = readr::col_character(),
            locus = readr::col_character(),
            transcript_ids = readr::col_character(),
            ensembl_gene_ids = readr::col_character(),
            ensembl_transcript_ids = readr::col_character()
          )
        )
    } else {
      warning("Require a pre-calculated annotation tibble in file: ",
              file_path)
    }
    rm(file_path, prefix_process_cuffdiff)

    return(gene_annotation_tibble)
  }

#' Read a comparison-specific isoform annotation tibble.
#'
#' @param genome_directory A \code{character} scalar with the genome directory
#'   path.
#' @param comparison_name A \code{character} scalar with the comparison name.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{tibble} with isoform annotation information.
#' @export
#'
#' @examples
#' comparison_name <- "global"
#' genome_directory <- "."
#' isoform_annotation_tibble <- bsfrc_read_isoform_annotation_tibble(
#'   genome_directory = genome_directory,
#'   comparison_name = comparison_name,
#'   verbose = FALSE)
bsfrc_read_isoform_annotation_tibble <- function(genome_directory,
                                                 comparison_name,
                                                 verbose = FALSE) {
  isoform_annotation_tibble <- NULL

  prefix_process_cuffdiff <-
    bsfrc_get_prefix_process_cuffdiff(comparison_name = comparison_name)

  file_path <- file.path(genome_directory,
                         prefix_process_cuffdiff,
                         paste(
                           paste(prefix_process_cuffdiff, "isoforms", "annotation", sep = "_"),
                           "tsv",
                           sep = "."
                         ))

  if (file.exists(file_path) &&
      file.info(file_path)$size > 0L) {
    if (verbose) {
      message("Loading a isoform annotation tibble ...")
    }

    isoform_annotation_tibble <-
      readr::read_tsv(file = file_path)
  } else {
    warning("Require a pre-calculated annotation tibble in file: ",
            file_path)
  }
  rm(file_path, prefix_process_cuffdiff)

  return(isoform_annotation_tibble)
}

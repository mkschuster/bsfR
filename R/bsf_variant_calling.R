#
# Library module of bsfR variant calling functions.
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

#' Get a variant calling analysis prefix.
#'
#' @return A \code{character} scalar with the variant calling analysis prefix.
#' @export
#'
#' @examples
#' \dontrun{
#' prefix_variant_calling <- bsfvc_get_prefix()
#' }
bsfvc_get_prefix <- function() {
  return(paste("variant",
               "calling",
               sep = "_"))
}

#' Import Ensembl annotation from a GTF file.
#'
#' A GTF file is imported into a GRanges object and optionally filtered for just
#' "basic" annotation, before optional flanks are added to both, start and end.
#' The GRanges object is then reduced to get a non-redundant set of transcribed
#' regions suitable for variant calling.
#'
#' Summary list:
#'   "transcribed_ranges" GRanges object of non-redundant transcribed genomic regions
#'   "exon_path" Ensembl GTF file path.
#'   "exon_number_raw" Total number of unfiltered Ensembl exon GRanges objects.
#'   "exon_width_raw" Total number of unfiltered Ensembl Exon bases.
#'   "exon_number" Total number of filtered Ensembl exon GRanges objects.
#'   "exon_width" Total number of filtered Ensembl Exon bases.
#'   "exon_flank_width" Total number of filtered Ensembl Exon bases including flanks.
#'   "transcribed_number" Total number of reduced Ensembl exon GRanges objects.
#'   "transcribed_width" Total number of reduced Ensembl Exon bases including flanks.
#'
#' @param exon_path A \code{character} scalar with the GTF file path.
#' @param genome_version A \code{character} scalar with the genome assembly
#'   version.
#' @param basic A \code{logical} scalar to filter for "basic" annotation.
#' @param exon_flanks A \code{integer} scalar to add exon flanks.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{list} with summary information.
#' @importFrom GenomicRanges reduce resize width
#' @importFrom S4Vectors mcols
#' @importFrom rtracklayer import
#' @export
#'
#' @examples
#' \dontrun{
#' ensembl_granges <- bsfvc_import_ensembl(
#'   exon_path = exon_path,
#'   genome_version = "hs37d5",
#'   basic = TRUE,
#'   exon_flanks = 0L,
#'   verbose = FALSE)
#' }
bsfvc_import_ensembl <-
  function(exon_path,
           genome_version = "hs37d5",
           basic = TRUE,
           exon_flanks = 0L,
           verbose = FALSE) {
    # Import Ensembl gene, transcript and exon annotation as GRanges vector object,
    # where only exon components are relevant for this analysis.
    summary_list <- list()
    summary_list$exon_path <- exon_path
    if (verbose) {
      message("Importing exon ranges: ", summary_list$exon_path)
    }

    exon_ranges <-
      rtracklayer::import(con = exon_path,
                          # format = "gtf",
                          genome = genome_version,
                          feature.type = "exon")

    summary_list$exon_number_raw <- length(x = exon_ranges)
    if (verbose) {
      message("Number of raw exon ranges: ", summary_list$exon_number_raw)
    }

    summary_list$exon_width_raw <-
      sum(GenomicRanges::width(x = exon_ranges))
    if (verbose) {
      message("Cumulative width of raw exon ranges: ",
              summary_list$exon_width_raw)
    }

    # Ensembl now annotates a "tag" in GFF files with value "basic" indicating
    # standard (basic) transcript models. Can the exon ranges be subset by such a
    # tag?
    if ("tag" %in% names(x = S4Vectors::mcols(x = exon_ranges)) &&
        basic) {
      if (verbose) {
        message("Filtering by GTF 'tag = \"basic\"' annotation")
      }
      # Use the %in% operator for character matching as it sets NA values to FALSE,
      # automatically.
      exon_ranges <-
        exon_ranges[S4Vectors::mcols(x = exon_ranges)$tag %in% "basic",]
    }

    summary_list$exon_number <- length(x = exon_ranges)
    if (verbose) {
      message("Number of exon ranges: ", summary_list$exon_number)
    }

    summary_list$exon_width <-
      sum(GenomicRanges::width(x = exon_ranges))
    if (verbose) {
      message("Cumulative width of exon ranges: ", summary_list$exon_width)
    }

    # Resize the GenomicRanges to apply flanking regions, by default 0L, to the
    # exon ranges.
    exon_ranges <-
      GenomicRanges::resize(
        x = exon_ranges,
        width = GenomicRanges::width(x = exon_ranges) + exon_flanks,
        fix = "end"
      )
    exon_ranges <-
      GenomicRanges::resize(
        x = exon_ranges,
        width = GenomicRanges::width(x = exon_ranges) + exon_flanks,
        fix = "start"
      )

    summary_list$exon_flank_width <-
      sum(GenomicRanges::width(x = exon_ranges))
    if (verbose) {
      message("Cumulative width of exon ranges with flanks: ",
              summary_list$exon_flank_width)
    }

    # Reduce the non-redundant Ensembl exons to their footprint on the genome to get
    # transcribed regions. The revmap column contains the mapping to original
    # exon_ranges components.
    if (verbose) {
      message("Reducing exon ranges to transcribed ranges.")
    }
    summary_list$transcribed_ranges <-
      GenomicRanges::reduce(
        x = exon_ranges,
        drop.empty.ranges = TRUE,
        with.revmap = TRUE,
        ignore.strand = TRUE
      )

    summary_list$transcribed_number <-
      length(x = summary_list$transcribed_ranges)
    if (verbose) {
      message("Number of transcribed ranges: ",
              summary_list$transcribed_number)
    }

    summary_list$transcribed_width <-
      sum(GenomicRanges::width(x = summary_list$transcribed_ranges))
    if (verbose) {
      message("Cumulative width of transcribed ranges: ",
              summary_list$transcribed_width)
    }

    return(summary_list)
  }

#' Import exome target annotation.
#'
#' Summary list:
#'   "target_path" Target region GTF file path
#'   "target_ranges" Target GRanges object
#'   "target_number_raw"
#' @param target_path A \code{character} scalar with the GTF file path.
#' @param genome_version A \code{character} scalar with the genome version.
#' @param verbose A \code{logical} scalar to emit messages.
#'
#' @return A \code{list} with summary information.
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges width
#' @export
#'
#' @examples
#' \dontrun{
#' target_granges <- bsfvc_import_targets(
#'   target_path = target_path,
#'   genome_version = "hs37d5",
#'   verbose = FALSE)
#' }
bsfvc_import_targets <-
  function(target_path,
           genome_version = NULL,
           verbose = FALSE) {
    # Read the file of targeted regions, if available. Although the GATK Callable
    # Loci analysis is generally only run on these target regions, this file
    # provides the target (probe) names of the enrichment design.
    summary_list <- list()
    summary_list$target_path <- target_path

    if (verbose) {
      message("Importing target range annotation: ",
              summary_list$target_path)
    }
    # The rtrackayer::import() function reads the genome version from the "db"
    # attribute of the BED "track" line.
    summary_list$target_ranges <-
      rtracklayer::import(con = summary_list$target_path,
                          # format = "gtf",
                          genome = genome_version)

    summary_list$target_number_raw <-
      length(x = summary_list$target_ranges)
    if (verbose) {
      message("Number of target ranges: ", summary_list$target_number_raw)
    }

    summary_list$target_width_raw <-
      sum(GenomicRanges::width(x = summary_list$target_ranges))
    if (verbose) {
      message("Cumulative width of target ranges: ",
              summary_list$target_width_raw)
    }

    return(summary_list)
  }

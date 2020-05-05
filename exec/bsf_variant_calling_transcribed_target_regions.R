#!/usr/bin/env Rscript
#
# BSF R script to constrain target enrichment regions to exon regions.
#
# All Ensembl "exon" features are imported from a GTF file and reduced into a
# (non-overlapping) set of transcribed regions of the genome. Optionally, target
# regions can be imported and overlapped with the transcribed regions above to
# export the minmal set of transcribed target regions in BED format.
#
# Without target regions, the set of transcribed regions is exported as a BED
# file.
#
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

suppressPackageStartupMessages(expr = library(package = "optparse"))

argument_list <-
  optparse::parse_args(object = optparse::OptionParser(
    option_list = list(
      optparse::make_option(
        opt_str = c("--verbose", "-v"),
        action = "store_true",
        default = TRUE,
        help = "Print extra output [default]",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = c("--quiet", "-q"),
        action = "store_false",
        default = FALSE,
        dest = "verbose",
        help = "Print little output",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = c("--exon-path"),
        dest = "exon_path",
        help = "File path to the Ensembl gene, transcript and exon annotation GTF",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--target-path"),
        dest = "target_path",
        help = "File path to the enrichment targets BED",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--output-path"),
        dest = "output_path",
        help = "File path for the transcribed target BED",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--genome-version"),
        dest = "genome_version",
        help = "Genome (assembly) version",
        type = "character"
      )
    )
  ))

suppressPackageStartupMessages(expr = library(package = "bsfR"))
suppressPackageStartupMessages(expr = library(package = "Biostrings"))

# Keep overall statistics of constrained GRanges in a list.
constrained_list <- list()

# Import Ensembl "exon" annotation as GTF.
ensembl_list <- bsfR::bsfvc_import_ensembl(
  exon_path = argument_list$exon_path,
  genome_version = argument_list$genome_version,
  basic = TRUE,
  exon_flanks = 0L,
  verbose = argument_list$verbose
)

# Read the file of target regions, if available. Although the GATK Callable
# Loci analysis is generally only run on these target regions, this file
# provides the target (probe) names of the enrichment design.

if (!is.null(x = argument_list$target_path)) {
  target_list <-
    bsfR::bsfvc_import_targets(
      target_path = argument_list$target_path,
      genome_version = argument_list$genome_version,
      verbose = argument_list$verbose
    )

  if (argument_list$verbose) {
    message("Overlapping target and transcribed GRanges.")
  }
  overlap_frame <-
    IRanges::mergeByOverlaps(query = target_list$target_ranges,
                             subject = ensembl_list$transcribed_ranges)
  # Adjust start and end to the minimally overlapping regions.
  constrained_list$constrained_ranges <-
    GRanges(
      seqnames = seqnames(x = overlap_frame$target_ranges),
      ranges = IRanges(
        start = pmax(
          start(x = overlap_frame$target_ranges),
          start(x = overlap_frame$transcribed_ranges)
        ),
        end = pmin(
          end(x = overlap_frame$target_ranges),
          end(x = overlap_frame$transcribed_ranges)
        )
      )
    )

  constrained_list$constrained_number <-
    length(x = constrained_list$constrained_ranges)

  constrained_list$constrained_width <-
    sum(GenomicRanges::width(x = constrained_list$constrained_ranges))

  rm(overlap_frame, target_list)
} else {
  # If target regions are not avaiable, all transcribed GRanges count.
  if (argument_list$verbose) {
    message("Not importing any target range annotation.")
  }
  constrained_list$constrained_ranges <-
    ensembl_list$transcribed_ranges
  constrained_list$constrained_number <-
    ensembl_list$transcribed_number
  constrained_list$constrained_width <-
    ensembl_list$transcribed_width
}
if (argument_list$verbose) {
  message("Number of transcribed target ranges: ",
          constrained_list$constrained_number)
  message(
    "Cumulative width of transcribed target ranges: ",
    constrained_list$constrained_width
  )
}

# Export to BED format.
export.bed(object = constrained_list$constrained_ranges,
           con = argument_list$output_path)

rm(constrained_list,
   ensembl_list,
   argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

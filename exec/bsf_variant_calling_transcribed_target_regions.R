#!/usr/bin/env Rscript
#
# BSF R script to constrain target enrichment regions to exon regions.
#
# All Ensembl "exon" features are imported from a GTF file and reduced into a
# (non-overlapping) set of transcribed regions of the genome. Optionally, target
# regions can be imported and overlapped with the transcribed regions above to
# export the minimal set of transcribed target regions in BED format.
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
      ),
      optparse::make_option(
        opt_str = c("--flanks"),
        default = 0L,
        dest = "flanks",
        help = "Exon and Target flanking regions [0L]",
        type = "integer"
      ),
      optparse::make_option(
        opt_str = c("--full-annotation"),
        action = "store_false",
        default = TRUE,
        dest = "basic",
        help = "Import full rather than basic Ensembl annotation",
        type = "logical"
      )
    )
  ))

suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
suppressPackageStartupMessages(expr = library(package = "bsfR"))
suppressPackageStartupMessages(expr = library(package = "Biostrings"))

summary_list <-
  bsfR::bsfvc_import_constrained_granges(
    exon_path = argument_list$exon_path,
    exon_flanks = argument_list$flanks,
    exon_basic = argument_list$basic,
    target_path = argument_list$target_path,
    target_flanks = argument_list$flanks,
    genome_version = argument_list$genome_version,
    verbose = argument_list$verbose
  )

# Export to BED format.
rtracklayer::export.bed(object = summary_list$constrained_granges,
                        con = argument_list$output_path)

rm(summary_list,
   argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessioninfo::session_info())

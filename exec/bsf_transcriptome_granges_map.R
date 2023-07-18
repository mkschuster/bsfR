#!/usr/bin/env Rscript
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


# BSF R Script to map between NCBI and UCSC sequence region names. A
# GenomicRanges::GRanges object can be adjusted to a compatible
# GenomeInfoDb::Seqinfo object.

# Option Parsing ----------------------------------------------------------


suppressPackageStartupMessages(expr = library(package = "optparse"))

argument_list <-
  optparse::parse_args(object = optparse::OptionParser(
    option_list = list(
      optparse::make_option(
        opt_str = "--verbose",
        action = "store_true",
        default = TRUE,
        help = "Print extra output [default]",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = "--quiet",
        action = "store_false",
        default = FALSE,
        dest = "verbose",
        help = "Print little output",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = "--granges-path",
        dest = "granges_path",
        help = "Source GenomicRanges::GRanges path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--seqinfo-path",
        dest = "seqinfo_path",
        help = "Target GenomeInfoDb::Seqinfo path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--ncbi-version",
        dest = "ncbi_version",
        help = "NCBI genome version",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--ucsc-version",
        dest = "ucsc_version",
        help = "UCSC genome version",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--ncbi-granges",
        action = "store_true",
        dest = "ncbi_to_ucsc",
        help = "NCBI to UCSC mapping",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = "--ucsc-granges",
        action = "store_false",
        dest = "ncbi_to_ucsc",
        help = "UCSC to NCBI mapping",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = "--output-directory",
        default = ".",
        dest = "output_directory",
        help = "Output directory path [.]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--output-path",
        dest = "output_path",
        help = "Output GRanges path []",
        type = "character"
      )
    )
  ))

if (is.null(x = argument_list$output_path)) {
  stop("Missing --output-path option")
}

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "readr"))
# Bioconductor
suppressPackageStartupMessages(expr = library(package = "BiocVersion"))
suppressPackageStartupMessages(expr = library(package = "GenomicRanges"))
suppressPackageStartupMessages(expr = library(package = "bsfR"))

source_granges <- readr::read_rds(file = argument_list$granges_path)

target_seqinfo <- readr::read_rds(file = argument_list$seqinfo_path)

genome_list <-
  bsfR::bsfg_get_genome_list(
    resource_directory = ".",
    ensembl_version = "109",
    ncbi_version = argument_list$ncbi_version,
    ucsc_version = argument_list$ucsc_version
  )

assembly_report_tibble <-
  bsfR::bsfg_get_assembly_report(genome_list = genome_list)

target_granges <-
  if (argument_list$ncbi_to_ucsc) {
    bsfR::bsfg_convert_seqlevels_ncbi_to_ucsc(
      ncbi_granges = source_granges,
      ucsc_seqinfo = target_seqinfo,
      assembly_report_tibble = assembly_report_tibble,
      verbose = argument_list$verbose
    )
  } else {
    bsfR::bsfg_convert_seqlevels_ucsc_to_ncbi(
      ucsc_granges = source_granges,
      ncbi_seqinfo = target_seqinfo,
      assembly_report_tibble = assembly_report_tibble,
      verbose = argument_list$verbose
    )
  }

readr::write_rds(x = target_granges,
                 file = argument_list$output_path,
                 compress = "gz")

rm(
  target_granges,
  assembly_report_tibble,
  genome_list,
  target_seqinfo,
  source_granges,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessioninfo::session_info())

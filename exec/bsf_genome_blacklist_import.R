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


# BSF R script to import black lists.
#
# - Import an (ENCODE) BED file as GenomicRanges::GRanges object.
#
# - Reduce the black list GenomicRanges::GRanges object.
# - Sort the black list GenomicRanges::GRanges object.
# - Serialise the black list GenomicRanges::GRanges object to disk.
#
# Output files:
# - <output_directory>/<blacklist_prefix>_blacklist_granges.rds

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
        opt_str = "--blacklist-path",
        dest = "blacklist_path",
        help = "Blacklist (BED) path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--seqinfo-path",
        dest = "seqinfo_path",
        help = "GenomeInfoDb::Seqinfo path (preferred over --genome-version)",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--genome-version",
        dest = "genome_version",
        help = "Genome version (--seqinfo-path is preferred)",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--blacklist-prefix",
        dest = "blacklist_prefix",
        help = "Blacklist file name prefix [basename(<blacklist_path>)]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--output-directory",
        default = ".",
        dest = "output_directory",
        help = "Output directory path [.]",
        type = "character"
      )
    )
  ))

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "readr"))
suppressPackageStartupMessages(expr = library(package = "stringr"))
# Bioconductor
suppressPackageStartupMessages(expr = library(package = "BiocVersion"))
suppressPackageStartupMessages(expr = library(package = "GenomeInfoDb"))
suppressPackageStartupMessages(expr = library(package = "rtracklayer"))

if (is.null(x = argument_list$blacklist_path)) {
  stop("Missing --blacklist-path option")
}

if (is.null(x = argument_list$blacklist_prefix)) {
  # Get the base name of the blacklist path and
  # split off a (BED, GTF, ...) extension.

  blacklist_name_character <-
    stringr::str_split_1(
      string = base::basename(path = argument_list$blacklist_path),
      pattern = stringr::fixed(pattern = ".")
    )

  blacklist_name_index <-
    base::max(1L, base::length(x = blacklist_name_character) - 1L)

  argument_list$blacklist_prefix <-
    paste(blacklist_name_character[1L:blacklist_name_index], collapse = ".")

  rm(blacklist_name_index, blacklist_name_character)
}

# Read or retrieve a GenomeInfoDb::Seqinfo object.

seqinfo_object <-
  if (is.null(x = argument_list$seqinfo_path)) {
    GenomeInfoDb::Seqinfo(genome = argument_list$genome_version)
  } else {
    readr::read_rds(file = argument_list$seqinfo_path)
  }

if (is.null(x = seqinfo_object)) {
  stop("Could not retrieve a GenomeInfoDb::Seqinfo object")
}

black_granges <-
  rtracklayer::import(con = argument_list$blacklist_path,
                      # format = "bed",
                      genome = seqinfo_object)

# Reduce and sort the imported blacklist GenomicRanges::GRanges.

black_granges <-
  GenomicRanges::reduce(
    x = black_granges,
    drop.empty.ranges = TRUE,
    with.revmap = TRUE,
    ignore.strand = TRUE
  )
black_granges <- GenomicRanges::sort(x = black_granges)

# Serialize the blacklist GenomicRanges::GRanges to disk.

readr::write_rds(
  x = black_granges,
  file = file.path(argument_list$output_directory, paste(
    paste(argument_list$blacklist_prefix, "blacklist", "granges", sep = "_"),
    "rds",
    sep = "."
  )),
  compress = "gz"
)

rm(black_granges,
   seqinfo_object,
   argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessioninfo::session_info())

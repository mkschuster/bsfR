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


# BSF R Script to merge two GenomicRanges::GRanges objects.
#
# Output file:
# - <output_directory>/<transcriptome_version>_granges.rds
# - <output_directory>/<transcriptome_version>.gtf

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
        opt_str = "--granges-path-1",
        dest = "granges_path_1",
        help = "GRanges path 1",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--granges-path-2",
        dest = "granges_path_2",
        help = "GRanges path 2",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--transcriptome-version",
        dest = "transcriptome_version",
        help = "Transcriptome version",
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

if (is.null(x = argument_list$granges_path_1)) {
  stop("Missing --granges-path-1 option")
}

if (is.null(x = argument_list$granges_path_2)) {
  stop("Missing --granges-path-2 option")
}

if (is.null(x = argument_list$transcriptome_version)) {
  stop("Missing --transcriptome-version option")
}

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "readr"))
# Bioconductor
suppressPackageStartupMessages(expr = library(package = "BiocVersion"))
suppressPackageStartupMessages(expr = library(package = "GenomicRanges"))
suppressPackageStartupMessages(expr = library(package = "rtracklayer"))

granges_object_1 <-
  readr::read_rds(file = argument_list$granges_path_1)
granges_object_2 <-
  readr::read_rds(file = argument_list$granges_path_2)

granges_object <- c(granges_object_1, granges_object_2)

readr::write_rds(
  x = granges_object,
  file = file.path(argument_list$output_directory,
                   paste(
                     paste(argument_list$transcriptome_version, "granges", sep = "_"),
                     "rds",
                     sep = "."
                   )),
  compress = "gz"
)

# Re-export the merged GenomicRanges::GRanges objects in GTF format.
# The rtracklayer::export() function checks for a single genome.
seqinfo_object <- GenomicRanges::seqinfo(x = granges_object)
GenomeInfoDb::genome(x = seqinfo_object) <-
  paste(unique(x = GenomeInfoDb::genome(x = seqinfo_object)), collapse = "_")
GenomicRanges::seqinfo(x = granges_object) <- seqinfo_object
rm(seqinfo_object)

rtracklayer::export(object = granges_object,
                    con = file.path(
                      argument_list$output_directory,
                      paste(argument_list$transcriptome_version,
                            "gtf",
                            sep = ".")
                    ))

rm(granges_object,
   granges_object_2,
   granges_object_1,
   argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessioninfo::session_info())

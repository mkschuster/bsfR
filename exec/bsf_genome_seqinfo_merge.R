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


# BSF R Script to merge two GenomeInfoDb::Seqinfo objects.
#
# Output files:
# - <output_directory>/<genome_version>_seqinfo.rds

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
        opt_str = "--seqinfo-path-1",
        dest = "seqinfo_path_1",
        help = "Seqinfo file path (genome_seqinfo.rds)",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--seqinfo-path-2",
        dest = "seqinfo_path_2",
        help = "Seqinfo file path (genome_seqinfo.rds)",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--genome-version",
        dest = "genome_version",
        help = "Genome assembly version",
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

if (is.null(x = argument_list$genome_version)) {
  stop("Missing --genome-version option")
}

if (is.null(x = argument_list$seqinfo_path_1)) {
  stop("Missing --seqinfo-path-1 option")
}

if (is.null(x = argument_list$seqinfo_path_2)) {
  stop("Missing --seqinfo-path-2 option")
}


# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "readr"))
# Bioconductor
suppressPackageStartupMessages(expr = library(package = "BiocVersion"))
suppressPackageStartupMessages(expr = library(package = "GenomeInfoDb"))

seqinfo_1 <- readr::read_rds(file = argument_list$seqinfo_path_1)
seqinfo_2 <- readr::read_rds(file = argument_list$seqinfo_path_2)

seqinfo_object <- GenomeInfoDb::merge(x = seqinfo_1, y = seqinfo_2)

readr::write_rds(
  x = seqinfo_object,
  file = file.path(argument_list$output_directory,
                   paste(
                     paste(argument_list$genome_version, "seqinfo", sep = "_"),
                     "rds",
                     sep = "."
                   )),
  compress = "gz"
)

rm(seqinfo_object, seqinfo_2, seqinfo_1, argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessioninfo::session_info())

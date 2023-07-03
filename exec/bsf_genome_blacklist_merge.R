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


# BSF R script to merge black list GenomicRanges::GRanges objects.
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
        opt_str = "--blacklist-path-1",
        dest = "blacklist_path_1",
        help = "Blacklist GRanges path 1",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--blacklist-path-2",
        dest = "blacklist_path_2",
        help = "Blacklist GRanges path 2",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--blacklist-prefix",
        dest = "blacklist_prefix",
        help = "Blacklist file name prefix",
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

if (is.null(x = argument_list$blacklist_path_1)) {
  stop("Missing --blacklist-path-1 option")
}

if (is.null(x = argument_list$blacklist_path_2)) {
  stop("Missing --blacklist-path-2 option")
}

if (is.null(x = argument_list$blacklist_prefix)) {
  stop("Missing --blacklist-prefix option")
}

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "readr"))
# Bioconductor
suppressPackageStartupMessages(expr = library(package = "BiocVersion"))
suppressPackageStartupMessages(expr = library(package = "GenomicRanges"))

granges_1 <- readr::read_rds(file = argument_list$blacklist_path_1)
granges_2 <- readr::read_rds(file = argument_list$blacklist_path_2)

black_granges <- c(granges_1, granges_2)

readr::write_rds(x = black_granges,
                 file = file.path(argument_list$output_directory,
                                  paste(
                                    paste(argument_list$blacklist_prefix,
                                          "blacklist",
                                          "granges",
                                          sep = "_"),
                                    "rds",
                                    sep = "."
                                  )))

rm(black_granges,
   granges_2,
   granges_1,
   argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessioninfo::session_info())

#!/usr/bin/env Rscript
#
# Copyright 2013 - 2022 Michael K. Schuster
#
# Biomedical Sequencing Facility (BSF), part of the genomics core facility
# of the Research Center for Molecular Medicine (CeMM) of the
# Austrian Academy of Sciences and the Medical University of Vienna (MUW).
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


# BSF R script to run a ChIPQC on top of a DiffBind analysis.

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
        opt_str = "--comparison",
        dest = "comparison",
        help = "Comparison name",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--factor",
        dest = "factor",
        help = "ChIP factor",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--threads",
        default = 1L,
        dest = "threads",
        help = "Number of parallel processing threads [1]",
        type = "integer"
      )
    )
  ))

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# Bioconductor
suppressPackageStartupMessages(expr = library(package = "BiocParallel"))
suppressPackageStartupMessages(expr = library(package = "BiocVersion"))
suppressPackageStartupMessages(expr = library(package = "ChIPQC"))
suppressPackageStartupMessages(expr = library(package = "DiffBind"))

# Set the number of parallel threads in the MulticoreParam instance.
BiocParallel::register(BPPARAM = BiocParallel::MulticoreParam(workers = argument_list$threads))


prefix <-
  paste("chipseq",
        "diff",
        "bind",
        argument_list$comparison,
        argument_list$factor,
        sep = "_")


# Initialise a DBA object -------------------------------------------------


diffbind_dba <- NULL
file_path <- file.path(prefix, paste0(prefix, '_DBA.RData'))
if (file.exists(file_path) &&
    file.info(file_path)$size > 0L) {
  message("Loading a DiffBind DBA object")
  diffbind_dba <-
    DiffBind::dba.load(dir = prefix, pre = paste0(prefix, "_"))
} else {
  stop("Require a pre-calculated DiffBind DBA object in file: ",
       file_path)
}
rm(file_path)

# Initialise a ChIPQCexperiment object ------------------------------------


chipqc_experiment <- ChIPQC::ChIPQC(experiment = diffbind_dba)

# By not setting the "chromosomes" option, only the first chromosome with peaks will be inspected.

ChIPQC::ChIPQCreport(
  object = chipqc_experiment,
  reportFolder = paste(
    "chipseq",
    "chipqc",
    argument_list$comparison,
    argument_list$factor,
    sep = "_"
  )
)

rm(chipqc_experiment, diffbind_dba, prefix, argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessioninfo::session_info())

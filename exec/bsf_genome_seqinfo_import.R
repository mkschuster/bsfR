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


# BSF R Script to import a GenomeInfoDb::Seqinfo object.
#
# If supported by the GenomeInfoDb package, the GenomeInfoDb::Seqinfo object can
# be imported automatically from the Ensembl or UCSC Genome Browser purely via
# the genome assembly version. Alternatively, a Samtools faidx file in
# conjunction with a comma-separated list of circular sequences can be
# specified.
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
        opt_str = "--genome-version",
        dest = "genome_version",
        help = "Genome assembly version",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--faidx-path",
        dest = "faidx_path",
        help = "Samtools faidx file path (genome.fa.fai)",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--circular-sequences",
        dest = "circular_sequences",
        help = "Comma-separated list of circular sequence names",
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

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "readr"))
suppressPackageStartupMessages(expr = library(package = "stringr"))
# Bioconductor
suppressPackageStartupMessages(expr = library(package = "BiocVersion"))
suppressPackageStartupMessages(expr = library(package = "GenomeInfoDb"))

seqinfo_object <- NULL
if (is.null(x = argument_list$faidx_path)) {
  # Try fetching the Seqinfo object via the assembly name from NCBI or UCSC.
  seqinfo_object <-
    GenomeInfoDb::Seqinfo(genome = argument_list$genome_version)
} else {
  # http://www.htslib.org/doc/faidx.html
  faidx_tibble <- readr::read_tsv(
    file = argument_list$faidx_path,
    col_names = c("NAME", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"),
    col_types = readr::cols(
      "NAME" = readr::col_character(),
      "LENGTH" = readr::col_integer(),
      "OFFSET" = readr::col_double(),
      "LINEBASES" = readr::col_integer(),
      "LINEWIDTH" = readr::col_integer()
    )
  )

  problem_tibble <- readr::problems(x = faidx_tibble)
  if (nrow(x = problem_tibble) > 0) {
    print(x = problem_tibble)
  }
  rm(problem_tibble)

  seqinfo_object <-
    GenomeInfoDb::Seqinfo(
      seqnames = faidx_tibble$NAME,
      seqlengths = faidx_tibble$LENGTH,
      isCircular = rep_len(x = FALSE, length.out = length(x = faidx_tibble$NAME)),
      genome = argument_list$genome_version
    )

  if (!is.null(x = argument_list$circular_sequences)) {
    # Split on "," characters.
    circular_character <- stringr::str_split_1(
      string = argument_list$circular_sequences,
      pattern = stringr::fixed(pattern = ",")
    )

    # Trim remaining white space around "," characters.
    circular_character <- stringr::str_trim(string = circular_character)

    # Since the %in% operator returns TRUE for circular sequences and FALSE for
    # all other sequences, the resulting logical vector can be assigned
    # directly.

    GenomeInfoDb::isCircular(x = seqinfo_object) <-
      GenomeInfoDb::seqnames(x = seqinfo_object) %in% circular_character

    rm(circular_character)
  }

  rm(faidx_tibble)
}

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

rm(seqinfo_object, argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessioninfo::session_info())

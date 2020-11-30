#!/usr/bin/env Rscript
#
# BSF R script to read Exonerate Verbose Useful Labelled Gapped Alignment Report
# (VULGAR) lines, convert them into GenomicRanges::GRanges objects and export
# them via rtracklayer::export() in BED, GFF or GTF format.
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
        opt_str = c("--source-path"),
        dest = "source_path",
        help = "Exonerate alignment file path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--target-path"),
        dest = "target_path",
        help = "BED file path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--target-format"),
        default = "bed",
        dest = "target_format",
        help = "Target file format for rtracklayer::export() [bed]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--chunk-size"),
        default = 1000L,
        dest = "chunk_size",
        help = "Number of lines to process per chunk [1000]",
        type = "integer"
      )
    )
  ))

suppressPackageStartupMessages(expr = library(package = "tidyverse"))
suppressPackageStartupMessages(expr = library(package = "rtracklayer"))

# Use a tibble to capture the following fields of a VULGAR line.
vulgar_tibble <- NULL

variable_names <-   c(
  "query_name",
  "query_start",
  "query_end",
  "query_strand",
  "target_name",
  "target_start",
  "target_end",
  "target_strand",
  "score",
  "cigar"
)

exonerate_connection <-
  file(description = argument_list$source_path, open = "rt")
while (TRUE) {
  exonerate_lines <-
    base::readLines(con = exonerate_connection, n = argument_list$chunk_size)
  if (length(x = exonerate_lines) == 0L) {
    rm(exonerate_lines)
    break()
  }

  # Subset the Exonerate alignment report lines to the relevant VULGAR lines,
  # which look like:
  #
  # vulgar: nCoV-2019_1_LEFT 0 24 + NC_045512.2 30 54 + 120 M 24 24
  # vulgar: nCoV-2019_1_RIGHT 0 25 + NC_045512.2 410 385 - 125 M 25 25
  vulgar_lines <-
    exonerate_lines[base::startsWith(x = exonerate_lines, prefix = "vulgar:")]
  rm(exonerate_lines)

  # Split the character vector on space characters into a character matrix of
  # eleven columns, so that the CIGAR string remains unsplit in the last column.
  # Exclude the first column (vulgar:) and assign column names to the remaining
  # character matrix. Coerce the character matrix into a tibble and bind it to
  # the summary VULGAR tibble by rows.
  character_matrix <- stringr::str_split_fixed(
    string = vulgar_lines,
    pattern = stringr::fixed(pattern = " "),
    n = 11L
  )[,-1L]
  colnames(x = character_matrix) <- variable_names

  vulgar_tibble <-
    dplyr::bind_rows(vulgar_tibble,
                     tibble::as_tibble(x = character_matrix))

  rm(character_matrix, vulgar_lines)
}
close(con = exonerate_connection)
rm(exonerate_connection)

# Coerce start and end coordinates into integer vectors and recode unknown
# strand information from Exonerate, denoted by "." to GenomicRanges::GRanges,
# denoted by "*".
vulgar_tibble <-
  dplyr::mutate(
    .data = vulgar_tibble,
    query_start = as.integer(x = .data$query_start),
    query_end = as.integer(x = .data$query_end),
    query_strand = dplyr::recode(.x = .data$query_strand, "." = "*"),
    target_start = as.integer(x = .data$target_start),
    target_end = as.integer(x = .data$target_end),
    target_strand = dplyr::recode(.x = .data$target_strand, "." = "*")
  )

# Start and end coordinates require assigning, depending on the strand. While
# Exonerate reports start > end coordinates on the reverse strand,
# GenomicRanges::GRanges always expects start < end. In addition, Exonerate uses
# in-between coordinates like BED, so that the start coordinate needs
# incrementing by one.

exonerate_granges <-
  GenomicRanges::GRanges(
    seqnames = vulgar_tibble$target_name,
    ranges = IRanges::IRanges(
      start = dplyr::if_else(
        condition = vulgar_tibble$target_strand == "+",
        true = vulgar_tibble$target_start,
        false = vulgar_tibble$target_end,
        missing = vulgar_tibble$target_start
      ) + 1L,
      end = dplyr::if_else(
        condition = vulgar_tibble$target_strand == "+",
        true = vulgar_tibble$target_end,
        false = vulgar_tibble$target_start,
        missing = vulgar_tibble$target_end
      ),
      names = vulgar_tibble$query_name
    ),
    strand = vulgar_tibble$target_strand
  )

rtracklayer::export(
  object = GenomicRanges::sort(x = exonerate_granges, ignore.strand = TRUE),
  con = argument_list$target_path,
  format = argument_list$target_format
)

rm(exonerate_granges,
   variable_names,
   vulgar_tibble,
   argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

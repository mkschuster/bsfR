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


# BSF R script to parse and summarise Tophat2 alignment summary files.
# The resulting data frame (rnaseq_tophat_alignment_summary.tsv),
# as well as plots of alignment rates per sample in PDF (rnaseq_tophat_alignment_summary.pdf)
# and PNG format (rnaseq_tophat_alignment_summary.png) are written into the current working
# directory.

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
        opt_str = "--genome-directory",
        default = ".",
        dest = "genome_directory",
        help = "Genome directory path [.]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--output-directory",
        default = ".",
        dest = "output_directory",
        help = "Output directory path [.]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--directory-pattern",
        default = "^rnaseq_tophat_(.*)$",
        dest = "directory_pattern",
        help = "Tophat2 directory name pattern [^rnaseq_tophat_(.*)$]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--plot-factor",
        default = 0.5,
        dest = "plot_factor",
        help = "Plot width increase per 24 samples [0.5]",
        type = "double"
      ),
      optparse::make_option(
        opt_str = "--plot-width",
        default = 7.0,
        dest = "plot_width",
        help = "Plot width in inches [7.0]",
        type = "double"
      ),
      optparse::make_option(
        opt_str = "--plot-height",
        default = 7.0,
        dest = "plot_height",
        help = "Plot height in inches [7.0]",
        type = "double"
      )
    )
  ))

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "ggplot2"))
suppressPackageStartupMessages(expr = library(package = "purrr"))
suppressPackageStartupMessages(expr = library(package = "readr"))
# BSF
suppressPackageStartupMessages(expr = library(package = "bsfR"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

#' @title Parse a Tophat2 alignment report.
#'
#' @description Parses a Tophat2 alignment report.
#'
#' @noRd
#' @param directory_path A \code{character} scalar of the Tophat2 alignment report
#'   directory path.
#'
#' @return: A named \code{list}.
#' \describe{
#' \item{input}{A \code{integer} scalar with the number of input reads.}
#' \item{mapped}{A \code{integer} scalar with the number of mapped reads.}
#' \item{multiple}{A \code{integer} scalar with the number of multiply mapped reads.}
#' \item{above}{A \code{integer} scalar with the number of multiply mapped reads above the threshold.}
#' \item{threshold}{A \code{integer} scalar with the threshold of multi-mapped reads.}
#' \item{total_rate}{A \code{double} scalar with the total mapping rate as percentage.}
#' \item{error}{A \code{character} scalar with an error message.}
#' }
parse_tophat2_report <- function(directory_path) {
  result_list <- list()

  # This is the layout of a Tophat2 alignment summary report file.
  #
  #   [1] "Reads:"
  #   [2] "          Input     :  21791622"
  #   [3] "           Mapped   :  21518402 (98.7% of input)"
  #   [4] "            of these:   2010462 ( 9.3%) have multiple alignments (8356 have >20)"
  #   [5] "98.7% overall read mapping rate."

  report_lines <- readr::read_lines(
    file = file.path(
      directory_path,
      "align_summary.txt"
    ),
    n_max = 100L
  )

  # Find the first line of the alignment report.
  line_number <-
    which(grepl(pattern = "Reads:", x = report_lines[line_number]))

  # Skip this sample report in case it does not have exactly one alignment report.
  if (length(x = line_number) != 1L) {
    message("    parsing failed!")
    rm(report_lines, line_number)
    result_list$error <- "parsing error"
    return(result_list)
  }

  line_number <- line_number + 1L

  # Parse the second line of "input" reads as integer.
  result_list$input <- as.integer(
    x = sub(
      pattern = "[[:space:]]+Input[[:space:]]+:[[:space:]]+([[:digit:]]+)",
      replacement = "\\1",
      x = report_lines[line_number]
    )
  )
  line_number <- line_number + 1L

  # Parse the third line of "mapped" reads as integer.
  result_list$mapped <- as.integer(
    x = sub(
      pattern = "[[:space:]]+Mapped[[:space:]]+:[[:space:]]+([[:digit:]]+) .*",
      replacement = "\\1",
      x = report_lines[line_number]
    )
  )
  line_number <- line_number + 1L

  # Get the number of "multiply" aligned reads from the fourth line.
  result_list$multiple <- as.integer(
    x = sub(
      pattern = ".+:[[:space:]]+([[:digit:]]+) .*",
      replacement = "\\1",
      x = report_lines[line_number]
    )
  )

  # Get the number of multiply aligned reads "above" the multiple alignment
  # threshold as integer.
  result_list$above <- as.integer(
    x = sub(
      pattern = ".+alignments \\(([[:digit:]]+) have.+",
      replacement = "\\1",
      x = report_lines[line_number]
    )
  )

  # Get the multiple alignment "threshold" as integer.
  result_list$threshold <- as.integer(x = sub(
    pattern = ".+ >([[:digit:]]+)\\)",
    replacement = "\\1",
    x = report_lines[line_number]
  ))
  line_number <- line_number + 1L

  result_list$total_rate <- as.double(
    x = sub(
      pattern = "([[:digit:].]+)% overall read mapping rate",
      replacement = "\\1",
      x = report_lines[line_number]
    )
  )

  rm(line_number, report_lines)

  return(result_list)
}

message("Processing Tophat2 alignment reports ...")

sample_tibble <- bsfR::bsfu_summarise_report_directories(
  directory_path = argument_list$genome_directory,
  directory_pattern = argument_list$directory_pattern,
  report_function = parse_tophat2_report
)

# Write the alignment summary tibble to disk and create plots.

readr::write_tsv(
  x = sample_tibble,
  file = file.path(
    argument_list$output_directory,
    "rnaseq_tophat_alignment_summary.tsv"
  )
)

# Alignment Summary Plot --------------------------------------------------


plot_paths <-
  file.path(argument_list$output_directory,
            paste(
              paste("rnaseq", "tophat", "alignment", "summary", sep = "_"),
              graphics_formats,
              sep = "."
            ))

ggplot_object <- ggplot2::ggplot(data = sample_tibble)

ggplot_object <-
  ggplot_object + ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = .data$mapped,
      y = .data$mapped / .data$input,
      colour = .data$sample_name
    )
  )

ggplot_object <-
  ggplot_object + ggplot2::labs(
    x = "Reads Number",
    y = "Reads Fraction",
    colour = "Sample",
    title = "Tophat2 Alignment Summary"
  )

# Reduce the label font size and the legend key size and allow a maximum of 24
# guide legend rows.
ggplot_object <-
  ggplot_object + ggplot2::guides(colour = ggplot2::guide_legend(
    keywidth = ggplot2::rel(x = 0.8),
    keyheight = ggplot2::rel(x = 0.8),
    nrow = 24L
  ))

ggplot_object <-
  ggplot_object + ggplot2::theme(legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7)))

# Scale the plot width with the number of samples, by adding a quarter of
# the original width for each 24 samples.
plot_width <-
  argument_list$plot_width + (ceiling(x = base::nrow(x = sample_tibble) / 24L) - 1L) * argument_list$plot_width * argument_list$plot_factor

for (plot_path in plot_paths) {
  ggplot2::ggsave(
    filename = plot_path,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height,
    limitsize = FALSE
  )
}

rm(
  plot_path,
  plot_width,
  ggplot_object,
  plot_paths,
  sample_tibble,
  parse_tophat2_report,
  graphics_formats,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessioninfo::session_info())

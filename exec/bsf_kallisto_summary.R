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


# BSF R script to summarise Kallisto alignment reports.
#
# This script reads sample-specific run_info.json files and collates the
# information in a summary tibble. The fraction of pseudo-aligned reads and
# unique pseudo-aligned reads is plotted separately, while a summary plot
# correlates the fractions and numbers of pseudo-aligned reads.

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
        default = "^kallisto_sample_(.*)$",
        dest = "directory_pattern",
        help = "Kallisto directory path pattern [^kallisto_sample_(.*)$]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--prefix",
        default = "kallisto_summary",
        dest = "prefix",
        help = "File name prefix",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--plot-width",
        default = 7.0,
        dest = "plot_width",
        help = "Plot width in inches [7.0]",
        type = "numeric"
      ),
      optparse::make_option(
        opt_str = "--plot-height",
        default = 7.0,
        dest = "plot_height",
        help = "Plot height in inches [7.0]",
        type = "numeric"
      )
    )
  ))

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "dplyr"))
suppressPackageStartupMessages(expr = library(package = "ggplot2"))
suppressPackageStartupMessages(expr = library(package = "purrr"))
suppressPackageStartupMessages(expr = library(package = "readr"))
suppressPackageStartupMessages(expr = library(package = "tidyr"))
# CRAN
suppressPackageStartupMessages(expr = library(package = "rjson"))
# BSF
suppressPackageStartupMessages(expr = library(package = "bsfR"))

# Save plots in the following formats.
graphics_formats <- c("pdf" = "pdf", "png" = "png")

#' @title Parse a Kallisto pseudo-alignment report.
#'
#' @description Parses a Kallisto pseudo-alignment report.
#'
#' @noRd
#' @param directory_path A \code{character} scalar of the Kallisto alignment report
#'   directory path.
#' @references argument_list$directory_pattern A \code{character} scalar of the
#'   directory pattern to match sample names from Kallisto output directories.
#'
#' @return: A named \code{list}.
parse_kallisto_report <- function(directory_path) {
  return(rjson::fromJSON(file = base::file.path(directory_path, "run_info.json")))
}

# Summarise per-sample run reports ----------------------------------------


message("Processing Kallisto pseudo-alignment reports ...")

sample_tibble <- bsfR::bsfu_summarise_report_directories(
  directory_path = argument_list$genome_directory,
  directory_pattern = argument_list$directory_pattern,
  report_function = parse_kallisto_report
)

readr::write_tsv(
  x = sample_tibble,
  file = file.path(argument_list$output_directory,
                   paste(
                     paste(argument_list$prefix,
                           "sample",
                           sep = "_"),
                     "tsv",
                     sep = "."
                   )),
  col_names = TRUE
)

# Calculate the number of unmapped reads for plotting.
sample_tibble <-
  dplyr::mutate(
    .data = sample_tibble,
    "n_unmapped" = .data$n_processed - .data$n_unique - .data$n_pseudoaligned
  )

# Scatter plot of read number versus alignment rate per sample ------------


message("Creating a scatter plot of read number versus alignment rate per sample")

ggplot_object <-
  ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = dplyr::select(
        .data = sample_tibble,
        "sample" = .data$sample_name,
        "unmapped" = .data$n_unmapped,
        "multi" = .data$n_pseudoaligned,
        "unique" = .data$n_unique,
        "processed" = .data$n_processed
      ),
      cols = c(.data$unmapped, .data$multi, .data$unique),
      names_to = "status",
      values_to = "number"
    )
  )

ggplot_object <-
  ggplot_object + ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = .data$number,
      y = .data$number / .data$processed,
      colour = .data$sample,
      shape = .data$status
    ),
    alpha = I(1 / 3)
  )

ggplot_object <-
  ggplot_object + ggplot2::labs(
    x = "Number",
    y = "Fraction",
    colour = "Sample",
    shape = "Mapping Status",
    title = "Kallisto Summary"
  )

ggplot_object <-
  ggplot_object + ggplot2::theme(legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7)))

# Scale the plot width with the number of samples, by adding a quarter of
# the original width for each 24 samples.
plot_width <-
  argument_list$plot_width + (ceiling(x = base::nrow(x = sample_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
  ggplot2::ggsave(
    filename = file.path(
      argument_list$output_directory,
      paste(
        paste(argument_list$prefix,
              "sample",
              sep = "_"),
        graphics_format,
        sep = "."
      )
    ),
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height,
    limitsize = FALSE
  )
}
rm(graphics_format, plot_width, ggplot_object)

# Column plot of read numbers per sample --------------------------------


message("Creating a column plot of read numbers per sample")

ggplot_object <-
  ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = dplyr::select(
        .data = sample_tibble,
        "sample" = .data$sample_name,
        "unmapped" = .data$n_unmapped,
        "multi" = .data$n_pseudoaligned,
        "unique" = .data$n_unique
      ),
      cols = c(.data$unmapped, .data$multi, .data$unique),
      names_to = "status",
      values_to = "number"
    )
  )

ggplot_object <-
  ggplot_object + ggplot2::geom_col(
    mapping = ggplot2::aes(
      x = .data$sample,
      y = .data$number,
      fill = .data$status
    ),
    alpha = I(1 / 3)
  )

ggplot_object <-
  ggplot_object + ggplot2::labs(x = "Sample",
                                y = "Reads Number",
                                fill = "Mapping Status",
                                title = "Kallisto Pseudo-Aligned Numbers per Sample")

# Reduce the label font size and the legend key size and allow a maximum of 24
# guide legend rows.
ggplot_object <-
  ggplot_object + ggplot2::guides(colour = ggplot2::guide_legend(
    keywidth = ggplot2::rel(x = 0.8),
    keyheight = ggplot2::rel(x = 0.8),
    nrow = 24L
  ))

ggplot_object <-
  ggplot_object + ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      size = ggplot2::rel(x = 0.7),
      hjust = 0.0,
      vjust = 0.5,
      angle = 90.0
    ),
    legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7))
  )

# Scale the plot width with the number of samples, by adding a quarter of
# the original width for each 24 samples.
plot_width <-
  argument_list$plot_width + (ceiling(x = base::nrow(x = sample_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
  ggplot2::ggsave(
    filename = file.path(
      argument_list$output_directory,
      paste(
        paste(argument_list$prefix,
              "pseudaligned",
              "number",
              sep = "_"),
        graphics_format,
        sep = "."
      )
    ),
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height,
    limitsize = FALSE
  )
}
rm(graphics_format, plot_width, ggplot_object)

# Column plot of read fractions per sample ------------------------------


message("Creating a column plot of read fractions per sample")

ggplot_object <-
  ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = dplyr::transmute(
        .data = sample_tibble,
        "sample" = .data$sample_name,
        "unmapped" = .data$n_unmapped / .data$n_processed,
        "multi" = .data$n_pseudoaligned / .data$n_processed,
        "unique" = .data$n_unique / .data$n_processed
      ),
      cols = c(.data$unmapped, .data$multi, .data$unique),
      names_to = "status",
      values_to = "fraction"
    )
  )

ggplot_object <-
  ggplot_object + ggplot2::geom_col(
    mapping = ggplot2::aes(
      x = .data$sample,
      y = .data$fraction,
      fill = .data$status
    ),
    alpha = I(1 / 3)
  )

ggplot_object <-
  ggplot_object + ggplot2::labs(x = "Sample",
                                y = "Reads Fraction",
                                fill = "Mapping Status",
                                title = "Kallisto Pseudo-Aligned Fractions per Sample")

# Reduce the label font size and the legend key size and allow a maximum of 24
# guide legend rows.
ggplot_object <-
  ggplot_object + ggplot2::guides(colour = ggplot2::guide_legend(
    keywidth = ggplot2::rel(x = 0.8),
    keyheight = ggplot2::rel(x = 0.8),
    nrow = 24L
  ))

ggplot_object <-
  ggplot_object + ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      size = ggplot2::rel(x = 0.7),
      hjust = 0.0,
      vjust = 0.5,
      angle = 90.0
    ),
    legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7))
  )

# Scale the plot width with the number of samples, by adding a quarter of
# the original width for each 24 samples.
plot_width <-
  argument_list$plot_width + (ceiling(x = base::nrow(x = sample_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
  ggplot2::ggsave(
    filename = file.path(
      argument_list$output_directory,
      paste(
        paste(argument_list$prefix,
              "pseudaligned",
              "fraction",
              sep = "_"),
        graphics_format,
        sep = "."
      )
    ),
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height,
    limitsize = FALSE
  )
}
rm(graphics_format, plot_width, ggplot_object)

rm(sample_tibble,
   parse_kallisto_report,
   graphics_formats,
   argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessioninfo::session_info())

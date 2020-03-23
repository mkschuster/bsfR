#!/usr/bin/env Rscript
#
# BSF R script to summarise Kallisto alignment reports.
#
# This script reads sample-specifc run_info.json files and collates the
# information in a summary tibble. The fraction of pseudo-aligned reads and
# unique pseudo-aligned reads is plotted separately, while a summary plot
# correlates the fractions and numbers of pseudo-aligned reads.
#
#
# Copyright 2013 - 2019 Michael K. Schuster
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

argument_list <- parse_args(object = OptionParser(
  option_list = list(
    make_option(
      opt_str = c("--verbose", "-v"),
      action = "store_true",
      default = TRUE,
      help = "Print extra output [default]",
      type = "logical"
    ),
    make_option(
      opt_str = c("--quiet", "-q"),
      action = "store_false",
      default = FALSE,
      dest = "verbose",
      help = "Print little output",
      type = "logical"
    ),
    make_option(
      opt_str = c("--pattern-file"),
      default = "^run_info\\.json$",
      dest = "pattern_file",
      help = "Kallisto report file name pattern",
      type = "character"
    ),
    make_option(
      opt_str = c("--pattern-sample"),
      default = "^kallisto_sample_(.*)$",
      dest = "pattern_sample",
      help = "Kallisto report sample name pattern",
      type = "character"
    ),
    make_option(
      opt_str = c("--prefix"),
      default = "kallisto_summary",
      dest = "prefix",
      help = "File name prefix",
      type = "character"
    ),
    make_option(
      opt_str = c("--plot-width"),
      default = 7.0,
      dest = "plot_width",
      help = "Plot width in inches [7.0]",
      type = "numeric"
    ),
    make_option(
      opt_str = c("--plot-height"),
      default = 7.0,
      dest = "plot_height",
      help = "Plot height in inches [7.0]",
      type = "numeric"
    )
  )
))

suppressPackageStartupMessages(expr = library(package = "tidyverse"))
suppressPackageStartupMessages(expr = library(package = "rjson"))

# Save plots in the following formats.
graphics_formats <- c("pdf" = "pdf", "png" = "png")

# Summarise per-sample run reports ----------------------------------------


summary_frame <- tibble::tibble(
  "file_name" =
    base::list.files(pattern = argument_list$pattern_file, recursive = TRUE),
  "sample_name" = gsub(
    pattern = argument_list$pattern_sample,
    replacement = "\\1",
    x = base::basename(path = base::dirname(path = .data$file_name))
  )
)
message("Processing Kallisto reports for number of samples: ",
        nrow(x = summary_frame))

kallisto_frame <- tibble::tibble()

for (i in seq_len(length.out = nrow(x = summary_frame))) {
  message("  ", summary_frame[i, "sample_name", drop = TRUE])

  json_list <-
    rjson::fromJSON(file = summary_frame[i, "file_name", drop = TRUE])

  kallisto_frame <-
    dplyr::bind_rows(.data = kallisto_frame, json_list)
  rm(json_list)
}
rm(i)

summary_frame <- dplyr::bind_cols(summary_frame, kallisto_frame)
rm(kallisto_frame)

readr::write_tsv(
  x = summary_frame,
  path = paste(paste(argument_list$prefix,
                     "sample",
                     sep = "_"),
               "tsv",
               sep = "."),
  col_names = TRUE
)

# Scatter plot of read number versus alignment rate per sample ------------


ggplot_object <- ggplot2::ggplot(data = summary_frame)
ggplot_object <-
  ggplot_object + ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = .data$n_pseudoaligned,
      y = .data$n_pseudoaligned / .data$n_processed,
      colour = .data$sample_name
    ),
    alpha = I(1 / 3)
  )
ggplot_object <-
  ggplot_object + ggplot2::labs(
    x = "Number",
    y = "Fraction",
    colour = "Sample",
    title = "Kallisto Summary"
  )
ggplot_object <-
  ggplot_object + ggplot2::theme(legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7)))

# Scale the plot width with the number of samples, by adding a quarter of
# the original width for each 24 samples.
plot_width <-
  argument_list$plot_width + (ceiling(x = nrow(x = summary_frame) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
  ggplot2::ggsave(
    filename = paste(
      paste(argument_list$prefix,
            "sample",
            sep = "_"),
      graphics_format,
      sep = "."
    ),
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height,
    limitsize = FALSE
  )
}
rm(graphics_format, plot_width, ggplot_object)

# Fraction of pseudoaligned reads per sample ------------------------------


ggplot_object <- ggplot2::ggplot(data = summary_frame)
ggplot_object <-
  ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(
    x = .data$sample_name,
    y = .data$n_pseudoaligned / .data$n_processed
  ))
ggplot_object <-
  ggplot_object + ggplot2::labs(x = "Sample",
                                y = "Fraction",
                                title = "Fraction of Pseudoaligned Reads")
ggplot_object <-
  ggplot_object + ggplot2::theme(axis.text.x = ggplot2::element_text(
    angle = 90,
    hjust = 0,
    size = ggplot2::rel(x = 0.7)
  ))

# Scale the plot width with the number of samples, by adding a quarter of
# the original width for each 24 samples.
plot_width <-
  argument_list$plot_width + (ceiling(x = nrow(x = summary_frame) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
  ggplot2::ggsave(
    filename = paste(
      paste(argument_list$prefix,
            "pseudoaligned",
            sep = "_"),
      graphics_format,
      sep = "."
    ),
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height,
    limitsize = FALSE
  )
}
rm(graphics_format, plot_width, ggplot_object)

# Fraction of unique reads per sample -------------------------------------


ggplot_object <- ggplot2::ggplot(data = summary_frame)
ggplot_object <-
  ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(
    x = .data$sample_name,
    y = .data$n_unique / .data$n_processed
  ))
ggplot_object <-
  ggplot_object + ggplot2::labs(x = "Sample",
                                y = "Fraction",
                                title = "Fraction of Uniquely Pseudoaligned Reads")
ggplot_object <-
  ggplot_object + ggplot2::theme(axis.text.x = ggplot2::element_text(
    angle = 90,
    hjust = 0,
    size = ggplot2::rel(x = 0.7)
  ))

# Scale the plot width with the number of samples, by adding a quarter of
# the original width for each 24 samples.
plot_width <-
  argument_list$plot_width + (ceiling(x = nrow(x = summary_frame) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
  ggplot2::ggsave(
    filename = paste(
      paste(argument_list$prefix,
            "uniquely_pseudoaligned",
            sep = "_"),
      graphics_format,
      sep = "."
    ),
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height,
    limitsize = FALSE
  )
}
rm(graphics_format, plot_width, ggplot_object)

rm(summary_frame,
   graphics_formats,
   argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

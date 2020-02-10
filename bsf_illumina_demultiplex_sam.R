#!/usr/bin/env Rscript
#
# BSF R script to aggregate and plot Picard tools IlluminaSamDemux and
# Illumina2bam tools BamIndexDecoder metrics files.
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

argument_list <-
  parse_args(object = OptionParser(
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
        opt_str = c("--directory-path"),
        dest = "directory_path",
        help = "Directory path of EXPERIMENT_LANE_metrics.tsv files",
        type = "character"
      ),
      make_option(
        opt_str = c("--file-pattern"),
        default = "_metrics.tsv$",
        dest = "file_pattern",
        help = "File pattern to capture EXPERIMENT_LANE_metrics.tsv files",
        type = "character"
      ),
      make_option(
        opt_str = c("--file-path"),
        dest = "file_path",
        help = "File path of a EXPERIMENT_LANE_metrics.tsv file",
        type = "character"
      ),
      make_option(
        opt_str = c("--plot-factor"),
        default = 0.5,
        dest = "plot_factor",
        help = "Plot width increase per 24 samples [0.5]",
        type = "numeric"
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

if (is.null(x = argument_list$directory_path) &
    is.null(x = argument_list$file_path)) {
  stop("Missing --directory-path or --file-path option")
}

message("Loading packages")
suppressPackageStartupMessages(expr = library(package = "tidyverse"))
suppressPackageStartupMessages(expr = library(package = "PicardReports"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

message("Reading metrics files obeying pattern: ",
        argument_list$file_pattern)
metrics_files <- NULL
if (is.null(x = argument_list$directory_path)) {
  argument_list$directory_path <-
    base::dirname(path = argument_list$file_path)
  metrics_files <-
    base::basename(path = argument_list$file_path)
} else {
  metrics_files <-
    base::dir(
      path = argument_list$directory_path,
      pattern = argument_list$file_pattern,
      recursive = TRUE
    )
}
message("Number of metrics files: ", length(x = metrics_files))

for (i in seq_along(along.with = metrics_files)) {
  message("Processing Picard report: ", metrics_files[i])
  picard_report <-
    PicardMetricsFromFilePath(file_path = file.path(argument_list$directory_path, metrics_files[i]))

  # PicardBarcodeMetrics reports

  if (is(object = picard_report, class2 = "PicardBarcodeMetrics")) {
    # Increase the plot width per 24 samples.
    plot_width <-
      argument_list$plot_width + (ceiling(x = nrow(x = metrics(object = picard_report)) / 24L) - 1L) * argument_list$plot_width * argument_list$plot_factor

    # Plot cluster *fractions* by sample.

    ggplot_object <- plot_fractions(object = picard_report)
    for (graphics_format in graphics_formats) {
      ggplot2::ggsave(
        filename = file.path(
          argument_list$directory_path,
          base::dirname(path = metrics_files[i]),
          # Sub-directory of the metrics report.
          paste(
            paste(
              PicardReports::name(object = picard_report),
              "metrics",
              "fraction",
              sep = "_"
            ),
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
    rm(graphics_format, ggplot_object)

    # Plot cluster *numbers* by sample.

    ggplot_object <- plot_numbers(object = picard_report)
    for (graphics_format in graphics_formats) {
      ggplot2::ggsave(
        filename = file.path(
          argument_list$directory_path,
          base::dirname(path = metrics_files[i]),
          # Sub-directory of the metrics report.
          paste(
            paste(
              PicardReports::name(object = picard_report),
              "metrics",
              "number",
              sep = "_"
            ),
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
    rm(graphics_format, ggplot_object)

    rm(plot_width)
  }

  if (is(object = picard_report, class2 = "PicardAlignmentSummaryMetrics")) {

    file_prefix <- base::sub(
      pattern = "^(.*)_metrics.tsv$",
      replacement = "\\1",
      x = base::basename(path = metrics_files[i])
    )

    plot_object <-
      plot_cluster_numbers(
        object = picard_report,
        name = file_prefix
      )

    for (graphics_format in graphics_formats) {
      ggplot2::ggsave(
        filename = file.path(
          argument_list$directory_path,
          base::dirname(path = metrics_files[i]),
          # Sub-directory of the metrics report.
          paste(
            paste(
              file_prefix,
              "metrics",
              "number",
              sep = "_"
            ),
            graphics_format,
            sep = "."
          )
        ),
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height,
        limitsize = FALSE
      )
    }
    rm(graphics_format, ggplot_object, file_prefix)
  }
  rm(picard_report)
}

rm(i,
   metrics_files,
   graphics_formats,
   argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

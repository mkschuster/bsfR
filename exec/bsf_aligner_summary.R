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


# BSF R script to summarise Picard Alignment Summary Metrics for each sample
# and plot at the read group or sample level.

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
        opt_str = "--prefix",
        dest = "prefix",
        help = "File name prefix",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--directory-path"),
        default = ".",
        dest = "directory_path",
        help = "Directory path [.]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--output-directory"),
        default = ".",
        dest = "output_directory",
        help = "Output directory path [.]",
        type = "character"
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
      ),
      optparse::make_option(
        opt_str = c("--plot-limit-png"),
        default = 150.0,
        dest = "plot_limit_png",
        help = "Maximum size of the PNG device in inches [150.0]",
        type = "double"
      )
    )
  ))

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "ggplot2"))
# BSF
suppressPackageStartupMessages(expr = library(package = "PicardReports"))


# Save plots in the following formats.
graphics_formats <- c("pdf" = "pdf", "png" = "png")

# Assign a file prefix.
prefix_summary <- paste(argument_list$prefix, "summary", sep = "_")
# Assign a Picard Alignment Summary Metrics (pasm) prefix.
prefix_pasm <- "pasm"
# Assign a Picard Duplication Summary Metrics (pdsm) prefix.
prefix_pdsm <- "pdsm"

# Picard Alignment Summary Metrics ----------------------------------------


# This block is essentially identical to the
# PicardReports/exec/picard_reports_alignment_summary.R script, but uses
# different file search patterns, as well as plot and TSV file paths.

picard_frame <-
  PicardReports::PicardFrameFromDirectoryPath(
    directory_path = argument_list$directory_path,
    file_pattern = paste0(
      "^",
      argument_list$prefix,
      "_sample_(.*)_alignment_summary_metrics.tsv$"
    )
  )

type_list <- list("READ_GROUP" = "read_group", "SAMPLE" = "sample")

for (type in names(x = type_list)) {
  # Write the metrics tibble in TSV format --------------------------------


  readr::write_tsv(
    x = PicardReports::metrics_tibble(object = picard_frame, type = type),
    file = file.path(
      argument_list$output_directory,
      paste(
        paste(prefix_summary,
              prefix_pasm,
              "metrics",
              type_list[[type]],
              sep = "_"),
        "tsv",
        sep = "."
      )
    )
  )

  # Plot the absolute number versus the fraction --------------------------


  ggplot_object <-
    PicardReports::plot_number_versus_fraction(object = picard_frame, type = type)

  plot_width <-
    argument_list$plot_width + (ceiling(x = nlevels(x = ggplot_object$data[, type, drop = TRUE]) / 24L) - 1L) * argument_list$plot_width * 0.75

  for (graphics_format in graphics_formats) {
    if (graphics_format == "png" &&
        plot_width > argument_list$plot_limit_png) {
      message("PNG plot exceeding maximum size: ",
              plot_width,
              " > ",
              argument_list$plot_limit_png)
      next
    }

    ggplot2::ggsave(
      filename = file.path(
        argument_list$output_directory,
        paste(
          paste(prefix_summary,
                prefix_pasm,
                "alignment",
                type_list[[type]],
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

  # Plot the absolute number of aligned pass-filter reads -----------------


  ggplot_object <-
    PicardReports::plot_aligned_pass_filter_numbers(object = picard_frame, type = type)

  plot_width <-
    argument_list$plot_width + (ceiling(x = nlevels(x = ggplot_object$data[, type, drop = TRUE]) / 24L) - 1L) * argument_list$plot_width * 0.75

  for (graphics_format in graphics_formats) {
    if (graphics_format == "png" &&
        plot_width > argument_list$plot_limit_png) {
      message("PNG plot exceeding maximum size: ",
              plot_width,
              " > ",
              argument_list$plot_limit_png)
      next
    }

    ggplot2::ggsave(
      filename = file.path(
        argument_list$output_directory,
        paste(
          paste(prefix_summary,
                prefix_pasm,
                "absolute",
                type_list[[type]],
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

  # Plot the fraction of aligned pass-filter reads ------------------------


  ggplot_object <-
    PicardReports::plot_aligned_pass_filter_fractions(object = picard_frame, type = type)

  plot_width <-
    argument_list$plot_width + (ceiling(x = nlevels(x = ggplot_object$data[, type, drop = TRUE]) / 24L) - 1L) * argument_list$plot_width * 0.75

  for (graphics_format in graphics_formats) {
    if (graphics_format == "png" &&
        plot_width > argument_list$plot_limit_png) {
      message("PNG plot exceeding maximum size: ",
              plot_width,
              " > ",
              argument_list$plot_limit_png)
      next
    }

    ggplot2::ggsave(
      filename = file.path(
        argument_list$output_directory,
        paste(
          paste(
            prefix_summary,
            prefix_pasm,
            "percentage",
            type_list[[type]],
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
  rm(graphics_format, plot_width, ggplot_object)

  # Plot the strand balance of aligned pass-filter reads ------------------


  ggplot_object <-
    PicardReports::plot_strand_balance(object = picard_frame, type = type)

  plot_width <-
    argument_list$plot_width + (ceiling(x = nlevels(x = ggplot_object$data[, type, drop = TRUE]) / 24L) - 1L) * argument_list$plot_width * 0.75

  for (graphics_format in graphics_formats) {
    if (graphics_format == "png" &&
        plot_width > argument_list$plot_limit_png) {
      message("PNG plot exceeding maximum size: ",
              plot_width,
              " > ",
              argument_list$plot_limit_png)
      next
    }

    ggplot2::ggsave(
      filename = file.path(
        argument_list$output_directory,
        paste(
          paste(
            prefix_summary,
            prefix_pasm,
            "strand_balance",
            type_list[[type]],
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
  rm(graphics_format, plot_width, ggplot_object)
}
rm(type, type_list, picard_frame)

# Picard Duplication Metrics ----------------------------------------------


# This block is essentially identical with
# PicardReports/exec/picard_reports_duplication.R, but uses different
# file search patterns, as well as plot and TSV file paths.

picard_frame <-
  PicardReports::PicardFrameFromDirectoryPath(
    directory_path = argument_list$directory_path,
    file_pattern = paste0(
      "^",
      argument_list$prefix,
      "_sample_(.*)_duplicate_metrics.tsv$"
    )
  )

# Write the summary tibble in TSV format ----------------------------------


readr::write_tsv(
  x = PicardReports::metrics_tibble(object = picard_frame),
  file = file.path(argument_list$output_directory,
                   paste(
                     paste(prefix_summary,
                           prefix_pdsm,
                           "metrics",
                           "sample",
                           sep = "_"),
                     "tsv",
                     sep = "."
                   ))
)

# Plot the duplication fractions per sample (or library) ------------------


ggplot_object <-
  PicardReports::plot_fractions(object = picard_frame)

plot_width <-
  argument_list$plot_width + (ceiling(x = nrow(x = ggplot_object$data) / 24L) - 1L) * argument_list$plot_width * 0.75

for (graphics_format in graphics_formats) {
  if (graphics_format == "png" &&
      plot_width > argument_list$plot_limit_png) {
    message("PNG plot exceeding maximum size: ",
            plot_width,
            " > ",
            argument_list$plot_limit_png)
    next
  }

  ggplot2::ggsave(
    filename = file.path(
      argument_list$output_directory,
      paste(
        paste(prefix_summary,
              prefix_pdsm,
              "percentage",
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

# Plot the duplication levels per sample (or library) ---------------------


ggplot_object <-
  PicardReports::plot_levels(object = picard_frame)

plot_width <-
  argument_list$plot_width + (ceiling(x = nrow(x = ggplot_object$data) / 24L) - 1L) * argument_list$plot_width * 0.75

for (graphics_format in graphics_formats) {
  if (graphics_format == "png" &&
      plot_width > argument_list$plot_limit_png) {
    message("PNG plot exceeding maximum size: ",
            plot_width,
            " > ",
            argument_list$plot_limit_png)
    next
  }

  ggplot2::ggsave(
    filename = file.path(
      argument_list$output_directory,
      paste(
        paste(prefix_summary,
              prefix_pdsm,
              "levels",
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

rm(
  prefix_pdsm,
  prefix_pasm,
  prefix_summary,
  argument_list,
  graphics_formats
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessioninfo::session_info())

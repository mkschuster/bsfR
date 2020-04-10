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
        opt_str = c("--pattern-file"),
        default = "^run_info\\.json$",
        dest = "pattern_file",
        help = "Kallisto report file name pattern",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--pattern-sample"),
        default = "^kallisto_sample_(.*)$",
        dest = "pattern_sample",
        help = "Kallisto report sample name pattern",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--prefix"),
        default = "kallisto_summary",
        dest = "prefix",
        help = "File name prefix",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--plot-width"),
        default = 7.0,
        dest = "plot_width",
        help = "Plot width in inches [7.0]",
        type = "numeric"
      ),
      optparse::make_option(
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


sample_tibble <- tibble::tibble(
  "file_name" =
    base::list.files(pattern = argument_list$pattern_file, recursive = TRUE),
  "sample_name" = gsub(
    pattern = argument_list$pattern_sample,
    replacement = "\\1",
    x = base::basename(path = base::dirname(path = .data$file_name))
  )
)
message("Processing Kallisto reports for number of samples: ",
        nrow(x = sample_tibble))

kallisto_tibble <- tibble::tibble()

for (i in seq_len(length.out = nrow(x = sample_tibble))) {
  message("  ", sample_tibble[i, "sample_name", drop = TRUE])

  json_list <-
    rjson::fromJSON(file = sample_tibble[i, "file_name", drop = TRUE])

  kallisto_tibble <-
    dplyr::bind_rows(.data = kallisto_tibble, json_list)
  rm(json_list)
}
rm(i)

sample_tibble <- dplyr::bind_cols(sample_tibble, kallisto_tibble)
rm(kallisto_tibble)

readr::write_tsv(
  x = sample_tibble,
  path = paste(paste(argument_list$prefix,
                     "sample",
                     sep = "_"),
               "tsv",
               sep = "."),
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
  argument_list$plot_width + (ceiling(x = nrow(x = sample_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25
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
      angle = 90,
      hjust = 0,
      size = ggplot2::rel(x = 0.7)
    ),
    legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7))
  )
# Scale the plot width with the number of samples, by adding a quarter of
# the original width for each 24 samples.
plot_width <-
  argument_list$plot_width + (ceiling(x = nrow(x = sample_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
  ggplot2::ggsave(
    filename = paste(
      paste(argument_list$prefix,
            "pseudaligned",
            "number",
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
      angle = 90,
      hjust = 0,
      size = ggplot2::rel(x = 0.7)
    ),
    legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7))
  )
# Scale the plot width with the number of samples, by adding a quarter of
# the original width for each 24 samples.
plot_width <-
  argument_list$plot_width + (ceiling(x = nrow(x = sample_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
  ggplot2::ggsave(
    filename = paste(
      paste(argument_list$prefix,
            "pseudaligned",
            "fraction",
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

rm(sample_tibble,
   graphics_formats,
   argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

#!/usr/bin/env Rscript
#
# BSF R script to create summary statistics of aligned versus counted reads.
#
# The script reads the STAR summary statistics of uniquely, multi and unmapped
# reads and joins the total_counts variable added to the column DataFrame of the
# RangedSummarizedExperiment object by the bsf_rnaseq_deseq_analysis.R script.
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
#

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
      opt_str = c("--design-name"),
      # default = "global",
      dest = "design_name",
      help = "Design name",
      type = "character"
    ),
    make_option(
      opt_str = c("--genome-directory"),
      default = ".",
      dest = "genome_directory",
      help = "Genome directory path [.]",
      type = "character"
    ),
    make_option(
      opt_str = c("--output-directory"),
      default = ".",
      dest = "output_directory",
      help = "Output directory path [.]",
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

if (is.null(x = argument_list$design_name)) {
  stop("Missing --design-name option")
}

suppressPackageStartupMessages(expr = library(package = "bsfR"))
suppressPackageStartupMessages(expr = library(package = "tidyverse"))
suppressPackageStartupMessages(expr = library(package = "DESeq2"))

# Save plots in the following formats.
graphics_formats <- c("pdf" = "pdf", "png" = "png")

prefix <-
  bsfR::bsfrd_get_prefix_deseq(design_name = argument_list$design_name)

output_directory <-
  file.path(argument_list$output_directory, prefix)
if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

# Read data ---------------------------------------------------------------


# Read the DESeqDataSet object to get the total counts.

ranged_summarized_experiment <-
  bsfR::bsfrd_read_summarized_experiment(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name
  )

# Extract just the sample name and the total_counts as the number of "counted" reads.
rse_count_tibble <- tibble::tibble(
  "sample" = as.character(x = SummarizedExperiment::colData(x = ranged_summarized_experiment)$sample),
  "counted" = SummarizedExperiment::colData(x = ranged_summarized_experiment)$total_counts
)

# Read the STAR aligner summary tibble.

star_summary_tibble <-
  readr::read_tsv(
    file = file.path(
      argument_list$genome_directory,
      "star_aligner_summary_table_sample.tsv"
    ),
    col_types = readr::cols(
      sample = readr::col_character(),
      reads_input = readr::col_double(),
      reads_unique = readr::col_double(),
      reads_multi = readr::col_double(),
      junctions_total = readr::col_double(),
      junctions_sjdb = readr::col_double(),
      junctions_gtag = readr::col_double(),
      junctions_gcag = readr::col_double(),
      junctions_atac = readr::col_double(),
      junctions_non_canonical = readr::col_double()
    )
  )

# Select just the count variables.
star_count_tibble <- dplyr::select(
  .data = star_summary_tibble,
  .data$sample,
  "total" = .data$reads_input,
  "unique" = .data$reads_unique,
  "multi" = .data$reads_multi
)
# Calculate the unmapped variable.
star_count_tibble <-
  dplyr::mutate(.data = star_count_tibble,
                "unmapped" = .data$total - .data$unique - .data$multi)

# Plot read counts --------------------------------------------------------


merged_count_tibble <-
  dplyr::inner_join(x = star_count_tibble,
                    y = rse_count_tibble,
                    by = "sample")

rm(star_count_tibble, rse_count_tibble)

ggplot_object <- ggplot2::ggplot(
  data = tidyr::gather(
    data = merged_count_tibble,
    key = "status",
    value = "number",
    .data$total,
    .data$unique,
    .data$multi,
    .data$unmapped,
    .data$counted
  )
)

ggplot_object <-
  ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(
    x = .data$sample,
    y = .data$number,
    colour = .data$status
  ))

ggplot_object <-
  ggplot_object + ggplot2::labs(
    x = "Sample",
    y = "Number",
    colour = "Status",
    title = "Number of Reads per Sample",
    subtitle = paste("Design:", argument_list$design_name)
  )

ggplot_object <-
  ggplot_object + ggplot2::theme(axis.text.x = ggplot2::element_text(
    angle = 90,
    hjust = 0,
    size = rel(x = 0.8)
  ))

for (graphics_format in graphics_formats) {
  ggplot2::ggsave(
    filename = file.path(output_directory,
                         paste(
                           paste(prefix,
                                 "alignment",
                                 "summary",
                                 "number",
                                 sep = "_"),
                           graphics_format,
                           sep = "."
                         )),
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
}
rm(graphics_format, ggplot_object)

# Plot read fractions -----------------------------------------------------


merged_fraction_tibble <- tibble::tibble(
  "sample" = merged_count_tibble$sample,
  "unique" = merged_count_tibble$unique / merged_count_tibble$total,
  "multi" = merged_count_tibble$multi / merged_count_tibble$total,
  "unmapped" = merged_count_tibble$unmapped / merged_count_tibble$total,
  "counted" = merged_count_tibble$counted / merged_count_tibble$total
)

ggplot_object <- ggplot2::ggplot(
  data = tidyr::gather(
    data = merged_fraction_tibble,
    key = "status",
    value = "fraction",
    .data$unique,
    .data$multi,
    .data$unmapped,
    .data$counted
  )
)

ggplot_object <-
  ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(
    x = .data$sample,
    y = .data$fraction,
    colour = .data$status
  ))

ggplot_object <-
  ggplot_object + ggplot2::labs(
    x = "Sample",
    y = "Fraction",
    colour = "Status",
    title = "Fraction of Reads per Sample",
    subtitle = paste("Design:", argument_list$design_name)
  )

ggplot_object <-
  ggplot_object + ggplot2::theme(axis.text.x = ggplot2::element_text(
    angle = 90,
    hjust = 0,
    size = rel(x = 0.8)
  ))

for (graphics_format in graphics_formats) {
  ggplot2::ggsave(
    filename = file.path(output_directory,
                         paste(
                           paste(prefix,
                                 "alignment",
                                 "summary",
                                 "fraction",
                                 sep = "_"),
                           graphics_format,
                           sep = "."
                         )),
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
}
rm(graphics_format, ggplot_object)

rm(
  merged_fraction_tibble,
  merged_count_tibble,
  star_summary_tibble,
  ranged_summarized_experiment,
  output_directory,
  prefix,
  graphics_formats,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessionInfo())
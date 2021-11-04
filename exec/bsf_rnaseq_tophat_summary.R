#!/usr/bin/env Rscript
#
# BSF R script to parse and summarise Tophat2 alignment summary files.
# The resulting data frame (rnaseq_tophat_alignment_summary.tsv),
# as well as plots of alignment rates per sample in PDF (rnaseq_tophat_alignment_summary.pdf)
# and PNG format (rnaseq_tophat_alignment_summary.png) are written into the current working
# directory.
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
        opt_str = c("--pattern-path"),
        default = "/rnaseq_tophat_.*$",
        dest = "pattern_path",
        help = "Tophat2 directory path pattern [/rnaseq_tophat_.*$]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--pattern-sample"),
        default = ".*/rnaseq_tophat_(.*)$",
        dest = "pattern_sample",
        help = "Tophat2 sample name pattern [.*/rnaseq_tophat_(.*)$]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--genome-directory"),
        default = ".",
        dest = "genome_directory",
        help = "Genome directory path [.]",
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
        opt_str = c("--plot-factor"),
        default = 0.5,
        dest = "plot_factor",
        help = "Plot width increase per 24 samples [0.5]",
        type = "numeric"
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

suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
suppressPackageStartupMessages(expr = library(package = "tidyverse"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

# Process all "rnaseq_tophat_*" directories in the current working directory.

message("Processing Tophat2 alignment report for sample:")
sample_tibble <- tibble::tibble(
  # R character vector of directory paths.
  "directory_path" = grep(
    pattern = argument_list$pattern_path,
    x = list.dirs(
      path = argument_list$genome_directory,
      full.names = TRUE,
      recursive = FALSE
    ),
    value = TRUE
  ),
  # R character vector of sample names.
  "sample_name" = base::gsub(
    pattern = argument_list$pattern_sample,
    replacement = "\\1",
    x = .data$directory_path
  ),
  # R integer vector of the number of input reads.
  "input" = NA_integer_,
  # R integer vector of the number of mapped reads.
  "mapped" = NA_integer_,
  # R integer vector of the number of multiply mapped reads.
  "multiple" = NA_integer_,
  # R integer vector of the number of multiply mapped reads above the threshold.
  "above" = NA_integer_,
  # R integer vector of the mapping threshold.
  "threshold" = NA_integer_
)

for (i in seq_len(length.out = base::nrow(x = sample_tibble))) {
  message("  ", sample_tibble$sample_name[i])

  # This is the layout of a Tophat align_summary.txt file.
  #
  #   [1] "Reads:"
  #   [2] "          Input     :  21791622"
  #   [3] "           Mapped   :  21518402 (98.7% of input)"
  #   [4] "            of these:   2010462 ( 9.3%) have multiple alignments (8356 have >20)"
  #   [5] "98.7% overall read mapping rate."

  align_summary <-
    readLines(
      con = file.path(
        argument_list$genome_directory,
        sample_tibble$directory_path[i],
        "align_summary.txt"
      ),
      n = 100L
    )

  # Parse the second line of "input" reads.
  sample_tibble$input[i] <- as.integer(
    x = sub(
      pattern = "[[:space:]]+Input[[:space:]]+:[[:space:]]+([[:digit:]]+)",
      replacement = "\\1",
      x = align_summary[2L]
    )
  )

  # Parse the third line of "mapped" reads.
  sample_tibble$mapped[i] <- as.integer(
    x = sub(
      pattern = "[[:space:]]+Mapped[[:space:]]+:[[:space:]]+([[:digit:]]+) .*",
      replacement = "\\1",
      x = align_summary[3L]
    )
  )

  # Get the number of "multiply" aligned reads from the fourth line.
  sample_tibble$multiple[i] <- as.integer(
    x = sub(
      pattern = ".+:[[:space:]]+([[:digit:]]+) .*",
      replacement = "\\1",
      x = align_summary[4L]
    )
  )

  # Get the number of multiply aligned reads "above" the multiple alignment threshold.
  sample_tibble$above[i] <- as.integer(
    x = sub(
      pattern = ".+alignments \\(([[:digit:]]+) have.+",
      replacement = "\\1",
      x = align_summary[4L]
    )
  )

  # Get the multiple alignment "threshold".
  sample_tibble$threshold[i] <- as.integer(x = sub(
    pattern = ".+ >([[:digit:]]+)\\)",
    replacement = "\\1",
    x = align_summary[4L]
  ))

  rm(align_summary)
}
rm(i)

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
    title = "TopHat Alignment Summary"
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
  graphics_formats,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessioninfo::session_info())

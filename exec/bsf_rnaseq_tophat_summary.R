#!/usr/bin/env Rscript
#
# BSF R script to parse and summarise Tophat alignment summary files.
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
      opt_str = c("--genome-version"),
      default = NULL,
      dest = "genome_version",
      help = "Genome version",
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

suppressPackageStartupMessages(expr = library(package = "tidyverse"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

#' Parse Tophat alignment summary (align_summary.txt) files and
#' return a data frame.
#'
#' @param summary_frame: Data frame with alignment summary statistics
#' @return: Data frame with alignment summary statistics

process_align_summary <- function(summary_frame) {
  if (is.null(x = summary_frame)) {
    stop("Missing summary_frame argument")
  }

  # Initialise further columns in the data frame.
  frame_length <- nrow(x = summary_frame)
  # R integer vector of the number of input reads.
  summary_frame$input = integer(length = frame_length)
  # R integer vector of the number of mapped reads.
  summary_frame$mapped = integer(length = frame_length)
  # R integer vector of the number of multiply mapped reads.
  summary_frame$multiple = integer(length = frame_length)
  # R integer vector of the number of multiply mapped reads above the threshold.
  summary_frame$above = integer(length = frame_length)
  # R integer vector of the mapping threshold.
  summary_frame$threshold = integer(length = frame_length)
  rm(frame_length)

  for (i in seq_len(length.out = nrow(x = summary_frame))) {
    prefix_tophat <-
      paste("rnaseq", "tophat", summary_frame[i, "sample"], sep = "_")
    file_path <- file.path(prefix_tophat, "align_summary.txt")

    if (!file.exists(file_path)) {
      warning("Missing Tophat alignment summary file ", file_path)
      rm(prefix_tophat, file_path)
      next
    }

    align_summary <- readLines(con = file_path)

    # This is the layout of a Tophat align_summary.txt file.
    #
    #   [1] "Reads:"
    #   [2] "          Input     :  21791622"
    #   [3] "           Mapped   :  21518402 (98.7% of input)"
    #   [4] "            of these:   2010462 ( 9.3%) have multiple alignments (8356 have >20)"
    #   [5] "98.7% overall read mapping rate."

    # Parse the second line of "input" reads.
    summary_frame[i, "input"] <- as.integer(
      x = sub(
        pattern = "[[:space:]]+Input[[:space:]]+:[[:space:]]+([[:digit:]]+)",
        replacement = "\\1",
        x = align_summary[2L]
      )
    )

    # Parse the third line of "mapped" reads.
    summary_frame[i, "mapped"] <- as.integer(
      x = sub(
        pattern = "[[:space:]]+Mapped[[:space:]]+:[[:space:]]+([[:digit:]]+) .*",
        replacement = "\\1",
        x = align_summary[3L]
      )
    )

    # Get the number of "multiply" aligned reads from the fourth line.
    summary_frame[i, "multiple"] <- as.integer(
      x = sub(
        pattern = ".+:[[:space:]]+([[:digit:]]+) .*",
        replacement = "\\1",
        x = align_summary[4L]
      )
    )

    # Get the number of multiply aligned reads "above" the multiple alignment threshold.
    summary_frame[i, "above"] <- as.integer(
      x = sub(
        pattern = ".+alignments \\(([[:digit:]]+) have.+",
        replacement = "\\1",
        x = align_summary[4L]
      )
    )

    # Get the multiple alignment "threshold".
    summary_frame[i, "threshold"] <- as.integer(x = sub(
      pattern = ".+ >([[:digit:]]+)\\)",
      replacement = "\\1",
      x = align_summary[4L]
    ))

    rm(align_summary, prefix_tophat)
  }
  rm(i)

  return(summary_frame)
}

# Process all "rnaseq_tophat_*" directories in the current working directory.
# Initialise a data frame with row names of all "rnaseq_tophat_*" directories
# via their common prefix and parse the sample name simply by removing the prefix.

summary_frame <- data.frame(
  row.names = sub(
    pattern = "^rnaseq_tophat_",
    replacement = "",
    x = grep(
      pattern = '^rnaseq_tophat_',
      x = list.dirs(full.names = FALSE, recursive = FALSE),
      value = TRUE
    )
  ),
  stringsAsFactors = FALSE
)

# Set the sample also explictly as a data.frame column, required for plotting.

summary_frame$sample <- row.names(x = summary_frame)
summary_frame <-
  process_align_summary(summary_frame = summary_frame)

# Write the alignment summary frame to disk and create plots.

file_path <- "rnaseq_tophat_alignment_summary.tsv"
write.table(
  x = summary_frame,
  file = file_path,
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)
rm(file_path)

# Alignment Summary Plot --------------------------------------------------


ggplot_object <- ggplot2::ggplot(data = summary_frame)
ggplot_object <-
  ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(
    x = .data$mapped,
    y = .data$mapped / .data$input,
    colour = .data$sample
  ))
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
  argument_list$plot_width + (ceiling(x = nrow(x = summary_frame) / 24L) - 1L) * argument_list$plot_width * argument_list$plot_factor

for (graphics_format in graphics_formats) {
  ggplot2::ggsave(
    filename = paste("rnaseq_tophat_alignment_summary", graphics_format, sep = "."),
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height,
    limitsize = FALSE
  )
}

rm(
  graphics_format,
  ggplot_object,
  plot_width,
  summary_frame,
  graphics_formats,
  argument_list,
  process_align_summary
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

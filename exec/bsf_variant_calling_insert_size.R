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


# BSF R script to compile and plot insert size information of a variant calling analysis.

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
        opt_str = "--file-path",
        dest = "file_path",
        help = "BAM file path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--chunk-size",
        default = 1000000L,
        dest = "chunk_size",
        help = "Chunk size, i.e. number of BAM lines to process at once [1,000,000]",
        type = "integer"
      ),
      optparse::make_option(
        opt_str = "--density-plot",
        action = "store_true",
        default = FALSE,
        dest = "density_plot",
        help = "Draw a density rather than a column plot",
        type = "logical"
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

if (is.null(x = argument_list$file_path)) {
  stop("Missing --file-path option")
}

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "ggplot2"))
suppressPackageStartupMessages(expr = library(package = "readr"))
suppressPackageStartupMessages(expr = library(package = "tibble"))
# Bioconductor
suppressPackageStartupMessages(expr = library(package = "BiocVersion"))
suppressPackageStartupMessages(expr = library(package = "Rsamtools"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

prefix <- "variant_calling_diagnose_sample"

sample_name <-
  base::gsub(
    pattern = "^variant_calling_process_sample_(.*)_realigned.bam$",
    replacement = "\\1",
    x = base::basename(path = argument_list$file_path)
  )

# Scan BAM files ----------------------------------------------------------


records_read_total <- 0L
records_read_chunk <- 0L
summary_vector <- integer()

bam_file <- BamFile(file = argument_list$file_path,
                    yieldSize = argument_list$chunk_size)
open(con = bam_file)
while (TRUE) {
  region_list <- Rsamtools::scanBam(
    file = bam_file,
    param = Rsamtools::ScanBamParam(
      flag = Rsamtools::scanBamFlag(
        isPaired = TRUE,
        isProperPair = TRUE,
        isSecondaryAlignment = FALSE,
        isNotPassingQualityControls = FALSE,
        isDuplicate = FALSE
      ),
      what = c("isize")
    )
  )
  # The Rsamtools::scanBam() function returns a list of lists.
  records_read_chunk <- 0L
  for (i in seq_along(along.with = region_list)) {
    # Add 1L to adjust to R vector indices starting at 1.
    chunk_vector <-
      tabulate(bin = abs(x = region_list[[i]]$isize + 1L))
    maximum_length <-
      max(length(x = summary_vector), length(x = chunk_vector))
    length(x = summary_vector) <- maximum_length
    length(x = chunk_vector) <- maximum_length
    summary_vector[is.na(x = summary_vector)] <- 0L
    chunk_vector[is.na(x = chunk_vector)] <- 0L
    summary_vector <- summary_vector + chunk_vector
    records_read_chunk <-
      records_read_chunk + length(x = region_list[[i]]$isize)
    rm(maximum_length, chunk_vector)
  }
  rm(i)
  records_read_total <- records_read_total + records_read_chunk
  # print(x = paste0("Records processed: ", records_read_total))
  rm(region_list)
  if (records_read_chunk == 0L) {
    break()
  }
}
close(con = bam_file)
rm(bam_file)

summary_tibble <-
  tibble::tibble(size = seq_along(along.with = summary_vector) - 1L,
                 # Subtract 1L to adjust to R vector indices starting at 1.
                 frequency = summary_vector)
rm(summary_vector)

readr::write_tsv(x = summary_tibble,
                 file = paste(paste(prefix,
                                    sample_name,
                                    "insert_size",
                                    sep = "_"),
                              "tsv",
                              sep = "."))

# Plot insert size distribution -------------------------------------------


plot_paths <- paste(paste(prefix,
                          sample_name,
                          "insert_size",
                          sep = "_"),
                    graphics_formats,
                    sep = ".")

ggplot_object <- ggplot2::ggplot(data = summary_tibble)
if (argument_list$density_plot) {
  ggplot_object <-
    ggplot_object + ggplot2::geom_density(mapping = ggplot2::aes(x = .data$size, y = .data$frequency),
                                          stat = "identity")
} else {
  ggplot_object <-
    ggplot_object + ggplot2::geom_col(mapping = ggplot2::aes(x = .data$size, y = .data$frequency))
}
ggplot_object <-
  ggplot_object + ggplot2::labs(
    x = "Insert Size",
    y = "Frequency",
    title = "Insert Size Density",
    subtitle = sample_name
  )
for (plot_path in plot_paths) {
  ggplot2::ggsave(
    filename = plot_path,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height,
    limitsize = FALSE
  )
}
rm(
  plot_path,
  ggplot_object,
  plot_paths,
  summary_tibble,
  records_read_chunk,
  records_read_total,
  sample_name,
  prefix,
  argument_list,
  graphics_formats
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessioninfo::session_info())

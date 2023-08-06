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


# BSF R script to summarise STAR aligner alignment reports.

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
        opt_str = "--directory-path",
        default = ".",
        dest = "directory_path",
        help = "STAR alignment report directory path [.]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--file-pattern",
        default = "^star_align_(.*)_Log\\.final\\.out$",
        dest = "file_pattern",
        help = "STAR alignment report sample name pattern [^star_align_(.*)_Log\\.final\\.out$]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--sample-table-path",
        dest = "sample_table_path",
        help = "Sample to Read Group TSV table",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--prefix",
        default = "star_summary",
        dest = "prefix",
        help = "File name prefix",
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
suppressPackageStartupMessages(expr = library(package = "dplyr"))
suppressPackageStartupMessages(expr = library(package = "ggplot2"))
suppressPackageStartupMessages(expr = library(package = "readr"))
suppressPackageStartupMessages(expr = library(package = "stringr"))
suppressPackageStartupMessages(expr = library(package = "tibble"))
suppressPackageStartupMessages(expr = library(package = "tidyr"))

# Save plots in the following formats.
graphics_formats <- c("pdf" = "pdf", "png" = "png")

# Parse STAR aligner log files --------------------------------------------


parse_report <-
  function(file_path) {
    # Empty lines and single column rows render parsing STAR Log.final.out files
    # more cumbersome than necessary. Use base::readLines() and
    # stringr::str_split_fixed() to split into a matrix.
    #
    # > stringr::str_split_fixed(string = star_lines, pattern = fixed(pattern = "\t"), n = 2)
    # [,1]                                                [,2]
    # [01,] "                                 Started job on |" "Apr 01 20:49:14"
    # [02,] "                             Started mapping on |" "Apr 01 20:50:12"
    # [03,] "                                    Finished on |" "Apr 01 20:51:57"
    # [04,] "       Mapping speed, Million of reads per hour |" "43.78"
    # [05,] ""                                                  ""
    # [06,] "                          Number of input reads |" "1276962"
    # [07,] "                      Average input read length |" "302"
    # [08,] "                                    UNIQUE READS:" ""
    # [09,] "                   Uniquely mapped reads number |" "1128810"
    # [10,] "                        Uniquely mapped reads % |" "88.40%"
    # [11,] "                          Average mapped length |" "297.00"
    # [12,] "                       Number of splices: Total |" "17699"
    # [13,] "            Number of splices: Annotated (sjdb) |" "8"
    # [14,] "                       Number of splices: GT/AG |" "475"
    # [15,] "                       Number of splices: GC/AG |" "16"
    # [16,] "                       Number of splices: AT/AC |" "0"
    # [17,] "               Number of splices: Non-canonical |" "17208"
    # [18,] "                      Mismatch rate per base, % |" "0.92%"
    # [19,] "                         Deletion rate per base |" "0.01%"
    # [20,] "                        Deletion average length |" "2.05"
    # [21,] "                        Insertion rate per base |" "0.00%"
    # [22,] "                       Insertion average length |" "1.07"
    # [23,] "                             MULTI-MAPPING READS:" ""
    # [24,] "        Number of reads mapped to multiple loci |" "3276"
    # [25,] "             % of reads mapped to multiple loci |" "0.26%"
    # [26,] "        Number of reads mapped to too many loci |" "1"
    # [27,] "             % of reads mapped to too many loci |" "0.00%"
    # [28,] "                                  UNMAPPED READS:" ""
    # [29,] "       % of reads unmapped: too many mismatches |" "0.00%"
    # [30,] "                 % of reads unmapped: too short |" "11.28%"
    # [31,] "                     % of reads unmapped: other |" "0.06%"
    # [32,] "                                  CHIMERIC READS:" ""
    # [33,] "                       Number of chimeric reads |" "0"
    # [34,] "                            % of chimeric reads |" "0.00%"

    star_character <- stringr::str_split_fixed(
      string = base::readLines(con = file_path),
      pattern = "\t",
      n = 2L
    )[c(1L:4L, 6L:7L, 9L:22L, 24L:27L, 29L:31L, 33L:34L), 2L]

    names(x = star_character) <- c(
      "started_job",
      "started_mapping",
      "finished",
      "mapping_speed",
      # [5,]
      "input_reads",
      "average_length",
      # [8,] UNIQUE READS
      "uniquely_mapped_reads",
      "uniquely_mapped_percentage",
      "average_mapped_length",
      "number_splice_total",
      "number_splice_sjdb",
      "number_splice_gtag",
      "number_splice_gcag",
      "number_splice_atac",
      "number_splice_non_canonical",
      "mismatch_rate",
      "deletion_rate",
      "deletion_average_length",
      "insertion_rate",
      "insertion_average_length",
      # [23,] MULTI-MAPPING READS
      "multi_mapped_number",
      "multi_mapped_percentage",
      "multi_unmapped_number",
      "multi_unmapped_percentage",
      # [28,] UNMAPPED READS
      "unmapped_mismatched_percentage",
      "unmapped_short_percentage",
      "unmapped_other_percentage",
      # [32,] CHIMERIC READS
      "chimeric_number",
      "chimeric_percentage"
    )

    return(star_character)
  }

read_group_tibble <-
  bsfR::bsfu_summarise_report_files(
    directory_path = argument_list$directory_path,
    file_pattern = argument_list$file_pattern,
    report_function = parse_report,
    variable_name = "read_group_name"
  )

# Convert types.
read_group_tibble <-
  readr::type_convert(
    df = read_group_tibble,
    col_types = readr::cols(
      "read_group_name" = readr::col_character(),
      # The following three dates lack the year, so are not terribly useful.
      "started_job" = readr::col_character(),
      "started_mapping" = readr::col_character(),
      "finished" = readr::col_character(),
      "mapping_speed" = readr::col_double(),
      # [5,]
      "input_reads" = readr::col_integer(),
      "average_length" = readr::col_double(),
      # [8,] UNIQUE READS
      "uniquely_mapped_reads" = readr::col_integer(),
      "uniquely_mapped_percentage" = readr::col_double(),
      "average_mapped_length" = readr::col_double(),
      "number_splice_total" = readr::col_integer(),
      "number_splice_sjdb" = readr::col_integer(),
      "number_splice_gtag" = readr::col_integer(),
      "number_splice_gcag" = readr::col_integer(),
      "number_splice_atac" = readr::col_integer(),
      "number_splice_non_canonical" = readr::col_integer(),
      "mismatch_rate" = readr::col_double(),
      "deletion_rate" = readr::col_double(),
      "deletion_average_length" = readr::col_double(),
      "insertion_rate" = readr::col_double(),
      "insertion_average_length" = readr::col_double(),
      # [23,] MULTI-MAPPING READS
      "multi_mapped_number" = readr::col_integer(),
      "multi_mapped_percentage" = readr::col_double(),
      "multi_unmapped_number" = readr::col_integer(),
      "multi_unmapped_percentage" = readr::col_double(),
      # [28,] UNMAPPED READS
      "unmapped_mismatched_percentage" = readr::col_double(),
      "unmapped_short_percentage" = readr::col_double(),
      "unmapped_other_percentage" = readr::col_double(),
      # [32,] CHIMERIC READS
      "chimeric_number" = readr::col_integer(),
      "chimeric_percentage" = readr::col_double(),
      "file_path" = readr::col_character(),
    )
  )

message("Writing read group-level summary table")
readr::write_tsv(
  x = read_group_tibble,
  file = file.path(argument_list$output_directory, paste(
    paste(argument_list$prefix, "table", "read_group", sep = "_"),
    "tsv",
    sep = "."
  )),
  col_names = TRUE
)

# Scatter plot of read number versus alignment rate per read group --------


message("Creating a scatter plot of read number versus alignment rate per read group")

ggplot_object <- ggplot2::ggplot(
  data = tidyr::pivot_longer(
    data = dplyr::transmute(
      .data = read_group_tibble,
      "read_group" = .data$read_group_name,
      "input" = .data$input_reads,
      "unique" = .data$uniquely_mapped_reads,
      "multi" = .data$multi_mapped_number,
      "unmapped" = .data$input_reads - .data$uniquely_mapped_reads - .data$multi_mapped_number
    ),
    cols = c("unmapped", "multi", "unique"),
    names_to = "status",
    values_to = "mapped"
  )
)

ggplot_object <-
  ggplot_object + ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = .data$mapped,
      y = .data$mapped / .data$input,
      colour = .data$read_group,
      shape = .data$status
    ),
    alpha = I(1 / 3)
  )

ggplot_object <-
  ggplot_object + ggplot2::labs(
    x = "Reads Number",
    y = "Reads Fraction",
    colour = "Read Group",
    shape = "Mapping Status",
    title = "STAR Alignment Summary per Read Group"
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

# Scale the plot width with the number of read groups, by adding a quarter of
# the original width for each 24 read groups.
# Because read group names are quite long, extend already for the first column.
plot_width <-
  argument_list$plot_width + (ceiling(x = base::nrow(x = read_group_tibble) / 24L) - 1L) * argument_list$plot_width * 0.75

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
        paste(argument_list$prefix,
              "alignment",
              "read_group",
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

# Column plot of read numbers per read group ------------------------------


message("Creating a column plot of read numbers per read group")

ggplot_object <- ggplot2::ggplot(
  data = tidyr::pivot_longer(
    data = dplyr::transmute(
      .data = read_group_tibble,
      "read_group" = .data$read_group_name,
      "unique" = .data$uniquely_mapped_reads,
      "multi" = .data$multi_mapped_number,
      "unmapped" =
        .data$input_reads - .data$uniquely_mapped_reads - .data$multi_mapped_number
    ),
    cols = c("unmapped", "multi", "unique"),
    names_to = "status",
    values_to = "number"
  )
)

ggplot_object <-
  ggplot_object + ggplot2::geom_col(
    mapping = ggplot2::aes(
      x = .data$read_group,
      y = .data$number,
      fill = .data$status
    ),
    alpha = I(1 / 3)
  )

ggplot_object <-
  ggplot_object + ggplot2::labs(x = "Read Group",
                                y = "Reads Number",
                                fill = "Mapping Status",
                                title = "STAR Aligner Mapped Numbers per Read Group")

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

# Scale the plot width with the number of read groups, by adding a quarter of
# the original width for each 24 read groups.
plot_width <-
  argument_list$plot_width + (ceiling(x = base::nrow(x = read_group_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25

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
        paste(argument_list$prefix,
              "mapped",
              "number",
              "read_group",
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

# Column plot of read fractions per read group ----------------------------


message("Creating a column plot of read fractions per read group")

ggplot_object <- ggplot2::ggplot(
  data = tidyr::pivot_longer(
    data = dplyr::transmute(
      .data = read_group_tibble,
      "read_group" = .data$read_group_name,
      "unique" = .data$uniquely_mapped_reads / .data$input_reads,
      "multi" = .data$multi_mapped_number / .data$input_reads,
      "unmapped" = (
        .data$input_reads - .data$uniquely_mapped_reads - .data$multi_mapped_number
      ) / .data$input_reads
    ),
    cols = c("unmapped", "multi", "unique"),
    names_to = "status",
    values_to = "fraction"
  )
)

ggplot_object <-
  ggplot_object + ggplot2::geom_col(
    mapping = ggplot2::aes(
      x = .data$read_group,
      y = .data$fraction,
      fill = .data$status
    ),
    alpha = I(1 / 3)
  )

ggplot_object <-
  ggplot_object + ggplot2::labs(x = "Read Group",
                                y = "Reads Fraction",
                                fill = "Mapping Status",
                                title = "STAR Aligner Mapped Fractions per Read Group")

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

# Scale the plot width with the number of read groups, by adding a quarter of
# the original width for each 24 read groups.
plot_width <-
  argument_list$plot_width + (ceiling(x = base::nrow(x = read_group_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25

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
          argument_list$prefix,
          "mapped",
          "fraction",
          "read_group",
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

# Column plot of splice junction numbers per read group -------------------


message("Creating a column plot of splice junction numbers per read group")

ggplot_object <- ggplot2::ggplot(
  data = tidyr::pivot_longer(
    data = dplyr::select(
      .data = read_group_tibble,
      "read_group" = "read_group_name",
      "gtag" = "number_splice_gtag",
      "gcag" = "number_splice_gcag",
      "atac" = "number_splice_atac",
      "non_canonical" = "number_splice_non_canonical"
    ),
    cols = c("non_canonical", "atac", "gcag", "gtag"),
    names_to = "splice_junction",
    values_to = "number"
  )
)

ggplot_object <-
  ggplot_object + ggplot2::geom_col(
    mapping = ggplot2::aes(
      x = .data$read_group,
      y = .data$number,
      fill = .data$splice_junction
    ),
    alpha = I(1 / 3)
  )

ggplot_object <-
  ggplot_object + ggplot2::labs(x = "Read Group",
                                y = "Splice Junction Number",
                                fill = "Splice Junction",
                                title = "STAR Aligner Splice Junction Numbers per Read Group")

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

# Scale the plot width with the number of read groups, by adding a quarter of
# the original width for each 24 read groups.
plot_width <-
  argument_list$plot_width + (ceiling(x = base::nrow(x = read_group_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25

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
          argument_list$prefix,
          "junction",
          "number",
          "read_group",
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

# Column plot of splice junction fractions per read group -----------------


message("Creating a column plot of splice junction fractions per read group")

ggplot_object <- ggplot2::ggplot(
  data = tidyr::pivot_longer(
    data = dplyr::transmute(
      .data = read_group_tibble,
      "read_group" = .data$read_group_name,
      "gtag" = .data$number_splice_gtag / .data$number_splice_total,
      "gcag" = .data$number_splice_gcag / .data$number_splice_total,
      "atac" = .data$number_splice_atac / .data$number_splice_total,
      "non_canonical" = .data$number_splice_non_canonical / .data$number_splice_total
    ),
    cols = c("non_canonical", "atac", "gcag", "gtag"),
    names_to = "junction",
    values_to = "fraction"
  )
)

ggplot_object <-
  ggplot_object + ggplot2::geom_col(
    mapping = ggplot2::aes(
      x = .data$read_group,
      y = .data$fraction,
      fill = .data$junction
    ),
    alpha = I(1 / 3)
  )

ggplot_object <-
  ggplot_object + ggplot2::labs(x = "Read Group",
                                y = "Splice Junction Fraction",
                                fill = "Splice Junction",
                                title = "STAR Aligner Splice Junction Fractions per Read Group")

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

# Scale the plot width with the number of read groups, by adding a quarter of
# the original width for each 24 samples.
plot_width <-
  argument_list$plot_width + (ceiling(x = base::nrow(x = read_group_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25

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
          argument_list$prefix,
          "junction",
          "fraction",
          "read_group",
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

# Integrate read group data on sample level -------------------------------


# Can the read group-level data be integrated on the sample level?

if (file.exists(argument_list$sample_table_path)) {
  message("Integrating on sample-level ...")
  sample_tibble <-
    dplyr::left_join(
      # Read the read group to sample mapping data frame provided by BSF Python.
      x = readr::read_tsv(
        file = argument_list$sample_table_path,
        col_names = TRUE,
        col_types = readr::cols(
          "sample" = readr::col_character(),
          "read_group" = readr::col_character()
        )
      ),
      # Create a smaller tibble for mapping and junction information and
      # include the read group for merging.
      y = dplyr::select(
        .data = read_group_tibble,
        "read_group" = "read_group_name",
        "reads_input" = "input_reads",
        "reads_unique" = "uniquely_mapped_reads",
        "reads_multi" = "multi_mapped_number",
        "junctions_total" = "number_splice_total",
        "junctions_sjdb" = "number_splice_sjdb",
        "junctions_gtag" = "number_splice_gtag",
        "junctions_gcag" = "number_splice_gcag",
        "junctions_atac" = "number_splice_atac",
        "junctions_non_canonical" = "number_splice_non_canonical",
      ),
      by = "read_group"
    )

  sample_tibble <-
    dplyr::group_by(.data = sample_tibble, .data$sample)

  sample_tibble <-
    dplyr::summarise(.data = sample_tibble, dplyr::across(
      .cols = c(
        "reads_input",
        "reads_unique",
        "reads_multi",
        "junctions_total",
        "junctions_sjdb",
        "junctions_gtag",
        "junctions_gcag",
        "junctions_atac",
        "junctions_non_canonical"
      ),
      .fns = ~ sum(.x)
    ))

  # Write the sample-level tibble with mapping and junction information to disk.
  message("Writing sample-level summary table")
  readr::write_tsv(
    x = sample_tibble,
    file = file.path(
      argument_list$output_directory,
      paste(
        paste(argument_list$prefix, "table", "sample", sep = "_"),
        "tsv",
        sep = "."
      )
    ),
    col_names = TRUE
  )

  # For plotting, add calculations of mapped and unmapped reads.
  sample_tibble <- dplyr::mutate(
    .data = sample_tibble,
    "reads_mapped" = .data$reads_unique + .data$reads_multi,
    "reads_unmapped" = .data$reads_input - .data$reads_mapped
  )

  # Scatter plot of read number versus alignment rate per sample ----------


  message("Creating a scatter plot of read number versus alignment rate per sample")

  ggplot_object <- ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = dplyr::select(
        .data = sample_tibble,
        "sample" = "sample",
        "input" = "reads_input",
        "multi" = "reads_multi",
        "unique" = "reads_unique",
        "unmapped" = "reads_unmapped"
      ),
      cols = c("unmapped", "multi", "unique"),
      names_to = "status",
      values_to = "mapped"
    )
  )

  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$mapped,
        y = .data$mapped / .data$input,
        colour = .data$sample,
        shape = .data$status
      ),
      alpha = I(1 / 3)
    )

  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Reads Number",
      y = "Reads Fraction",
      colour = "Sample",
      shape = "Mapping Status",
      title = "STAR Alignment Summary per Sample"
    )

  # Reduce the label font size and the legend key size and allow a maximum of 24
  # guide legend rows.
  ggplot_object <-
    ggplot_object + ggplot2::guides(
      colour = ggplot2::guide_legend(
        keywidth = ggplot2::rel(x = 0.8),
        keyheight = ggplot2::rel(x = 0.8),
        nrow = 24L
      )
    )

  ggplot_object <-
    ggplot_object + ggplot2::theme(legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7)))

  # Scale the plot width with the number of read groups, by adding a quarter of
  # the original width for each 24 read groups.
  plot_width <-
    argument_list$plot_width + (ceiling(x = base::nrow(x = sample_tibble) / 24L) - 1L) * argument_list$plot_width * 0.33

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
          paste(argument_list$prefix,
                "alignment",
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

  ggplot_object <- ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = dplyr::select(
        .data = sample_tibble,
        "sample" = "sample",
        "unmapped" = "reads_unmapped",
        "multi" = "reads_multi",
        "unique" = "reads_unique"
      ),
      cols = c("unmapped", "multi", "unique"),
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
    ggplot_object + ggplot2::labs(
      x = "Sample",
      y = "Reads Number",
      fill = "Mapping Status",
      title = "STAR Aligner Mapped Numbers per Sample"
    )

  # Reduce the label font size and the legend key size and allow a maximum of 24
  # guide legend rows.
  ggplot_object <-
    ggplot_object + ggplot2::guides(
      colour = ggplot2::guide_legend(
        keywidth = ggplot2::rel(x = 0.8),
        keyheight = ggplot2::rel(x = 0.8),
        nrow = 24L
      )
    )

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
          paste(argument_list$prefix,
                "mapped",
                "number",
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

  # Column plot of read fractions per sample ------------------------------


  message("Creating a column plot of read fractions per sample")

  ggplot_object <- ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = dplyr::transmute(
        .data = sample_tibble,
        "sample" = .data$sample,
        "unmapped" = .data$reads_unmapped / .data$reads_input,
        "multi" = .data$reads_multi / .data$reads_input,
        "unique" = .data$reads_unique / .data$reads_input
      ),
      cols = c("unmapped", "multi", "unique"),
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
    ggplot_object + ggplot2::labs(
      x = "Sample",
      y = "Reads Fraction",
      fill = "Mapping Status",
      title = "STAR Aligner Mapped Fractions per Sample"
    )

  # Reduce the label font size and the legend key size and allow a maximum of 24
  # guide legend rows.
  ggplot_object <-
    ggplot_object + ggplot2::guides(
      colour = ggplot2::guide_legend(
        keywidth = ggplot2::rel(x = 0.8),
        keyheight = ggplot2::rel(x = 0.8),
        nrow = 24L
      )
    )

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
          paste(argument_list$prefix,
                "mapped",
                "fraction",
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

  # Column plot of splice junction numbers per sample ---------------------


  message("Creating a column plot of splice junction numbers per sample")

  ggplot_object <- ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = dplyr::select(
        .data = sample_tibble,
        "sample" = "sample",
        "gtag" = "junctions_gtag",
        "gcag" = "junctions_gcag",
        "atac" = "junctions_atac",
        "non_canonical" = "junctions_non_canonical"
      ),
      cols = c("non_canonical", "atac", "gcag", "gtag"),
      names_to = "junction",
      values_to = "number"
    )
  )

  ggplot_object <-
    ggplot_object + ggplot2::geom_col(
      mapping = ggplot2::aes(
        x = .data$sample,
        y = .data$number,
        fill = .data$junction
      ),
      alpha = I(1 / 3)
    )

  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Sample",
      y = "Splice Junction Number",
      fill = "Splice Junction",
      title = "STAR Aligner Splice Junction Numbers per Sample"
    )

  # Reduce the label font size and the legend key size and allow a maximum of 24
  # guide legend rows.
  ggplot_object <-
    ggplot_object + ggplot2::guides(
      colour = ggplot2::guide_legend(
        keywidth = ggplot2::rel(x = 0.8),
        keyheight = ggplot2::rel(x = 0.8),
        nrow = 24L
      )
    )

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
          paste(argument_list$prefix,
                "junction",
                "number",
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

  # Column plot of splice junction fractions per sample -------------------


  message("Creating a column plot of splice junction fractions per sample")

  ggplot_object <- ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = dplyr::transmute(
        .data = sample_tibble,
        "sample" = .data$sample,
        "gtag" = .data$junctions_gtag / .data$junctions_total,
        "gcag" = .data$junctions_gcag / .data$junctions_total,
        "atac" = .data$junctions_atac / .data$junctions_total,
        "non_canonical" = .data$junctions_non_canonical / .data$junctions_total
      ),
      cols = c("non_canonical", "atac", "gcag", "gtag"),
      names_to = "junction",
      values_to = "fraction"
    )
  )

  ggplot_object <-
    ggplot_object + ggplot2::geom_col(
      mapping = ggplot2::aes(
        x = .data$sample,
        y = .data$fraction,
        fill = .data$junction
      ),
      alpha = I(1 / 3)
    )

  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Sample",
      y = "Splice Junction Fraction",
      fill = "Splice Junction",
      title = "STAR Aligner Splice Junction Fractions per Sample"
    )

  # Reduce the label font size and the legend key size and allow a maximum of 24
  # guide legend rows.
  ggplot_object <-
    ggplot_object + ggplot2::guides(
      colour = ggplot2::guide_legend(
        keywidth = ggplot2::rel(x = 0.8),
        keyheight = ggplot2::rel(x = 0.8),
        nrow = 24L
      )
    )

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
            argument_list$prefix,
            "junction",
            "fraction",
            "sample",
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

  rm(sample_tibble)
}

rm(
  read_group_tibble,
  parse_report,
  graphics_formats,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessioninfo::session_info())

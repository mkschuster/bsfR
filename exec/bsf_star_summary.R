#!/usr/bin/env Rscript
#
# BSF R script to summarise STAR aligner alignment reports.
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
        default = "^star_align_.*_Log\\.final\\.out$",
        dest = "pattern_file",
        help = "STAR alignment report file name pattern [^star_align_.*_Log\\.final\\.out$]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--pattern-sample"),
        default = "^star_align_(.*)_Log\\.final\\.out$",
        dest = "pattern_sample",
        help = "STAR alignment report sample name pattern [^star_align_(.*)_Log\\.final\\.out$]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--prefix"),
        default = "star_summary",
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

suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
suppressPackageStartupMessages(expr = library(package = "tidyverse"))

# Save plots in the following formats.
graphics_formats <- c("pdf" = "pdf", "png" = "png")
# Maximum size for the PNG device in inches.
graphics_maximum_size_png <- 100.0

# Parse STAR aligner log files --------------------------------------------


read_group_tibble <- tibble::tibble(
  "file_name" =
    base::list.files(pattern = argument_list$pattern_file, recursive = TRUE),
  "read_group_name" = base::gsub(
    pattern = argument_list$pattern_sample,
    replacement = "\\1",
    x = base::basename(path = .data$file_name)
  )
)
message(
  "Processing STAR alignment reports for number of read groups: ",
  nrow(x = read_group_tibble)
)

# Empty lines and single column rows render parsing STAR Log.final.out files
# more cumbersome than necessary. Use base::readLines() and
# stringr::str_split_fixed() to split into a matrix.
#
# > stringr::str_split_fixed(string = star_lines, pattern = fixed(pattern = "\t"), n = 2)
# [,1]                                                [,2]
# [1,] "                                 Started job on |" "Apr 01 20:49:14"
# [2,] "                             Started mapping on |" "Apr 01 20:50:12"
# [3,] "                                    Finished on |" "Apr 01 20:51:57"
# [4,] "       Mapping speed, Million of reads per hour |" "43.78"
# [5,] ""                                                  ""
# [6,] "                          Number of input reads |" "1276962"
# [7,] "                      Average input read length |" "302"
# [8,] "                                    UNIQUE READS:" ""
# [9,] "                   Uniquely mapped reads number |" "1128810"
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

variable_names <- c(
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

star_tibble <- tibble::tibble()
for (i in seq_len(length.out = nrow(x = read_group_tibble))) {
  message("  ", read_group_tibble$read_group_name[i])

  star_character <- stringr::str_split_fixed(
    string = base::readLines(con = read_group_tibble$file_name[i]),
    pattern = "\t",
    n = 2L
  )[c(1L:4L, 6L:7L, 9L:22L, 24L:27L, 29L:31L, 33L:34L), 2L]
  names(x = star_character) <- variable_names
  star_tibble <- dplyr::bind_rows(star_tibble, star_character)
  rm(star_character)
}
rm(i)

read_group_tibble <-
  dplyr::bind_cols(read_group_tibble, star_tibble)
rm(star_tibble, variable_names)

# Convert factor to character to integer vectors.
read_group_tibble <- dplyr::mutate(
  .data = read_group_tibble,
  "input_reads" = as.integer(x = .data$input_reads),
  "uniquely_mapped_reads" = as.integer(x = .data$uniquely_mapped_reads),
  "number_splice_total" = as.integer(x = .data$number_splice_total),
  "number_splice_sjdb" = as.integer(x = .data$number_splice_sjdb),
  "number_splice_gtag" = as.integer(x = .data$number_splice_gtag),
  "number_splice_gcag" = as.integer(x = .data$number_splice_gcag),
  "number_splice_atac" = as.integer(x = .data$number_splice_atac),
  "number_splice_non_canonical" = as.integer(x = .data$number_splice_non_canonical),
  "multi_mapped_number" = as.integer(x = .data$multi_mapped_number),
  "chimeric_number" = as.integer(x = .data$chimeric_number)
)

message("Writing read group-level summary table")
readr::write_tsv(
  x = read_group_tibble,
  file = paste(
    paste(argument_list$prefix, "table", "read_group", sep = "_"),
    "tsv",
    sep = "."
  ),
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
    cols = c(.data$unmapped, .data$multi, .data$unique),
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
  argument_list$plot_width + (ceiling(x = nrow(x = read_group_tibble) / 24L) - 1L) * argument_list$plot_width * 0.75
for (graphics_format in graphics_formats) {
  if (graphics_format == "png" &&
      plot_width > graphics_maximum_size_png) {
    message("PNG plot exceeding maximum size: ",
            plot_width,
            " > ",
            graphics_maximum_size_png)
    next
  }
  ggplot2::ggsave(
    filename = paste(
      paste(argument_list$prefix,
            "alignment",
            "read_group",
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
    cols = c(.data$unmapped, .data$multi, .data$unique),
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
  argument_list$plot_width + (ceiling(x = nrow(x = read_group_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
  if (graphics_format == "png" &&
      plot_width > graphics_maximum_size_png) {
    message("PNG plot exceeding maximum size: ",
            plot_width,
            " > ",
            graphics_maximum_size_png)
    next
  }
  ggplot2::ggsave(
    filename = paste(
      paste(argument_list$prefix,
            "mapped",
            "number",
            "read_group",
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
    cols = c(.data$unmapped, .data$multi, .data$unique),
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
  argument_list$plot_width + (ceiling(x = nrow(x = read_group_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
  if (graphics_format == "png" &&
      plot_width > graphics_maximum_size_png) {
    message("PNG plot exceeding maximum size: ",
            plot_width,
            " > ",
            graphics_maximum_size_png)
    next
  }
  ggplot2::ggsave(
    filename = paste(
      paste(
        argument_list$prefix,
        "mapped",
        "fraction",
        "read_group",
        sep = "_"
      ),
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

# Column plot of splice junction numbers per read group -------------------


message("Creating a column plot of splice junction numbers per read group")

ggplot_object <- ggplot2::ggplot(
  data = tidyr::pivot_longer(
    data = dplyr::select(
      .data = read_group_tibble,
      "read_group" = .data$read_group_name,
      "gtag" = .data$number_splice_gtag,
      "gcag" = .data$number_splice_gcag,
      "atac" = .data$number_splice_atac,
      "non_canonical" = .data$number_splice_non_canonical
    ),
    cols = c(.data$non_canonical, .data$atac, .data$gcag, .data$gtag),
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
  argument_list$plot_width + (ceiling(x = nrow(x = read_group_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
  if (graphics_format == "png" &&
      plot_width > graphics_maximum_size_png) {
    message("PNG plot exceeding maximum size: ",
            plot_width,
            " > ",
            graphics_maximum_size_png)
    next
  }
  ggplot2::ggsave(
    filename = paste(
      paste(
        argument_list$prefix,
        "junction",
        "number",
        "read_group",
        sep = "_"
      ),
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
    cols = c(.data$non_canonical, .data$atac, .data$gcag, .data$gtag),
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
  argument_list$plot_width + (ceiling(x = nrow(x = read_group_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
  if (graphics_format == "png" &&
      plot_width > graphics_maximum_size_png) {
    message("PNG plot exceeding maximum size: ",
            plot_width,
            " > ",
            graphics_maximum_size_png)
    next
  }
  ggplot2::ggsave(
    filename = paste(
      paste(
        argument_list$prefix,
        "junction",
        "fraction",
        "read_group",
        sep = "_"
      ),
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

# Integrate read group data on sample level -------------------------------


# Can the read group-level data be integrated on the sample level?

file_path <-
  paste0(argument_list$prefix, "_read_group_to_sample.tsv")
if (file.exists(file_path)) {
  message("Integrating on sample-level ...")
  sample_tibble <-
    dplyr::left_join(
      # Read the read group to sample mapping data frame provided by BSF Python.
      x = readr::read_tsv(
        file = file_path,
        col_names = TRUE,
        col_types = readr::cols(
          sample = readr::col_character(),
          read_group = readr::col_character()
        )
      ),
      # Create a smaller tibble for mapping and junction information and
      # include the read group for merging.
      y = dplyr::select(
        .data = read_group_tibble,
        "read_group" = .data$read_group_name,
        "reads_input" = .data$input_reads,
        "reads_unique" = .data$uniquely_mapped_reads,
        "reads_multi" = .data$multi_mapped_number,
        "junctions_total" = .data$number_splice_total,
        "junctions_sjdb" = .data$number_splice_sjdb,
        "junctions_gtag" = .data$number_splice_gtag,
        "junctions_gcag" = .data$number_splice_gcag,
        "junctions_atac" = .data$number_splice_atac,
        "junctions_non_canonical" = .data$number_splice_non_canonical,
      ),
      by = "read_group"
    )
  sample_tibble <-
    dplyr::group_by(.data = sample_tibble, .data$sample)
  sample_tibble <- dplyr::summarise_at(
    .tbl = sample_tibble,
    .vars = c(
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
    .funs = ~ sum(.)
  )

  # Write the sample-level tibble with mapping and junction information to disk.
  message("Writing sample-level summary table")
  readr::write_tsv(
    x = sample_tibble,
    file = paste(
      paste(argument_list$prefix, "table", "sample", sep = "_"),
      "tsv",
      sep = "."
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
        "sample" = .data$sample,
        "input" = .data$reads_input,
        "multi" = .data$reads_multi,
        "unique" = .data$reads_unique,
        "unmapped" = .data$reads_unmapped
      ),
      cols = c(.data$unmapped, .data$multi, .data$unique),
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
    argument_list$plot_width + (ceiling(x = nrow(x = sample_tibble) / 24L) - 1L) * argument_list$plot_width * 0.33
  for (graphics_format in graphics_formats) {
    if (graphics_format == "png" &&
        plot_width > graphics_maximum_size_png) {
      message("PNG plot exceeding maximum size: ",
              plot_width,
              " > ",
              graphics_maximum_size_png)
      next
    }
    ggplot2::ggsave(
      filename = paste(
        paste(argument_list$prefix,
              "alignment",
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

  ggplot_object <- ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = dplyr::select(
        .data = sample_tibble,
        "sample" = .data$sample,
        "unmapped" = .data$reads_unmapped,
        "multi" = .data$reads_multi,
        "unique" = .data$reads_unique
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
    argument_list$plot_width + (ceiling(x = nrow(x = sample_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25
  for (graphics_format in graphics_formats) {
    if (graphics_format == "png" &&
        plot_width > graphics_maximum_size_png) {
      message("PNG plot exceeding maximum size: ",
              plot_width,
              " > ",
              graphics_maximum_size_png)
      next
    }
    ggplot2::ggsave(
      filename = paste(
        paste(argument_list$prefix,
              "mapped",
              "number",
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
    argument_list$plot_width + (ceiling(x = nrow(x = sample_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25
  for (graphics_format in graphics_formats) {
    if (graphics_format == "png" &&
        plot_width > graphics_maximum_size_png) {
      message("PNG plot exceeding maximum size: ",
              plot_width,
              " > ",
              graphics_maximum_size_png)
      next
    }
    ggplot2::ggsave(
      filename = paste(
        paste(argument_list$prefix,
              "mapped",
              "fraction",
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

  # Column plot of splice junction numbers per sample ---------------------


  message("Creating a column plot of splice junction numbers per sample")

  ggplot_object <- ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = dplyr::select(
        .data = sample_tibble,
        "sample" = .data$sample,
        "gtag" = .data$junctions_gtag,
        "gcag" = .data$junctions_gcag,
        "atac" = .data$junctions_atac,
        "non_canonical" = .data$junctions_non_canonical
      ),
      cols = c(.data$non_canonical, .data$atac, .data$gcag, .data$gtag),
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
    argument_list$plot_width + (ceiling(x = nrow(x = sample_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25
  for (graphics_format in graphics_formats) {
    if (graphics_format == "png" &&
        plot_width > graphics_maximum_size_png) {
      message("PNG plot exceeding maximum size: ",
              plot_width,
              " > ",
              graphics_maximum_size_png)
      next
    }
    ggplot2::ggsave(
      filename = paste(
        paste(argument_list$prefix,
              "junction",
              "number",
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
      cols = c(.data$non_canonical, .data$atac, .data$gcag, .data$gtag),
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
    argument_list$plot_width + (ceiling(x = nrow(x = sample_tibble) / 24L) - 1L) * argument_list$plot_width * 0.25
  for (graphics_format in graphics_formats) {
    if (graphics_format == "png" &&
        plot_width > graphics_maximum_size_png) {
      message("PNG plot exceeding maximum size: ",
              plot_width,
              " > ",
              graphics_maximum_size_png)
      next
    }
    ggplot2::ggsave(
      filename = paste(
        paste(
          argument_list$prefix,
          "junction",
          "fraction",
          "sample",
          sep = "_"
        ),
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

  rm(sample_tibble)
}
rm(file_path)

rm(read_group_tibble,
   graphics_maximum_size_png,
   graphics_formats,
   argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessioninfo::session_info())

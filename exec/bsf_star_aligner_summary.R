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
        default = "^star_aligner_align_.*_Log\\.final\\.out$",
        dest = "pattern_file",
        help = "STAR alignment report file name pattern",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--pattern-sample"),
        default = "^star_aligner_align_(.*)_Log\\.final\\.out$",
        dest = "pattern_sample",
        help = "STAR alignment report sample name pattern",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--prefix"),
        default = "star_aligner_summary",
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

# Save plots in the following formats.
graphics_formats <- c("pdf" = "pdf", "png" = "png")

# Parse STAR aligner log files --------------------------------------------


summary_frame <- data.frame(
  file_name =
    base::list.files(pattern = argument_list$pattern_file, recursive = TRUE),
  stringsAsFactors = FALSE
)
message("Processing STAR alignment reports for number of read groups: ",
        nrow(x = summary_frame))

star_frame <- NULL

for (i in seq_len(length.out = nrow(x = summary_frame))) {
  summary_frame[i, "read_group_name"] <-
    gsub(
      pattern = argument_list$pattern_sample,
      replacement = "\\1",
      x = base::basename(path = summary_frame[i, "file_name", drop = TRUE])
    )
  message("  ", summary_frame[i, "read_group_name", drop = TRUE])

  # This is the layout of a STAR alignment summary file.
  #
  # Started job on |       May 27 15:24:19
  # Started mapping on |       May 27 15:32:57
  # Finished on |       May 27 15:40:35
  # Mapping speed, Million of reads per hour |       134.93
  #
  # Number of input reads |       17166324
  # Average input read length |       50
  # UNIQUE READS:
  #   Uniquely mapped reads number |       13192124
  # Uniquely mapped reads % |       76.85%
  #   Average mapped length |       49.86
  # Number of splices: Total |       1602861
  # Number of splices: Annotated (sjdb) |       1584501
  # Number of splices: GT/AG |       1591365
  # Number of splices: GC/AG |       9392
  # Number of splices: AT/AC |       1173
  # Number of splices: Non-canonical |       931
  # Mismatch rate per base, % |       0.19%
  #   Deletion rate per base |       0.00%
  # Deletion average length |       1.44
  # Insertion rate per base |       0.00%
  # Insertion average length |       1.28
  # MULTI-MAPPING READS:
  #   Number of reads mapped to multiple loci |       3849655
  # % of reads mapped to multiple loci |       22.43%
  #   Number of reads mapped to too many loci |       17781
  # % of reads mapped to too many loci |       0.10%
  #   UNMAPPED READS:
  #   % of reads unmapped: too many mismatches |       0.00%
  #   % of reads unmapped: too short |       0.44%
  #   % of reads unmapped: other |       0.18%
  #   CHIMERIC READS:
  #   Number of chimeric reads |       0
  # % of chimeric reads |       0.00%

  temporary_frame <- as.data.frame(t(
    x = read.table(
      file = summary_frame[i, "file_name", drop = TRUE],
      sep = "|",
      fill = TRUE,
      strip.white = TRUE,
      stringsAsFactors = FALSE
    )
  ))
  if (is.null(x = star_frame)) {
    # The first temporary frame has the column headers as the first line.
    # Unfortunately they are too complex to be useful.
    # Names are assigned manually below.
    star_frame <- temporary_frame[2L, , drop = FALSE]
  } else {
    star_frame <-
      rbind(star_frame, temporary_frame[2L, , drop = FALSE], stringsAsFactors = FALSE)
  }
  rm(temporary_frame)
}
names(x = star_frame) <-
  c(
    "started_job",
    "started_mapping",
    "finished",
    "mapping_speed",
    "input_reads",
    "average_length",
    "unique_reads",
    # placeholder V7
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
    "average_deletion_length",
    "insertion_rate",
    "average_insertion_length",
    "multi_mapping_reads",
    # placeholder V22
    "multi_mapped_number",
    "multi_mapped_percentage",
    "multi_unmapped_percentage",
    "unmapped_reads",
    # placeholder V26
    "unmapped_mismatched_percentage",
    "unmapped_short_percentage",
    "unmapped_other_percentage",
    "chimeric_reads",
    # placeholder V30
    "chimeric_number",
    "chimeric_percentage"
  )
# Remove placeholder columns of no value.
star_frame <-
  star_frame[, c(1L:6L, 8L:21L, 23L:25L, 27L:29L, 31L:32L), drop = FALSE]
summary_frame <-
  cbind(summary_frame, star_frame, stringsAsFactors = FALSE)
rm(i, star_frame)

# Convert factor to character to integer vectors.
summary_frame$input_reads <-
  as.integer(x = as.character(x = summary_frame$input_reads))
summary_frame$uniquely_mapped_reads <-
  as.integer(x = as.character(x = summary_frame$uniquely_mapped_reads))
summary_frame$number_splice_total <-
  as.integer(x = as.character(x = summary_frame$number_splice_total))
summary_frame$number_splice_sjdb <-
  as.integer(x = as.character(x = summary_frame$number_splice_sjdb))
summary_frame$number_splice_gtag <-
  as.integer(x = as.character(x = summary_frame$number_splice_gtag))
summary_frame$number_splice_gcag <-
  as.integer(x = as.character(x = summary_frame$number_splice_gcag))
summary_frame$number_splice_atac <-
  as.integer(x = as.character(x = summary_frame$number_splice_atac))
summary_frame$number_splice_non_canonical <-
  as.integer(x = as.character(x = summary_frame$number_splice_non_canonical))
summary_frame$multi_mapped_number <-
  as.integer(x = as.character(x = summary_frame$multi_mapped_number))
summary_frame$chimeric_number <-
  as.integer(x = as.character(x = summary_frame$chimeric_number))

message("Writing read group-level summary table")
write.table(
  x = summary_frame,
  file = paste0(argument_list$prefix, "_table_read_group.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

# Scatter plot of read number versus alignment rate per read group --------


message("Creating a scatter plot of read number versus alignment rate per read group")

plotting_frame <-
  tidyr::pivot_longer(
    data = tibble::tibble(
      "read_group" = summary_frame$read_group_name,
      "input" = summary_frame$input_reads,
      "unique" = summary_frame$uniquely_mapped_reads,
      "multi" = summary_frame$multi_mapped_number,
      "unmapped" = summary_frame$input_reads - summary_frame$uniquely_mapped_reads - summary_frame$multi_mapped_number,
      stringsAsFactors = FALSE
    ),
    cols = c(.data$unmapped, .data$multi, .data$unique),
    names_to = "status",
    values_to = "mapped"
  )

ggplot_object <- ggplot2::ggplot(data = plotting_frame)
ggplot_object <-
  ggplot_object + ggplot2::geom_point(
    mapping = ggplot2::aes(
      x = .data$mapped,
      y = .data$mapped / .data$input,
      colour = .data$read_group,
      shape = .data$status
    )
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
  argument_list$plot_width + (ceiling(x = nrow(x = summary_frame) / 24L) - 0L) * argument_list$plot_width * 1.0
for (graphics_format in graphics_formats) {
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
rm(graphics_format, plot_width, ggplot_object, plotting_frame)

# Column plot of read numbers per read group ------------------------------


message("Creating a column plot of read numbers per read group")

plotting_frame <-
  tidyr::pivot_longer(
    data = tibble::tibble(
      "read_group" = summary_frame$read_group_name,
      "unique" = summary_frame$uniquely_mapped_reads,
      "multi" = summary_frame$multi_mapped_number,
      "unmapped" =
        summary_frame$input_reads - summary_frame$uniquely_mapped_reads - summary_frame$multi_mapped_number
    ),
    cols = c(.data$unmapped, .data$multi, .data$unique),
    names_to = "status",
    values_to = "number"
  )

ggplot_object <- ggplot2::ggplot(data = plotting_frame)
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
      angle = 90,
      hjust = 0,
      size = ggplot2::rel(x = 0.7)
    ),
    legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7))
  )
# Scale the plot width with the number of read groups, by adding a quarter of
# the original width for each 24 read groups.
plot_width <-
  argument_list$plot_width + (ceiling(x = nrow(x = summary_frame) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
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
rm(graphics_format, plot_width, ggplot_object, plotting_frame)

# Column plot of read fractions per read group ----------------------------


message("Creating a column plot of read fractions per read group")
plotting_frame <-
  tidyr::pivot_longer(
    data = tibble::tibble(
      "read_group" = summary_frame$read_group_name,
      "unique" = summary_frame$uniquely_mapped_reads / summary_frame$input_reads,
      "multi" = summary_frame$multi_mapped_number / summary_frame$input_reads,
      "unmapped" = (
        summary_frame$input_reads - summary_frame$uniquely_mapped_reads - summary_frame$multi_mapped_number
      ) / summary_frame$input_reads,
      stringsAsFactors = FALSE
    ),
    cols = c(.data$unmapped, .data$multi, .data$unique),
    names_to = "status",
    values_to = "fraction"
  )

ggplot_object <- ggplot2::ggplot(data = plotting_frame)
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
      angle = 90,
      hjust = 0,
      size = ggplot2::rel(x = 0.7)
    ),
    legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7))
  )
# Scale the plot width with the number of read groups, by adding a quarter of
# the original width for each 24 read groups.
plot_width <-
  argument_list$plot_width + (ceiling(x = nrow(x = summary_frame) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
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
rm(graphics_format, plot_width, ggplot_object, plotting_frame)

# Column plot of splice junction numbers per read group -------------------


message("Creating a column plot of splice junction numbers per read group")

plotting_frame <-
  tidyr::pivot_longer(
    data = tibble::tibble(
      "read_group" = summary_frame$read_group_name,
      "gtag" = summary_frame$number_splice_gtag,
      "gcag" = summary_frame$number_splice_gcag,
      "atac" = summary_frame$number_splice_atac,
      "non_canonical" = summary_frame$number_splice_non_canonical,
      stringsAsFactors = FALSE
    ),
    cols = c(.data$non_canonical, .data$atac, .data$gcag, .data$gtag),
    names_to = "splice_junction",
    values_to = "number"
  )

ggplot_object <- ggplot2::ggplot(data = plotting_frame)
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
      angle = 90,
      hjust = 0,
      size = ggplot2::rel(x = 0.7)
    ),
    legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7))
  )
# Scale the plot width with the number of read groups, by adding a quarter of
# the original width for each 24 read groups.
plot_width <-
  argument_list$plot_width + (ceiling(x = nrow(x = summary_frame) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
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
rm(graphics_format, plot_width, ggplot_object, plotting_frame)

# Column plot of splice junction fractions per read group -----------------


message("Creating a column plot of splice junction fractions per read group")

plotting_frame <-
  tidyr::pivot_longer(
    data = tibble::tibble(
      "read_group" = summary_frame$read_group_name,
      "gtag" = summary_frame$number_splice_gtag / summary_frame$number_splice_total,
      "gcag" = summary_frame$number_splice_gcag / summary_frame$number_splice_total,
      "atac" = summary_frame$number_splice_atac / summary_frame$number_splice_total,
      "non_canonical" = summary_frame$number_splice_non_canonical / summary_frame$number_splice_total
    ),
    cols = c(.data$non_canonical, .data$atac, .data$gcag, .data$gtag),
    names_to = "junction",
    values_to = "fraction"
  )

ggplot_object <- ggplot2::ggplot(data = plotting_frame)
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
      angle = 90,
      hjust = 0,
      size = ggplot2::rel(x = 0.7)
    ),
    legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7))
  )
# Scale the plot width with the number of read groups, by adding a quarter of
# the original width for each 24 samples.
plot_width <-
  argument_list$plot_width + (ceiling(x = nrow(x = summary_frame) / 24L) - 1L) * argument_list$plot_width * 0.25
for (graphics_format in graphics_formats) {
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
rm(graphics_format, plot_width, ggplot_object, plotting_frame)

# Integrate read group data on sample level -------------------------------


# Can the read group-level data be integrated on the sample level?

file_path <-
  paste0(argument_list$prefix, "_read_group_to_sample.tsv")
if (file.exists(file_path)) {
  message("Integrating on sample-level ...")
  merged_frame <-
    merge.data.frame(
      # Read the read group to sample mapping data frame provided by BSF Python.
      x = read.table(
        file = file_path,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
      ),
      # Create a smaller data frame for mapping and junction information and include the read group for merging.
      y = data.frame(
        "read_group" = summary_frame$read_group_name,
        "reads_input" = summary_frame$input_reads,
        "reads_unique" = summary_frame$uniquely_mapped_reads,
        "reads_multi" = summary_frame$multi_mapped_number,
        # The reads_mapped and reads_unmapped columns are just linear combinations.
        # Add them in for plotting later.
        # "reads_mapped" = summary_frame$uniquely_mapped_reads + summary_frame$multi_mapped_number,
        # "reads_unmapped" =
        #   summary_frame$input_reads - summary_frame$uniquely_mapped_reads - summary_frame$multi_mapped_number,
        "junctions_total" = summary_frame$number_splice_total,
        "junctions_sjdb" = summary_frame$number_splice_sjdb,
        "junctions_gtag" = summary_frame$number_splice_gtag,
        "junctions_gcag" = summary_frame$number_splice_gcag,
        "junctions_atac" = summary_frame$number_splice_atac,
        "junctions_non_canonical" = summary_frame$number_splice_non_canonical,
        stringsAsFactors = FALSE
      ),
      by = "read_group"
    )
  # Re-structure the plotting frame for aggregation.
  # Keep only the numeric vectors and set the read group names
  # as the row names. Then aggregate by a list of sample names.
  # seq.int(from = 3L, to = ncol(x = merged_frame))
  column_names <- names(x = merged_frame)
  plotting_frame <-
    merged_frame[, column_names != "read_group" &
                   column_names != "sample", drop = FALSE]
  rm(column_names)
  row.names(x = plotting_frame) <- merged_frame$read_group
  aggregate_frame <-
    aggregate.data.frame(x = plotting_frame,
                         by = list(merged_frame$sample),
                         FUN = sum)
  # Rename the "Group.1" column as the result of the aggreagtion by sample list into "sample".
  names(x = aggregate_frame)[names(x = aggregate_frame) == "Group.1"] <-
    "sample"
  # Write the aggregated frame with mapping and junction information to disk.
  message("Writing sample-level summary table")
  write.table(
    x = aggregate_frame,
    file = paste0(argument_list$prefix, "_table_sample.tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  # For plotting, add in mapped und unmapped reads.
  aggregate_frame$reads_mapped <-
    aggregate_frame$reads_unique + aggregate_frame$reads_multi
  aggregate_frame$reads_unmapped <-
    aggregate_frame$reads_input - aggregate_frame$reads_mapped

  # Scatter plot of read number versus alignment rate per sample ----------


  message("Creating a scatter plot of read number versus alignment rate per sample")
  plotting_frame <-
    tidyr::pivot_longer(
      data = tibble::tibble(
        "sample" = aggregate_frame$sample,
        "input" = aggregate_frame$reads_input,
        "multi" = aggregate_frame$reads_multi,
        "unique" = aggregate_frame$reads_unique,
        "unmapped" = aggregate_frame$reads_unmapped,
        stringsAsFactors = FALSE
      ),
      cols = c(.data$unmapped, .data$multi, .data$unique),
      names_to = "status",
      values_to = "mapped"
    )

  ggplot_object <- ggplot2::ggplot(data = plotting_frame)
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$mapped,
        y = .data$mapped / .data$input,
        colour = .data$sample,
        shape = .data$status
      )
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
    argument_list$plot_width + (ceiling(x = nrow(x = aggregate_frame) / 24L) - 1L) * argument_list$plot_width * 0.25
  for (graphics_format in graphics_formats) {
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
  rm(graphics_format, plot_width, ggplot_object, plotting_frame)

  # Column plot of read numbers per sample --------------------------------


  message("Creating a column plot of read numbers per sample")
  plotting_frame <-
    tidyr::pivot_longer(
      data = tibble::tibble(
        "sample" = aggregate_frame$sample,
        "unmapped" = aggregate_frame$reads_unmapped,
        "multi" = aggregate_frame$reads_multi,
        "unique" = aggregate_frame$reads_unique,
        stringsAsFactors = FALSE
      ),
      cols = c(.data$unmapped, .data$multi, .data$unique),
      names_to = "status",
      values_to = "number"
    )

  ggplot_object <- ggplot2::ggplot(data = plotting_frame)
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
        angle = 90,
        hjust = 0,
        size = ggplot2::rel(x = 0.7)
      ),
      legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7))
    )
  # Scale the plot width with the number of samples, by adding a quarter of
  # the original width for each 24 samples.
  plot_width <-
    argument_list$plot_width + (ceiling(x = nrow(x = aggregate_frame) / 24L) - 1L) * argument_list$plot_width * 0.25
  for (graphics_format in graphics_formats) {
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
  rm(graphics_format, plot_width, ggplot_object, plotting_frame)

  # Column plot of read fractions per sample ------------------------------


  message("Creating a column plot of read fractions per sample")

  plotting_frame <-
    tidyr::pivot_longer(
      data = tibble::tibble(
        "sample" = aggregate_frame$sample,
        "unmapped" = aggregate_frame$reads_unmapped / aggregate_frame$reads_input,
        "multi" = aggregate_frame$reads_multi / aggregate_frame$reads_input,
        "unique" = aggregate_frame$reads_unique / aggregate_frame$reads_input
      ),
      cols = c(.data$unmapped, .data$multi, .data$unique),
      names_to = "status",
      values_to = "fraction"
    )

  ggplot_object <- ggplot2::ggplot(data = plotting_frame)
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
        angle = 90,
        hjust = 0,
        size = ggplot2::rel(x = 0.7)
      ),
      legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7))
    )
  # Scale the plot width with the number of samples, by adding a quarter of
  # the original width for each 24 samples.
  plot_width <-
    argument_list$plot_width + (ceiling(x = nrow(x = aggregate_frame) / 24L) - 1L) * argument_list$plot_width * 0.25
  for (graphics_format in graphics_formats) {
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
  rm(graphics_format, plot_width, ggplot_object, plotting_frame)

  # Column plot of splice junction numbers per sample ---------------------


  message("Creating a column plot of splice junction numbers per sample")

  plotting_frame <-
    tidyr::pivot_longer(
      data = tibble::tibble(
        "sample" = aggregate_frame$sample,
        "gtag" = aggregate_frame$junctions_gtag,
        "gcag" = aggregate_frame$junctions_gcag,
        "atac" = aggregate_frame$junctions_atac,
        "non_canonical" = aggregate_frame$junctions_non_canonical
      ),
      cols = c(.data$non_canonical, .data$atac, .data$gcag, .data$gtag),
      names_to = "junction",
      values_to = "number"
    )

  ggplot_object <- ggplot2::ggplot(data = plotting_frame)
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
        angle = 90,
        hjust = 0,
        size = ggplot2::rel(x = 0.7)
      ),
      legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7))
    )
  # Scale the plot width with the number of samples, by adding a quarter of
  # the original width for each 24 samples.
  plot_width <-
    argument_list$plot_width + (ceiling(x = nrow(x = aggregate_frame) / 24L) - 1L) * argument_list$plot_width * 0.25
  for (graphics_format in graphics_formats) {
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
  rm(graphics_format, plot_width, ggplot_object, plotting_frame)

  # Column plot of splice junction fractions per sample -------------------


  message("Creating a column plot of splice junction fractions per sample")

  plotting_frame <-
    tidyr::pivot_longer(
      data = tibble::tibble(
        "sample" = aggregate_frame$sample,
        "gtag" = aggregate_frame$junctions_gtag / aggregate_frame$junctions_total,
        "gcag" = aggregate_frame$junctions_gcag / aggregate_frame$junctions_total,
        "atac" = aggregate_frame$junctions_atac / aggregate_frame$junctions_total,
        "non_canonical" = aggregate_frame$junctions_non_canonical / aggregate_frame$junctions_total
      ),
      cols = c(.data$non_canonical, .data$atac, .data$gcag, .data$gtag),
      names_to = "junction",
      values_to = "fraction"
    )

  ggplot_object <- ggplot2::ggplot(data = plotting_frame)
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
        angle = 90,
        hjust = 0,
        size = ggplot2::rel(x = 0.7)
      ),
      legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7))
    )
  # Scale the plot width with the number of samples, by adding a quarter of
  # the original width for each 24 samples.
  plot_width <-
    argument_list$plot_width + (ceiling(x = nrow(x = aggregate_frame) / 24L) - 1L) * argument_list$plot_width * 0.25
  for (graphics_format in graphics_formats) {
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
  rm(graphics_format, plot_width, ggplot_object, plotting_frame)

  rm(aggregate_frame, merged_frame)
}
rm(file_path)

rm(summary_frame,
   graphics_formats,
   argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

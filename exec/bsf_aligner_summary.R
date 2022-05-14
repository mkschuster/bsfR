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
        opt_str = "--plot-width",
        default = 7.0,
        dest = "plot_width",
        help = "Plot width in inches [7.0]",
        type = "numeric"
      ),
      optparse::make_option(
        opt_str = "--plot-height",
        default = 7.0,
        dest = "plot_height",
        help = "Plot height in inches [7.0]",
        type = "numeric"
      )
    )
  ))

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "ggplot2"))
suppressPackageStartupMessages(expr = library(package = "tidyr"))


# Save plots in the following formats.
graphics_formats <- c("pdf" = "pdf", "png" = "png")
# Maximum size for the PNG device in inches.
graphics_maximum_size_png <- 200.0

# Assign a file prefix.
prefix_summary <- paste(argument_list$prefix, "summary", sep = "_")
# Assign a Picard Alignment Summary Metrics (pasm) prefix.
prefix_pasm <- "pasm"
# Assign a Picard Duplication Summary Metrics (pdsm) prefix.
prefix_pdsm <- "pdsm"

# Picard Alignment Summary Metrics ----------------------------------------


# Process Picard Alignment Summary Metrics reports.

message("Processing Picard Alignment Summary Metrics report for sample:")
combined_metrics_sample <- NULL
combined_metrics_read_group <- NULL

file_names <-
  base::list.files(
    pattern = paste0(
      "^",
      argument_list$prefix,
      "_sample_.*_alignment_summary_metrics.tsv$"
    ),
    recursive = TRUE
  )
for (file_name in file_names) {
  sample_name <-
    base::gsub(
      pattern = paste0(
        "^",
        argument_list$prefix,
        "_sample_(.*?)_alignment_summary_metrics.tsv$"
      ),
      replacement = "\\1",
      x = base::basename(path = file_name)
    )
  message("  ", sample_name)
  # Since the Picard AlignmentSummaryMetrics uses a hash character (#) in the
  # read group component to separate platform unit and sample name, the Picard
  # reports need special parsing.
  # Find the ## METRICS CLASS line and parse without allowing further comments.
  metrics_lines <- readLines(con = file_name)
  metrics_line <-
    which(x = grepl(pattern = "## METRICS CLASS", x = metrics_lines))
  picard_metrics_total <-
    utils::read.table(
      file = file_name,
      header = TRUE,
      sep = "\t",
      skip = metrics_line[1L],
      fill = TRUE,
      comment.char = "",
      stringsAsFactors = FALSE
    )
  rm(metrics_line, metrics_lines)
  # To support numeric sample names the utils::read.table(stringsAsFactors = FALSE) is turned off.
  # Convert SAMPLE, LIBRARY and READ_GROUP into character vectors.
  picard_metrics_total$SAMPLE <-
    as.character(x = picard_metrics_total$SAMPLE)
  picard_metrics_total$LIBRARY <-
    as.character(x = picard_metrics_total$LIBRARY)
  picard_metrics_total$READ_GROUP <-
    as.character(x = picard_metrics_total$READ_GROUP)

  # The Picard Alignment Metrics report has changed format through versions.
  # Columns PF_READS_IMPROPER_PAIRS and PCT_PF_READS_IMPROPER_PAIRS were added at a later stage.
  if (!"PF_READS_IMPROPER_PAIRS" %in% names(x = picard_metrics_total)) {
    picard_metrics_total$PF_READS_IMPROPER_PAIRS <- 0L
  }

  if (!"PCT_PF_READS_IMPROPER_PAIRS" %in% names(x = picard_metrics_total)) {
    picard_metrics_total$PCT_PF_READS_IMPROPER_PAIRS <- 0.0
  }

  # Select only rows showing the SAMPLE summary, i.e. showing SAMPLE, but no LIBRARY and READ_GROUP information.
  picard_metrics_sample <-
    picard_metrics_total[(!is.na(x = picard_metrics_total$SAMPLE)) &
                           (picard_metrics_total$SAMPLE != "") &
                           (picard_metrics_total$LIBRARY == "") &
                           (picard_metrics_total$READ_GROUP == ""),]
  combined_metrics_sample <-
    base::rbind(combined_metrics_sample, picard_metrics_sample)
  rm(picard_metrics_sample)

  # Select only rows showing READ_GROUP summary, i.e. showing READ_GROUP information.
  picard_metrics_read_group <-
    picard_metrics_total[(picard_metrics_total$READ_GROUP != ""),]
  combined_metrics_read_group <-
    base::rbind(combined_metrics_read_group, picard_metrics_read_group)
  rm(picard_metrics_read_group)

  rm(sample_name, picard_metrics_total)
}
rm(file_name, file_names)

if (!is.null(x = combined_metrics_sample)) {
  # Order the data frame by SAMPLE.
  combined_metrics_sample <-
    combined_metrics_sample[order(combined_metrics_sample$SAMPLE),]
  # Manually convert CATEGORY and SAMPLE columns into factors, which are handy for plotting.
  combined_metrics_sample$CATEGORY <-
    as.factor(x = combined_metrics_sample$CATEGORY)
  combined_metrics_sample$SAMPLE <-
    as.factor(x = combined_metrics_sample$SAMPLE)
  # Add an additional LABEL factor column defined as a concatenation of SAMPLE and CATEGORY.
  combined_metrics_sample$LABEL <-
    as.factor(x = paste(
      combined_metrics_sample$SAMPLE,
      combined_metrics_sample$CATEGORY,
      sep =
        "_"
    ))
  utils::write.table(
    x = combined_metrics_sample,
    file = paste(prefix_summary, prefix_pasm, "metrics_sample.tsv", sep = "_"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  # Order the data frame by READ_GROUP
  combined_metrics_read_group <-
    combined_metrics_read_group[order(combined_metrics_read_group$READ_GROUP),]
  # Manually convert CATEGORY and READ_GROUP columns into factors, which are handy for plotting.
  combined_metrics_read_group$CATEGORY <-
    as.factor(x = combined_metrics_read_group$CATEGORY)
  combined_metrics_read_group$READ_GROUP <-
    as.factor(x = combined_metrics_read_group$READ_GROUP)
  # Add an additional LABEL factor column defined as a concatenation of READ_GROUP and CATEGORY.
  combined_metrics_read_group$LABEL <-
    as.factor(
      x = paste(
        combined_metrics_read_group$READ_GROUP,
        combined_metrics_read_group$CATEGORY,
        sep =
          "_"
      )
    )
  utils::write.table(
    x = combined_metrics_read_group,
    file = paste(prefix_summary, prefix_pasm, "metrics_read_group.tsv", sep = "_"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  # Plot the absolute number versus the fraction per sample ---------------


  message("Plotting the absolute number versus the fraction per sample")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_sample[, c("CATEGORY", "SAMPLE", "PF_READS_ALIGNED", "PF_READS"), drop = FALSE])
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$PF_READS_ALIGNED,
        y = .data$PF_READS_ALIGNED / .data$PF_READS,
        colour = .data$SAMPLE,
        shape = .data$CATEGORY
      ),
      alpha = I(1 / 3)
    )
  ggplot_object <- ggplot_object + ggplot2::labs(
    x = "Reads Number",
    y = "Reads Fraction",
    colour = "Sample",
    shape = "Category",
    title = "Alignment Summary per Sample"
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
  # Adjust the plot width according to batches of 24 samples or read groups.
  plot_width <-
    argument_list$plot_width + (ceiling(x = nlevels(x = combined_metrics_sample$SAMPLE) / 24L) - 1L) * argument_list$plot_width * 0.33
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
        paste(prefix_summary,
              prefix_pasm,
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

  # Plot the absolute number versus the fraction per read group -----------


  message("Plotting the absolute number versus the fraction per read group")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_read_group[, c("CATEGORY", "READ_GROUP", "PF_READS_ALIGNED", "PF_READS"), drop = FALSE])
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$PF_READS_ALIGNED,
        y = .data$PF_READS_ALIGNED / .data$PF_READS,
        colour = .data$READ_GROUP,
        shape = .data$CATEGORY
      ),
      alpha = I(1 / 3)
    )
  ggplot_object <- ggplot_object + ggplot2::labs(
    x = "Reads Number",
    y = "Reads Fraction",
    colour = "Read Group",
    shape = "Category",
    title = "Alignment Summary per Read Group"
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
  # Adjust the plot width according to batches of 24 samples or read groups.
  plot_width <-
    argument_list$plot_width + (ceiling(x = nlevels(x = combined_metrics_read_group$READ_GROUP) / 24L) - 1L) * argument_list$plot_width * 0.75
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
          prefix_summary,
          prefix_pasm,
          "alignment",
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

  # Adjust the plot width according to batches of 24 samples or read groups.
  plot_width_sample <-
    argument_list$plot_width + (ceiling(x = nlevels(x = combined_metrics_sample$SAMPLE) / 24L) - 1L) * argument_list$plot_width * 0.25
  plot_width_read_group <-
    argument_list$plot_width + (ceiling(x = nlevels(x = combined_metrics_read_group$READ_GROUP) / 24L) - 1L) * argument_list$plot_width * 0.25

  # Plot the absolute number of aligned pass-filter reads per sample ------


  message("Plotting the absolute number of aligned pass-filter reads per sample")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_sample[, c("CATEGORY", "SAMPLE", "PF_READS_ALIGNED"), drop = FALSE])
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$SAMPLE,
        y = .data$PF_READS_ALIGNED,
        colour = .data$CATEGORY
      ),
      alpha = I(1 / 3)
    )
  ggplot_object <- ggplot_object + ggplot2::labs(
    x = "Sample",
    y = "Reads Number",
    colour = "Category",
    title = "Aligned Pass-Filter Reads per Sample"
  )
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = ggplot2::rel(x = 0.8),
        hjust = 0.0,
        vjust = 0.5,
        angle = 90.0
      )
    )
  for (graphics_format in graphics_formats) {
    if (graphics_format == "png" &&
        plot_width_sample > graphics_maximum_size_png) {
      message(
        "PNG plot exceeding maximum size: ",
        plot_width_sample,
        " > ",
        graphics_maximum_size_png
      )
      next
    }
    ggplot2::ggsave(
      filename = paste(
        paste(prefix_summary,
              prefix_pasm,
              "absolute",
              "sample",
              sep = "_"),
        graphics_format,
        sep = "."
      ),
      plot = ggplot_object,
      width = plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot the absolute number of aligned pass-filter reads per read group ----


  message("Plotting the absolute number of aligned pass-filter reads per read group")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_read_group[, c("CATEGORY", "READ_GROUP", "PF_READS_ALIGNED"), drop = FALSE])
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$READ_GROUP,
        y = .data$PF_READS_ALIGNED,
        colour = .data$CATEGORY
      ),
      alpha = I(1 / 3)
    )
  ggplot_object <- ggplot_object + ggplot2::labs(
    x = "Read Group",
    y = "Reads Number",
    colour = "Category",
    title = "Aligned Pass-Filter Reads per Read Group"
  )
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = ggplot2::rel(x = 0.8),
        hjust = 0.0,
        vjust = 0.5,
        angle = 90.0
      )
    )
  for (graphics_format in graphics_formats) {
    if (graphics_format == "png" &&
        plot_width_read_group > graphics_maximum_size_png) {
      message(
        "PNG plot exceeding maximum size: ",
        plot_width_read_group,
        " > ",
        graphics_maximum_size_png
      )
      next
    }
    ggplot2::ggsave(
      filename = paste(
        paste(
          prefix_summary,
          prefix_pasm,
          "absolute",
          "read_group",
          sep = "_"
        ),
        graphics_format,
        sep = "."
      ),
      plot = ggplot_object,
      width = plot_width_read_group,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot the percentage of aligned pass-filter reads per sample -----------


  message("Plotting the percentage of aligned pass-filter reads per sample")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_sample[, c("CATEGORY", "SAMPLE", "PCT_PF_READS_ALIGNED"), drop = FALSE])
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$SAMPLE,
        y = .data$PCT_PF_READS_ALIGNED,
        colour = .data$CATEGORY
      ),
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Sample",
      y = "Reads Fraction",
      colour = "Category",
      title = "Aligned Pass-Filter Reads per Sample"
    )
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = ggplot2::rel(x = 0.8),
        hjust = 0.0,
        vjust = 0.5,
        angle = 90.0
      )
    )
  for (graphics_format in graphics_formats) {
    if (graphics_format == "png" &&
        plot_width_sample > graphics_maximum_size_png) {
      message(
        "PNG plot exceeding maximum size: ",
        plot_width_sample,
        " > ",
        graphics_maximum_size_png
      )
      next
    }
    ggplot2::ggsave(
      filename = paste(
        paste(prefix_summary,
              prefix_pasm,
              "percentage",
              "sample",
              sep = "_"),
        graphics_format,
        sep = "."
      ),
      plot = ggplot_object,
      width = plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot the percentage of aligned pass-filter reads per read group -------


  message("Plotting the percentage of aligned pass-filter reads per read group")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_read_group[, c("CATEGORY", "READ_GROUP", "PCT_PF_READS_ALIGNED"), drop = FALSE])
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$READ_GROUP,
        y = .data$PCT_PF_READS_ALIGNED,
        colour = .data$CATEGORY
      ),
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Read Group",
      y = "Reads Fraction",
      colour = "Category",
      title = "Aligned Pass-Filter Reads per Read Group"
    )
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = ggplot2::rel(x = 0.8),
        hjust = 0.0,
        vjust = 0.5,
        angle = 90.0
      )
    )
  for (graphics_format in graphics_formats) {
    if (graphics_format == "png" &&
        plot_width_read_group > graphics_maximum_size_png) {
      message(
        "PNG plot exceeding maximum size: ",
        plot_width_read_group,
        " > ",
        graphics_maximum_size_png
      )
      next
    }
    ggplot2::ggsave(
      filename = paste(
        paste(
          prefix_summary,
          prefix_pasm,
          "percentage",
          "read_group",
          sep = "_"
        ),
        graphics_format,
        sep = "."
      ),
      plot = ggplot_object,
      width = plot_width_read_group,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot the strand balance of aligned pass-filter reads per sample -------


  message("Plotting the strand balance of aligned pass-filter reads per sample")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_sample[, c("CATEGORY", "SAMPLE", "STRAND_BALANCE"), drop = FALSE])
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$SAMPLE,
        y = .data$STRAND_BALANCE,
        colour = .data$CATEGORY
      ),
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Sample",
      y = "Reads Fraction",
      colour = "Category",
      title = "Strand Balance of Aligned Pass-Filter Reads per Sample"
    )
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = ggplot2::rel(x = 0.8),
        hjust = 0.0,
        vjust = 0.5,
        angle = 90.0
      )
    )
  for (graphics_format in graphics_formats) {
    if (graphics_format == "png" &&
        plot_width_sample > graphics_maximum_size_png) {
      message(
        "PNG plot exceeding maximum size: ",
        plot_width_sample,
        " > ",
        graphics_maximum_size_png
      )
      next
    }
    ggplot2::ggsave(
      filename = paste(
        paste(
          prefix_summary,
          prefix_pasm,
          "strand_balance",
          "sample",
          sep = "_"
        ),
        graphics_format,
        sep = "."
      ),
      plot = ggplot_object,
      width = plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot the strand balance of aligned pass-filter reads per read group ----


  message("Plotting the strand balance of aligned pass-filter reads per read group")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_read_group[, c("CATEGORY", "READ_GROUP", "STRAND_BALANCE"), drop = FALSE])
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$READ_GROUP,
        y = .data$STRAND_BALANCE,
        colour = .data$CATEGORY
      ),
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Read Group",
      y = "Reads Fraction",
      colour = "Category",
      title = "Strand Balance of Aligned Pass-Filter Reads per Read Group"
    )
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = ggplot2::rel(x = 0.8),
        hjust = 0.0,
        vjust = 0.5,
        angle = 90.0
      )
    )
  for (graphics_format in graphics_formats) {
    if (graphics_format == "png" &&
        plot_width_read_group > graphics_maximum_size_png) {
      message(
        "PNG plot exceeding maximum size: ",
        plot_width_read_group,
        " > ",
        graphics_maximum_size_png
      )
      next
    }
    ggplot2::ggsave(
      filename = paste(
        paste(
          prefix_summary,
          prefix_pasm,
          "strand_balance",
          "read_group",
          sep = "_"
        ),
        graphics_format,
        sep = "."
      ),
      plot = ggplot_object,
      width = plot_width_read_group,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  rm(plot_width_read_group, plot_width_sample)
}
rm(combined_metrics_read_group, combined_metrics_sample)

# Picard Duplication Metrics ----------------------------------------------


# Process Picard Duplication Metrics reports.

message("Processing Picard Duplication Metrics reports for sample:")
combined_metrics_sample <- NULL

file_names <-
  base::list.files(
    pattern = paste0(
      "^",
      argument_list$prefix,
      "_sample_.*_duplicate_metrics.tsv$"
    ),
    recursive = TRUE
  )
for (file_name in file_names) {
  sample_name <-
    base::gsub(
      pattern = paste0(
        "^",
        argument_list$prefix,
        "_sample_(.*?)_duplicate_metrics.tsv$"
      ),
      replacement = "\\1",
      x = base::basename(path = file_name)
    )
  message("  ", sample_name)

  # Picard Tools added a histogram section that needs excluding from parsing.
  # Find the lines starting with "## METRICS CLASS" and "## HISTOGRAM" and read that many lines.
  metrics_lines <- readLines(con = file_name)
  metrics_line <-
    which(x = grepl(pattern = "## METRICS CLASS", x = metrics_lines))
  histogram_line <-
    which(x = grepl(pattern = "## HISTOGRAM", x = metrics_lines))
  if (length(x = histogram_line)) {
    # Set the number of rows to read excluding 3 more lines,
    # the "## HISTOGRAM" line, the blank line and the header line.
    number_read <- histogram_line[1L] - metrics_line[1L] - 3L
    number_skip <- metrics_line[1L]
  } else {
    number_read <- -1L
    number_skip <- metrics_line[1L]
  }
  picard_metrics_sample <-
    utils::read.table(
      file = file_name,
      header = TRUE,
      sep = "\t",
      nrows = number_read,
      skip = number_skip,
      fill = TRUE,
      comment.char = "",
      stringsAsFactors = FALSE
    )
  rm(number_read,
     number_skip,
     histogram_line,
     metrics_line,
     metrics_lines)

  # Add the sample name, which is not part of the Picard report.
  picard_metrics_sample$SAMPLE <- as.character(x = sample_name)

  # The Picard Duplication Metrics report has changed format through versions.
  # Column SECONDARY_OR_SUPPLEMENTARY_RDS was added at a later stage.
  if (!"SECONDARY_OR_SUPPLEMENTARY_RDS" %in% picard_metrics_sample) {
    picard_metrics_sample$SECONDARY_OR_SUPPLEMENTARY_RDS <- 0L
  }

  combined_metrics_sample <-
    base::rbind(combined_metrics_sample, picard_metrics_sample)

  rm(sample_name, picard_metrics_sample)
}
rm(file_name, file_names)

if (!is.null(x = combined_metrics_sample)) {
  # Order the sample frame by SAMPLE.
  combined_metrics_sample <-
    combined_metrics_sample[order(combined_metrics_sample$SAMPLE),]
  # Convert the SAMPLE column into factors, which come more handy for plotting.
  combined_metrics_sample$SAMPLE <-
    as.factor(x = combined_metrics_sample$SAMPLE)
  # Add additional percentages into the table.
  combined_metrics_sample$PERCENT_UNPAIRED_READ_DUPLICATION <-
    combined_metrics_sample$UNPAIRED_READ_DUPLICATES / combined_metrics_sample$UNPAIRED_READS_EXAMINED
  combined_metrics_sample$PERCENT_READ_PAIR_DUPLICATION <-
    combined_metrics_sample$READ_PAIR_DUPLICATES / combined_metrics_sample$READ_PAIRS_EXAMINED
  combined_metrics_sample$PERCENT_READ_PAIR_OPTICAL_DUPLICATION <-
    combined_metrics_sample$READ_PAIR_OPTICAL_DUPLICATES / combined_metrics_sample$READ_PAIRS_EXAMINED
  utils::write.table(
    x = combined_metrics_sample,
    file = paste(prefix_summary, prefix_pdsm, "metrics_sample.tsv", sep = "_"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  # Adjust the plot width according to batches of 24 samples or read groups.
  plot_width <-
    argument_list$plot_width + (ceiling(x = nlevels(x = combined_metrics_sample$SAMPLE) / 24L) - 1L) * argument_list$plot_width * 0.3

  # Plot Percent Duplication per Sample -----------------------------------


  message("Plotting the percent duplication per sample")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_sample)
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(x = .data$SAMPLE, y = .data$PERCENT_DUPLICATION),
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(x = "Sample", y = "Duplication Fraction", title = "Duplication Fraction per Sample")
  ggplot_object <-
    ggplot_object + ggplot2::guides(colour = ggplot2::guide_legend(nrow = 24L))
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = ggplot2::rel(x = 0.8),
        hjust = 0.0,
        vjust = 0.5,
        angle = 90.0
      )
    )
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
        paste(prefix_summary,
              prefix_pdsm,
              "percentage",
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
  rm(graphics_format, ggplot_object)

  # Plot Duplication Classes per Sample -----------------------------------


  # Plot PERCENT_UNPAIRED_READ_DUPLICATION, PERCENT_READ_PAIR_DUPLICATION,
  # PERCENT_READ_PAIR_OPTICAL_DUPLICATION and PERCENT_DUPLICATION per sample.

  message("Plotting the duplication levels per sample")

  ggplot_object <- ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = combined_metrics_sample,
      cols = c(
        .data$PERCENT_UNPAIRED_READ_DUPLICATION,
        .data$PERCENT_READ_PAIR_DUPLICATION,
        .data$PERCENT_READ_PAIR_OPTICAL_DUPLICATION,
        .data$PERCENT_DUPLICATION
      ),
      names_to = "DUPLICATION",
      values_to = "fraction"
    )
  )
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$SAMPLE,
        y = .data$fraction,
        colour = .data$DUPLICATION
      ),
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Sample",
      y = "Duplication Fraction",
      colour = "Duplication",
      title = "Duplication Levels per Sample"
    )
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = ggplot2::rel(x = 0.8),
        hjust = 0.0,
        vjust = 0.5,
        angle = 90.0
      )
    )
  for (graphics_format in graphics_formats) {
    ggplot2::ggsave(
      filename = paste(
        paste(prefix_summary,
              prefix_pdsm,
              "levels",
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

  rm(graphics_format, ggplot_object)

  rm(plot_width)
}
rm(combined_metrics_sample)

rm(
  prefix_pdsm,
  prefix_pasm,
  prefix_summary,
  argument_list,
  graphics_maximum_size_png,
  graphics_formats
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessioninfo::session_info())

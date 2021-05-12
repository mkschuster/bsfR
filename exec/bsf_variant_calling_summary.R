#!/usr/bin/env Rscript
#
# BSF R script to summarise a variant calling analysis. Picard Duplication
# Metrics, Picard Alignment Summary Metrics and Picard Hybrid Selection Metrics
# reports are read for each sample and plotted at the read group or sample
# level.
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
        opt_str = c("--prefix"),
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
# Maximum size for the PNG device in inches.
graphics_maximum_size_png <- 100.0

# Assign a file prefix.
prefix_summary <- "variant_calling_summary"
if (is.null(x = argument_list$prefix)) {
  # If a prefix was not provided, try to get it from a cohort-level file.
  file_names <-
    base::list.files(pattern = "^variant_calling_process_cohort_.*_genotyped_raw_snp_raw_indel.vcf.gz$")
  for (file_name in file_names) {
    cohort_name <-
      gsub(pattern = "^variant_calling_process_cohort_(.*?)_genotyped_raw_snp_raw_indel.vcf.gz$",
           replacement = "\\1",
           x = file_name)
    message("Cohort name: ", cohort_name)
    prefix_summary <-
      paste("variant_calling_summary", cohort_name, sep = "_")
    rm(cohort_name)
  }
  rm(file_names)
} else {
  prefix_summary <- argument_list$prefix
}

# Picard Duplication Metrics ----------------------------------------------


# Process Picard Duplication Metrics reports.

message("Processing Picard Duplication Metrics reports for sample:")
combined_metrics_sample <- NULL

file_names <-
  base::list.files(pattern = "^variant_calling_process_sample_.*_duplicate_metrics.tsv$")
for (file_name in file_names) {
  sample_name <-
    gsub(pattern = "^variant_calling_process_sample_(.*?)_duplicate_metrics.tsv$",
         replacement = "\\1",
         x = file_name)
  message("  ", sample_name)

  # Picard Tools added a histogram section that needs excluding from parsing.
  # Find the lines starting with "## METRICS CLASS" and "## HISTOGRAM" and read
  # that many lines.
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
  if (is.null(x = picard_metrics_sample$SECONDARY_OR_SUPPLEMENTARY_RDS)) {
    picard_metrics_sample$SECONDARY_OR_SUPPLEMENTARY_RDS <- 0L
  }

  if (is.null(x = combined_metrics_sample)) {
    combined_metrics_sample <- picard_metrics_sample
  } else {
    combined_metrics_sample <-
      rbind(combined_metrics_sample, picard_metrics_sample)
  }
  rm(picard_metrics_sample)
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
    file = paste(prefix_summary, "duplication_metrics_sample.tsv", sep = "_"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  # Adjust the plot width according to batches of 24 samples or read groups.
  plot_width_sample <-
    argument_list$plot_width + (ceiling(x = nlevels(x = combined_metrics_sample$SAMPLE) / 24L) - 1L) * argument_list$plot_width * 0.3
  # message("Plot width sample: ", plot_width_sample)

  # Plot Duplication Fraction per Sample ----------------------------------


  message("Plotting the duplication fraction per sample")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_sample)
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(x = .data$SAMPLE, y = .data$PERCENT_DUPLICATION))
  ggplot_object <-
    ggplot_object + ggplot2::labs(x = "Sample", y = "Duplication Fraction", title = "Duplication Fraction per Sample")
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
        prefix_summary,
        paste("duplication_percentage_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_sample > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_sample,
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
    ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(
      x = .data$SAMPLE,
      y = .data$fraction,
      colour = .data$DUPLICATION
    ))
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Sample",
      y = "Fraction",
      colour = "Duplication Level",
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
        prefix_summary,
        paste("duplication_levels_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_sample > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }

  rm(graphics_format, ggplot_object)

  rm(plot_width_sample)
}
rm(combined_metrics_sample)

# Picard Alignment Summary Metrics ----------------------------------------


# Process Picard Alignment Summary Metrics reports.

message("Processing Picard Alignment Summary Metrics reports for sample:")
combined_metrics_sample <- NULL
combined_metrics_read_group <- NULL

file_names <-
  base::list.files(pattern = "^variant_calling_process_sample_.*_alignment_summary_metrics.tsv$")
for (file_name in file_names) {
  sample_name <-
    gsub(pattern = "^variant_calling_process_sample_(.*?)_alignment_summary_metrics.tsv$",
         replacement = "\\1",
         x = file_name)
  message("  ", sample_name)
  # Since the Illumina2bam tools BamIndexDecoder uses a hash character (#) in
  # the read group component to separate platform unit and sample name, the
  # Picard reports need special parsing. Find the ## METRICS CLASS line and
  # parse without allowing further comments.
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
  # To support numeric sample names the read.table(stringsAsFactors = FALSE) is
  # turned off. Convert SAMPLE, LIBRARY and READ_GROUP into character vectors.
  picard_metrics_total$SAMPLE <-
    as.character(x = picard_metrics_total$SAMPLE)
  picard_metrics_total$LIBRARY <-
    as.character(x = picard_metrics_total$LIBRARY)
  picard_metrics_total$READ_GROUP <-
    as.character(x = picard_metrics_total$READ_GROUP)

  # The Picard Alignment Metrics report has changed format through versions.
  # Columns PF_READS_IMPROPER_PAIRS and PCT_PF_READS_IMPROPER_PAIRS were added
  # at a later stage.
  if (is.null(x = picard_metrics_total$PF_READS_IMPROPER_PAIRS)) {
    picard_metrics_total$PF_READS_IMPROPER_PAIRS <- 0L
  }
  if (is.null(x = picard_metrics_total$PCT_PF_READS_IMPROPER_PAIRS)) {
    picard_metrics_total$PCT_PF_READS_IMPROPER_PAIRS <- 0.0
  }

  # Select only rows showing the SAMPLE summary, i.e. showing SAMPLE, but no
  # LIBRARY and READ_GROUP information.
  picard_metrics_sample <-
    picard_metrics_total[(!is.na(x = picard_metrics_total$SAMPLE)) &
                           (picard_metrics_total$SAMPLE != "") &
                           (picard_metrics_total$LIBRARY == "") &
                           (picard_metrics_total$READ_GROUP == ""), ]
  if (is.null(x = combined_metrics_sample)) {
    combined_metrics_sample <- picard_metrics_sample
  } else {
    combined_metrics_sample <-
      rbind(combined_metrics_sample, picard_metrics_sample)
  }
  rm(picard_metrics_sample)

  # Select only rows showing READ_GROUP summary, i.e. showing READ_GROUP
  # information.
  picard_metrics_read_group <-
    picard_metrics_total[(picard_metrics_total$READ_GROUP != ""), ]
  if (is.null(x = combined_metrics_read_group)) {
    combined_metrics_read_group <- picard_metrics_read_group
  } else {
    combined_metrics_read_group <-
      rbind(combined_metrics_read_group, picard_metrics_read_group)
  }
  rm(picard_metrics_read_group)

  rm(sample_name, picard_metrics_total)
}
rm(file_name, file_names)

if (!is.null(x = combined_metrics_sample)) {
  # Order the data frame by SAMPLE.
  combined_metrics_sample <-
    combined_metrics_sample[order(combined_metrics_sample$SAMPLE),]
  # Manually convert CATEGORY and SAMPLE columns into factors, which are handy
  # for plotting.
  combined_metrics_sample$CATEGORY <-
    as.factor(x = combined_metrics_sample$CATEGORY)
  combined_metrics_sample$SAMPLE <-
    as.factor(x = combined_metrics_sample$SAMPLE)
  # Add an additional LABEL factor column defined as a concatenation of SAMPLE
  # and CATEGORY.
  combined_metrics_sample$LABEL <-
    as.factor(x = paste(
      combined_metrics_sample$SAMPLE,
      combined_metrics_sample$CATEGORY,
      sep =
        "_"
    ))
  utils::write.table(
    x = combined_metrics_sample,
    file = paste(prefix_summary, "alignment_metrics_sample.tsv", sep = "_"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  # Order the data frame by READ_GROUP
  combined_metrics_read_group <-
    combined_metrics_read_group[order(combined_metrics_read_group$READ_GROUP),]
  # Manually convert CATEGORY and READ_GROUP columns into factors, which are
  # handy for plotting.
  combined_metrics_read_group$CATEGORY <-
    as.factor(x = combined_metrics_read_group$CATEGORY)
  combined_metrics_read_group$READ_GROUP <-
    as.factor(x = combined_metrics_read_group$READ_GROUP)
  # Add an additional LABEL factor column defined as a concatenation of
  # READ_GROUP and CATEGORY.
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
    file = paste(prefix_summary, "alignment_metrics_read_group.tsv", sep = "_"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  # Adjust the plot width according to batches of 24 samples or read groups.
  plot_width_sample <-
    argument_list$plot_width + (ceiling(x = nlevels(x = combined_metrics_sample$SAMPLE) / 24L) - 1L) * argument_list$plot_width * 0.25
  plot_width_read_group <-
    argument_list$plot_width + (ceiling(x = nlevels(x = combined_metrics_read_group$READ_GROUP) / 24L) - 1L) * argument_list$plot_width * 0.35
  # message("Plot width sample: ", plot_width_sample)
  # message("Plot width read group: ", plot_width_read_group)

  # Plot the aligned pass-filter reads number per sample ------------------


  message("Plotting the aligned pass-filter reads number per sample")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_sample[, c("CATEGORY", "SAMPLE", "PF_READS_ALIGNED"), drop = FALSE])
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(
      x = .data$SAMPLE,
      y = .data$PF_READS_ALIGNED,
      colour = .data$CATEGORY
    ))
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Sample",
      y = "Number",
      colour = "Category",
      title = "Aligned Pass-Filter Reads Number per Sample"
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
        prefix_summary,
        paste("alignment_absolute_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_sample > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot the aligned pass-filter reads number per read group --------------


  message("Plotting the aligned pass-filter reads number per read group")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_read_group[, c("CATEGORY", "READ_GROUP", "PF_READS_ALIGNED"), drop = FALSE])
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$READ_GROUP,
        y = .data$PF_READS_ALIGNED,
        colour = .data$CATEGORY
      )
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Read Group",
      y = "Number",
      colour = "Category",
      title = "Aligned Pass-Filter Reads Number per Read Group"
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
        prefix_summary,
        paste("alignment_absolute_read_group", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_read_group > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_read_group,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot the aligned pass-filter reads fraction per sample ----------------


  message("Plotting the aligned pass-filter reads fraction per sample")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_sample[, c("CATEGORY", "SAMPLE", "PCT_PF_READS_ALIGNED"), drop = FALSE])
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$SAMPLE,
        y = .data$PCT_PF_READS_ALIGNED,
        colour = .data$CATEGORY
      )
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Sample",
      y = "Fraction",
      colour = "Category",
      title = "Aligned Pass-Filter Reads Fraction per Sample"
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
        prefix_summary,
        paste("alignment_percentage_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_sample > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot the aligned pass-filter reads fraction per read group ------------


  message("Plotting the aligned pass-filter reads fraction per read group")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_read_group[, c("CATEGORY", "READ_GROUP", "PCT_PF_READS_ALIGNED"), drop = FALSE])
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$READ_GROUP,
        y = .data$PCT_PF_READS_ALIGNED,
        colour = .data$CATEGORY
      )
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Read Group",
      y = "Fraction",
      colour = "Category",
      title = "Aligned Pass-Filter Reads Fraction per Read Group"
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
        prefix_summary,
        paste("alignment_percentage_read_group", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_read_group > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_read_group,
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
    ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(
      x = .data$SAMPLE,
      y = .data$STRAND_BALANCE,
      colour = .data$CATEGORY
    ))
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Sample",
      y = "Fraction",
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
    ggplot2::ggsave(
      filename = paste(
        prefix_summary,
        paste("alignment_strand_balance_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_sample > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_sample,
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
      )
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Read Group",
      y = "Fraction",
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
    ggplot2::ggsave(
      filename = paste(
        prefix_summary,
        paste(
          "alignment_strand_balance_read_group",
          graphics_format,
          sep = "."
        ),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_read_group > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_read_group,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  rm(plot_width_read_group, plot_width_sample)
}
rm(combined_metrics_read_group, combined_metrics_sample)

# Picard Hybrid Selection Metrics -----------------------------------------


# Process Picard Hybrid Selection Metrics reports.

message("Processing Picard Hybrid Selection Metrics reports for sample:")
combined_metrics_sample <- NULL
combined_metrics_read_group <- NULL

file_names <-
  base::list.files(pattern = "^variant_calling_diagnose_sample_.*_hybrid_selection_metrics.tsv$")
for (file_name in file_names) {
  sample_name <-
    gsub(pattern = "^variant_calling_diagnose_sample_(.*?)_hybrid_selection_metrics.tsv$",
         replacement = "\\1",
         x = file_name)
  message("  ", sample_name)
  # Picard Tools added a histogram section that needs excluding from parsing.
  # Find the lines starting with "## METRICS CLASS" and "## HISTOGRAM" and read
  # that many lines.
  metrics_lines <- readLines(con = file_name)
  metrics_line <-
    which(x = grepl(pattern = "## METRICS CLASS", x = metrics_lines))
  histogram_line <-
    which(x = grepl(pattern = "## HISTOGRAM", x = metrics_lines))
  if (length(x = histogram_line)) {
    number_read <- histogram_line[1L] - metrics_line[1L] - 3L
    number_skip <- metrics_line[1L]
  } else {
    number_read <- -1L
    number_skip <- metrics_line[1L]
  }
  picard_metrics_total <-
    utils::read.table(
      file = file_name,
      header = TRUE,
      sep = "\t",
      # Set the number of rows excluding 3 more lines,
      # the "## HISTOGRAM" line, the blank line and the header line.
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
  # To support numeric sample names the read.table(stringsAsFactors = FALSE) is
  # turned off. Convert SAMPLE, LIBRARY and READ_GROUP into character vectors.
  picard_metrics_total$SAMPLE <-
    as.character(x = picard_metrics_total$SAMPLE)
  picard_metrics_total$LIBRARY <-
    as.character(x = picard_metrics_total$LIBRARY)
  picard_metrics_total$READ_GROUP <-
    as.character(x = picard_metrics_total$READ_GROUP)

  # The Picard Hybrid Selection Metrics report has changed format through versions.
  # Column PCT_TARGET_BASES_1X was added at a later stage.
  if (!"PCT_TARGET_BASES_1X" %in% names(x = picard_metrics_total)) {
    picard_metrics_total$PCT_TARGET_BASES_1X <- 0.0
  }

  if (!"MAX_TARGET_COVERAGE" %in% names(x = picard_metrics_total)) {
    picard_metrics_total$MAX_TARGET_COVERAGE <- 0L
  }

  if (!"PCT_EXC_ADAPTER" %in% names(x = picard_metrics_total)) {
    picard_metrics_total$PCT_EXC_ADAPTER <- 0.0
  }

  if (!"PF_BASES" %in% names(x = picard_metrics_total)) {
    picard_metrics_total$PF_BASES <- 0L
  }

  # Select only rows showing the SAMPLE summary, i.e. showing SAMPLE, but no
  # LIBRARY and READ_GROUP information.
  picard_metrics_sample <-
    picard_metrics_total[(!is.na(x = picard_metrics_total$SAMPLE)) &
                           (picard_metrics_total$SAMPLE != "") &
                           (picard_metrics_total$LIBRARY == "") &
                           (picard_metrics_total$READ_GROUP == ""), ]
  if (is.null(x = combined_metrics_sample)) {
    combined_metrics_sample <- picard_metrics_sample
  } else {
    combined_metrics_sample <-
      rbind(combined_metrics_sample, picard_metrics_sample)
  }
  rm(picard_metrics_sample)

  # Select only rows showing READ_GROUP summary, i.e. showing READ_GROUP
  # information.
  picard_metrics_read_group <-
    picard_metrics_total[(picard_metrics_total$READ_GROUP != ""), ]
  if (is.null(x = combined_metrics_read_group)) {
    combined_metrics_read_group <- picard_metrics_read_group
  } else {
    combined_metrics_read_group <-
      rbind(combined_metrics_read_group, picard_metrics_read_group)
  }
  rm(picard_metrics_read_group)

  rm(sample_name, picard_metrics_total)
}
rm(file_name, file_names)

# The Picard Hybrid Selection Metrics is currently optional.

if (!is.null(x = combined_metrics_sample)) {
  # Sort the data frame by SAMPLE.
  combined_metrics_sample <-
    combined_metrics_sample[order(combined_metrics_sample$SAMPLE), ]
  # Manually convert BAIT_SET and SAMPLE columns into factors, which are handy
  # for plotting.
  combined_metrics_sample$BAIT_SET <-
    as.factor(x = combined_metrics_sample$BAIT_SET)
  combined_metrics_sample$SAMPLE <-
    as.factor(x = combined_metrics_sample$SAMPLE)
  utils::write.table(
    x = combined_metrics_sample,
    file = paste(prefix_summary, "hybrid_metrics_sample.tsv", sep = "_"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  # Sort the data frame by READ_GROUP.
  combined_metrics_read_group <-
    combined_metrics_read_group[order(combined_metrics_read_group$READ_GROUP), ]
  # Manually convert BAIT_SET and READ_GROUP columns into factors, which are
  # handy for plotting.
  combined_metrics_read_group$BAIT_SET <-
    as.factor(x = combined_metrics_read_group$BAIT_SET)
  combined_metrics_read_group$READ_GROUP <-
    as.factor(x = combined_metrics_read_group$READ_GROUP)
  utils::write.table(
    x = combined_metrics_read_group,
    file = paste(prefix_summary, "hybrid_metrics_read_group.tsv", sep = "_"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  # Adjust the plot width according to batches of 24 samples or read groups.
  plot_width_sample <- argument_list$plot_width + (ceiling(x = (
    nlevels(x = combined_metrics_sample$SAMPLE) / 24L
  )) - 1L) * argument_list$plot_width * 0.25
  plot_width_read_group <- argument_list$plot_width + (ceiling(x = (
    nlevels(x = combined_metrics_read_group$READ_GROUP) / 24L
  )) - 1L) * argument_list$plot_width * 0.25
  # message("Plot width sample: ", plot_width_sample)
  # message("Plot width read group: ", plot_width_read_group)

  # Plot the percentage of unique pass-filter reads per sample ------------


  message("Plotting the percentage of unique pass-filter reads per sample")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_sample)
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(x = .data$SAMPLE, y = .data$PCT_PF_UQ_READS))
  ggplot_object <-
    ggplot_object + ggplot2::labs(x = "Sample" , y = "Fraction PF Unique", title = "Unique Pass-Filter Reads per Sample")
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
        prefix_summary,
        paste("hybrid_unique_percentage_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_sample > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot the percentage of unique pass-filter reads per read group --------


  message("Plotting the percentage of unique pass-filter reads per read group")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_read_group)
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$READ_GROUP,
        y = .data$PCT_PF_UQ_READS,
        shape = .data$BAIT_SET
      )
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Read Group",
      y = "PF Unique Fraction",
      shape = "Bait Set",
      title = "Unique Pass-Filter Reads per Read Group"
    )
  # For more than six shapes (scale_shape()), a manual scale
  # (scale_shape_manual()) needs setting up.
  # https://ggplot2.tidyverse.org/reference/scale_shape.html
  ggplot_object <-
    ggplot_object + ggplot2::scale_shape_manual(values = seq_len(length.out = nlevels(x = combined_metrics_read_group$BAIT_SET)))
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = ggplot2::rel(x = 0.5),
        hjust = 0.0,
        vjust = 0.5,
        angle = 90.0
      )
    )
  for (graphics_format in graphics_formats) {
    ggplot2::ggsave(
      filename = paste(
        prefix_summary,
        paste(
          "hybrid_unique_percentage_read_group",
          graphics_format,
          sep = "."
        ),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_read_group > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_read_group,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot the mean target coverage per sample ------------------------------


  message("Plotting the mean target coverage per sample")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_sample)
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(x = .data$SAMPLE, y = .data$MEAN_TARGET_COVERAGE))
  ggplot_object <-
    ggplot_object + ggplot2::labs(x = "Sample", y = "Mean Target Coverage", title = "Mean Target Coverage per Sample")
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
        prefix_summary,
        paste("hybrid_target_coverage_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_sample > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot the mean target coverage per read group --------------------------


  message("Plotting the mean target coverage per read group")
  ggplot_object <-
    ggplot2::ggplot(data = combined_metrics_read_group)
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$READ_GROUP,
        y = .data$MEAN_TARGET_COVERAGE,
        shape = .data$BAIT_SET
      )
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Read Group",
      y = "Mean Target Coverage",
      shape = "Bait Set",
      title = "Mean Target Coverage per Read Group"
    )
  # For more than six shapes (scale_shape()), a manual scale
  # (scale_shape_manual()) needs setting up.
  # https://ggplot2.tidyverse.org/reference/scale_shape.html
  ggplot_object <-
    ggplot_object + ggplot2::scale_shape_manual(values = seq_len(length.out = nlevels(x = combined_metrics_read_group$BAIT_SET)))
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = ggplot2::rel(x = 0.5),
        hjust = 0.0,
        vjust = 0.5,
        angle = 90.0
      )
    )
  for (graphics_format in graphics_formats) {
    ggplot2::ggsave(
      filename = paste(
        prefix_summary,
        paste(
          "hybrid_target_coverage_read_group",
          graphics_format,
          sep = "."
        ),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_read_group > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_read_group,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot the percentage of exluded bases per sample -----------------------


  message("Plotting the percentage of excluded bases per sample")

  ggplot_object <- ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = combined_metrics_sample,
      cols = c(
        .data$PCT_EXC_DUPE,
        .data$PCT_EXC_MAPQ,
        .data$PCT_EXC_BASEQ,
        .data$PCT_EXC_OVERLAP,
        .data$PCT_EXC_OFF_TARGET
      ),
      names_to = "EXCLUDED",
      values_to = "fraction"
    )
  )
  ggplot_object <-
    ggplot_object + ggplot2::geom_col(
      mapping = ggplot2::aes(
        x = .data$SAMPLE,
        y = .data$fraction,
        fill = .data$EXCLUDED
      ),
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Sample",
      y = "Fraction",
      fill = "Excluded",
      title = "Excluded Bases per Sample"
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
        prefix_summary,
        paste("hybrid_excluded_bases_sample",
              graphics_format,
              sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_sample > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot the percentage of excluded bases per read group ------------------


  message("Plotting the percentage of excluded bases per read group")

  ggplot_object <- ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = combined_metrics_read_group,
      cols = c(
        .data$PCT_EXC_DUPE,
        .data$PCT_EXC_MAPQ,
        .data$PCT_EXC_BASEQ,
        .data$PCT_EXC_OVERLAP,
        .data$PCT_EXC_OFF_TARGET
      ),
      names_to = "EXCLUDED",
      values_to = "fraction"
    )
  )
  ggplot_object <-
    ggplot_object + ggplot2::geom_col(
      mapping = ggplot2::aes(
        x = .data$READ_GROUP,
        y = .data$fraction,
        fill = .data$EXCLUDED,
        colour = .data$BAIT_SET
      ),
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Read Group",
      y = "Fraction",
      fill = "Excluded",
      colour = "Bait Set",
      title = "Excluded Bases per Read Group"
    )
  ggplot_object <-
    ggplot_object + ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1L),
      colour = ggplot2::guide_legend(order = 2L)
    )
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = ggplot2::rel(x = 0.5),
        hjust = 0.0,
        vjust = 0.5,
        angle = 90.0
      )
    )
  for (graphics_format in graphics_formats) {
    ggplot2::ggsave(
      filename = paste(
        prefix_summary,
        paste(
          "hybrid_excluded_bases_read_group",
          graphics_format,
          sep = "."
        ),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_sample > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot target coverage levels per sample --------------------------------


  # Plot PCT_TARGET_BASES_1X, PCT_TARGET_BASES_2X, PCT_TARGET_BASES_10X,
  # PCT_TARGET_BASES_20X, PCT_TARGET_BASES_30X, PCT_TARGET_BASES_40X,
  # PCT_TARGET_BASES_50X, PCT_TARGET_BASES_100X per sample.
  message("Plotting the coverage levels per sample")

  ggplot_object <- ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = combined_metrics_sample,
      cols = c(
        .data$PCT_TARGET_BASES_1X,
        .data$PCT_TARGET_BASES_2X,
        .data$PCT_TARGET_BASES_10X,
        .data$PCT_TARGET_BASES_20X,
        .data$PCT_TARGET_BASES_30X,
        .data$PCT_TARGET_BASES_40X,
        .data$PCT_TARGET_BASES_50X,
        .data$PCT_TARGET_BASES_100X
      ),
      names_to = "COVERAGE",
      values_to = "fraction"
    )
  )
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(
      x = .data$SAMPLE,
      y = .data$fraction,
      colour = .data$COVERAGE
    ))
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Sample",
      y = "Fraction",
      colour = "Coverage Level",
      title = "Coverage Levels per Sample"
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
        prefix_summary,
        paste(
          "hybrid_target_coverage_levels_sample",
          graphics_format,
          sep = "."
        ),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_sample > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot target coverage levels per read group ----------------------------


  # Plot PCT_TARGET_BASES_1X, PCT_TARGET_BASES_2X, PCT_TARGET_BASES_10X,
  # PCT_TARGET_BASES_20X, PCT_TARGET_BASES_30X, PCT_TARGET_BASES_40X,
  # PCT_TARGET_BASES_50X, PCT_TARGET_BASES_100X per read group.
  message("Plotting the coverage levels per read group")

  ggplot_object <- ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = combined_metrics_read_group,
      cols = c(
        .data$PCT_TARGET_BASES_1X,
        .data$PCT_TARGET_BASES_2X,
        .data$PCT_TARGET_BASES_10X,
        .data$PCT_TARGET_BASES_20X,
        .data$PCT_TARGET_BASES_30X,
        .data$PCT_TARGET_BASES_40X,
        .data$PCT_TARGET_BASES_50X,
        .data$PCT_TARGET_BASES_100X
      ),
      names_to = "COVERAGE",
      values_to = "fraction"
    )
  )
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$READ_GROUP,
        y = .data$fraction,
        colour = .data$COVERAGE,
        shape = .data$BAIT_SET
      )
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Read Group",
      y = "Fraction",
      colour = "Coverage",
      shape = "Bait Set",
      title = "Coverage Levels per Read Group"
    )
  # For more than six shapes (scale_shape()), a manual scale
  # (scale_shape_manual()) needs setting up.
  # https://ggplot2.tidyverse.org/reference/scale_shape.html
  ggplot_object <-
    ggplot_object + ggplot2::scale_shape_manual(values = seq_len(length.out = nlevels(x = combined_metrics_read_group$BAIT_SET)))
  ggplot_object <-
    ggplot_object + ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1L),
      shape = ggplot2::guide_legend(order = 2L)
    )
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = ggplot2::rel(x = 0.5),
        hjust = 0.0,
        vjust = 0.5,
        angle = 90.0
      )
    )
  for (graphics_format in graphics_formats) {
    ggplot2::ggsave(
      filename = paste(
        prefix_summary,
        paste(
          "hybrid_target_coverage_levels_read_group",
          graphics_format,
          sep = "."
        ),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_sample > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot the nominal coverage per sample ----------------------------------


  # Plot the nominal coverage (i.e. PF_BASES_ALIGNED / TARGET_TERRITORY) per
  # sample.

  message("Plotting the nominal coverage per sample")
  plotting_frame <-
    combined_metrics_sample[, c("SAMPLE",
                                "READ_GROUP",
                                "BAIT_SET",
                                "PF_BASES_ALIGNED",
                                "TARGET_TERRITORY")]
  plotting_frame$NOMINAL_COVERAGE <-
    plotting_frame$PF_BASES_ALIGNED / plotting_frame$TARGET_TERRITORY

  ggplot_object <- ggplot2::ggplot(data = plotting_frame)
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(x = .data$SAMPLE,
                                                               y = .data$NOMINAL_COVERAGE))
  ggplot_object <-
    ggplot_object + ggplot2::labs(x = "Sample", y = "Nominal Coverage", title = "Nominal Coverage per Sample")
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
        prefix_summary,
        paste("hybrid_nominal_coverage_sample",
              graphics_format,
              sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_sample > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object, plotting_frame)

  # Plot the nominal coverage per read group ------------------------------


  # Plot the nominal coverage (i.e. PF_BASES_ALIGNED / TARGET_TERRITORY) per
  # read group.

  message("Plotting the nominal coverage per read group")
  plotting_frame <-
    combined_metrics_read_group[, c("SAMPLE",
                                    "READ_GROUP",
                                    "BAIT_SET",
                                    "PF_BASES_ALIGNED",
                                    "TARGET_TERRITORY")]
  plotting_frame$NOMINAL_COVERAGE <-
    plotting_frame$PF_BASES_ALIGNED / plotting_frame$TARGET_TERRITORY

  ggplot_object <- ggplot2::ggplot(data = plotting_frame)
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$READ_GROUP,
        y = .data$NOMINAL_COVERAGE,
        shape = .data$BAIT_SET
      )
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Read Group",
      y = "Nominal Coverage",
      shape = "Bait Set",
      title = "Nominal Coverage per Read Group"
    )
  # For more than six shapes (scale_shape()), a manual scale
  # (scale_shape_manual()) needs setting up.
  # https://ggplot2.tidyverse.org/reference/scale_shape.html
  ggplot_object <-
    ggplot_object + ggplot2::scale_shape_manual(values = seq_len(length.out = nlevels(x = combined_metrics_read_group$BAIT_SET)))
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = ggplot2::rel(x = 0.5),
        hjust = 0.0,
        vjust = 0.5,
        angle = 90.0
      )
    )
  for (graphics_format in graphics_formats) {
    ggplot2::ggsave(
      filename = paste(
        prefix_summary,
        paste(
          "hybrid_nominal_coverage_read_group",
          graphics_format,
          sep = "."
        ),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_read_group > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_read_group,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object, plotting_frame)

  rm(plot_width_read_group, plot_width_sample)
}
rm(combined_metrics_read_group, combined_metrics_sample)

# Non-callable Summary Reports --------------------------------------------


# Process non-callable summary reports.

message("Processing non-callable summary reports for sample:")

# Initialise a data frame with all possible columns and rows at once.
combined_metrics_sample <- data.frame(
  row.names = base::list.files(pattern = "^variant_calling_diagnose_sample_.*_non_callable_summary.tsv$")
)
combined_metrics_sample$file_name <-
  row.names(x = combined_metrics_sample)
combined_metrics_sample$exon_path <-
  character(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$exon_number_raw <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$exon_width_raw <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$exon_number <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$exon_width <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$exon_width_flank <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$transcribed_number <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$transcribed_width <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$target_path <-
  character(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$target_number_raw <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$target_width_raw <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$target_width_flank <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$constrained_number <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$constrained_width <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$callable_loci_path <-
  character(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$sample_name <-
  character(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$non_callable_number_raw <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$non_callable_width_raw <-
  integer(length = nrow(x = combined_metrics_sample))
for (level in c(
  "REF_N",
  # "PASS" according to documentation, but should be "CALLABLE" in practice.
  "NO_COVERAGE",
  "LOW_COVERAGE",
  "EXCESSIVE_COVERAGE",
  "POOR_MAPPING_QUALITY"
)) {
  combined_metrics_sample[, paste("non_callable_number_constrained", level, sep = ".")] <-
    integer(length = nrow(x = combined_metrics_sample))
  combined_metrics_sample[, paste("non_callable_width_constrained", level, sep = ".")] <-
    integer(length = nrow(x = combined_metrics_sample))
}
rm(level)

for (i in seq_len(length.out = nrow(x = combined_metrics_sample))) {
  sample_name <-
    gsub(pattern = "^variant_calling_diagnose_sample_(.*?)_non_callable_summary.tsv$",
         replacement = "\\1",
         x = combined_metrics_sample$file_name[i])
  message("  ", sample_name)
  non_callable_metrics_sample <-
    utils::read.table(
      file = combined_metrics_sample$file_name[i],
      header = TRUE,
      colClasses = c(
        "exon_path" = "character",
        "exon_number_raw" = "integer",
        "exon_width_raw" = "integer",
        "exon_number" = "integer",
        "exon_width" = "integer",
        "exon_width_flank" = "integer",
        "transcribed_number" = "integer",
        "transcribed_width" = "integer",
        "target_path" = "character",
        "target_number_raw" = "integer",
        "target_width_raw" = "integer",
        "target_width_flank" = "integer",
        "constrained_number" = "integer",
        "constrained_width" = "integer",
        "callable_loci_path" = "character",
        "sample_name" = "character",
        "non_callable_number_raw" = "integer",
        "non_callable_width_raw" = "integer",
        "non_callable_number_constrained.TOTAL" = "integer",
        "non_callable_width_constrained.TOTAL" = "integer",
        "non_callable_number_constrained.REF_N" = "integer",
        "non_callable_width_constrained.REF_N" = "integer",
        # "non_callable_number_constrained.CALLABLE" = "integer",
        # "non_callable_width_constrained.CALLABLE" = "integer",
        "non_callable_number_constrained.NO_COVERAGE" = "integer",
        "non_callable_width_constrained.NO_COVERAGE" = "integer",
        "non_callable_number_constrained.LOW_COVERAGE" = "integer",
        "non_callable_width_constrained.LOW_COVERAGE" = "integer",
        "non_callable_number_constrained.EXCESSIVE_COVERAGE" = "integer",
        "non_callable_width_constrained.EXCESSIVE_COVERAGE" = "integer",
        "non_callable_number_constrained.POOR_MAPPING_QUALITY" = "integer",
        "non_callable_width_constrained.POOR_MAPPING_QUALITY" = "integer"
      ),
      fill = TRUE
    )
  for (column_name in names(x = combined_metrics_sample)) {
    if (column_name %in% names(x = non_callable_metrics_sample)) {
      combined_metrics_sample[i, column_name] <-
        non_callable_metrics_sample[[1, column_name]]
    } else {
      # With the exception of the "file_name" component, set all components
      # undefined in the sample-specific data frame to NA in the combined data
      # frame.
      if (column_name != "file_name") {
        combined_metrics_sample[i, column_name] <- NA
      }
    }
  }
  rm(column_name, non_callable_metrics_sample, sample_name)
}
rm(i)

if (nrow(x = combined_metrics_sample) > 0L) {
  # Sort the data frame by sample_name.
  combined_metrics_sample <-
    combined_metrics_sample[order(combined_metrics_sample$sample_name), ]
  # Convert the sample_name column into factors, which come more handy for
  # plotting.
  combined_metrics_sample$sample_name <-
    as.factor(x = combined_metrics_sample$sample_name)

  # Adjust the plot width according to batches of 24 samples or read groups.
  plot_width_sample <- argument_list$plot_width + (ceiling(x = (
    nlevels(x = combined_metrics_sample$sample_name) / 24L
  )) - 1L) * argument_list$plot_width * 0.25
  # message("Plot width sample: ", plot_width_sample)

  utils::write.table(
    x = combined_metrics_sample,
    file = paste(prefix_summary, "non_callable_metrics_sample.tsv", sep = "_"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  # Plot the number of non-callable loci per sample -----------------------


  message("Plotting the number of non-callable loci per sample")
  plotting_frame <- data.frame(
    sample_name = combined_metrics_sample$sample_name,
    constrained_number = combined_metrics_sample$constrained_number
  )
  # Only columns that begin with "^non_callable_number_constrained\." are
  # required, the remainder is the mapping status.
  column_names <- names(x = combined_metrics_sample)
  mapping_status <-
    gsub(pattern = "^non_callable_number_constrained\\.(.*)$",
         replacement = "\\1",
         x = column_names)
  # Extract only those columns which had a match and set the mapping status as
  # their new name.
  for (i in which(x = grepl(pattern = "^non_callable_number_constrained\\.", x = column_names))) {
    plotting_frame[, mapping_status[i]] <-
      combined_metrics_sample[, column_names[i]]
  }
  rm(i, mapping_status, column_names)

  # Now, pivot the data frame, but keep sample_name and constrained_width
  # as identifiers.
  plotting_frame <- tidyr::pivot_longer(
    data = plotting_frame,
    cols = -c(.data$sample_name, .data$constrained_number),
    names_to = "mapping_status",
    values_to = "number"
  )
  # Calculate the fractions on the basis of the constrained target number.
  plotting_frame$fraction <-
    plotting_frame$number / plotting_frame$constrained_number
  # For the moment, remove lines with "TOTAL".
  plotting_frame <-
    plotting_frame[plotting_frame$mapping_status != "TOTAL", ]

  ggplot_object <- ggplot2::ggplot(data = plotting_frame)
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$sample_name,
        y = .data$number,
        colour = .data$mapping_status
      )
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Sample",
      y = "Number",
      colour = "Non-Callable",
      title = "Number of Non-Callable Loci per Sample"
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
        prefix_summary,
        paste("non_callable_absolute_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_sample > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object, plotting_frame)

  # Plot the fraction of non-callable loci per sample ---------------------


  message("Plotting the fraction of non-callable loci per sample")
  # Reorganise the combined_metrics_sample data frame for plotting non-callable
  # target widths.
  plotting_frame <- data.frame(
    sample_name = combined_metrics_sample$sample_name,
    constrained_width = combined_metrics_sample$constrained_width
  )
  # Only columns that begin with "^non_callable_width_constrained\." are
  # required, the remainder is the mapping status.
  column_names <- names(x = combined_metrics_sample)
  mapping_status <-
    gsub(pattern = "^non_callable_width_constrained\\.(.*)$",
         replacement = "\\1",
         x = column_names)
  # Extract only those columns which had a match and set the mapping status as
  # their new name.
  for (i in which(x = grepl(pattern = "^non_callable_width_constrained\\.", x = column_names))) {
    plotting_frame[, mapping_status[i]] <-
      combined_metrics_sample[, column_names[i]]
  }
  rm(i, mapping_status, column_names)

  # Now, pivot the data frame, but keep sample_name and constrained_width
  # as identifiers.
  plotting_frame <- tidyr::pivot_longer(
    data = plotting_frame,
    cols = -c(.data$sample_name, .data$constrained_width),
    names_to = "mapping_status",
    values_to = "width"
  )
  # Calculate the fractions on the basis of the constrained target width.
  plotting_frame$fraction <-
    plotting_frame$width / plotting_frame$constrained_width
  # For the moment, remove lines with "TOTAL".
  plotting_frame <-
    plotting_frame[plotting_frame$mapping_status != "TOTAL", ]

  ggplot_object <- ggplot2::ggplot(data = plotting_frame)
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$sample_name,
        y = .data$fraction,
        colour = .data$mapping_status
      )
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Sample",
      y = "Fraction",
      colour = "Non-Callable",
      title = "Fraction of Non-Callable Loci per Sample"
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
        prefix_summary,
        paste("non_callable_percentage_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width_sample > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object, plotting_frame)

  rm(plot_width_sample)
}
rm(combined_metrics_sample)

rm(prefix_summary,
   argument_list,
   graphics_maximum_size_png,
   graphics_formats)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

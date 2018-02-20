#! /usr/bin/env Rscript
#
# BSF R script to summarise a variant caling analysis. Picard Duplication
# Metrics, Picard Alignment Summary Metrics and Picard Hybrid Selection Metrics
# reports are read for each sample and plotted at the read group or sample level.
#
#
# Copyright 2013 - 2017 Michael K. Schuster
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
      opt_str = c("--prefix"),
      dest = "prefix",
      help = "File name prefix",
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

suppressPackageStartupMessages(expr = library(package = "ggplot2"))
suppressPackageStartupMessages(expr = library(package = "reshape2"))

# Save plots in the following formats.
graphics_formats <- c("pdf", "png")

# Assign a file prefix.
prefix_summary <- "variant_calling_summary"
if (is.null(x = argument_list$prefix)) {
  # If a prefix was not provided, try to get it from a cohort-level file.
  file_names <-
    list.files(pattern = "^variant_calling_process_cohort_.*_genotyped_raw_snp_raw_indel.vcf.gz$")
  for (file_name in file_names) {
    cohort_name <-
      gsub(pattern = "^variant_calling_process_cohort_(.*?)_genotyped_raw_snp_raw_indel.vcf.gz$",
           replacement = "\\1",
           x = file_name)
    message(paste0("Cohort name: ", cohort_name))
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
  list.files(pattern = "^variant_calling_process_sample_.*_duplicate_metrics.tsv$")
for (file_name in file_names) {
  sample_name <-
    gsub(pattern = "^variant_calling_process_sample_(.*?)_duplicate_metrics.tsv$",
         replacement = "\\1",
         x = file_name)
  message(paste0("  ", sample_name))
  
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
    read.table(
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
  write.table(
    x = combined_metrics_sample,
    file = paste(prefix_summary, "duplication_metrics_sample.tsv", sep = "_"),
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
    ggplot(data = combined_metrics_sample)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Percent Duplication per Sample")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = SAMPLE, y = PERCENT_DUPLICATION))
  ggplot_object <-
    ggplot_object + guides(colour = guide_legend(nrow = 24L))
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      size = rel(x = 0.8)
    ))
  for (graphics_format in graphics_formats) {
    ggsave(
      filename = paste(
        prefix_summary,
        paste("duplication_percentage_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = plot_width,
      height = argument_list$plot_height
    )
  }
  rm(graphics_format, ggplot_object)
  
  # Plot Duplication Classes per Sample -----------------------------------
  
  
  # Plot PERCENT_UNPAIRED_READ_DUPLICATION, PERCENT_READ_PAIR_DUPLICATION,
  # PERCENT_READ_PAIR_OPTICAL_DUPLICATION and PERCENT_DUPLICATION per sample.
  
  message("Plotting the duplication levels per sample")
  plotting_frame <- melt(
    data = combined_metrics_sample,
    id.vars = c("SAMPLE"),
    measure.vars = c(
      "PERCENT_UNPAIRED_READ_DUPLICATION",
      "PERCENT_READ_PAIR_DUPLICATION",
      "PERCENT_READ_PAIR_OPTICAL_DUPLICATION",
      "PERCENT_DUPLICATION"
    ),
    variable.name = "DUPLICATION",
    value.name = "fraction"
  )
  
  ggplot_object <- ggplot(data = plotting_frame)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Duplication Levels per Sample")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = SAMPLE,
                                             y = fraction,
                                             colour = DUPLICATION))
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      size = rel(x = 0.8)
    ))
  for (graphics_format in graphics_formats) {
    ggsave(
      filename = paste(
        prefix_summary,
        paste("duplication_levels_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = plot_width,
      height = argument_list$plot_height
    )
  }
  
  rm(graphics_format, ggplot_object, plotting_frame)
  
  rm(plot_width)
}
rm(combined_metrics_sample)

# Picard Alignment Summary Metrics ----------------------------------------


# Process Picard Alignment Summary Metrics reports.

message("Processing Picard Alignment Summary Metrics reports for sample:")
combined_metrics_sample <- NULL
combined_metrics_read_group <- NULL

file_names <-
  list.files(pattern = "^variant_calling_process_sample_.*_alignment_summary_metrics.tsv$")
for (file_name in file_names) {
  sample_name <-
    gsub(pattern = "^variant_calling_process_sample_(.*?)_alignment_summary_metrics.tsv$",
         replacement = "\\1",
         x = file_name)
  message(paste0("  ", sample_name))
  # Since the Illumina2bam tools BamIndexDecoder uses a hash character (#) in the read group component
  # to separate platform unit and sample name, the Picard reports need special parsing.
  # Find the ## METRICS CLASS line and parse without allowing further comments.
  metrics_lines <- readLines(con = file_name)
  metrics_line <-
    which(x = grepl(pattern = "## METRICS CLASS", x = metrics_lines))
  picard_metrics_total <-
    read.table(
      file = file_name,
      header = TRUE,
      sep = "\t",
      skip = metrics_line[1L],
      fill = TRUE,
      comment.char = "",
      stringsAsFactors = FALSE
    )
  rm(metrics_line, metrics_lines)
  # To support numeric sample names the read.table(stringsAsFactors = FALSE) is turned off.
  # Convert SAMPLE, LIBRARY and READ_GROUP into character vectors.
  picard_metrics_total$SAMPLE <-
    as.character(x = picard_metrics_total$SAMPLE)
  picard_metrics_total$LIBRARY <-
    as.character(x = picard_metrics_total$LIBRARY)
  picard_metrics_total$READ_GROUP <-
    as.character(x = picard_metrics_total$READ_GROUP)
  
  # The Picard Alignment Metrics report has changed format through versions.
  # Columns PF_READS_IMPROPER_PAIRS and PCT_PF_READS_IMPROPER_PAIRS were added at a later stage.
  if (is.null(x = picard_metrics_total$PF_READS_IMPROPER_PAIRS)) {
    picard_metrics_total$PF_READS_IMPROPER_PAIRS <- 0L
  }
  if (is.null(x = picard_metrics_total$PCT_PF_READS_IMPROPER_PAIRS)) {
    picard_metrics_total$PCT_PF_READS_IMPROPER_PAIRS <- 0.0
  }
  
  # Select only rows showing the SAMPLE summary, i.e. showing SAMPLE, but no LIBRARY and READ_GROUP information.
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
  
  # Select only rows showing READ_GROUP summary, i.e. showing READ_GROUP information.
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
  write.table(
    x = combined_metrics_sample,
    file = paste(prefix_summary, "alignment_metrics_sample.tsv", sep = "_"),
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
  write.table(
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
  
  # Plot the absolute number of aligned pass-filter reads per sample ------
  
  
  message("Plotting the absolute number of aligned pass-filter reads per sample")
  ggplot_object <-
    ggplot(data = combined_metrics_sample)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Aligned Pass-Filter Reads per Sample")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = CATEGORY, y = PF_READS, colour = SAMPLE))
  ggplot_object <-
    ggplot_object + guides(colour = guide_legend(nrow = 24L))
  for (graphics_format in graphics_formats) {
    ggsave(
      filename = paste(
        prefix_summary,
        paste("alignment_absolute_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = plot_width_sample,
      height = argument_list$plot_height
    )
  }
  rm(graphics_format, ggplot_object)
  
  # Plot the absolute number of aligned pass-filter reads per read group ----
  
  
  message("Plotting the absolute number of aligned pass-filter reads per read group")
  ggplot_object <-
    ggplot(data = combined_metrics_read_group)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Aligned Pass-Filter Reads per Read Group")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = CATEGORY, y = PF_READS, colour = READ_GROUP))
  ggplot_object <-
    ggplot_object + guides(colour = guide_legend(nrow = 24L))
  ggplot_object <-
    ggplot_object + theme(legend.text = element_text(size = rel(x = 0.5)))
  for (graphics_format in graphics_formats) {
    ggsave(
      filename = paste(
        prefix_summary,
        paste("alignment_absolute_read_group", graphics_format, sep = "."),
        sep = "_"
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
    ggplot(data = combined_metrics_sample)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Aligned Pass-Filter Reads per Sample")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = CATEGORY, y = PCT_PF_READS_ALIGNED, colour = SAMPLE))
  ggplot_object <-
    ggplot_object + guides(colour = guide_legend(nrow = 24L))
  for (graphics_format in graphics_formats) {
    ggsave(
      filename = paste(
        prefix_summary,
        paste("alignment_percentage_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = plot_width_sample,
      height = argument_list$plot_height
    )
  }
  rm(graphics_format, ggplot_object)
  
  # Plot the percentage of aligned pass-filter reads per read group -------
  
  
  message("Plotting the percentage of aligned pass-filter reads per read group")
  ggplot_object <-
    ggplot(data = combined_metrics_read_group)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Aligned Pass-Filter Reads per Read Group")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = CATEGORY, y = PCT_PF_READS_ALIGNED, colour = READ_GROUP))
  ggplot_object <-
    ggplot_object + guides(colour = guide_legend(nrow = 24L))
  ggplot_object <-
    ggplot_object + theme(legend.text = element_text(size = rel(x = 0.5)))
  for (graphics_format in graphics_formats) {
    ggsave(
      filename = paste(
        prefix_summary,
        paste("alignment_percentage_read_group", graphics_format, sep = "."),
        sep = "_"
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

# Picard Hybrid Selection Metrics -----------------------------------------


# Process Picard Hybrid Selection Metrics reports.

message("Processing Picard Hybrid Selection Metrics reports for sample:")
combined_metrics_sample <- NULL
combined_metrics_read_group <- NULL

file_names <-
  list.files(pattern = "^variant_calling_diagnose_sample_.*_hybrid_selection_metrics.tsv$")
for (file_name in file_names) {
  sample_name <-
    gsub(pattern = "^variant_calling_diagnose_sample_(.*?)_hybrid_selection_metrics.tsv$",
         replacement = "\\1",
         x = file_name)
  message(paste0("  ", sample_name))
  # Picard Tools added a histogram section that needs excluding from parsing.
  # Find the lines starting with "## METRICS CLASS" and "## HISTOGRAM" and read that many lines.
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
    read.table(
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
  # To support numeric sample names the read.table(stringsAsFactors = FALSE) is turned off.
  # Convert SAMPLE, LIBRARY and READ_GROUP into character vectors.
  picard_metrics_total$SAMPLE <-
    as.character(x = picard_metrics_total$SAMPLE)
  picard_metrics_total$LIBRARY <-
    as.character(x = picard_metrics_total$LIBRARY)
  picard_metrics_total$READ_GROUP <-
    as.character(x = picard_metrics_total$READ_GROUP)
  
  # The Picard Hybrid Selection Metrics report has changed format through versions.
  # Column PCT_TARGET_BASES_1X was added at a later stage.
  if (is.null(x = picard_metrics_total$PCT_TARGET_BASES_1X)) {
    picard_metrics_total$PCT_TARGET_BASES_1X <- 0.0
  }
  
  if (is.null(x = picard_metrics_total$MAX_TARGET_COVERAGE)) {
    picard_metrics_total$MAX_TARGET_COVERAGE <- 0L
  }
  
  # Select only rows showing the SAMPLE summary, i.e. showing SAMPLE, but no LIBRARY and READ_GROUP information.
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
  
  # Select only rows showing READ_GROUP summary, i.e. showing READ_GROUP information.
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
  # Manually convert BAIT_SET and SAMPLE columns into factors, which are handy for plotting.
  combined_metrics_sample$BAIT_SET <-
    as.factor(x = combined_metrics_sample$BAIT_SET)
  combined_metrics_sample$SAMPLE <-
    as.factor(x = combined_metrics_sample$SAMPLE)
  write.table(
    x = combined_metrics_sample,
    file = paste(prefix_summary, "hybrid_metrics_sample.tsv", sep = "_"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  
  # Sort the data frame by READ_GROUP.
  combined_metrics_read_group <-
    combined_metrics_read_group[order(combined_metrics_read_group$READ_GROUP), ]
  # Manually convert BAIT_SET and READ_GROUP columns into factors, which are handy for plotting.
  combined_metrics_read_group$BAIT_SET <-
    as.factor(x = combined_metrics_read_group$BAIT_SET)
  combined_metrics_read_group$READ_GROUP <-
    as.factor(x = combined_metrics_read_group$READ_GROUP)
  write.table(
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
  
  # Plot the percentage of unique pass-filter reads per sample ------------
  
  
  message("Plotting the percentage of unique pass-filter reads per sample")
  ggplot_object <-
    ggplot(data = combined_metrics_sample)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Unique Pass-Filter Reads per Sample")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = SAMPLE, y = PCT_PF_UQ_READS))
  ggplot_object <-
    ggplot_object + guides(colour = guide_legend(nrow = 24L))
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      size = rel(x = 0.8)
    ))
  for (graphics_format in graphics_formats) {
    ggsave(
      filename = paste(
        prefix_summary,
        paste("hybrid_unique_percentage_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = plot_width_sample,
      height = argument_list$plot_height
    )
  }
  rm(graphics_format, ggplot_object)
  
  # Plot the percentage of unique pass-filter reads per read group --------
  
  
  message("Plotting the percentage of unique pass-filter reads per read group")
  ggplot_object <-
    ggplot(data = combined_metrics_read_group)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Unique Pass-Filter Reads per Read Group")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = READ_GROUP, y = PCT_PF_UQ_READS, shape = BAIT_SET))
  # ggplot2 only adds six shapes automatically. Since there may be more, add them manually.
  ggplot_object <-
    ggplot_object + scale_shape_manual(values = seq_len(length.out = nlevels(x = combined_metrics_read_group$BAIT_SET)))
  ggplot_object <-
    ggplot_object + guides(colour = guide_legend(nrow = 24L))
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      size = rel(x = 0.5)
    ))
  for (graphics_format in graphics_formats) {
    ggsave(
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
      width = plot_width_read_group,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)
  
  # Plot the mean target coverage per sample ------------------------------
  
  
  message("Plotting the mean target coverage per sample")
  ggplot_object <-
    ggplot(data = combined_metrics_sample)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Mean Target Coverage per Sample")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = SAMPLE, y = MEAN_TARGET_COVERAGE))
  ggplot_object <-
    ggplot_object + guides(colour = guide_legend(nrow = 24L))
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      size = rel(x = 0.8)
    ))
  for (graphics_format in graphics_formats) {
    ggsave(
      filename = paste(
        prefix_summary,
        paste("hybrid_target_coverage_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = plot_width_sample,
      height = argument_list$plot_height
    )
  }
  rm(graphics_format, ggplot_object)
  
  # Plot the mean target coverage per read group --------------------------
  
  
  message("Plotting the mean target coverage per read group")
  ggplot_object <-
    ggplot(data = combined_metrics_read_group)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Mean Target Coverage per Read Group")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = READ_GROUP, y = MEAN_TARGET_COVERAGE, shape = BAIT_SET))
  # ggplot2 only adds six shapes automatically. Since there may be more, add them manually.
  ggplot_object <-
    ggplot_object + scale_shape_manual(values = seq_len(length.out = nlevels(x = combined_metrics_read_group$BAIT_SET)))
  ggplot_object <-
    ggplot_object + guides(colour = guide_legend(nrow = 24L))
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      size = rel(x = 0.5)
    ))
  for (graphics_format in graphics_formats) {
    ggsave(
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
      width = plot_width_read_group,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)
  
  # Plot the percentage of exluded bases per sample -----------------------
  
  
  message("Plotting the percentage of excluded bases per sample")
  plotting_frame <- melt(
    data = combined_metrics_sample,
    id.vars = c("SAMPLE", "BAIT_SET"),
    measure.vars = c(
      "PCT_EXC_DUPE",
      "PCT_EXC_MAPQ",
      "PCT_EXC_BASEQ",
      "PCT_EXC_OVERLAP",
      "PCT_EXC_OFF_TARGET"
    ),
    variable.name = "EXCLUDED",
    value.name = "fraction"
  )
  
  ggplot_object <- ggplot(data = plotting_frame)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Excluded Bases per Sample")
  # ggplot_object <-
  #   ggplot_object + geom_point(mapping = aes(x = SAMPLE,
  #                                            y = fraction,
  #                                            colour = EXCLUDED))
  ggplot_object <-
    ggplot_object + geom_col(mapping = aes(x = SAMPLE, y = fraction, fill = EXCLUDED),
                             alpha = I(1 / 3))
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      size = rel(x = 0.8)
    ))
  for (graphics_format in graphics_formats) {
    ggsave(
      filename = paste(
        prefix_summary,
        paste("hybrid_excluded_bases_sample",
              graphics_format,
              sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = plot_width_sample,
      height = argument_list$plot_height
    )
  }
  rm(graphics_format, ggplot_object, plotting_frame)
  
  # Plot the percentage of exluded bases per read group -------------------
  
  
  message("Plotting the percentage of excluded bases per read group")
  plotting_frame <- melt(
    data = combined_metrics_read_group,
    id.vars = c("READ_GROUP", "BAIT_SET"),
    measure.vars = c(
      "PCT_EXC_DUPE",
      "PCT_EXC_MAPQ",
      "PCT_EXC_BASEQ",
      "PCT_EXC_OVERLAP",
      "PCT_EXC_OFF_TARGET"
    ),
    variable.name = "EXCLUDED",
    value.name = "fraction"
  )
  
  ggplot_object <- ggplot(data = plotting_frame)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Excluded Bases per Read Group")
  # ggplot_object <-
  #   ggplot_object + geom_point(mapping = aes(
  #     x = READ_GROUP,
  #     y = fraction,
  #     colour = EXCLUDED,
  #     shape = BAIT_SET
  #   ))
  ggplot_object <-
    ggplot_object + geom_col(
      mapping = aes(
        x = READ_GROUP,
        y = fraction,
        fill = EXCLUDED,
        colour = BAIT_SET
      ),
      alpha = I(1 / 3)
    )
  # ggplot2 only adds six shapes automatically. Since there may be more, add them manually.
  ggplot_object <-
    ggplot_object + scale_shape_manual(values = seq_len(length.out = nlevels(x = combined_metrics_read_group$BAIT_SET)))
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      size = rel(x = 0.5)
    ))
  for (graphics_format in graphics_formats) {
    ggsave(
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
      width = plot_width_sample,
      height = argument_list$plot_height
    )
  }
  rm(graphics_format, ggplot_object, plotting_frame)
  
  # Plot target coverage levels per sample --------------------------------
  
  
  # Plot PCT_TARGET_BASES_1X, PCT_TARGET_BASES_2X, PCT_TARGET_BASES_10X, PCT_TARGET_BASES_20X,
  # PCT_TARGET_BASES_30X, PCT_TARGET_BASES_40X, PCT_TARGET_BASES_50X, PCT_TARGET_BASES_100X per sample.
  message("Plotting the coverage levels per sample")
  plotting_frame <- melt(
    data = combined_metrics_sample,
    id.vars = c("SAMPLE", "BAIT_SET"),
    measure.vars = c(
      "PCT_TARGET_BASES_1X",
      "PCT_TARGET_BASES_2X",
      "PCT_TARGET_BASES_10X",
      "PCT_TARGET_BASES_20X",
      "PCT_TARGET_BASES_30X",
      "PCT_TARGET_BASES_40X",
      "PCT_TARGET_BASES_50X",
      "PCT_TARGET_BASES_100X"
    ),
    variable.name = "COVERAGE",
    value.name = "fraction"
  )
  
  ggplot_object <- ggplot(data = plotting_frame)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Coverage Levels per Sample")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = SAMPLE,
                                             y = fraction,
                                             colour = COVERAGE))
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      size = rel(x = 0.8)
    ))
  for (graphics_format in graphics_formats) {
    ggsave(
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
      width = plot_width_sample,
      height = argument_list$plot_height
    )
  }
  rm(graphics_format, ggplot_object, plotting_frame)
  
  # Plot target coverage levels per read group ----------------------------
  
  
  # Plot PCT_TARGET_BASES_1X, PCT_TARGET_BASES_2X, PCT_TARGET_BASES_10X, PCT_TARGET_BASES_20X,
  # PCT_TARGET_BASES_30X, PCT_TARGET_BASES_40X, PCT_TARGET_BASES_50X, PCT_TARGET_BASES_100X per read group.
  message("Plotting the coverage levels per read group")
  plotting_frame <- melt(
    data = combined_metrics_read_group,
    id.vars = c("READ_GROUP", "BAIT_SET"),
    measure.vars = c(
      "PCT_TARGET_BASES_1X",
      "PCT_TARGET_BASES_2X",
      "PCT_TARGET_BASES_10X",
      "PCT_TARGET_BASES_20X",
      "PCT_TARGET_BASES_30X",
      "PCT_TARGET_BASES_40X",
      "PCT_TARGET_BASES_50X",
      "PCT_TARGET_BASES_100X"
    ),
    variable.name = "COVERAGE",
    value.name = "fraction"
  )
  
  ggplot_object <- ggplot(data = plotting_frame)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Coverage Levels per Read Group")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(
      x = READ_GROUP,
      y = fraction,
      colour = COVERAGE,
      shape = BAIT_SET
    ))
  # ggplot2 only adds six shapes automatically. Since there may be more, add them manually.
  ggplot_object <-
    ggplot_object + scale_shape_manual(values = seq_len(length.out = nlevels(x = combined_metrics_read_group$BAIT_SET)))
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      size = rel(x = 0.5)
    ))
  for (graphics_format in graphics_formats) {
    ggsave(
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
      width = plot_width_sample,
      height = argument_list$plot_height
    )
  }
  rm(graphics_format, ggplot_object, plotting_frame)
  
  # Plot the nominal coverage per sample ----------------------------------
  
  
  # Plot the nominal coverage (i.e. PF_BASES_ALIGNED / TARGET_TERRITORY) per sample.
  
  message("Plotting the nominal coverage per sample")
  plotting_frame <-
    combined_metrics_sample[, c("SAMPLE",
                                "READ_GROUP",
                                "BAIT_SET",
                                "PF_BASES_ALIGNED",
                                "TARGET_TERRITORY")]
  plotting_frame$NOMINAL_COVERAGE <-
    plotting_frame$PF_BASES_ALIGNED / plotting_frame$TARGET_TERRITORY
  
  ggplot_object <- ggplot(data = plotting_frame)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Nominal Coverage per Sample")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = SAMPLE,
                                             y = NOMINAL_COVERAGE))
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      size = rel(x = 0.8)
    ))
  for (graphics_format in graphics_formats) {
    ggsave(
      filename = paste(
        prefix_summary,
        paste("hybrid_nominal_coverage_sample",
              graphics_format,
              sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object, plotting_frame)
  
  # Plot the nominal coverage per read group ------------------------------
  
  
  # Plot the nominal coverage (i.e. PF_BASES_ALIGNED / TARGET_TERRITORY) per read group.
  
  message("Plotting the nominal coverage per read group")
  plotting_frame <-
    combined_metrics_read_group[, c("SAMPLE",
                                    "READ_GROUP",
                                    "BAIT_SET",
                                    "PF_BASES_ALIGNED",
                                    "TARGET_TERRITORY")]
  plotting_frame$NOMINAL_COVERAGE <-
    plotting_frame$PF_BASES_ALIGNED / plotting_frame$TARGET_TERRITORY
  
  ggplot_object <- ggplot(data = plotting_frame)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Nominal Coverage per Read Group")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = READ_GROUP,
                                             y = NOMINAL_COVERAGE,
                                             shape = BAIT_SET))
  # ggplot2 only adds six shapes automatically. Since there may be more, add them manually.
  ggplot_object <-
    ggplot_object + scale_shape_manual(values = seq_len(length.out = nlevels(x = combined_metrics_read_group$BAIT_SET)))
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      size = rel(x = 0.5)
    ))
  for (graphics_format in graphics_formats) {
    ggsave(
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
      width = plot_width_read_group,
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
  row.names = list.files(pattern = "^variant_calling_diagnose_sample_.*_non_callable_summary.tsv$")
)
combined_metrics_sample$file_name <-
  row.names(x = combined_metrics_sample)
combined_metrics_sample$exon_path <-
  character(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$exon_number <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$exon_width <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$exon_flank_width <-
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
combined_metrics_sample$target_number_constrained <-
  integer(length = nrow(x = combined_metrics_sample))
combined_metrics_sample$target_width_constrained <-
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
  # "PASS" according to documentation, but should be "CALLABLE" in practice
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
         x = combined_metrics_sample[i, "file_name"])
  message(paste0("  ", sample_name))
  non_callable_metrics_sample <-
    read.table(
      file = combined_metrics_sample[i, "file_name"],
      header = TRUE,
      colClasses = c(
        "exon_path" = "character",
        # Raw exons before filtering processed_transcripts, etc. via 'tag "basic";'
        # "exon_number_raw" = "integer",
        # "exon_width_raw" = "integer",
        "exon_number" = "integer",
        "exon_width" = "integer",
        "exon_flank_width" = "integer",
        "transcribed_number" = "integer",
        "transcribed_width" = "integer",
        "target_path" = "character",
        "target_number_raw" = "integer",
        "target_width_raw" = "integer",
        "target_number_constrained" = "integer",
        "target_width_constrained" = "integer",
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
      # With the exception of the "file_name" component, set all components undefined in the
      # sample-specific data frame to NA in the combined data frame.
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
  # Convert the sample_name column into factors, which come more handy for plotting.
  combined_metrics_sample$sample_name <-
    as.factor(x = combined_metrics_sample$sample_name)
  
  # Adjust the plot width according to batches of 24 samples or read groups.
  plot_width_sample <- argument_list$plot_width + (ceiling(x = (
    nlevels(x = combined_metrics_sample$sample_name) / 24L
  )) - 1L) * argument_list$plot_width * 0.25
  
  write.table(
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
    target_number_constrained = combined_metrics_sample$target_number_constrained
  )
  # Only columns that begin with "^non_callable_number_constrained\." are required, the remainder is the mapping status.
  column_names <- names(x = combined_metrics_sample)
  mapping_status <-
    gsub(pattern = "^non_callable_number_constrained\\.(.*)$",
         replacement = "\\1",
         x = column_names)
  # Extract only those columns which had a match and set the mapping status as their new name.
  for (i in which(x = grepl(pattern = "^non_callable_number_constrained\\.", x = column_names))) {
    plotting_frame[, mapping_status[i]] <-
      combined_metrics_sample[, column_names[i]]
  }
  rm(i, mapping_status, column_names)
  
  # Now, melt the data frame, but keep sample_name and target_width_constrained as identifiers.
  plotting_frame <- melt(
    data = plotting_frame,
    id.vars = c("sample_name", "target_number_constrained"),
    variable.name = "mapping_status",
    value.name = "number"
  )
  # Calculate the fractions on the basis of the constrained target number.
  plotting_frame[, "fraction"] <-
    plotting_frame$number / plotting_frame$target_number_constrained
  # For the moment remove lines with "TOTAL".
  plotting_frame <-
    plotting_frame[plotting_frame$mapping_status != "TOTAL", ]
  
  ggplot_object <- ggplot(data = plotting_frame)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Number of Non-Callable Loci per Sample")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = sample_name, y = number, colour = mapping_status))
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      size = rel(x = 0.8)
    ))
  for (graphics_format in graphics_formats) {
    ggsave(
      filename = paste(
        prefix_summary,
        paste("non_callable_absolute_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object, plotting_frame)
  
  # Plot the fraction of non-callable loci per sample ---------------------
  
  
  message("Plotting the fraction of non-callable loci per sample")
  # Reorganise the combined_metrics_sample data frame for plotting non-callable target widths.
  plotting_frame <- data.frame(
    sample_name = combined_metrics_sample$sample_name,
    target_width_constrained = combined_metrics_sample$target_width_constrained
  )
  # Only columns that begin with "^non_callable_width_constrained\." are required, the remainder is the mapping status.
  column_names <- names(x = combined_metrics_sample)
  mapping_status <-
    gsub(pattern = "^non_callable_width_constrained\\.(.*)$",
         replacement = "\\1",
         x = column_names)
  # Extract only those columns which had a match and set the mapping status as their new name.
  for (i in which(x = grepl(pattern = "^non_callable_width_constrained\\.", x = column_names))) {
    plotting_frame[, mapping_status[i]] <-
      combined_metrics_sample[, column_names[i]]
  }
  rm(i, mapping_status, column_names)
  
  # Now, melt the data frame, but keep sample_name and target_width_constrained as identifiers.
  plotting_frame <- melt(
    data = plotting_frame,
    id.vars = c("sample_name", "target_width_constrained"),
    variable.name = "mapping_status",
    value.name = "width"
  )
  # Calculate the fractions on the basis of the constrained target width.
  plotting_frame[, "fraction"] <-
    plotting_frame$width / plotting_frame$target_width_constrained
  # For the moment remove lines with "TOTAL".
  plotting_frame <-
    plotting_frame[plotting_frame$mapping_status != "TOTAL", ]
  
  ggplot_object <- ggplot(data = plotting_frame)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Fraction of Non-Callable Loci per Sample")
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = sample_name, y = fraction, colour = mapping_status))
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      size = rel(x = 0.8)
    ))
  for (graphics_format in graphics_formats) {
    ggsave(
      filename = paste(
        prefix_summary,
        paste("non_callable_percentage_sample", graphics_format, sep = "."),
        sep = "_"
      ),
      plot = ggplot_object,
      width = plot_width_sample,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object, plotting_frame)
  
  rm(plot_width_sample)
}
rm(combined_metrics_sample)

rm(prefix_summary, argument_list, graphics_formats)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

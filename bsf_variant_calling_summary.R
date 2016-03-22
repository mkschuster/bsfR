#! /usr/bin/env Rscript
#
# BSF R script to summarise a variant caling analysis. Picard Duplication 
# Metrics, Picard Alignment Summary Metrics and Picard Hybrid Selection Metrics 
# reports are read for each sample and plotted at the aliquot or sample level.
#
#
# Copyright 2013 -2015 Michael K. Schuster
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
suppressPackageStartupMessages(expr = library(package = "ggplot2"))

# Get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults.

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
      opt_str = c("--plot-width"),
      default = 7,
      dest = "plot_width",
      help = "Plot width in inches",
      type = "numeric"
    ),
    make_option(
      opt_str = c("--plot-height"),
      default = 7,
      dest = "plot_height",
      help = "Plot height in inches",
      type = "numeric"
    )
  )
))

# Process Picard Duplication Metrics reports.

combined_metrics_sample <- NULL

file_names <-
  list.files(pattern = "variant_calling_process_sample_.*_duplicate_metrics.tsv")
for (file_name in file_names) {
  sample_name <-
    gsub(pattern = "variant_calling_process_sample_(.*?)_duplicate_metrics.tsv",
         replacement = "\\1",
         x = file_name)
  picard_metrics_sample <-
    read.table(
      file = file_name,
      header = TRUE,
      sep = "\t",
      nrows = 1,
      skip = 6
    )
  # Add the sample name, which is not part of the Picard report.
  picard_metrics_sample$SAMPLE <- sample_name
  
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
  # Adjust the plot width according to batches of 24 samples or aliquots.
  plot_width <-
    argument_list$plot_width + (ceiling(x = nlevels(x = combined_metrics_sample$SAMPLE) / 24) - 1) * argument_list$plot_width * 0.3
  
  # Plot the absolute number of aligned pass-filter reads per sample.
  plot_object <-
    ggplot(data = combined_metrics_sample)
  plot_object <-
    plot_object + ggtitle(label = "Percent Duplication per Sample")
  plot_object <-
    plot_object + geom_point(mapping = aes(x = SAMPLE, y = PERCENT_DUPLICATION))
  plot_object <-
    plot_object + guides(color = guide_legend(nrow = 24))
  ggsave(
    filename = "variant_calling_summary_duplication_sample.png",
    plot = plot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  ggsave(
    filename = "variant_calling_summary_duplication_sample.pdf",
    plot = plot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  rm(plot_object)
  
  rm(plot_width)
}
rm(combined_metrics_sample)

# Process Picard Alignment Summary Metrics reports.

combined_metrics_sample <- NULL
combined_metrics_aliquot <- NULL

file_names <-
  list.files(pattern = "variant_calling_process_sample_.*_alignment_summary_metrics.tsv")
for (file_name in file_names) {
  sample_name <-
    gsub(pattern = "variant_calling_process_sample_(.*?)_alignment_summary_metrics.tsv",
         replacement = "\\1",
         x = file_name)
  picard_metrics_total <-
    read.table(file = file_name,
               header = TRUE,
               sep = "\t")
  
  # Modify the row names so that the names do not clash.
  # row.names(picard_metrics_total) <- paste(sample_name, row.names(x = picard_metrics_total), sep = "_")
  
  # Select only rows showing the SAMPLE summary, i.e. showing SAMPLE, but no LIBRARY and READ_GROUP information.
  picard_metrics_sample <-
    picard_metrics_total[(picard_metrics_total$SAMPLE != "") &
                           (picard_metrics_total$LIBRARY == "") &
                           (picard_metrics_total$READ_GROUP == ""),]
  if (is.null(x = combined_metrics_sample)) {
    combined_metrics_sample <- picard_metrics_sample
  } else {
    combined_metrics_sample <-
      rbind(combined_metrics_sample, picard_metrics_sample)
  }
  rm(picard_metrics_sample)
  
  # Select only rows showing READ_GROUP summary, i.e. showing READ_GROUP information.
  picard_metrics_aliquot <-
    picard_metrics_total[(picard_metrics_total$READ_GROUP != ""),]
  if (is.null(x = combined_metrics_aliquot)) {
    combined_metrics_aliquot <- picard_metrics_aliquot
  } else {
    combined_metrics_aliquot <-
      rbind(combined_metrics_aliquot, picard_metrics_aliquot)
  }
  rm(picard_metrics_aliquot)
  
  rm(sample_name, picard_metrics_total)
}
rm(file_name, file_names)

if (!is.null(x = combined_metrics_sample)) {
  # Add an additional LABEL factor column defined as a concatenation of SAMPLE or READ_GROUP and CATEGORY.
  combined_metrics_sample$LABEL <-
    as.factor(x = paste(
      combined_metrics_sample$SAMPLE,
      combined_metrics_sample$CATEGORY,
      sep =
        "_"
    ))
  combined_metrics_aliquot$LABEL <-
    as.factor(
      x = paste(
        combined_metrics_aliquot$READ_GROUP,
        combined_metrics_aliquot$CATEGORY,
        sep =
          "_"
      )
    )
  # Add an additional ALIQUOT factor column defined as a concatenation of SAMPLE and READ_GROUP.
  combined_metrics_aliquot$ALIQUOT <-
    as.factor(
      x = paste(
        combined_metrics_aliquot$SAMPLE,
        combined_metrics_aliquot$READ_GROUP,
        sep =
          "_"
      )
    )
  
  # Adjust the plot width according to batches of 24 samples or aliquots.
  plot_width_aliquot <-
    argument_list$plot_width + (ceiling(x = nlevels(x = combined_metrics_aliquot$ALIQUOT) / 24) - 1) * argument_list$plot_width * 0.35
  plot_width_sample <-
    argument_list$plot_width + (ceiling(x = nlevels(x = combined_metrics_sample$SAMPLE) / 24) - 1) * argument_list$plot_width * 0.25
  
  # Plot the absolute number of aligned pass-filter reads per sample.
  plot_object <-
    ggplot(data = combined_metrics_sample)
  plot_object <-
    plot_object + ggtitle(label = "Aligned Pass-Filter Reads per Sample")
  plot_object <-
    plot_object + geom_point(mapping = aes(x = CATEGORY, y = PF_READS, color = SAMPLE))
  plot_object <-
    plot_object + guides(color = guide_legend(nrow = 24))
  ggsave(
    filename = "variant_calling_summary_alignment_reads_sample.png",
    plot = plot_object,
    width = plot_width_sample,
    height = argument_list$plot_height
  )
  ggsave(
    filename = "variant_calling_summary_alignment_reads_sample.pdf",
    plot = plot_object,
    width = plot_width_sample,
    height = argument_list$plot_height
  )
  rm(plot_object)
  
  # Plot the absolute number of aligned pass-filter reads per aliquot.
  plot_object <-
    ggplot(data = combined_metrics_aliquot)
  plot_object <-
    plot_object + ggtitle(label = "Aligned Pass-Filter Reads per Aliquot")
  plot_object <-
    plot_object + geom_point(mapping = aes(x = CATEGORY, y = PF_READS, color = ALIQUOT))
  plot_object <-
    plot_object + guides(color = guide_legend(nrow = 24))
  plot_object <-
    plot_object + theme(legend.text = element_text(size = rel(x = 0.5)))
  ggsave(
    filename = "variant_calling_summary_alignment_reads_read_group.png",
    plot = plot_object,
    width = plot_width_aliquot,
    height = argument_list$plot_height,
    limitsize = FALSE
  )
  ggsave(
    filename = "variant_calling_summary_alignment_reads_read_group.pdf",
    plot = plot_object,
    width = plot_width_aliquot,
    height = argument_list$plot_height,
    limitsize = FALSE
  )
  rm(plot_object)
  
  # Plot the percentage of aligned pass-filter reads per sample.
  plot_object <-
    ggplot(data = combined_metrics_sample)
  plot_object <-
    plot_object + ggtitle(label = "Aligned Pass-Filter Reads per Sample")
  plot_object <-
    plot_object + geom_point(mapping = aes(x = CATEGORY, y = PCT_PF_READS_ALIGNED, color = SAMPLE))
  plot_object <-
    plot_object + guides(color = guide_legend(nrow = 24))
  ggsave(
    filename = "variant_calling_summary_alignment_percentage_sample.png",
    plot = plot_object,
    width = plot_width_sample,
    height = argument_list$plot_height
  )
  ggsave(
    filename = "variant_calling_summary_alignment_percentage_sample.pdf",
    plot = plot_object,
    width = plot_width_sample,
    height = argument_list$plot_height
  )
  rm(plot_object)
  
  # Plot the percentage of aligned pass-filter reads per aliquot.
  plot_object <-
    ggplot(data = combined_metrics_aliquot)
  plot_object <-
    plot_object + ggtitle(label = "Aligned Pass-Filter Reads per Aliquot")
  plot_object <-
    plot_object + geom_point(mapping = aes(x = CATEGORY, y = PCT_PF_READS_ALIGNED, color = ALIQUOT))
  plot_object <-
    plot_object + guides(color = guide_legend(nrow = 24))
  plot_object <-
    plot_object + theme(legend.text = element_text(size = rel(x = 0.5)))
  ggsave(
    filename = "variant_calling_summary_alignment_percentage_read_group.png",
    plot = plot_object,
    width = plot_width_aliquot,
    height = argument_list$plot_height,
    limitsize = FALSE
  )
  ggsave(
    filename = "variant_calling_summary_alignment_percentage_read_group.pdf",
    plot = plot_object,
    width = plot_width_aliquot,
    height = argument_list$plot_height,
    limitsize = FALSE
  )
  rm(plot_object)
  
  rm(plot_width_aliquot, plot_width_sample)
}
rm(combined_metrics_aliquot, combined_metrics_sample)

# Process Picard Hybrid Selection Metrics reports.

combined_metrics_sample <- NULL
combined_metrics_aliquot <- NULL

file_names <-
  list.files(pattern = "variant_calling_diagnose_sample_.*_hybrid_selection_metrics.tsv")
for (file_name in file_names) {
  sample_name <-
    gsub(pattern = "variant_calling_diagnose_sample_(.*?)_hybrid_selection_metrics.tsv",
         replacement = "\\1",
         x = file_name)
  # message(paste("Processing hybrid selection metrics file", file_name, "for sample", sample_name))
  picard_metrics_total <-
    read.table(file = file_name,
               header = TRUE,
               sep = "\t")
  
  # Modify the row names so that the names do not clash.
  # row.names(picard_metrics_total) <- paste(sample_name, row.names(x = picard_metrics_total), sep = "_")
  
  # Select only rows showing the SAMPLE summary, i.e. showing SAMPLE, but no LIBRARY and READ_GROUP information.
  picard_metrics_sample <-
    picard_metrics_total[(picard_metrics_total$SAMPLE != "") &
                           (picard_metrics_total$LIBRARY == "") &
                           (picard_metrics_total$READ_GROUP == ""),]
  if (is.null(x = combined_metrics_sample)) {
    combined_metrics_sample <- picard_metrics_sample
  } else {
    combined_metrics_sample <-
      rbind(combined_metrics_sample, picard_metrics_sample)
  }
  rm(picard_metrics_sample)
  
  # Select only rows showing READ_GROUP summary, i.e. showing READ_GROUP information.
  picard_metrics_aliquot <-
    picard_metrics_total[(picard_metrics_total$READ_GROUP != ""),]
  if (is.null(x = combined_metrics_aliquot)) {
    combined_metrics_aliquot <- picard_metrics_aliquot
  } else {
    combined_metrics_aliquot <-
      rbind(combined_metrics_aliquot, picard_metrics_aliquot)
  }
  rm(picard_metrics_aliquot)
  
  rm(sample_name, picard_metrics_total)
}
rm(file_name, file_names)

# The Picard Hybrid Selection Metrics is currently optional.

if (!is.null(x = combined_metrics_sample)) {
  # Add an additional factor column ALIQUOT defined as a concatenation of SAMPLE and READ_GROUP.
  combined_metrics_aliquot$ALIQUOT <-
    as.factor(
      x = paste(
        combined_metrics_aliquot$SAMPLE,
        combined_metrics_aliquot$READ_GROUP,
        sep =
          "_"
      )
    )
  
  # Adjust the plot width according to batches of 24 samples or aliquots.
  plot_width_aliquot <- argument_list$plot_width + (ceiling(x = (
    nlevels(x = combined_metrics_aliquot$ALIQUOT) / 24
  )) - 1) * argument_list$plot_width * 0.25
  plot_width_sample <- argument_list$plot_width + (ceiling(x = (
    nlevels(x = combined_metrics_sample$SAMPLE) / 24
  )) - 1) * argument_list$plot_width * 0.25
  
  # Plot the percentage of unique pass-filter reads per sample.
  plot_object <-
    ggplot(data = combined_metrics_sample)
  plot_object <-
    plot_object + ggtitle(label = "Unique Pass-Filter Reads per Sample")
  plot_object <-
    plot_object + geom_point(mapping = aes(x = SAMPLE, y = PCT_PF_UQ_READS))
  plot_object <-
    plot_object + guides(color = guide_legend(nrow = 24))
  plot_object <-
    plot_object + theme(axis.text.x = element_text(
      angle = -90,
      hjust = 0,
      size = rel(x = 0.8)
    ))
  ggsave(
    filename = "variant_calling_summary_hybrid_unique_percentage_sample.png",
    plot = plot_object,
    width = plot_width_sample,
    height = argument_list$plot_height
  )
  ggsave(
    filename = "variant_calling_summary_hybrid_unique_percentage_sample.pdf",
    plot = plot_object,
    width = plot_width_sample,
    height = argument_list$plot_height
  )
  rm(plot_object)
  
  # Plot the percentage of unique pass-filter reads per aliquot.
  plot_object <-
    ggplot(data = combined_metrics_aliquot)
  plot_object <-
    plot_object + ggtitle(label = "Unique Pass-Filter Reads per Aliquot")
  plot_object <-
    plot_object + geom_point(mapping = aes(x = ALIQUOT, y = PCT_PF_UQ_READS))
  plot_object <-
    plot_object + guides(color = guide_legend(nrow = 24))
  plot_object <-
    plot_object + theme(axis.text.x = element_text(
      angle = -90,
      hjust = 0,
      size = rel(x = 0.8)
    ))
  ggsave(
    filename = "variant_calling_summary_hybrid_unique_percentage_read_group.png",
    plot = plot_object,
    width = plot_width_aliquot,
    height = argument_list$plot_height
  )
  ggsave(
    filename = "variant_calling_summary_hybrid_unique_percentage_read_group.pdf",
    plot = plot_object,
    width = plot_width_aliquot,
    height = argument_list$plot_height
  )
  rm(plot_object)
  
  # Plot the mean target coverage per sample.
  plot_object <-
    ggplot(data = combined_metrics_sample)
  plot_object <-
    plot_object + ggtitle(label = "Mean Target Coverage per Sample")
  plot_object <-
    plot_object + geom_point(mapping = aes(x = SAMPLE, y = MEAN_TARGET_COVERAGE))
  plot_object <-
    plot_object + guides(color = guide_legend(nrow = 24))
  plot_object <-
    plot_object + theme(axis.text.x = element_text(
      angle = -90,
      hjust = 0,
      size = rel(x = 0.8)
    ))
  ggsave(
    filename = "variant_calling_summary_hybrid_target_coverage_sample.png",
    plot = plot_object,
    width = plot_width_sample,
    height = argument_list$plot_height
  )
  ggsave(
    filename = "variant_calling_summary_hybrid_target_coverage_sample.pdf",
    plot = plot_object,
    width = plot_width_sample,
    height = argument_list$plot_height
  )
  rm(plot_object)
  
  # Plot the mean target coverage per aliquot.
  plot_object <-
    ggplot(data = combined_metrics_aliquot)
  plot_object <-
    plot_object + ggtitle(label = "Mean Target Coverage per Aliquot")
  plot_object <-
    plot_object + geom_point(mapping = aes(x = ALIQUOT, y = MEAN_TARGET_COVERAGE))
  plot_object <-
    plot_object + guides(color = guide_legend(nrow = 24))
  plot_object <-
    plot_object + theme(axis.text.x = element_text(
      angle = -90,
      hjust = 0,
      size = rel(x = 0.8)
    ))
  ggsave(
    filename = "variant_calling_summary_hybrid_target_coverage_read_group.png",
    plot = plot_object,
    width = plot_width_aliquot,
    height = argument_list$plot_height
  )
  ggsave(
    filename = "variant_calling_summary_hybrid_target_coverage_read_group.pdf",
    plot = plot_object,
    width = plot_width_aliquot,
    height = argument_list$plot_height
  )
  rm(plot_object)
  
  rm(plot_width_aliquot, plot_width_sample)
}
rm(combined_metrics_aliquot, combined_metrics_sample)
rm(argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

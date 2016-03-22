#! /usr/bin/env Rscript
#
# Copyright 2013 -2015 Michael K. Schuster
#
# Biomedical Sequencing Facility (BSF), part of the genomics core facility
# of the Research Center for Molecular Medicine (CeMM) of the
# Austrian Academy of Sciences and the Medical University of Vienna (MUW).
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
suppressPackageStartupMessages(expr = library(package = "DiffBind"))

# Get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,

argument_list <- parse_args(object = OptionParser(
  option_list = list(
    make_option(
      opt_str = c("-v", "--verbose"),
      action = "store_true",
      default = TRUE,
      help = "Print extra output [default]",
      type = "logical"
    ),
    make_option(
      opt_str = c("-q", "--quietly"),
      action = "store_false",
      default = FALSE,
      dest = "verbose",
      help = "Print little output",
      type = "logical"
    ),
    make_option(
      opt_str = c("-f", "--factor"),
      dest = "factor",
      help = "ChIP factor",
      type = "character"
    ),
    make_option(
      opt_str = c("-s", "--sample-annotation"),
      dest = "sample_annotation",
      help = "Sample annotation sheet",
      type = "character"
    ),
    make_option(
      opt_str = c("-g", "--genome-directory"),
      dest = "genome_directory",
      help = "Genome-specific analysis directory",
      type = "character"
    )
  )
))

prefix <-
  paste("chipseq", "diffbind", argument_list$factor, sep = "_")

output_directory <-
  file.path(argument_list$genome_directory, prefix)

if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

setwd(dir = output_directory)

message("Create a DiffBind dba object.")
diffbind_dba <-
  dba(sampleSheet = argument_list$sample_annotation,
      bCorPlot = FALSE)

message("Plotting correlation heatmap based on peak caller score data.")
file_name <-
  paste(prefix, "correlation_peak_caller_score.png", sep = "_")
png(filename = file.path(output_directory, file_name))
return_value <- dba.plotHeatmap(DBA = diffbind_dba, margin = 25)
active_device <- dev.off()

message("Counting reads.")
diffbind_dba <- dba.count(DBA = diffbind_dba, bCorPlot = FALSE)

message("Plotting correlation heatmap based on read counts.")
file_name <- paste(prefix, "correlation_read_counts.png", sep = "_")
png(filename = file.path(output_directory, file_name))
return_value <- dba.plotHeatmap(DBA = diffbind_dba, margin = 25)
active_device <- dev.off()

message("Establishing contrast by tissue.")
# The categories default to DBA_TISSUE, DBA_FACTOR, DBA_CONDITION and DBA_TREATMENT.
diffbind_dba <- dba.contrast(DBA = diffbind_dba, minMembers = 2)
# Check if setting contrasts was successful. It may not be if less than two replicates were available.
if (is.null(x = diffbind_dba$contrasts)) {
  # Set the mask manually. For the moment this only works for two samples.
  if (nrow(x = diffbind_dba$samples) == 2) {
    message("In lack of replicates, setting contrasts on the basis of the first two conditions.")
    diffbind_conditions <-
      unique(x = diffbind_dba$class[DBA_CONDITION, ])
    diffbind_dba <- dba.contrast(
      DBA = diffbind_dba,
      group1 = dba.mask(
        DBA = diffbind_dba,
        attribute = DBA_CONDITION,
        value = diffbind_conditions[1]
      ),
      group2 = dba.mask(
        DBA = diffbind_dba,
        attribute = DBA_CONDITION,
        value = diffbind_conditions[2]
      ),
      name1 = diffbind_conditions[1],
      name2 = diffbind_conditions[2]
    )
    rm(diffbind_conditions)
  }
}

message("Running differential binding affinity analysis.")
diffbind_dba <- dba.analyze(DBA = diffbind_dba, bCorPlot = FALSE)

message("Plotting correlation heatmap based on differential binding affinity analysis.")
file_name <- paste(prefix, "correlation_analysis.png", sep = "_")
png(filename = file.path(output_directory, file_name))
return_value <- dba.plotHeatmap(DBA = diffbind_dba, margin = 25)
active_device <- dev.off()

process_per_contrast <- function(contrast, group1, group2) {
  # Process per row of a contrasts data frame obtained via dba.show()
  # contrast the row.names() string of the data frame indicating the contrast number
  # group1 Group1 value
  # group2 Group2 value
  
  message(
    sprintf(
      "Write differentially bound sites for factor %s and contrast %s versus %s to disk.",
      argument_list$factor,
      group1,
      group2
    )
  )
  # TODO: The report function is quite peculiar in that it insists on a prefix DBA_
  # for the file name.
  diffbind_report <- dba.report(
    DBA = diffbind_dba,
    contrast = as.integer(contrast),
    bNormalized = TRUE,
    bCalled = TRUE,
    bCounts = TRUE,
    bCalledDetail = TRUE,
    file = sprintf("%s_report_%s__%s", argument_list$factor, group1, group2),
    DataType = DBA_DATA_FRAME
  )
  
  message(
    sprintf(
      "Creating MA plot for factor %s and contrast %s versus %s.",
      argument_list$factor,
      group1,
      group2
    )
  )
  file_name <-
    sprintf("%s_ma_plot_%s__%s.png", prefix, group1, group2)
  png(filename = file.path(output_directory, file_name))
  dba.plotMA(
    DBA = diffbind_dba,
    bNormalized = TRUE,
    bXY = FALSE,
    contrast = as.integer(contrast)
  )
  active_device <- dev.off()
  
  message(
    sprintf(
      "Creating scatter plot for factor %s and contrast %s versus %s.",
      argument_list$factor,
      group1,
      group2
    )
  )
  file_name <-
    sprintf("%s_scatter_plot_%s__%s.png", prefix, group1, group2)
  png(filename = file.path(output_directory, file_name))
  dba.plotMA(
    DBA = diffbind_dba,
    bNormalized = TRUE,
    bXY = TRUE,
    contrast = as.integer(contrast)
  )
  active_device <- dev.off()
  
  message(
    sprintf(
      "Creating PCA plot for factor %s and contrast %s versus %s.",
      argument_list$factor,
      group1,
      group2
    )
  )
  file_name <-
    sprintf("%s_pca_plot_%s__%s.png", prefix, group1, group2)
  png(filename = file.path(output_directory, file_name))
  dba.plotPCA(DBA = diffbind_dba, attributes = DBA_CONDITION)
  active_device <- dev.off()
  
  message(
    sprintf(
      "Creating Box plot for factor %s and contrast %s versus %s.",
      argument_list$factor,
      group1,
      group2
    )
  )
  file_name <-
    sprintf("%s_box_plot_%s__%s.png", prefix, group1, group2)
  png(filename = file.path(output_directory, file_name))
  diffbind_pvals <-
    dba.plotBox(DBA = diffbind_dba, bNormalized = TRUE)
  active_device <- dev.off()
  
  return(TRUE)
}

# Get a data frame with all contrasts to apply the above function to each row.
contrasts <- dba.show(DBA = diffbind_dba, bContrasts = TRUE)

# Replace '!' characters with 'not_'.
contrasts$Group1 <-
  gsub(pattern = "!",
       replacement = "not_",
       x = contrasts$Group1)
contrasts$Group2 <-
  gsub(pattern = "!",
       replacement = "not_",
       x = contrasts$Group2)

# Write the contrasts to disk.
write.csv(
  x = contrasts,
  file = paste(prefix, "contrasts.csv", sep = "_"),
  row.names = FALSE
)

# Apply the function to each row and discard the return value.
return_value <-
  mapply(FUN = process_per_contrast,
         row.names(contrasts),
         contrasts$Group1,
         contrasts$Group2)
rm(return_value)
rm(file_name)

# dba.overlap(DBA = diffbind_dba, mode = DBA_OLAP_RATE)

# message("Creating a Venn diagram")
# file_name <- paste(prefix, "box_plot.png", sep = "_")
# png(filename = file.path(output_directory, file_name))
# file_name <- dba.plotVenn(DBA = diffbind_dba)
# active_device <- dev.off()

rm(argument_list)

message("Save workspace image for manual post-processing.")
save.image()

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

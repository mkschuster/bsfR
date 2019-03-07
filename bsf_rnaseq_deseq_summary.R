#!/usr/bin/env Rscript
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --distribution=block
#SBATCH --mem=2048
#SBATCH --time=12:00:00
#SBATCH --partition=shortq
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --error=.bsf_rnaseq_deseq_summary_%j.err
#SBATCH --output=.bsf_rnaseq_deseq_summary_%j.out
#
# BSF R script to collate DESeq2 analysis summary tables.
#
#
# Copyright 2013 - 2019 Michael K. Schuster
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
      opt_str = c("--genome-directory"),
      default = ".",
      dest = "genome_directory",
      help = "Genome directory path",
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

suppressPackageStartupMessages(expr = library(package = "reshape2"))  # For reshape2::dcast()

summary_frame <- NULL
directory_names <-
  grep(
    pattern = "^rnaseq_deseq_.*",
    x = list.dirs(
      path = argument_list$genome_directory,
      full.names = FALSE,
      recursive = FALSE
    ),
    value = TRUE
  )
for (directory_name in directory_names) {
  design_name <-
    gsub(pattern = "^rnaseq_deseq_(.*?)$",
         replacement = "\\1",
         x = directory_name)
  file_names <-
    list.files(
      path = file.path(argument_list$genome_directory, directory_name),
      pattern = paste("^rnaseq_deseq", design_name, "contrasts_summary.tsv$", sep = "_"),
      full.names = TRUE
    )
  for (file_name in file_names) {
    design_table <-
      read.table(
        file = file_name,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
      )
    if (is.null(x = summary_frame)) {
      summary_frame <- design_table
    } else {
      summary_frame <- rbind(summary_frame, design_table)
    }
    rm(design_table)
  }
  rm(file_name, file_names, design_name)
}
rm(directory_name, directory_names)

# Replace all instances of NA in the Numerator and Denominator with "1".
summary_frame[is.na(x = summary_frame$Numerator), "Numerator"] <-
  "1"
summary_frame[is.na(x = summary_frame$Denominator), "Denominator"] <-
  "1"
# Replace all empty characters in the Numerator and Denominator with "1".
summary_frame[summary_frame$Denominator == "", "Denominator"] <- "1"

# Remove the "Label" variable, before casting.

subset_frame <- summary_frame[,!(names(x = summary_frame) %in% c("Label")), drop = FALSE]
casted_frame <-
  dcast(
    data = subset_frame,
    formula = Numerator + Denominator ~ Design,
    value.var = "Significant"
  )

write.table(
  x = casted_frame,
  file = "rnaseq_deseq_summary_design.tsv",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

# Find duplicated rows in "Design" and "Label" variables.

casted_frame <- summary_frame[, c("Design", "Label"), drop = FALSE]
subset_frame <- summary_frame[duplicated(x = casted_frame) | duplicated(x = casted_frame, fromLast = TRUE), ]
if (nrow(x = subset_frame)) {
  print(x = "Duplicated rows.")
  print(x = subset_frame)
}

# Remove the "Numerator" and "Denominator" variables, before casting.
subset_frame <- summary_frame[,!(names(x = summary_frame) %in% c("Numerator", "Denominator")), drop = FALSE]
casted_frame <-
  dcast(
    data = subset_frame,
    formula = Label ~ Design,
    value.var = "Significant"
  )

write.table(
  x = casted_frame,
  file = "rnaseq_deseq_summary_labels.tsv",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

if (FALSE) {
  write.table(
    x = summary_frame[,!(names(x = summary_frame) %in% c("Numerator", "Denominator")), drop = FALSE],
    file = "rnaseq_deseq_summary_simple.tsv",
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
}

rm(casted_frame, subset_frame, summary_frame, argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

#! /usr/bin/env Rscript
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --distribution=block
#SBATCH --mem=65536
#SBATCH --time=12:00:00
#SBATCH --partition=shortq
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --error=.bsf_rnaseq_deseq_results_%j.err
#SBATCH --output=.bsf_rnaseq_deseq_results_%j.out
#
# BSF R script to extract results of a DESeq2 analysis.
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
      opt_str = c("--design-name"),
      # default = "global",
      dest = "design_name",
      help = "design name",
      type = "character"
    ),
    make_option(
      opt_str = c("--lfc-threshold"),
      default = 0.0,
      dest = "lfc_threshold",
      help = "Log-fold change threshold [0.0]",
      type = "numeric"
    ),
    make_option(
      opt_str = c("--padj-threshold"),
      default = 0.1,
      dest = "padj_threshold",
      help = "Adjusted p-value threshold [0.1]",
      type = "numeric"
    ),
    make_option(
      opt_str = c("--threads"),
      default = 1L,
      dest = "threads",
      help = "Number of parallel processing threads",
      type = "integer"
    )
  )
))

# Check the input.

if (is.null(x = argument_list$design_name)) {
  stop("Missing --design-name option")
}

suppressPackageStartupMessages(expr = library(package = "BiocParallel"))
suppressPackageStartupMessages(expr = library(package = "DESeq2"))

message(paste0("Processing design '", argument_list$design_name, "'"))

# Set the number of parallel threads in the MulticoreParam instance.

BiocParallel::register(BPPARAM = MulticoreParam(workers = argument_list$threads))

# The working directory is the analyis genome directory.
# Create a new sub-directory for results if it does not exist.

prefix <-
  paste("rnaseq",
        "deseq",
        argument_list$design_name,
        sep = "_")

output_directory <- prefix
if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

# Annotation data.frame ---------------------------------------------------


# Load a pre-calculated annotation frame.
annotation_frame <- NULL

file_path <-
  file.path(output_directory,
            paste(prefix, "annotation.tsv", sep = "_"))
if (file.exists(file_path) &&
    file.info(file_path)$size > 0L) {
  message("Loading annotation frame")
  annotation_frame <-
    read.table(
      file = file_path,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
} else {
  stop(paste0("Requiring an annotation frame in file: ", file_path))
}
rm(file_path)

# DESeqDataSet ------------------------------------------------------------


# Load a pre-calculated DESeqDataSet object.
deseq_data_set <- NULL

file_path <-
  file.path(output_directory,
            paste0(prefix, "_deseq_data_set.Rdata"))
if (file.exists(file_path) &&
    file.info(file_path)$size > 0L) {
  message("Loading a DESeqDataSet object")
  load(file = file_path)
} else {
  stop(paste0(
    "Require a pre-calculated DESeqDataSet object in file: ",
    file_path
  ))
}
rm(file_path)


# Contrasts data.frame ----------------------------------------------------


# Read a data frame of contrasts with variables "Design", "Numerator", "Denominator" and "Label".
contrast_frame <-
  read.table(
    file = file.path(output_directory, paste(prefix, "contrasts.tsv", sep = "_")),
    header = TRUE,
    sep = "\t",
    colClasses = c(
      "Design" = "character",
      "Numerator" = "character",
      "Denominator" = "character",
      "Label" = "character"
    ),
    stringsAsFactors = FALSE
  )
# Subset to the selected design.
contrast_frame <-
  contrast_frame[contrast_frame$Design == argument_list$design_name, ]
if (nrow(contrast_frame) == 0L) {
  stop("No design remaining after selection for design name.")
}
contrast_frame$Significant <- 0L

for (i in seq_len(length.out = nrow(x = contrast_frame))) {
  # The "contrast" option of the DESeq results() function expects a list of numerator and denominator.
  contrast_list <-
    list(numerator = unlist(x = strsplit(x = contrast_frame[i, "Numerator"], split = ",")),
         denominator = unlist(x = strsplit(x = contrast_frame[i, "Denominator"], split = ",")))
  contrast_character <-
    paste(paste(contrast_list$numerator, collapse = "_"),
          "against",
          if (length(x = contrast_list$denominator) > 0L) {
            paste(contrast_list$denominator, collapse = "_")
          } else {
            "intercept"
          },
          sep = "_")
  
  message(paste0("Creating DESeqResults for ", contrast_character))
  
  deseq_results <-
    results(
      object = deseq_data_set,
      contrast = contrast_list,
      lfcThreshold = argument_list$lfc_threshold,
      alpha = argument_list$padj_threshold,
      format = "DataFrame",
      tidy = FALSE,
      # If tidy is TRUE, a classical data.frame is returned.
      parallel = TRUE
    )
  
  # print(x = summary(object = deseq_results))
  
  # MA Plot ---------------------------------------------------------------
  
  
  # Create a MA plot.
  file_path <-
    file.path(output_directory,
              paste(
                paste(prefix,
                      "contrast",
                      contrast_character,
                      "ma",
                      sep = "_"),
                "pdf",
                sep = "."
              ))
  pdf(file = file_path)
  plotMA(object = deseq_results)
  return_value <- dev.off()
  rm(return_value)
  
  file_path <-
    file.path(output_directory,
              paste(
                paste(prefix,
                      "contrast",
                      contrast_character,
                      "ma",
                      sep = "_"),
                "png",
                sep = "."
              ))
  png(file = file_path)
  plotMA(object = deseq_results)
  return_value <- dev.off()
  rm(return_value)
  
  # DESeqResults DataFrame ------------------------------------------------
  
  
  # Re-adjust the DESeqResults DataFrame for merging with the annotation DataFrame.
  deseq_results_frame <-
    DataFrame(
      gene_id = rownames(x = deseq_results),
      deseq_results[, c("baseMean",
                        "log2FoldChange",
                        "lfcSE",
                        "stat",
                        "pvalue",
                        "padj")],
      significant = factor(x = "no", levels = c("no", "yes"))
    )
  deseq_results_frame[!is.na(x = deseq_results_frame$padj) &
                        deseq_results_frame$padj <= argument_list$padj_threshold, "significant"] <-
    "yes"
  
  deseq_merge <-
    merge(x = annotation_frame, y = deseq_results_frame, by = "gene_id")
  
  file_path <-
    file.path(output_directory,
              paste(
                paste(prefix,
                      "contrast",
                      contrast_character,
                      "genes",
                      sep = "_"),
                "tsv",
                sep = "."
              ))
  
  write.table(
    x = deseq_merge,
    file = file_path,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
  )
  
  # Significant DESeqResults DataFrame ------------------------------------
  
  
  deseq_merge_significant <-
    subset(x = deseq_merge, padj <= argument_list$padj_threshold)
  
  # Record the number of significant genes.
  contrast_frame[i, "Significant"] <-
    nrow(x = deseq_merge_significant)
  
  file_path <-
    file.path(output_directory,
              paste(
                paste(prefix,
                      "contrast",
                      contrast_character,
                      "significant",
                      sep = "_"),
                "tsv",
                sep = "."
              ))
  
  write.table(
    x = deseq_merge_significant,
    file = file_path,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE
  )
  
  rm(
    deseq_merge_significant,
    deseq_merge,
    deseq_results_frame,
    deseq_results,
    file_path,
    contrast_character,
    contrast_list
  )
}

# Write summary frame -----------------------------------------------------


write.table(
  x = contrast_frame,
  file = file.path(output_directory, paste(prefix, "summary.tsv", sep = "_")),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)

rm(
  i,
  contrast_frame,
  deseq_data_set,
  annotation_frame,
  output_directory,
  prefix,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessionInfo())

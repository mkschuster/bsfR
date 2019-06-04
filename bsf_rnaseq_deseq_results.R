#!/usr/bin/env Rscript
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

# Check the input.

if (is.null(x = argument_list$design_name)) {
  stop("Missing --design-name option")
}

suppressPackageStartupMessages(expr = library(package = "BiocParallel"))
suppressPackageStartupMessages(expr = library(package = "DESeq2"))
suppressPackageStartupMessages(expr = library(package = "EnhancedVolcano"))

# Save plots in the following formats.

graphics_formats <- c("pdf", "png")

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
  contrast_frame[contrast_frame$Design == argument_list$design_name,]
if (nrow(contrast_frame) == 0L) {
  stop("No design remaining after selection for design name.")
}
contrast_frame$Significant <- 0L

for (i in seq_len(length.out = nrow(x = contrast_frame))) {
  # The "contrast" option of the DESeq results() function expects a list of numerator and denominator.
  contrast_list <-
    list("numerator" = unlist(x = strsplit(x = contrast_frame[i, "Numerator"], split = ",")),
         "denominator" = unlist(x = strsplit(x = contrast_frame[i, "Denominator"], split = ",")))
  contrast_character <-
    paste(paste(contrast_list$numerator, collapse = "_"),
          "against",
          if (length(x = contrast_list$denominator) > 0L) {
            paste(contrast_list$denominator, collapse = "_")
          } else {
            "intercept"
          },
          sep = "_")

  # Check for the significant genes table and if it exist already,
  # read it to get the number of significant genes for the summary data frame.

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

  if (file.exists(file_path) &&
      file.info(file_path)$size > 0L) {
    message(paste0("Skipping DESeqResults for ", contrast_character))

    deseq_merge_significant <-
      read.table(
        file = file_path,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
      )

    # Record the number of significant genes.
    contrast_frame[i, "Significant"] <-
      nrow(x = deseq_merge_significant)

    rm(deseq_merge_significant)
  } else {
    message(paste0("Creating DESeqResults for ", contrast_character))

    deseq_results_default <-
      DESeq2::results(
        object = deseq_data_set,
        contrast = contrast_list,
        lfcThreshold = argument_list$lfc_threshold,
        alpha = argument_list$padj_threshold,
        format = "DataFrame",
        tidy = FALSE,
        # If tidy is TRUE, a classical data.frame is returned.
        parallel = TRUE
      )

    # print(x = summary(object = deseq_results_default))

    # Run lfcShrink() if possible.
    deseq_results_shrunk <- NULL
    if (any(attr(x = deseq_data_set, which = "full_rank"))) {
      # The original DESEqDataSet object is full rank, so that log2-fold changes can be shrunk.
      # The any() function returns FALSE for NULL values.
      message(paste0("Shrinking log2-fold changes for ", contrast_character))
      deseq_results_shrunk <- DESeq2::lfcShrink(
        dds = deseq_data_set,
        contrast = contrast_list,
        res = deseq_results_default,
        type = "ashr",
        lfcThreshold = argument_list$lfc_threshold,
        format = "DataFrame",
        parallel = TRUE
      )
    }

    # MA Plot ---------------------------------------------------------------


    # Create a MA plot.
    message(paste0("Creating a MA plot for ", contrast_character))
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
    grDevices::pdf(
      file = file_path,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    if (is.null(x = deseq_results_shrunk)) {
      DESeq2::plotMA(object = deseq_results_default)
    } else {
      DESeq2::plotMA(object = deseq_results_shrunk)
    }
    base::invisible(x = grDevices::dev.off())

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
    grDevices::png(
      file = file_path,
      width = argument_list$plot_width,
      height = argument_list$plot_height,
      units = "in",
      res = 300L
    )
    if (is.null(x = deseq_results_shrunk)) {
      DESeq2::plotMA(object = deseq_results_default)
    } else {
      DESeq2::plotMA(object = deseq_results_shrunk)
    }
    base::invisible(x = grDevices::dev.off())

    # Enhanced Volcano plot -----------------------------------------------


    # To annotate gene symbols rather than Ensembl gene identifiers,
    # the DESeqResults DataFrame needs merging with the annotation frame.
    message(paste0("Creating an EnhancedVolcano plot for ", contrast_character))
    deseq_results_frame <-
      if (is.null(x = deseq_results_shrunk))
        DataFrame(
          gene_id = rownames(x = deseq_results_default),
          deseq_results_default[, c("baseMean",
                                    "log2FoldChange",
                                    "lfcSE",
                                    "stat",
                                    "pvalue",
                                    "padj")],
          significant = factor(x = "no", levels = c("no", "yes"))
        )
    else
      DataFrame(
        gene_id = rownames(x = deseq_results_shrunk),
        deseq_results_shrunk[, c("baseMean",
                                 "log2FoldChange",
                                 "lfcSE",
                                 "pvalue", # no "stat" column
                                 "padj")],
        significant = factor(x = "no", levels = c("no", "yes"))
      )

    deseq_results_merged <-
      merge(x = annotation_frame, y = deseq_results_frame, by = "gene_id")

    ggplot_object <- EnhancedVolcano::EnhancedVolcano(
      toptable = deseq_results_merged,
      lab = deseq_results_merged$gene_name,
      x = "log2FoldChange",
      y = "pvalue",
      selectLab = c()
    )
    rm(deseq_results_merged, deseq_results_frame)

    for (graphics_format in graphics_formats) {
      ggplot2::ggsave(
        filename = file.path(
          output_directory,
          paste(
            paste(prefix,
                  "contrast",
                  contrast_character,
                  "volcano",
                  sep = "_"),
            graphics_format,
            sep = "."
          )
        ),
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height
      )
    }
    rm(graphics_format,
       ggplot_object)

    # DESeqResults DataFrame ----------------------------------------------


    # Adjust the DESeqResults DataFrame for merging with the annotation DataFrame,
    # by setting the rownames() as "gene_id" variable.
    deseq_results_frame <-
      DataFrame(
        gene_id = rownames(x = deseq_results_default),
        deseq_results_default[, c("baseMean",
                                  "log2FoldChange",
                                  "lfcSE",
                                  "stat",
                                  "pvalue",
                                  "padj")],
        significant = factor(x = "no", levels = c("no", "yes"))
      )

    # Assign the "significant" factor on the basis of the adjusted p-value threshold.
    deseq_results_frame[!is.na(x = deseq_results_frame$padj) &
                          deseq_results_frame$padj <= argument_list$padj_threshold, "significant"] <-
      "yes"

    # Calculate ranks for ...
    # (1) the effect size (log2FoldChange), ...
    if (is.null(x = deseq_results_shrunk)) {
      deseq_results_frame$rank_log2_fold_change <-
        base::rank(
          x = -abs(x = deseq_results_frame$log2FoldChange),
          ties.method = c("min")
        )
    } else {
      deseq_results_frame$rank_log2_fold_change <-
        base::rank(
          x = -abs(x = deseq_results_shrunk$log2FoldChange),
          ties.method = c("min")
        )
    }

    # (2) the absolute level (baseMean) and ...
    deseq_results_frame$rank_base_mean <-
      base::rank(x = -deseq_results_frame$baseMean,
                 ties.method = c("min"))

    # (3) the statistical significance (padj).
    deseq_results_frame$rank_padj <-
      base::rank(x = deseq_results_frame$padj, ties.method = c("min"))

    # Calculate the maximum of the three ranks.
    deseq_results_frame$max_rank <-
      base::pmax(
        deseq_results_frame$rank_log2_fold_change,
        deseq_results_frame$rank_base_mean,
        deseq_results_frame$rank_padj
      )

    deseq_merge_complete <-
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
      x = deseq_merge_complete,
      file = file_path,
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE
    )

    # Significant DESeqResults DataFrame ------------------------------------


    deseq_merge_significant <-
      subset(x = deseq_merge_complete, padj <= argument_list$padj_threshold)

    # Record the number of significant genes.
    contrast_frame[i, "Significant"] <-
      nrow(x = deseq_merge_significant)

    file_path <-
      file.path(output_directory,
                paste(
                  paste(
                    prefix,
                    "contrast",
                    contrast_character,
                    "significant",
                    sep = "_"
                  ),
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
      deseq_merge_complete,
      deseq_results_frame,
      deseq_results_shrunk,
      deseq_results_default,
      contrast_character,
      contrast_list
    )
  }
  rm(file_path)
}

# Write summary frame -----------------------------------------------------


write.table(
  x = contrast_frame,
  file = file.path(
    output_directory,
    paste(prefix, "contrasts", "summary.tsv", sep = "_")
  ),
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
  graphics_formats,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessionInfo())

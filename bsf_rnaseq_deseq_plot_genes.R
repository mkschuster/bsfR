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
#SBATCH --error=.bsf_rnaseq_deseq_plot_genes_%j.err
#SBATCH --output=.bsf_rnaseq_deseq_plot_genes_%j.out
#
# BSF R script to plot (transformed) counts of individual genes of a DESeq2 analysis.
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
      help = "Design name",
      type = "character"
    ),
    make_option(
      opt_str = c("--genes-path"),
      dest = "genes_path",
      help = "File path of a table of 'gene_id' and 'gene_name' columns to plot",
      type = "character"
    ),
    make_option(
      opt_str = c("--groups"),
      dest = "groups",
      help = "Groups or factors",
      type = "character"
    ),
    make_option(
      opt_str = c("--normalised"),
      action = "store_true",
      default = TRUE,
      dest = "normalised",
      help = "Normalised gene counts [TRUE]",
      type = "logical"
    ),
    make_option(
      opt_str = c("--non-normalised"),
      action = "store_false",
      default = FALSE,
      dest = "normalised",
      help = "Non-normalised gene counts [FALSE]",
      type = "logical"
    ),
    make_option(
      opt_str = c("--maximum-number"),
      default = 25L,
      dest = "maximum_number",
      help = "Maximum number of genes to plot [25]",
      type = "integer"
    ),
    make_option(
      opt_str = c("--threads"),
      default = 1L,
      dest = "threads",
      help = "Number of parallel processing threads [1]",
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

if (is.null(x = argument_list$groups)) {
  stop("Missing --groups option")
}
if (is.null(x = argument_list$genes_path)) {
  stop("Missing --genes_path option")
}
if (is.null(x = argument_list$design_name)) {
  stop("Missing --design-name option")
}

suppressPackageStartupMessages(expr = library(package = "DESeq2"))  # for DESeq2::DESeqDataSet() and DESeq2::DESeq()
suppressPackageStartupMessages(expr = library(package = "ggplot2"))
suppressPackageStartupMessages(expr = library(package = "stringi"))  # For stringi::stri_split_fixed()

# Save plots in the following formats.

graphics_formats <- c("pdf", "png")

message(paste0("Processing design '", argument_list$design_name, "'"))

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

# Read a data frame of genes to plot.
genes_frame <-
  read.table(
    file = argument_list$genes_path,
    header = TRUE,
    sep = "\t",
    colClasses = c("gene_id" = "character",
                   "gene_name" = "character"),
    stringsAsFactors = FALSE
  )

if (!is.null(x = genes_frame$padj) &
    nrow(x = genes_frame) > argument_list$maximum_number) {
  message(paste0(
    "Plotting only ",
    argument_list$maximum_number,
    " out of ",
    nrow(x = genes_frame),
    " genes."
  ))
  # Order by adjusted p-value.
  genes_frame <- genes_frame[order(genes_frame$padj),]
  genes_frame <-
    genes_frame[seq_len(length.out = argument_list$maximum_number),]
}

for (i in seq_len(length.out = nrow(x = genes_frame))) {
  for (graphics_format in graphics_formats) {
    # Create a (transformed) counts plot.
    file_path <-
      file.path(output_directory,
                paste(
                  paste(prefix,
                        "gene",
                        genes_frame[i, "gene_id"],
                        genes_frame[i, "gene_name"],
                        sep = "_"),
                  graphics_format,
                  sep = "."
                ))
    if (file.exists(file_path) &&
        file.info(file_path)$size > 0L) {
      message(paste("Skipping plot", genes_frame[i, "gene_id"], genes_frame[i, "gene_name"], sep = " "))
    } else {
      message(paste("Creating plot", genes_frame[i, "gene_id"], genes_frame[i, "gene_name"], sep = " "))
      count_frame <- DESeq2::plotCounts(
        dds = deseq_data_set,
        gene = genes_frame[i, "gene_id"],
        intgroup = stri_split(str = argument_list$groups, fixed = ",")[[1L]],
        normalized = argument_list$normalised,
        xlab = paste(genes_frame[i, "gene_name"], "(", genes_frame[i, "gene_id"], ")", sep = " "),
        returnData = TRUE
      )
      # Create a now covariates variable pasting all values selected by the --groups option.
      count_frame$covariates <-
        factor(x = apply(
          X = subset(
            x = count_frame,
            select = -count,
            drop = FALSE
          ),
          MARGIN = 1,
          FUN = paste,
          collapse = ":"
        ))
      ggplot_object <- ggplot2::ggplot(data = count_frame)
      ggplot_object <-
        ggplot_object + ggplot2::geom_point(mapping = aes(x = covariates, y = count),
                                            alpha = I(1 / 3))
      ggplot_object <-
        ggplot_object + ggplot2::theme(axis.text.x = element_text(
          size = 8.0,
          hjust = 1.0,
          vjust = 0.5,
          angle = 90.0
        ))
      ggplot_object <-
        ggplot_object + ggplot2::ylab(label = if (argument_list$normalised) {
          "normalised counts"
        } else {
          "counts"
        })
      ggplot_object <-
        ggplot_object + ggplot2::ggtitle(label = paste("Gene", "Counts", genes_frame[i, "gene_id"], genes_frame[i, "gene_name"], sep = " "))
      ggplot2::ggsave(
        filename = file_path,
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height
      )
    }
    rm(file_path)
  }
  rm(graphics_format)
}

rm(
  i,
  genes_frame,
  deseq_data_set,
  argument_list,
  output_directory,
  prefix,
  graphics_formats
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessionInfo())

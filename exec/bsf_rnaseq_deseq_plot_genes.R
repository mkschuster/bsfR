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
        opt_str = c("--design-name"),
        # default = "global",
        dest = "design_name",
        help = "Design name",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--genes-path"),
        dest = "genes_path",
        help = "File path of a table of 'gene_id' and 'gene_name' columns to plot",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--groups"),
        dest = "groups",
        help = "Groups or factors",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--normalised"),
        action = "store_true",
        default = TRUE,
        dest = "normalised",
        help = "Normalised gene counts [TRUE]",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = c("--non-normalised"),
        action = "store_false",
        default = FALSE,
        dest = "normalised",
        help = "Non-normalised gene counts [FALSE]",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = c("--maximum-number"),
        default = 25L,
        dest = "maximum_number",
        help = "Maximum number of genes to plot [25]",
        type = "integer"
      ),
      optparse::make_option(
        opt_str = c("--genome-directory"),
        default = ".",
        dest = "genome_directory",
        help = "Genome directory path [.]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--output-directory"),
        default = ".",
        dest = "output_directory",
        help = "Output directory path [.]",
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

suppressPackageStartupMessages(expr = library(package = "bsfR"))
suppressPackageStartupMessages(expr = library(package = "DESeq2"))
suppressPackageStartupMessages(expr = library(package = "tidyverse"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

message("Processing design '", argument_list$design_name, "'")

prefix_deseq <-
  bsfR::bsfrd_get_prefix_deseq(design_name = argument_list$design_name)

# The working directory is the analysis genome directory.
# Create a new sub-directory for results if it does not exist.

output_directory <-
  file.path(argument_list$output_directory, prefix_deseq)
if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

# Load a pre-calculated DESeqDataSet object.
deseq_data_set <-
  bsfR::bsfrd_read_deseq_data_set(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name
  )

# Read a tibble of genes to plot.
genes_tibble <-
  bsfR::bsfrd_read_gene_set_tibble(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name,
    gene_set_path = argument_list$genes_path
  )

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

if ("padj" %in% names(x = genes_frame) &&
    nrow(x = genes_frame) > argument_list$maximum_number) {
  message(
    "Plotting only ",
    argument_list$maximum_number,
    " out of ",
    nrow(x = genes_frame),
    " genes."
  )
  # Order by adjusted p-value.
  genes_frame <-
    genes_frame[order(genes_frame$padj), , drop = FALSE]
  genes_frame <-
    genes_frame[seq_len(length.out = argument_list$maximum_number), , drop = FALSE]
}

for (i in seq_len(length.out = nrow(x = genes_frame))) {
  for (graphics_format in graphics_formats) {
    # Create a (transformed) counts plot.
    file_path <-
      file.path(output_directory,
                paste(
                  paste(prefix_deseq,
                        "gene",
                        genes_frame[i, "gene_id", drop = TRUE],
                        genes_frame[i, "gene_name", drop = TRUE],
                        sep = "_"),
                  graphics_format,
                  sep = "."
                ))
    if (file.exists(file_path) &&
        file.info(file_path)$size > 0L) {
      message("Skipping plot", genes_frame[i, "gene_id", drop = TRUE], genes_frame[i, "gene_name", drop = TRUE], sep = " ")
    } else {
      message("Creating plot", genes_frame[i, "gene_id", drop = TRUE], genes_frame[i, "gene_name", drop = TRUE], sep = " ")
      count_frame <- DESeq2::plotCounts(
        dds = deseq_data_set,
        gene = genes_frame[i, "gene_id", drop = TRUE],
        intgroup = stringr::str_split(
          string = argument_list$groups,
          pattern = stringr::fixed(pattern = ",")
        )[[1L]],
        normalized = argument_list$normalised,
        xlab = paste(genes_frame[i, "gene_name", drop = TRUE], "(", genes_frame[i, "gene_id", drop = TRUE], ")", sep = " "),
        returnData = TRUE
      )
      # Create a new covariates variable pasting all values selected by the --groups option.
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
        ggplot_object + ggplot2::geom_point(
          mapping = ggplot2::aes(x = .data$covariates, y = .data$count),
          alpha = I(1 / 3)
        )
      ggplot_object <-
        ggplot_object + ggplot2::labs(
          x = "Covariates",
          y = if (argument_list$normalised) {
            "Normalised Counts"
          } else {
            "Counts"
          },
          title = paste("Gene", "Counts", genes_frame[i, "gene_id", drop = TRUE], genes_frame[i, "gene_name", drop = TRUE], sep = " ")
        )
      ggplot_object <-
        ggplot_object + ggplot2::theme(axis.text.x = ggplot2::element_text(
          size = 8.0,
          hjust = 1.0,
          vjust = 0.5,
          angle = 90.0
        ))
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
  prefix_deseq,
  graphics_formats
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessionInfo())

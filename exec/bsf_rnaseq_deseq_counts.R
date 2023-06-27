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


# BSF R script to plot (transformed) counts of individual genes of a DESeq2 analysis.

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
        opt_str = "--design-name",
        # default = "global",
        dest = "design_name",
        help = "Design name",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--gene-path",
        dest = "gene_path",
        help = "Gene set file path for counts plotting [NULL]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--groups",
        dest = "groups",
        help = "Groups or factors",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--normalised",
        action = "store_true",
        default = TRUE,
        dest = "normalised",
        help = "Normalised gene counts [TRUE]",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = "--non-normalised",
        action = "store_false",
        default = FALSE,
        dest = "normalised",
        help = "Non-normalised gene counts [FALSE]",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = "--genome-directory",
        default = ".",
        dest = "genome_directory",
        help = "Genome directory path [.]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--output-directory",
        default = ".",
        dest = "output_directory",
        help = "Output directory path [.]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--plot-width",
        default = 7.0,
        dest = "plot_width",
        help = "Plot width in inches [7.0]",
        type = "double"
      ),
      optparse::make_option(
        opt_str = "--plot-height",
        default = 7.0,
        dest = "plot_height",
        help = "Plot height in inches [7.0]",
        type = "double"
      )
    )
  ))

if (is.null(x = argument_list$groups)) {
  stop("Missing --groups option")
}
if (is.null(x = argument_list$gene_path)) {
  stop("Missing --gene-path option")
}
if (is.null(x = argument_list$design_name)) {
  stop("Missing --design-name option")
}

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "dplyr"))
suppressPackageStartupMessages(expr = library(package = "ggplot2"))
suppressPackageStartupMessages(expr = library(package = "stringr"))
# Bioconductor
suppressPackageStartupMessages(expr = library(package = "BiocVersion"))
suppressPackageStartupMessages(expr = library(package = "DESeq2"))
# BSF
suppressPackageStartupMessages(expr = library(package = "bsfR"))

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
    design_name = argument_list$design_name,
    verbose = argument_list$verbose
  )

# Read a tibble of genes to plot.
genes_tibble <-
  bsfR::bsfrd_read_gene_set_tibble(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name,
    gene_set_path = argument_list$gene_path,
    verbose = argument_list$verbose
  )

# If the tibble exists, test for NA values in the gene_id variable.
missing_tibble <-
  dplyr::filter(.data = genes_tibble, is.na(x = .data$gene_id))
if (base::nrow(x = missing_tibble) > 0L) {
  print(x = "The following gene_name values could not be resolved into gene_id values:")
  print(x = missing_tibble)
}
rm(missing_tibble)

# Finally, filter out all observations with NA values in the gene_id variable.
genes_tibble <-
  dplyr::filter(.data = genes_tibble,!is.na(x = .data$gene_id))

for (i in seq_len(length.out = base::nrow(x = genes_tibble))) {
  for (graphics_format in graphics_formats) {
    # Create a (transformed) counts plot.
    file_path <-
      file.path(output_directory,
                paste(
                  paste(
                    prefix_deseq,
                    "gene",
                    genes_tibble$gene_id[i],
                    genes_tibble$gene_name[i],
                    sep = "_"
                  ),
                  graphics_format,
                  sep = "."
                ))
    if (file.exists(file_path) &&
        file.info(file_path)$size > 0L) {
      message("Skipping plot",
              genes_tibble$gene_id[i],
              genes_tibble$gene_name[i],
              sep = " ")
    } else {
      message("Creating plot",
              genes_tibble$gene_id[i],
              genes_tibble$gene_name[i],
              sep = " ")
      count_frame <- DESeq2::plotCounts(
        dds = deseq_data_set,
        gene = genes_tibble$gene_id[i],
        intgroup = stringr::str_split(
          string = argument_list$groups,
          pattern = stringr::fixed(pattern = ",")
        )[[1L]],
        normalized = argument_list$normalised,
        xlab = paste(
          genes_tibble$gene_name[i],
          "(",
          genes_tibble$gene_id[i],
          ")",
          sep = " "
        ),
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
          title = paste(
            "Gene",
            "Counts",
            genes_tibble$gene_id[i],
            genes_tibble$gene_name[i],
            sep = " "
          )
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
      ggplot2::ggsave(
        filename = file_path,
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height
      )
      rm(ggplot_object, count_frame)
    }
    rm(file_path)
  }
  rm(graphics_format)
}

rm(
  i,
  genes_tibble,
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

print(x = sessioninfo::session_info())

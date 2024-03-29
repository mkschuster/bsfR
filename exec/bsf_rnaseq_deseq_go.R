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


# BSF R script to annotate DESeq2 results with the Gene Ontology.

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
        opt_str = "--annotation-dbi",
        dest = "annotation_dbi",
        help = "Organism-specific AnnotationDbi package [NULL]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--padj-threshold",
        default = 0.1,
        dest = "padj_threshold",
        help = "Threshold for the adjusted p-value [0.1]",
        type = "double"
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

if (is.null(x = argument_list$design_name)) {
  stop("Missing --design-name option")
}

if (is.null(x = argument_list$annotation_dbi)) {
  stop("Missing --annotation-dbi option")
}

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "ggplot2"))
# CRAN
# suppressPackageStartupMessages(expr = library(package = "Nozzle.R1"))
# Bioconductor
suppressPackageStartupMessages(expr = library(package = "topGO"))
# BSF
suppressPackageStartupMessages(expr = library(package = "bsfR"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

go_names <-
  c("BP" = "Biological Process", "CC" = "Cellular Component", "MF" = "Molecular Function")

prefix_go <-
  bsfR::bsfrd_get_prefix_go(design_name = argument_list$design_name)

output_directory <-
  file.path(argument_list$output_directory, prefix_go)
if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

gene_selection_function <- function(x) {
  # Select genes by their adjusted p-value.
  return(x <= argument_list$padj_threshold)
}

contrast_tibble <-
  bsfR::bsfrd_read_contrast_tibble(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name,
    summary = TRUE,
    verbose = argument_list$verbose
  )

for (contrast_index in seq_len(length.out = base::nrow(x = contrast_tibble))) {
  contrast_character <-
    bsfR::bsfrd_get_contrast_character(contrast_tibble = contrast_tibble, index = contrast_index)

  # Read the annotated results tibble with all genes for this contrast.
  deseq_results_tibble <-
    bsfR::bsfrd_read_result_tibble(
      genome_directory = argument_list$genome_directory,
      design_name = argument_list$design_name,
      contrast_tibble = contrast_tibble,
      index = contrast_index,
      verbose = argument_list$verbose
    )

  if (is.null(x = deseq_results_tibble)) {
    rm(deseq_results_tibble, contrast_character)
    next()
  }

  # Reset the row names from the gene_id variable.
  # NOTE: The as.data.frame.tbl_df() function does not use the row.names option.
  deseq_results_frame <- base::as.data.frame(x = deseq_results_tibble)
  base::row.names(x = deseq_results_frame) <- deseq_results_tibble$gene_id

  # Filter in a data frame to have access to padj and SE.
  # Importantly, NA values need removing.
  message("Number of genes before NA removal: ",
          base::nrow(x = deseq_results_frame))

  deseq_go_frame <-
    deseq_results_frame[!is.na(x = deseq_results_frame$padj), , drop = FALSE]

  message("Number of genes after NA removal: ",
          base::nrow(x = deseq_go_frame))

  # Introduce a go_status factor variable for plotting.
  deseq_go_frame$go_status <-
    factor(
      x = character(length = base::nrow(x = deseq_go_frame)),
      levels = c("Used", "Not annotated", "Filtered")
    )

  # topGO needs a named double vector of adjusted p-values.
  all_genes_double <- deseq_go_frame$padj
  base::names(x = all_genes_double) <- deseq_go_frame$gene_id

  for (sub_go in c("BP", "CC", "MF")) {
    message("Create a topGOdata object for ontology ", sub_go)
    topgo_data <-
      new(
        Class = "topGOdata",
        # Options of the initialize() method of the topGOdata class.
        ontology = sub_go,
        allGenes = all_genes_double,
        geneSelectionFun = gene_selection_function,
        description = contrast_tibble$Label[contrast_index],
        # expressionMatrix = ,
        # phenotype =
        # nodeSize =
        annotationFun = topGO::annFUN.org,
        # annFUN.org(whichOnto, feasibleGenes = NULL, mapping, ID = "entrez")
        # The options "whichOnto" and "feasibleGenes" are filled in by the initialize() method.
        mapping = argument_list$annotation_dbi,
        ID = "ensembl"
      )
    print(x = topgo_data)

    if (TRUE) {
      # Plot the GO status for each gene in the context of the log2-fold change standard error.
      # Reset the "go_status" factor for each sub-ontology.
      deseq_go_frame$go_status <- "Used"
      deseq_go_frame$go_status[deseq_go_frame$gene_id %in% topGO::genes(object = topgo_data)] <-
        "Not annotated"
      deseq_go_frame$go_status[!gene_selection_function(x = deseq_go_frame$padj)] <-
        "Filtered"

      # gene_group_table <- table(deseq_go_frame$go_status)

      # Plot the groups.
      ggplot_object <- ggplot2::ggplot(data = deseq_go_frame)
      ggplot_object <-
        ggplot_object + ggplot2::geom_point(
          mapping = ggplot2::aes(
            x = .data$padj,
            y = .data$lfcSE,
            colour = .data$go_status
          ),
          alpha = alpha = I(1 / 3)
        )
      ggplot_object <-
        ggplot_object + ggplot2::labs(
          x = "Adjusted p-value",
          y = "Log2-fold change standard error",
          colour = "GO Status",
          title = contrast_tibble$Label[contrast_index],
          subtitle = paste("Ontology:", go_names[sub_go])
        )
      ggplot2::ggsave(
        filename = file.path(output_directory,
                             paste(
                               paste(
                                 prefix_go,
                                 "contrast",
                                 contrast_character,
                                 "go",
                                 sub_go,
                                 sep = "_"
                               ),
                               "pdf",
                               sep = "."
                             )),
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height
      )
      rm(ggplot_object)
    }

    message("Running topGO::runTest() for ontology ", sub_go)
    topgo_result <-
      topGO::runTest(object = topgo_data,
                     algorithm = "weight01",
                     statistic = "ks")
    print(x = topgo_result)

    result_frame <-
      topGO::GenTable(
        object = topgo_data,
        weight_ks = topgo_result,
        orderBy = "weight_ks",
        topNodes = 50L,
        numChar = 80L
      )

    utils::write.table(
      x = result_frame,
      file = file.path(output_directory,
                       paste(
                         paste(
                           prefix_go,
                           "contrast",
                           contrast_character,
                           "go",
                           sub_go,
                           sep = "_"
                         ),
                         "tsv",
                         sep = "."
                       )),
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )
    rm(result_frame, topgo_result, topgo_data)
  }
  rm(
    sub_go,
    all_genes_double,
    deseq_go_frame,
    deseq_results_frame,
    deseq_results_tibble,
    contrast_character
  )
}

rm(
  contrast_index,
  gene_selection_function,
  contrast_tibble,
  output_directory,
  prefix_go,
  graphics_formats,
  go_names,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessioninfo::session_info())

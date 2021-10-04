#!/usr/bin/env Rscript
#
# BSF R script to annotate DESeq2 results with the Enrichr tool.
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
        opt_str = c("--padj-threshold"),
        default = 0.1,
        dest = "padj_threshold",
        help = "Threshold for the adjusted p-value [0.1]",
        type = "numeric"
      ),
      optparse::make_option(
        opt_str = c("--l2fc-threshold"),
        default = 1.0,
        dest = "l2fc_threshold",
        help = "Threshold for the log2(fold-change) [1.0]",
        type = "numeric"
      ),
      optparse::make_option(
        opt_str = c("--maximum-terms"),
        default = 20L,
        dest = "maximum_terms",
        help = "Maximum number of terms [20]",
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

if (is.null(x = argument_list$design_name)) {
  stop("Missing --design-name option")
}

suppressPackageStartupMessages(expr = library(package = "enrichR"))
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
suppressPackageStartupMessages(expr = library(package = "Nozzle.R1"))
suppressPackageStartupMessages(expr = library(package = "tidyverse"))

# TODO: The list of Enrichr databases should be configurable.

enrichr_databases <- c(
  "BioCarta_2016",
  "GO_Biological_Process_2018",
  "GO_Cellular_Component_2018",
  "GO_Molecular_Function_2018",
  "KEGG_2019_Human",
  "KEGG_2019_Mouse",
  "MGI_Mammalian_Phenotype_Level_4_2019",
  "Mouse_Gene_Atlas",
  "Reactome_2016",
  "WikiPathways_2019_Human",
  "WikiPathways_2019_Mouse"
)

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

prefix_deseq <-
  bsfR::bsfrd_get_prefix_deseq(design_name = argument_list$design_name)

prefix_enrichr <-
  bsfR::bsfrd_get_prefix_enrichr(design_name = argument_list$design_name)

output_directory <-
  file.path(argument_list$output_directory, prefix_enrichr)
if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

#' Load a DESeq2 results tibble for a specific contrast.
#'
#' @param contrast_character A \code{character} scalar with the contrast.
#'
#' @return A \code{tibble} of DESeq2 results.
#' @seealso bsfrd_read_result_tibble
#' @noRd
#'
#' @examples
load_deseq_result_tibble <- function(contrast_character) {
  deseq_results_tibble <-
    bsfR::bsfrd_read_result_tibble(
      genome_directory = argument_list$genome_directory,
      design_name = argument_list$design_name,
      contrast_character = contrast_character,
      verbose = argument_list$verbose
    )

  if (is.null(x = deseq_results_tibble)) {
    stop("No DESeqResults frame for contrast ", contrast_character)
  }

  # Importantly, NA values need removing.
  message("Number of genes in total: ",
          nrow(x = deseq_results_tibble))

  deseq_results_tibble <-
    dplyr::filter(.data = deseq_results_tibble,!is.na(x = .data$padj))

  message("Number of genes after NA removal: ",
          nrow(x = deseq_results_tibble))

  deseq_results_tibble <-
    dplyr::filter(.data = deseq_results_tibble, .data$padj <= argument_list$padj_threshold)

  message("Number of genes after applying the padj threshold: ",
          nrow(x = deseq_results_tibble))

  deseq_results_tibble <-
    dplyr::filter(.data = deseq_results_tibble,
                  abs(x = .data$log2FoldChange) >= argument_list$l2fc_threshold)

  message("Number of genes after applying the l2fc threshold: ",
          nrow(x = deseq_results_tibble))

  return(deseq_results_tibble)
}


#' Load an Enrichr result tibble for a specific contrast and database.
#'
#' @param contrast_character A \code{character} scalar with the contrast.
#' @param enrichr_database A \code{character} scalar with the Enrichr database.
#'
#' @return A named \code{list} of Enrichr result \code{tibble} objects.
#'   \code{up}: \code{tibble} with results of of up-regulated genes.
#'   \code{down}: \code{tibble} with results of down-regulated genes.
#' @noRd
#'
#' @examples
load_enrichr_results <-
  function(contrast_character, enrichr_database) {
    result_list <- list()
    directions <- c("up", "down")

    file_paths <- file.path(output_directory,
                            paste(
                              paste(
                                prefix_enrichr,
                                contrast_character,
                                enrichr_database,
                                directions,
                                sep = "_"
                              ),
                              "tsv",
                              sep = "."
                            ))

    if (all(file.exists(file_paths)) &&
        all(file.info(file_paths)$size > 0L)) {
      for (direction_index in seq_along(along.with = directions)) {
        message("Loading: ",
                contrast_character,
                " ",
                enrichr_database,
                " ",
                directions[direction_index])

        result_list[[directions[direction_index]]] <-
          readr::read_tsv(
            file = file_paths[direction_index],
            col_types = readr::cols(
              Term = readr::col_character(),
              Overlap = readr::col_character(),
              P.value = readr::col_double(),
              Adjusted.P.value = readr::col_double(),
              Old.P.value = readr::col_double(),
              Old.Adjusted.P.value = readr::col_double(),
              Odds.Ratio = readr::col_double(),
              Combined.Score = readr::col_double(),
              Genes = readr::col_character()
            )
          )
      }
      rm(direction_index)
    } else {
      message("Loading: ",
              contrast_character,
              " ",
              enrichr_database)

      deseq_results_tibble <-
        load_deseq_result_tibble(contrast_character = contrast_character)

      for (direction_index in seq_along(along.with = directions)) {
        # Split the results into up- and down-regulated genes.
        message("Processing: ",
                contrast_character,
                " ",
                enrichr_database,
                " ",
                directions[direction_index])

        enrichr_tibble <-
          if (directions[direction_index] == "up") {
            dplyr::filter(.data = deseq_results_tibble, .data$log2FoldChange > 0.0)
          } else {
            dplyr::filter(.data = deseq_results_tibble, .data$log2FoldChange < 0.0)
          }

        readr::write_tsv(x = enrichr_tibble,
                         file = file.path(output_directory,
                                          paste(
                                            paste(
                                              prefix_enrichr,
                                              contrast_character,
                                              "genes",
                                              directions[direction_index],
                                              sep = "_"
                                            ),
                                            "tsv",
                                            sep = "."
                                          )))

        # Since Enrichr needs valid gene symbols in the gene_name variable,
        # post-filter after writing the TSV file.
        enrichr_result_list <- NULL

        if ("gene_name" %in% base::names(x = enrichr_tibble)) {
          enrichr_tibble <-
            dplyr::filter(.data = enrichr_tibble,!is.na(x = .data$gene_name))

          if (nrow(x = enrichr_tibble) > 0L) {
            enrichr_result_list <-
              enrichR::enrichr(genes = enrichr_tibble$gene_name,
                               databases = enrichr_databases)
          }
        }

        # Save the Enrichr results for each database, so that it is
        # automatically available for the next query.
        for (edb in enrichr_databases) {
          enrichr_result_tibble <-
            if (is.null(x = enrichr_result_list) ||
                nrow(x = enrichr_result_list[[edb]]) == 0L) {
              # Initialise an empty Enrichr results tibble if no genes were filtered.
              tibble::tibble(
                Term = character(),
                Overlap = character(),
                P.value = double(),
                Adjusted.P.value = double(),
                Old.P.value = double(),
                Old.Adjusted.P.value = double(),
                Odds.Ratio = double(),
                Combined.Score = double(),
                Genes = character()
              )
            } else {
              # Use readr::type_convert() since the Type variable is sometimes
              # read as logical rather than character.
              tibble::as_tibble(x = readr::type_convert(
                df = enrichr_result_list[[edb]],
                col_types = readr::cols(
                  Term = readr::col_character(),
                  Overlap = readr::col_character(),
                  P.value = readr::col_double(),
                  Adjusted.P.value = readr::col_double(),
                  Old.P.value = readr::col_double(),
                  Old.Adjusted.P.value = readr::col_double(),
                  Odds.Ratio = readr::col_double(),
                  Combined.Score = readr::col_double(),
                  Genes = readr::col_character()
                )
              ))
            }

          readr::write_tsv(x = enrichr_result_tibble,
                           file = file.path(output_directory,
                                            paste(
                                              paste(
                                                prefix_enrichr,
                                                contrast_character,
                                                edb,
                                                directions[direction_index],
                                                sep = "_"
                                              ),
                                              "tsv",
                                              sep = "."
                                            )))

          if (edb == enrichr_database) {
            result_list[[directions[direction_index]]] <- enrichr_result_tibble
          }
          rm(enrichr_result_tibble)
        }
        rm(edb, enrichr_tibble, enrichr_result_list)
      }
      rm(direction_index)
    }
    rm(file_paths, directions)

    return(result_list)
  }

# Contrasts Tibble --------------------------------------------------------


# Read a contrast tibble with variables "Design", "Numerator", "Denominator" and
# "Label".
contrast_tibble <-
  bsfR::bsfrd_read_contrast_tibble(
    genome_directory = argument_list$genome_directory,
    design_name = argument_list$design_name,
    summary = TRUE,
    verbose = argument_list$verbose
  )

# Create a "Contrasts" report section
nozzle_section_contrasts <-
  Nozzle.R1::newSection("Contrasts", class = SECTION.CLASS.RESULTS)

nozzle_section_contrasts <-
  Nozzle.R1::addTo(parent = nozzle_section_contrasts, Nozzle.R1::newTable(table = base::as.data.frame(x = contrast_tibble)))

nozzle_section_enrichr <-
  Nozzle.R1::newSection("Enrichr Reports", class = SECTION.CLASS.RESULTS)

for (contrast_index in seq_len(length.out = nrow(x = contrast_tibble))) {
  contrast_character <-
    bsfR::bsfrd_get_contrast_character(contrast_tibble = contrast_tibble, index = contrast_index)

  nozzle_section_contrast <-
    Nozzle.R1::newSubSection("Contrast ", contrast_tibble$Label[contrast_index])

  nozzle_section_contrast <-
    Nozzle.R1::addTo(
      parent = nozzle_section_contrast,
      newParagraph(
        "Lists of up-regulated (",
        Nozzle.R1::asLink(url = paste(
          paste(prefix_enrichr,
                contrast_character,
                "genes",
                "up",
                sep = "_"),
          "tsv",
          sep = "."
        ), "TSV"),
        ") and down-regulated (",
        Nozzle.R1::asLink(url = paste(
          paste(prefix_enrichr,
                contrast_character,
                "genes",
                "down",
                sep = "_"),
          "tsv",
          sep = "."
        ), "TSV"),
        ") genes."
      )
    )

  for (enrichr_index in seq_len(length.out = length(x = enrichr_databases))) {
    result_list <-
      load_enrichr_results(contrast_character = contrast_character,
                           enrichr_database = enrichr_databases[enrichr_index])

    result_tibble_up <- result_list$up
    if (nrow(x = result_tibble_up) > 0L) {
      result_tibble_up <-
        result_tibble_up[order(result_tibble_up$Combined.Score, decreasing = TRUE), , drop = FALSE]
      result_tibble_up <-
        result_tibble_up[seq_len(length.out = min(nrow(x = result_tibble_up), argument_list$maximum_terms)), , drop = FALSE]
      result_tibble_up$Direction <- "up"
    } else {
      result_tibble_up$Direction <- character(length = 0L)
    }

    result_tibble_down <- result_list$down
    if (nrow(x = result_tibble_down) > 0L) {
      result_tibble_down <-
        result_tibble_down[order(result_tibble_down$Combined.Score, decreasing = TRUE), , drop = FALSE]
      result_tibble_down <-
        result_tibble_down[seq_len(length.out = min(
          nrow(x = result_tibble_down),
          argument_list$maximum_terms
        )), , drop = FALSE]
      result_tibble_down$Direction <- "down"
    } else {
      result_tibble_down$Direction <- character(length = 0L)
    }

    result_tibble <-
      dplyr::bind_rows(result_tibble_up, result_tibble_down)
    rm(result_tibble_up, result_tibble_down)

    if (nrow(x = result_tibble) == 0L) {
      rm(result_tibble, result_list)
      next()
    }

    # Order from lowest (down-regulated) to highest (up-regulated) combined
    # score, which means that the up-regulated genes appear on top.
    # FIXME: Keep the order of scores.

    # result_tibble <- result_tibble[order(result_tibble$Combined.Score, decreasing = FALSE), , drop = FALSE]

    result_tibble$Term <-
      factor(x = result_tibble$Term,
             levels = unique(x = result_tibble$Term))

    ggplot_object <-
      ggplot2::ggplot(data = result_tibble)

    ggplot_object <- ggplot_object +
      ggplot2::geom_bar(
        mapping = ggplot2::aes(
          x = .data$Term,
          y = .data$Combined.Score,
          fill = .data$Direction
        ),
        stat = "identity",
        width = 0.5
        # position = "dodge"
      )

    ggplot_object <- ggplot_object +
      ggplot2::labs(
        x = "Term",
        y = "Combined Enrichr Score",
        fill = "Direction",
        title = "Enrichr Analysis",
        subtitle = enrichr_databases[enrichr_index]
      )

    # ggplot_object <- ggplot_object + ggplot2::facet_grid(rows = ggplot2::vars(Direction))

    ggplot_object <- ggplot_object +
      ggplot2::scale_fill_manual(
        name = "Expression",
        labels = c("Down regulated", "Up regulated"),
        values = c("down" = "#00ba38", "up" = "#f8766d")
      )

    ggplot_object <- ggplot_object + ggplot2::coord_flip()

    ggplot_object <-
      ggplot_object + ggplot2::theme(axis.text.y = ggplot2::element_text(size = ggplot2::rel(x = 0.5)))

    file_path_character <-
      paste(
        paste(
          prefix_enrichr,
          contrast_character,
          enrichr_databases[enrichr_index],
          sep = "_"
        ),
        graphics_formats,
        sep = "."
      )

    for (file_path in file_path_character) {
      ggplot2::ggsave(
        filename = file.path(output_directory, file_path),
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height
      )
    }

    nozzle_section_contrast <-
      Nozzle.R1::addTo(
        parent = nozzle_section_contrast,
        Nozzle.R1::newFigure(
          file = file_path_character[2L],
          "Enrichr combined score ",
          Nozzle.R1::asStrong(enrichr_databases[enrichr_index]),
          " for contrast ",
          Nozzle.R1::asStrong(contrast_tibble$Label[contrast_index]),
          fileHighRes = file_path_character[1L]
        )
      )

    nozzle_section_contrast <-
      Nozzle.R1::addTo(
        parent = nozzle_section_contrast,
        newParagraph(
          "Full Enrichr results ",
          Nozzle.R1::asStrong(enrichr_databases[enrichr_index]),
          " for up-regulated (",
          Nozzle.R1::asLink(url = paste(
            paste(
              prefix_enrichr,
              contrast_character,
              enrichr_databases[enrichr_index],
              "up",
              sep = "_"
            ),
            "tsv",
            sep = "."
          ), "TSV"),
          ") and down-regulated (",
          Nozzle.R1::asLink(url = paste(
            paste(
              prefix_enrichr,
              contrast_character,
              enrichr_databases[enrichr_index],
              "down",
              sep = "_"
            ),
            "tsv",
            sep = "."
          ), "TSV"),
          ") genes."
        )
      )

    rm(file_path,
       file_path_character,
       ggplot_object,
       result_tibble,
       result_list)
  }

  nozzle_section_enrichr <-
    Nozzle.R1::addTo(parent = nozzle_section_enrichr, nozzle_section_contrast)

  rm(enrichr_index, nozzle_section_contrast, contrast_character)
}
rm(contrast_index)

nozzle_report <-
  Nozzle.R1::newCustomReport("Enrichr Report", version = 0)

nozzle_report <-
  Nozzle.R1::addTo(parent = nozzle_report, nozzle_section_contrasts)

nozzle_report <-
  Nozzle.R1::addTo(parent = nozzle_report, nozzle_section_enrichr)

Nozzle.R1::writeReport(report = nozzle_report,
                       filename = file.path(output_directory,
                                            paste(prefix_enrichr,
                                                  "report",
                                                  sep = "_")))

rm(
  nozzle_report,
  nozzle_section_enrichr,
  nozzle_section_contrasts,
  contrast_tibble,
  enrichr_databases,
  output_directory,
  load_enrichr_results,
  load_deseq_result_tibble,
  prefix_enrichr,
  prefix_deseq,
  graphics_formats,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessioninfo::session_info())

#!/usr/bin/env Rscript
#
# This script uses the Bioconductor ChIPpeakAnno package to annotate the
# consensus peak set, as well as contrast peak sets of a differential binding
# analysis carried out by the Bioconductor DiffBind package via the
# bsf_chipseq_diff_bind.R script. The annotation is based on a Biocoductor TxDb
# object, as well as its corresponding GTF file that provides additional
# annotation, such as gene symbols and gene biotypes.
#
#
# Copyright 2013 - 2020 Michael K. Schuster
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

argument_list <-
  optparse::parse_args(object = optparse::OptionParser(
    option_list = list(
      optparse::make_option(
        opt_str = c("-v", "--verbose"),
        action = "store_true",
        default = TRUE,
        help = "Print extra output [default]",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = c("-q", "--quietly"),
        action = "store_false",
        default = FALSE,
        dest = "verbose",
        help = "Print little output",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = c("--comparison"),
        dest = "comparison",
        help = "Comparison name",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--factor"),
        dest = "factor",
        help = "ChIP factor",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--fdr-threshold"),
        default = 0.05,
        dest = "fdr_threshold",
        help = "FDR threshold [0.05]",
        type = "numeric"
      ),
      optparse::make_option(
        opt_str = c("--gtf-reference"),
        default = NULL,
        dest = "gtf_reference",
        help = "Reference transcriptome GTF file path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--txdb-path"),
        dest = "txdb_path",
        help = "Reference transcriptome Bioconductor TxDb file path",
        type = "character"
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

# Start of main script ----------------------------------------------------


suppressPackageStartupMessages(expr = library(package = "AnnotationDbi"))
suppressPackageStartupMessages(expr = library(package = "BiocVersion"))
suppressPackageStartupMessages(expr = library(package = "Biostrings"))
suppressPackageStartupMessages(expr = library(package = "ChIPpeakAnno"))
suppressPackageStartupMessages(expr = library(package = "DiffBind"))
suppressPackageStartupMessages(expr = library(package = "GenomicFeatures"))
suppressPackageStartupMessages(expr = library(package = "rtracklayer"))
suppressPackageStartupMessages(expr = library(package = "tidyverse"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

# Define the precedence of regions for the
# ChIPpeakAnno::assignChromosomeRegion() function so that each peak is assigned
# only once and all regions sum up to 100%.
region_precedence <- c("Promoters",
                       "immediateDownstream",
                       "fiveUTRs",
                       "threeUTRs",
                       "Exons",
                       "Introns")

prefix <-
  paste("chipseq",
        "diff",
        "bind",
        argument_list$comparison,
        argument_list$factor,
        sep = "_")

output_directory <-
  file.path(argument_list$output_directory, prefix)

if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

# Load a TxDb object ------------------------------------------------------


message("Loading a TxDb object")
txdb_object <- AnnotationDbi::loadDb(file = argument_list$txdb_path)

# Initialise a DBA object -------------------------------------------------


diffbind_dba <- NULL
file_path <-
  file.path(argument_list$genome_directory,
            prefix,
            paste0(prefix, '_DBA.RData'))
if (file.exists(file_path) &&
    file.info(file_path)$size > 0L) {
  message("Loading a DiffBind DBA object")
  diffbind_dba <-
    DiffBind::dba.load(
      dir = file.path(argument_list$genome_directory, prefix),
      pre = paste0(prefix, "_")
    )
} else {
  stop("Require a pre-calculated DiffBind DBA object in file: ",
       file_path)
}
rm(file_path)

# Get the consensus peak set ----------------------------------------------


message("Loading a DiffBind peakset (GRanges) object")
diffbind_peakset_granges <-
  DiffBind::dba.peakset(
    DBA = diffbind_dba,
    bRetrieve = TRUE,
    DataType = DiffBind::DBA_DATA_GRANGES
  )

# Assign chromosome regions -----------------------------------------------


message("Assigning chromosome regions: ",
        length(diffbind_peakset_granges))
chromosome_region_list <-
  ChIPpeakAnno::assignChromosomeRegion(peaks.RD = diffbind_peakset_granges,
                                       precedence = region_precedence,
                                       TxDb = txdb_object)

message("Plotting chromosome regions")
plot_paths <- file.path(output_directory,
                        paste(
                          paste(prefix,
                                "peak",
                                "set",
                                "regions",
                                sep = "_"),
                          graphics_formats,
                          sep = "."
                        ))

ggplot_object <-
  ggplot2::ggplot(data = tibble::tibble(
    name = dimnames(x = chromosome_region_list$percentage)$subjectHits,
    percentage = as.numeric(x = chromosome_region_list$percentage)
  ))

ggplot_object <-
  ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(x = .data$name, y = .data$percentage))

ggplot_object <-
  ggplot_object + ggplot2::labs(
    x = "Name",
    y = "Percentage",
    title = "Chromosome Regions",
    subtitle = "Combined Peak Set"
  )

ggplot_object <-
  ggplot_object + ggplot2::theme(axis.text.x = ggplot2::element_text(
    size = ggplot2::rel(x = 0.8),
    hjust = 0.0,
    vjust = 0.5,
    angle = 90.0
  ))

for (plot_path in plot_paths) {
  ggplot2::ggsave(
    filename = plot_path,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height,
    limitsize = FALSE
  )
}

# Write the chromosome region fractions to a TSV file.
readr::write_tsv(
  x = ggplot_object$data,
  file = file.path(output_directory,
                   paste(
                     paste(prefix,
                           "peak",
                           "set",
                           "regions",
                           sep = "_"),
                     "tsv",
                     sep = "."
                   )),
  col_names = TRUE
)
rm(plot_path, ggplot_object, plot_paths, chromosome_region_list)

# Initialise reference annotation -----------------------------------------


message("Reading reference GTF gene features")
gene_granges <-
  rtracklayer::import(
    con = argument_list$gtf_reference,
    format = "GTF",
    # genome = ,
    feature.type = "gene"
  )
message("Creating annotation frame")
gene_frame <-
  S4Vectors::mcols(x = gene_granges)[, c("gene_id",
                                         "gene_version",
                                         "gene_name",
                                         "gene_biotype",
                                         "gene_source")]
# Add the gene location as an Ensembl-like location, lacking the coordinate
# system name and version.
gene_frame$gene_location <-
  methods::as(object = gene_granges, Class = "character")
rm(gene_granges)

utils::write.table(
  x = gene_frame,
  file = file.path(output_directory,
                   paste(
                     paste(prefix, "gene", "set", sep = "_"), "tsv", sep = "."
                   )),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE
)

# Annotate the consensus peak set -----------------------------------------


message("Preparing annotation")
# Annotation based on a TxDb object.
# The MultiDb::OrganimsDb object would be used to provide gene symbols, but a
# matching version with regards to the assembly and the annotation is not always
# available. Hence, information is joined from the meta annotation columns of
# the GTF-based GenomicRanges::GRanges object.

annotation_granges <-
  ChIPpeakAnno::toGRanges(data = txdb_object, feature = "gene")

message("Annotating the consensus peak set: ",
        length(x = diffbind_peakset_granges))
annotated_granges <-
  ChIPpeakAnno::annotatePeakInBatch(myPeakList = diffbind_peakset_granges,
                                    AnnotationData = annotation_granges)

annotated_frame <-
  BiocGenerics::as.data.frame(x = annotated_granges, stringsAsFactors = FALSE)
rm(annotated_granges)

# Merge with the GTF meta annotation columns.
merged_frame <-
  base::merge.data.frame(x = annotated_frame,
                         y = gene_frame,
                         by.x = "feature",
                         by.y = "gene_id")
rm(annotated_frame)

# Order by the numeric "peak" variable.
merged_frame <-
  merged_frame[BiocGenerics::order(as.numeric(x = merged_frame$peak)),]

utils::write.table(
  x = merged_frame,
  file = file.path(output_directory,
                   paste(
                     paste(prefix,
                           "peak",
                           "set",
                           "genes",
                           "complete",
                           sep = "_"),
                     "tsv",
                     sep = "."
                   )),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)
rm(merged_frame)
# print(x = "Consensus peak set annotation warnings:")
# print(x = base::warnings())

#' Process a contrasts data frame obtained via DiffBind::dba.show() per row.
#'
#' @param contrast contrast frame row name indicating the contrast number
#' @param group1 A \code{character} scalar of contrast group 1
#' @param group2 A \code{character} scalar of contrast group 2
#' @param db_number An \code{integer} scalar with the number of differentially
#'   bound sites
#'
#' @return TRUE
#' @export
#'
#' @examples
process_per_contrast <-
  function(contrast, group1, group2, db_number) {
    message(
      sprintf(
        fmt = "Processing %s %s contrast %s versus %s",
        argument_list$comparison,
        argument_list$factor,
        group1,
        group2
      )
    )

    report_granges <- NULL
    base::tryCatch(
      expr = {
        report_granges <- DiffBind::dba.report(
          DBA = diffbind_dba,
          contrast = as.integer(x = contrast),
          th = 1.0,
          bNormalized = TRUE,
          bCalled = TRUE,
          bCounts = TRUE,
          bCalledDetail = TRUE,
          DataType = DiffBind::DBA_DATA_GRANGES
        )
      },
      error = function(cond) {
        message(
          "DiffBind::dba.report failed ",
          sprintf(
            "for comparison %s, factor %s and contrast %s versus %s ",
            argument_list$comparison,
            argument_list$factor,
            group1,
            group2
          ),
          "with message:\n",
          cond,
          appendLF = TRUE
        )
      }
    )

    message(
      sprintf(
        fmt = "  Annotating %s %s contrast %s versus %s peak set: %d",
        argument_list$comparison,
        argument_list$factor,
        group1,
        group2,
        length(x = report_granges)
      )
    )

    annotated_granges <-
      ChIPpeakAnno::annotatePeakInBatch(myPeakList = report_granges,
                                        AnnotationData = annotation_granges)
    # print(x = "Contrast peak set annotation warnings:")
    # print(x = base::warnings())

    annotated_frame <-
      BiocGenerics::as.data.frame(x = annotated_granges, stringsAsFactors = FALSE)
    rm(annotated_granges)

    merged_frame <-
      base::merge.data.frame(
        x = annotated_frame,
        y = gene_frame,
        by.x = "feature",
        by.y = "gene_id"
      )
    rm(annotated_frame)

    # Calculate ranks for ...
    # (1) the effect size (Fold), ...
    merged_frame$rank_fold_change <-
      base::rank(x = -abs(x = merged_frame$Fold),
                 ties.method = c("min"))

    # (2) the absolute level (Conc) and ...
    merged_frame$rank_conc <-
      base::rank(x = -merged_frame$Conc,
                 ties.method = c("min"))

    # (3) the statistical significance (padj).
    merged_frame$rank_fdr <-
      base::rank(x = merged_frame$FDR, ties.method = c("min"))

    # Order by the numeric "peak" variable.
    merged_frame <-
      merged_frame[BiocGenerics::order(as.numeric(x = merged_frame$peak)),]

    utils::write.table(
      x = merged_frame,
      file = file.path(
        output_directory,
        sprintf(
          "%s_peaks_%s__%s_genes_complete.tsv",
          prefix,
          group1,
          group2
        )
      ),
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )

    # Filter by the FDR threshold value.
    merged_frame <-
      merged_frame[merged_frame$FDR <= argument_list$fdr_threshold,]

    utils::write.table(
      x = merged_frame,
      file = file.path(
        output_directory,
        sprintf(
          "%s_peaks_%s__%s_genes_significant.tsv",
          prefix,
          group1,
          group2
        )
      ),
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )

    rm(merged_frame)

    significant_granges <-
      report_granges[report_granges$FDR <= argument_list$fdr_threshold,]
    message(sprintf(fmt = "  Assigning chromosome regions for significant peak set: %d",
                    length(x = significant_granges)))

    chromosome_region_list <-
      ChIPpeakAnno::assignChromosomeRegion(peaks.RD = significant_granges,
                                           precedence = region_precedence,
                                           TxDb = txdb_object)

    if (is.null(x = dimnames(x = chromosome_region_list$percentage)$subjectHits)) {
      warning("  Skipping an empty chromosome regions plot")
    } else {
      message("  Creating a chromosome regions plot")
      plot_paths <- file.path(
        output_directory,
        sprintf(
          "%s_peaks_%s__%s_regions.%s",
          prefix,
          group1,
          group2,
          graphics_formats
        )
      )

      ggplot_object <-
        ggplot2::ggplot(data = tibble::tibble(
          name = dimnames(x = chromosome_region_list$percentage)$subjectHits,
          percentage = as.numeric(x = chromosome_region_list$percentage)
        ))

      ggplot_object <-
        ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(x = .data$name, y = .data$percentage))

      ggplot_object <-
        ggplot_object + ggplot2::labs(
          x = "Name",
          y = "Percentage",
          title = "Chromosome Regions",
          subtitle = paste("Peak Set", group1, group2, sep = " ")
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

      for (plot_path in plot_paths) {
        ggplot2::ggsave(
          filename = plot_path,
          plot = ggplot_object,
          width = argument_list$plot_width,
          height = argument_list$plot_height,
          limitsize = FALSE
        )
      }
      # Write the fractions to a TSV file.
      readr::write_tsv(
        x = ggplot_object$data,
        file = file.path(
          output_directory,
          sprintf("%s_peaks_%s__%s_regions.tsv",
                  prefix,
                  group1,
                  group2)
        ),
        col_names = TRUE
      )
      rm(plot_path,
         ggplot_object,
         plot_paths)
    }
    rm(chromosome_region_list,
       significant_granges,
       report_granges)
  }

# Get a data frame with all contrasts to apply the above function to each row.
contrast_frame <-
  DiffBind::dba.show(DBA = diffbind_dba, bContrasts = TRUE)

# Replace '!' characters with 'not_'.

# DiffBind3 seems to use "Group", while DiffBind2 used "Group1".
group1_name <-
  if ("Group" %in% names(x = contrast_frame)) {
    "Group"
  } else {
    "Group1"
  }

contrast_frame[, group1_name] <-
  gsub(pattern = "!",
       replacement = "not_",
       x = contrast_frame[, group1_name])

contrast_frame$Group2 <-
  gsub(pattern = "!",
       replacement = "not_",
       x = contrast_frame$Group2)

return_value <-
  mapply(
    FUN = process_per_contrast,
    row.names(contrast_frame),
    contrast_frame[, group1_name],
    contrast_frame$Group2,
    # Since column 5 (DB.DESeq2) is a factor, it needs converting into a
    # character, before converting into an integer.
    as.integer(x = as.character(x = contrast_frame[, 5L]))
  )
rm(return_value, contrast_frame, group1_name, process_per_contrast)

rm(
  annotation_granges,
  gene_frame,
  diffbind_peakset_granges,
  diffbind_dba,
  txdb_object,
  output_directory,
  prefix,
  region_precedence,
  graphics_formats,
  argument_list
)

print(x = "Final warnings:")
print(x = base::warnings())

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

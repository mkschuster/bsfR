#!/usr/bin/env Rscript
#
# The bsf_chipseq_annotate.R script annotates the consensus peak set and
# contrast peak sets of a DiffBind analysis run via bsf_chipseq_diff_bind.R. The
# annotation is based on a Biocoductor TxDb object as well as GTF file that
# provides addiitonal annotation, such as gene symbols and gene biotypes.
#
#
# Copyright 2013 - 2019 Michael K. Schuster
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
      opt_str = c("--comparison"),
      dest = "comparison",
      help = "Comparison name",
      type = "character"
    ),
    make_option(
      opt_str = c("--factor"),
      dest = "factor",
      help = "ChIP factor",
      type = "character"
    ),
    make_option(
      opt_str = c("--fdr-threshold"),
      default = 0.05,
      dest = "fdr_threshold",
      help = "FDR threshold [0.05]",
      type = "numeric"
    ),
    make_option(
      opt_str = c("--gtf-reference"),
      default = NULL,
      dest = "gtf_reference",
      help = "Reference transcriptome GTF file path",
      type = "character"
    ),
    make_option(
      opt_str = c("--txdb-path"),
      dest = "txdb_path",
      help = "Reference transcriptome Bioconductor TxDb file path",
      type = "character"
    ),
    make_option(
      opt_str = c("--genome-directory"),
      default = ".",
      dest = "genome_directory",
      help = "Genome directory path [.]",
      type = "character"
    ),
    make_option(
      opt_str = c("--output-directory"),
      default = ".",
      dest = "output_directory",
      help = "Output directory path [.]",
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

# Start of main script ----------------------------------------------------


suppressPackageStartupMessages(expr = library(package = "AnnotationDbi"))
suppressPackageStartupMessages(expr = library(package = "ChIPpeakAnno"))
suppressPackageStartupMessages(expr = library(package = "DiffBind"))
suppressPackageStartupMessages(expr = library(package = "GenomicFeatures"))
suppressPackageStartupMessages(expr = library(package = "ggplot2"))
suppressPackageStartupMessages(expr = library(package = "tibble"))

# Save plots in the following formats.
graphics_formats <- c("pdf" = "pdf", "png" = "png")

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
  stop(paste0(
    "Require a pre-calculated DiffBind DBA object in file: ",
    file_path
  ))
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
  ChIPpeakAnno::assignChromosomeRegion(peaks.RD = diffbind_peakset_granges, TxDb = txdb_object)

message("Plotting chromosome regions")
plot_paths <- file.path(output_directory,
                        paste(
                          paste(prefix,
                                "peak",
                                "set",
                                "chromosome",
                                "regions",
                                sep = "_"),
                          graphics_formats,
                          sep = "."
                        ))

ggplot_object <-
  ggplot2::ggplot(data = tibble::tibble(
    name = attr(x = chromosome_region_list$percentage, which = "dimnames")$subjectHits,
    percentage = as.numeric(x = chromosome_region_list$percentage[])
  ))

ggplot_object <-
  ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(x = name, y = percentage))

ggplot_object <-
  ggplot_object + ggplot2::labs(
    x = "Name",
    y = "Percentage",
    title = "Chromosome Regions",
    subtitle = "Combined Peak Set"
  )

ggplot_object <-
  ggplot_object + ggplot2::theme(axis.text.x = ggplot2::element_text(
    size = 8.0,
    hjust = 1.0,
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
# Add the gene location as an Ensembl-like location, lacking the coordinate system name and version.
gene_frame$gene_location <-
  as(object = gene_granges, Class = "character")
rm(gene_granges)

write.table(
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
# matching version with regards to the assembly and the annotaiton is not always
# available. Hence, information is joined from the meta annotation columns of
# the GTF-based GRanges object.

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
  merged_frame[BiocGenerics::order(as.numeric(x = merged_frame$peak)), ]

write.table(
  x = merged_frame,
  file = file.path(output_directory,
                   paste(
                     paste(prefix,
                           "peak",
                           "set",
                           "genes",
                           sep = "_"),
                     "tsv",
                     sep = "."
                   )),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)
rm(merged_frame)
print(x = "Consensus peak set annotation warnings:")
print(x = base::warnings())

#' Process a contrasts data frame obtained via DiffBind::dba.show() per row.
#'
#' @param contrast contrast frame row name indicating the contrast number
#' @param group1 A \code{character} scalar of contrast group 1
#' @param group2 A \code{character} scalar of contrast group 2
#' @param db_number An \code{integer} scalar with the number of differntually bound sites
#'
#' @return TRUE
#' @export
#'
#' @examples
process_per_contrast <-
  function(contrast, group1, group2, db_number) {
    report_granges <- NULL
    tryCatch(
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
        fmt = "Annotating %s %s contrast %s versus %s peak set: %d",
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

    # Order by the numeric "peak" variable.
    merged_frame <-
      merged_frame[BiocGenerics::order(as.numeric(x = merged_frame$peak)), ]

    utils::write.table(
      x = merged_frame,
      file = file.path(
        output_directory,
        sprintf("%s_peaks_%s__%s_genes.tsv",
                prefix,
                group1,
                group2)
      ),
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )

    rm(merged_frame)

    print(x = "Contrast peak set annotation warnings:")
    print(x = base::warnings())
  }

# Get a data frame with all contrasts to apply the above function to each row.
contrast_frame <-
  DiffBind::dba.show(DBA = diffbind_dba, bContrasts = TRUE)

# Replace '!' characters with 'not_'.
contrast_frame$Group1 <-
  gsub(pattern = "!",
       replacement = "not_",
       x = contrast_frame$Group1)
contrast_frame$Group2 <-
  gsub(pattern = "!",
       replacement = "not_",
       x = contrast_frame$Group2)

# print(x = "Contrast frame")
# print(x = str(object = contrast_frame))
# print(x = row.names(x = contrast_frame))

return_value <-
  mapply(
    FUN = process_per_contrast,
    row.names(contrast_frame),
    contrast_frame$Group1,
    contrast_frame$Group2,
    # Since column 5 (DB.DESeq2) is a factor, it needs converting into a character,
    # before converting into an integer.
    as.integer(x = as.character(x = contrast_frame[, 5L]))
  )
rm(return_value, contrast_frame, process_per_contrast)

rm(
  annotation_granges,
  gene_frame,
  diffbind_peakset_granges,
  diffbind_dba,
  txdb_object,
  output_directory,
  prefix,
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

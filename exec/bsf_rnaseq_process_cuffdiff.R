#!/usr/bin/env Rscript
#
# BSF R script to post-processes Cuffdiff output utilising the cummeRbund
# package. This script creates a cummeRbund SQLite database, writes a set of QC
# plots in PDF and PNG format, splits the large and unwieldy differential data
# tables into pairwise comparisons and creates symbolic links to the original
# Cuffdiff differential tables.
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
      opt_str = c("--comparison-name"),
      dest = "comparison_name",
      help = "Comparison name",
      type = "character"
    ),
    make_option(
      opt_str = c("--gtf-assembly"),
      default = NULL,
      dest = "gtf_assembly",
      help = "GTF file specifying an assembled and merged transcriptome",
      type = "character"
    ),
    make_option(
      opt_str = c("--gtf-reference"),
      default = NULL,
      dest = "gtf_reference",
      help = "GTF file specifying a reference transcriptome",
      type = "character"
    ),
    make_option(
      opt_str = c("--genome-version"),
      default = NULL,
      dest = "genome_version",
      help = "Genome version",
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

# Validate the argument_list.

if (is.null(x = argument_list$comparison_name)) {
  stop("Missing --comparison_name option")
}

if (is.null(x = argument_list$gtf_assembly)) {
  stop("Missing --gtf_assembly option")
}

if (is.null(x = argument_list$gtf_reference)) {
  stop("Missing --gtf_reference option")
}

if (is.null(x = argument_list$genome_version)) {
  stop("Missing --genome_version option")
}

suppressPackageStartupMessages(expr = library(package = "cummeRbund"))
suppressPackageStartupMessages(expr = library(package = "rtracklayer"))

#' Convert character vector into a comma-separated value (CSV).
#' Apply unique() and sort(), before collapsing a character vector
#' into a single element, comma-separated value.
#'
#' @param x: Character vector
#' @return: Single element character vector of comma-sepatated values

character_to_csv <- function(x) {
  return(paste(sort(x = unique(x = x)), collapse = ","))
}

# Define CuffDiff and output directory names relative to the
# working directory.

cuffdiff_directory <-
  paste("rnaseq", "cuffdiff", argument_list$comparison_name, sep = "_")

output_directory <-
  paste("rnaseq",
        "process",
        "cuffdiff",
        argument_list$comparison_name,
        sep = "_")

# To avoid name clashes when downloading files, use the output directory
# name also as a prefix for all files therein.

prefix <- output_directory

# Create a cummeRbund database --------------------------------------------


# Read and index Cuffdiff output and create a CuffSet object.
# Load also the GTF file assembled by Cuffmerge, which fulfills foreign
# key constraints between the "features", "genes" and "isoforms" SQL tables.
# The CuffSet object has slots "genes", "isoforms", "TSS" and "CDS" that each
# are instances of the CuffData class. By default, CuffData accessor methods
# applied to a CuffSet class will operate on the "genes" slot.

message("Creating or loading a cummeRbund database")
cuff_set <-
  cummeRbund::readCufflinks(
    dir = cuffdiff_directory,
    gtfFile = argument_list$gtf_assembly,
    genome = argument_list$genome_version
  )

# Create a new sub-directory for plots if it does not exist.

if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

# Process Cuffdiff run information ----------------------------------------


frame_path <-
  file.path(output_directory, paste0(prefix, "_run_information.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0L) {
  message("Skipping a run information table")
} else {
  message("Creating a run information table")
  write.table(
    x = cummeRbund::runInfo(object = cuff_set),
    file = frame_path,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
}
rm(frame_path)

# Process Cuffdiff sample information -------------------------------------


sample_frame <- cummeRbund::samples(object = cuff_set)
# Get the number of samples.
sample_number <- nrow(x = sample_frame)
# Write the sample_frame table.
frame_path <-
  file.path(output_directory, paste0(prefix, "_samples.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0L) {
  message("Skipping a sample table")
} else {
  message("Creating a sample table")
  write.table(
    x = sample_frame,
    file = frame_path,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
}
rm(frame_path)
# Define an array of all possible pairwise sample comparisons or sample pairs.
sample_pairs <- combn(x = sample_frame$sample_name, m = 2L)
# Write a table of sample pair information by transposing the sample pairs
# array. This table allows the Python web code to link in pairwise plots.
frame_path <-
  file.path(output_directory, paste0(prefix, "_sample_pairs.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0L) {
  message("Skipping a sample pairs table")
} else {
  message("Creating a sample pairs table")
  write.table(
    x = aperm(a = sample_pairs),
    file = frame_path,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
}
rm(frame_path, sample_frame)

# Process Cuffdiff replicate information ----------------------------------


replicate_frame <- cummeRbund::replicates(object = cuff_set)
# Get the number of replicates.
replicate_number <- nrow(x = replicate_frame)
# Some plots require replicates. Check whether at least one row has a
# replicate column value greater than 0.
have_replicates <-
  (nrow(x = replicate_frame[replicate_frame$replicate > 0L,]) > 0L)
# Write the replicate_frame.
frame_path <-
  file.path(output_directory, paste0(prefix, "_replicates.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0L) {
  message("Skipping a replicate table")
} else {
  message("Creating a replicate table")
  write.table(
    x = replicate_frame,
    file = frame_path,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
}
rm(frame_path)
# Plot the log10(internal_scale) of the replicate_frame to
# visualise outliers.
plot_path_pdf <-
  file.path(output_directory, paste0(prefix, "_replicate_scale.pdf"))
plot_path_png <-
  file.path(output_directory, paste0(prefix, "_replicate_scale.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a library scale plot on replicates")
} else {
  message("Creating a library scale plot on replicates")
  ggplot_object <- ggplot2::ggplot(data = replicate_frame)
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(x = .data$rep_name, y = log10(.data$internal_scale)))
  ggplot_object <-
    ggplot_object + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -90, hjust = 0))
  plot_width <-
    argument_list$plot_width + (ceiling(x = replicate_number / 24L) - 1L) * argument_list$plot_width * 0.25
  ggplot2::ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  ggplot2::ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object, plot_width)
}
rm(plot_path_pdf, plot_path_png, replicate_frame)

# Starting QC plotting ----------------------------------------------------


message("Starting QC plotting")

# Dispersion Plot on Genes ------------------------------------------------


plot_path_pdf <-
  file.path(output_directory, paste0(prefix, "_genes_dispersion.pdf"))
plot_path_png <-
  file.path(output_directory, paste0(prefix, "_genes_dispersion.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a Dispersion Plot on Genes")
} else {
  message("Creating a Dispersion Plot on Genes")
  ggplot_object <-
    cummeRbund::dispersionPlot(object = cummeRbund::genes(object = cuff_set))
  ggplot2::ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggplot2::ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

# Dispersion Plot on Isoforms ---------------------------------------------


plot_path_pdf <-
  file.path(output_directory, paste0(prefix, "_isoforms_dispersion.pdf"))
plot_path_png <-
  file.path(output_directory, paste0(prefix, "_isoforms_dispersion.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a Dispersion Plot on Isoforms")
} else {
  message("Creating a Dispersion Plot on Isoforms")
  base::tryCatch(
    expr = {
      ggplot_object <-
        cummeRbund::dispersionPlot(object = cummeRbund::isoforms(object = cuff_set))
      ggplot2::ggsave(
        filename = plot_path_pdf,
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height
      )
      ggplot2::ggsave(
        filename = plot_path_png,
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height
      )
      rm(ggplot_object)
    },
    error = function(cond) {
      message("Dispersion Plot on Isoforms failed with message:")
      message(cond, appendLF = TRUE)
      return(NULL)
    }
  )
}
rm(plot_path_pdf, plot_path_png)

# Squared Coefficient of Variation (SCV) Plot on Genes --------------------


plot_path_pdf <-
  file.path(output_directory, paste0(prefix, "_genes_scv.pdf"))
plot_path_png <-
  file.path(output_directory, paste0(prefix, "_genes_scv.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a SCV Plot on Genes")
} else {
  # The plot requires replicates.
  if (have_replicates) {
    message("Creating a SCV Plot on Genes")
    ggplot_object <-
      cummeRbund::fpkmSCVPlot(object = cummeRbund::genes(object = cuff_set))
    ggplot2::ggsave(
      filename = plot_path_pdf,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ggplot2::ggsave(
      filename = plot_path_png,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    rm(ggplot_object)
  } else {
    message("Skipping a SCV Plot on Genes in lack of replicates")
  }
}
rm(plot_path_pdf, plot_path_png)

# Squared Coefficient of Variation (SCV) Plot on Isoforms -----------------


plot_path_pdf <-
  file.path(output_directory, paste0(prefix, "_isoforms_scv.pdf"))
plot_path_png <-
  file.path(output_directory, paste0(prefix, "_isoforms_scv.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a SCV Plot on Isoforms")
} else {
  # The plot requires replicates.
  if (have_replicates) {
    message("Creating a SCV Plot on Isoforms")
    ggplot_object <-
      cummeRbund::fpkmSCVPlot(object = cummeRbund::isoforms(object = cuff_set))
    ggplot2::ggsave(
      filename = plot_path_pdf,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ggplot2::ggsave(
      filename = plot_path_png,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    rm(ggplot_object)
  } else {
    message("Skipping a SCV Plot on Isoforms in lack of replicates")
  }
}
rm(plot_path_pdf, plot_path_png)

# Density Plot on Genes without replicates --------------------------------


plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_genes_density_wo_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_genes_density_wo_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a Density Plot on Genes without replicates")
} else {
  message("Creating a Density Plot on Genes without replicates")
  ggplot_object <-
    cummeRbund::csDensity(object = cummeRbund::genes(object = cuff_set),
                          replicates = FALSE)
  ggplot2::ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggplot2::ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

# Density Plot on Genes with replicates -----------------------------------


plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_genes_density_w_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_genes_density_w_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a Density Plot on Genes with replicates")
} else {
  message("Creating a Density Plot on Genes with replicates")
  ggplot_object <-
    cummeRbund::csDensity(object = cummeRbund::genes(object = cuff_set),
                          replicates = TRUE)
  ggplot2::ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggplot2::ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

# Density Plot on Isoforms without replicates -----------------------------


plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_density_wo_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_density_wo_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a Density Plot on Isoforms without replicates")
} else {
  message("Creating a Density Plot on Isoforms without replicates")
  ggplot_object <-
    cummeRbund::csDensity(object = cummeRbund::isoforms(object = cuff_set),
                          replicates = FALSE)
  ggplot2::ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggplot2::ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

# Density Plot on Isoforms with replicates --------------------------------


plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_density_w_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_density_w_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a Density Plot on Isoforms with replicates")
} else {
  message("Creating a Density Plot on Isoforms with replicates")
  ggplot_object <-
    cummeRbund::csDensity(object = cummeRbund::isoforms(object = cuff_set),
                          replicates = TRUE)
  ggplot2::ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggplot2::ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

# Box Plot on Genes with replicates ---------------------------------------


plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_genes_box_w_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_genes_box_w_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a Box Plot on Genes with replicates")
} else {
  message("Creating a Box Plot on Genes with replicates")
  # By default, the csBoxplot function adds a pseudocount of 1e-04 for log10
  # transformed fpkms. In case of a large number of missing values in
  # shallowly sequenced samples, this affects the plot.
  # Hence reproduce the plot here.
  # ggplot_object <-
  #   cummeRbund::csBoxplot(object = cummeRbund::genes(object = cuff_set),
  #                         replicates = TRUE)
  rep_fpkm_genes <-
    cummeRbund::repFpkm(object = cummeRbund::genes(object = cuff_set))
  # Rename the "rep_name" column into "condition".
  colnames(x = rep_fpkm_genes)[colnames(x = rep_fpkm_genes) == "rep_name"] <-
    "condition"
  ggplot_object <- ggplot2::ggplot(data = rep_fpkm_genes)
  ggplot_object <-
    ggplot_object + ggplot2::geom_boxplot(
      mapping = ggplot2::aes(
        x = .data$condition,
        y = log10(.data$fpkm),
        fill = .data$condition
      ),
      size = 0.3,
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = -90, hjust = 0),
      legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.8))
    )
  ggplot_object <-
    ggplot_object + ggplot2::scale_fill_hue(l = 50, h.start = 200)
  # Arrange a maximum of 24 replicates in each guide column.
  # ggplot_object <-
  #   ggplot_object + ggplot2::guides(fill = ggplot2::guide_legend(nrow = 24))
  # Use the base plot_witdh and add a quarter of the width for each additional
  # guide legend column.
  plot_width <-
    argument_list$plot_width + (ceiling(x = replicate_number / 24L) - 1L) * argument_list$plot_width * 0.25
  ggplot2::ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  ggplot2::ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object, plot_width, rep_fpkm_genes)
}
rm(plot_path_pdf, plot_path_png)

# Box Plot on Genes without replicates ------------------------------------


plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_genes_box_wo_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_genes_box_wo_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a Box Plot on Genes without replicates")
} else {
  message("Creating a Box Plot on Genes without replicates")
  # ggplot_object <-
  #   cummeRbund::csBoxplot(object = cummeRbund::genes(object = cuff_set),
  #                         replicates = FALSE)
  fpkm_genes <-
    cummeRbund::fpkm(object = cummeRbund::genes(object = cuff_set))
  # Rename the "sample_name" column into "condition".
  colnames(fpkm_genes)[colnames(fpkm_genes) == "sample_name"] <-
    "condition"
  ggplot_object <- ggplot2::ggplot(data = fpkm_genes)
  ggplot_object <-
    ggplot_object + ggplot2::geom_boxplot(
      mapping = ggplot2::aes(
        x = .data$condition,
        y = log10(.data$fpkm),
        fill = .data$condition
      ),
      size = 0.3,
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = -90, hjust = 0),
      legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.8))
    )
  ggplot_object <-
    ggplot_object + ggplot2::scale_fill_hue(l = 50, h.start = 200)
  # Arrange a maximum of 24 samples in each guide column.
  # ggplot_object <-
  #   ggplot_object + ggplot2::guides(fill = ggplot2::guide_legend(nrow = 24))
  # Use the base plot_witdh and add a quarter of the width for each additional
  # guide legend column.
  plot_width <-
    argument_list$plot_width + (ceiling(x = sample_number / 24L) - 1L) * argument_list$plot_width * 0.25
  ggplot2::ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  ggplot2::ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object, plot_width, fpkm_genes)
}
rm(plot_path_pdf, plot_path_png)

# Box Plot on Isoforms with replicates ------------------------------------


plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_box_w_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_box_w_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a Box Plot on Isoforms with replicates")
} else {
  message("Creating a Box Plot on Isoforms with replicates")
  # By default, the csBoxplot function adds a pseudocount of 1e-04 for log10
  # transformed fpkms. In case of a large number of missing values in
  # shallowly sequenced samples, this affects the plot.
  # Hence reproduce the plot here.
  # ggplot_object <-
  #   cummeRbund::csBoxplot(object = cummeRbund::isoforms(object = cuff_set),
  #                         replicates = TRUE)
  rep_fpkm_isoforms <-
    cummeRbund::repFpkm(object = cummeRbund::isoforms(object = cuff_set))
  # Rename the "rep_name" column into "condition".
  colnames(rep_fpkm_isoforms)[colnames(rep_fpkm_isoforms) == "rep_name"] <-
    "condition"
  ggplot_object <- ggplot2::ggplot(data = rep_fpkm_isoforms)
  ggplot_object <-
    ggplot_object + ggplot2::geom_boxplot(
      mapping = ggplot2::aes(
        x = .data$condition,
        y = log10(.data$fpkm),
        fill = .data$condition
      ),
      size = 0.3,
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = -90, hjust = 0),
      legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.8))
    )
  ggplot_object <-
    ggplot_object + ggplot2::scale_fill_hue(l = 50, h.start = 200)
  # Arrange a maximum of 24 replicates in each guide column.
  # ggplot_object <-
  #   ggplot_object + ggplot2::guides(fill = ggplot2::guide_legend(nrow = 24))
  # Use the base plot_witdh and add a quarter of the width for each additional
  # guide legend column.
  plot_width <-
    argument_list$plot_width + (ceiling(x = replicate_number / 24L) - 1L) * argument_list$plot_width * 0.25
  ggplot2::ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  ggplot2::ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object, plot_width, rep_fpkm_isoforms)
}
rm(plot_path_pdf, plot_path_png)

# Box Plot on Isoforms without replicates ---------------------------------


plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_box_wo_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_box_wo_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a Box Plot on Isoforms without replicates")
} else {
  message("Creating a Box Plot on Isoforms without replicates")
  # ggplot_object <-
  #   cummeRbund::csBoxplot(object = cummeRbund::isoforms(object = cuff_set),
  #                         replicates = FALSE)
  fpkm_isoforms <-
    cummeRbund::fpkm(object = cummeRbund::isoforms(object = cuff_set))
  # Rename the "sample_name" column into "condition".
  colnames(fpkm_isoforms)[colnames(fpkm_isoforms) == "sample_name"] <-
    "condition"
  ggplot_object <- ggplot2::ggplot(data = fpkm_isoforms)
  ggplot_object <-
    ggplot_object + ggplot2::geom_boxplot(
      mapping = ggplot2::aes(
        x = .data$condition,
        y = log10(.data$fpkm),
        fill = .data$condition
      ),
      size = 0.3,
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = -90, hjust = 0),
      legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.8))
    )
  ggplot_object <-
    ggplot_object + ggplot2::scale_fill_hue(l = 50, h.start = 200)
  # Arrange a maximum of 24 replicates in each guide column.
  # ggplot_object <-
  #   ggplot_object + ggplot2::guides(fill = ggplot2::guide_legend(nrow = 24))
  # Use the base plot_witdh and add a quarter of the with for each additional
  # guide legend column.
  plot_width <-
    argument_list$plot_width + (ceiling(x = sample_number / 24L) - 1L) * argument_list$plot_width * 0.25
  ggplot2::ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  ggplot2::ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object, plot_width, fpkm_isoforms)
}
rm(plot_path_pdf, plot_path_png)

# Scatter Matrix Plot on Genes --------------------------------------------


# Only include the plot for less than or equal to 20 samples.
if (sample_number <= 20L) {
  plot_path_pdf <-
    file.path(output_directory,
              paste0(prefix, "_genes_scatter_matrix.pdf"))
  plot_path_png <-
    file.path(output_directory,
              paste0(prefix, "_genes_scatter_matrix.png"))
  if (file.exists(plot_path_pdf) &&
      (file.info(plot_path_pdf)$size > 0L) &&
      file.exists(plot_path_png) &&
      (file.info(plot_path_png)$size > 0L)) {
    message("Skipping a Scatter Matrix Plot on Genes")
  } else {
    message("Creating a Scatter Matrix Plot on Genes")
    ggplot_object <-
      cummeRbund::csScatterMatrix(object = cummeRbund::genes(object = cuff_set))
    ggplot2::ggsave(
      filename = plot_path_pdf,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ggplot2::ggsave(
      filename = plot_path_png,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    rm(ggplot_object)
  }
  rm(plot_path_pdf, plot_path_png)
}

# Scatter Matrix Plot on Isoforms -----------------------------------------


if (sample_number <= 20L) {
  plot_path_pdf <-
    file.path(output_directory,
              paste0(prefix, "_isoforms_scatter_matrix.pdf"))
  plot_path_png <-
    file.path(output_directory,
              paste0(prefix, "_isoforms_scatter_matrix.png"))
  if (file.exists(plot_path_pdf) &&
      (file.info(plot_path_pdf)$size > 0L) &&
      file.exists(plot_path_png) &&
      (file.info(plot_path_png)$size > 0L)) {
    message("Skipping a Scatter Matrix Plot on Isoforms")
  } else {
    message("Creating a Scatter Matrix Plot on Isoforms")
    ggplot_object <-
      cummeRbund::csScatterMatrix(object = cummeRbund::isoforms(object = cuff_set))
    ggplot2::ggsave(
      filename = plot_path_pdf,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ggplot2::ggsave(
      filename = plot_path_png,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    rm(ggplot_object)
  }
  rm(plot_path_pdf, plot_path_png)
}

# Scatter Plot on Genes for each sample pair ------------------------------


for (i in seq_along(along.with = sample_pairs[1L,])) {
  plot_path_pdf <-
    file.path(
      output_directory,
      paste(
        prefix,
        sample_pairs[1L, i],
        sample_pairs[2L, i],
        "genes_scatter.pdf",
        sep = "_"
      )
    )
  plot_path_png <-
    file.path(
      output_directory,
      paste(
        prefix,
        sample_pairs[1L, i],
        sample_pairs[2L, i],
        "genes_scatter.png",
        sep = "_"
      )
    )
  if (file.exists(plot_path_pdf) &&
      (file.info(plot_path_pdf)$size > 0L) &&
      file.exists(plot_path_png) &&
      (file.info(plot_path_png)$size > 0L)) {
    message("Skipping a Scatter Plot on Genes for ",
            sample_pairs[1L, i],
            " versus ",
            sample_pairs[2L, i])
  } else {
    message("Creating a Scatter Plot on Genes for ",
            sample_pairs[1L, i],
            " versus ",
            sample_pairs[2L, i])

    # Unfortunately, the standard cummeRbund csScatter() function
    # does not allow colouring by significance.
    # ggplot_object <-
    #   cummeRbund::csScatter(
    #     object = cummeRbund::genes(object = cuff_set),
    #     x = sample_pairs[1L, i],
    #     y = sample_pairs[2L, i],
    #     colorByStatus = TRUE
    #   )

    # Re-implement scatter plots here.
    diff_data_genes <-
      cummeRbund::diffData(
        object = cummeRbund::genes(object = cuff_set),
        x = sample_pairs[1L, i],
        y = sample_pairs[2L, i],
        features = FALSE
      )
    ggplot_object <-
      ggplot2::ggplot(data = diff_data_genes,
                      mapping = ggplot2::aes(x = .data$value_1, y = .data$value_2))
    ggplot_object <- ggplot_object + ggplot2::theme_light()
    ggplot_object <-
      ggplot_object + ggplot2::labs(x = sample_pairs[1L, i], y = sample_pairs[2L, i])
    if (TRUE) {
      # Plot the non-significant genes with ggplot2::geom_hex(),
      # which performs much better with typical numbers of genes.
      ggplot_object <-
        ggplot_object + ggplot2::geom_hex(
          data = diff_data_genes[diff_data_genes$significant == "no",],
          alpha = I(1 / 3),
          show.legend = TRUE,
          bins = 50
        )
      # Manually set scale colours.
      ggplot_object <-
        ggplot_object +
        ggplot2::scale_fill_continuous(low = "#e6f0ff", high = "#0066ff")
      # Plot the significant genes with ggplot2::geom_point() in red.
      ggplot_object <-
        ggplot_object + ggplot2::geom_point(
          data = diff_data_genes[diff_data_genes$significant == "yes",],
          colour = "red",
          size = 1.2,
          alpha = I(1 / 3)
        )
    } else {
      # Plot significant and non-significant genes with ggplot2::geom_point().
      ggplot_object <-
        ggplot_object + ggplot2::geom_point(
          mapping = ggplot2::aes(colour = .data$significant),
          size = 1.2,
          alpha = I(1 / 3)
        )
      # Manually set scale colours.
      ggplot_object <-
        ggplot_object + ggplot2::scale_colour_manual(values = c("black", "red"))
    }
    ggplot_object <-
      ggplot_object + ggplot2::geom_abline(intercept = 0,
                                           slope = 1,
                                           linetype = 2)
    ggplot_object <-
      ggplot_object + ggplot2::geom_rug(size = 0.8, alpha = 0.01)
    # ggplot_object <-
    #   ggplot_object + ggplot2::stat_smooth(method = "lm", fill = "blue", alpha = 0.2)
    ggplot_object <-
      ggplot_object + ggplot2::scale_y_log10() + ggplot2::scale_x_log10()

    # Annotate the plot with the (Pearson) correlation coefficient and the
    # number of significantly up and downregulated genes.
    # For defining the data range for label placement,
    # eliminate rows with value 0.0.
    range_value_1 <-
      range(diff_data_genes[diff_data_genes$value_1 > 0.0,]$value_1)
    range_value_2 <-
      range(diff_data_genes[diff_data_genes$value_2 > 0.0,]$value_2)

    # Calculate the (Pearson) correlation coefficient and place it on the plot
    # in the lower right corner at 95% x and 5% y.
    ggplot_object <-
      ggplot_object + ggplot2::annotate(
        geom = "text",
        x = 10L ^ (log10(x = range_value_1[1L]) + (
          log10(x = range_value_1[2L]) - log10(x = range_value_1[1L])
        ) * 0.95),
        y = 10L ^ (log10(x = range_value_2[1L]) + (
          log10(x = range_value_2[2L]) - log10(x = range_value_2[1L])
        ) * 0.05),
        label = paste0("r = ", round(
          x = cor(x = diff_data_genes$value_1, y = diff_data_genes$value_2),
          digits = 3L
        ))
      )

    # Calculate the numbers of up and down regulated genes and place them
    # on the plot in the upper left corner at 5%x and 95% y.
    ggplot_object <-
      ggplot_object + ggplot2::annotate(
        geom = "text",
        x = 10L ^ (log10(x = range_value_1[1L]) + (
          log10(x = range_value_1[2L]) - log10(x = range_value_1[1L])
        ) * 0.05),
        y = 10L ^ (log10(x = range_value_2[1L]) + (
          log10(x = range_value_2[2L]) - log10(x = range_value_2[1L])
        ) * 0.95),
        label = paste0(
          "Up: ",
          nrow(x = diff_data_genes[diff_data_genes$log2_fold_change > 0.0 &
                                     diff_data_genes$significant == "yes",]),
          "\n",
          "Down: ",
          nrow(x = diff_data_genes[diff_data_genes$log2_fold_change < 0.0 &
                                     diff_data_genes$significant == "yes",])
        ),
        colour = "red",
        hjust = 0
      )

    rm(range_value_1, range_value_2)

    ggplot2::ggsave(
      filename = plot_path_pdf,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ggplot2::ggsave(
      filename = plot_path_png,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    rm(ggplot_object, diff_data_genes)
  }
  rm(plot_path_pdf, plot_path_png)
}

# Dendrogram Plot on Genes ------------------------------------------------


# The csDendro function returns a dendrogram object that cannot be saved with
# the ggplot2::ggsave() function.
plot_path_pdf <-
  file.path(output_directory, paste0(prefix, "_genes_dendrogram.pdf"))
if (file.exists(plot_path_pdf) &&
    file.info(plot_path_pdf)$size > 0L) {
  message("Skipping a Dendrogram Plot on Genes [PDF]")
} else {
  message("Creating a Dendrogram Plot on Genes [PDF]")
  grDevices::pdf(
    file = plot_path_pdf,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  cummeRbund::csDendro(object = cummeRbund::genes(object = cuff_set))
  base::invisible(x = grDevices::dev.off())
}
rm(plot_path_pdf)

plot_path_png <-
  file.path(output_directory, paste0(prefix, "_genes_dendrogram.png"))
if (file.exists(plot_path_png) &&
    file.info(plot_path_png)$size > 0L) {
  message("Skipping a Dendrogram Plot on Genes [PNG]")
} else {
  message("Creating a Dendrogram Plot on Genes [PNG]")
  grDevices::png(
    filename = plot_path_png,
    width = argument_list$plot_width,
    height = argument_list$plot_height,
    units = "in",
    res = 300L
  )
  cummeRbund::csDendro(object = cummeRbund::genes(object = cuff_set))
  base::invisible(x = grDevices::dev.off())
}
rm(plot_path_png)

# MA Plot on Genes for each sample pair based on FPKM values --------------


for (i in seq_along(along.with = sample_pairs[1L,])) {
  plot_path_pdf <-
    file.path(
      output_directory,
      paste(prefix, sample_pairs[1L, i], sample_pairs[2L, i], "genes_ma.pdf", sep = "_")
    )
  plot_path_png <-
    file.path(
      output_directory,
      paste(prefix, sample_pairs[1L, i], sample_pairs[2L, i], "genes_ma.png", sep = "_")
    )
  if (file.exists(plot_path_pdf) &&
      (file.info(plot_path_pdf)$size > 0L) &&
      file.exists(plot_path_png) &&
      (file.info(plot_path_png)$size > 0L)) {
    message("Skipping a MAplot on Genes for ",
            sample_pairs[1L, i],
            " versus ",
            sample_pairs[2L, i])
  } else {
    message("Creating a MAplot on Genes for ",
            sample_pairs[1L, i],
            " versus ",
            sample_pairs[2L, i])
    ggplot_object <-
      cummeRbund::MAplot(
        object = cummeRbund::genes(object = cuff_set),
        x = sample_pairs[1L, i],
        y = sample_pairs[2L, i]
      )
    ggplot2::ggsave(
      filename = plot_path_pdf,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ggplot2::ggsave(
      filename = plot_path_png,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    rm(ggplot_object)
  }
  rm(plot_path_pdf, plot_path_png)
}

# TODO: Create a MAplot on genes for each sample pair based on count data.

# Volcano Matrix Plot on Genes --------------------------------------------


plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_genes_volcano_matrix.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_genes_volcano_matrix.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a Volcano Matrix Plot on Genes")
} else {
  message("Creating a Volcano Matrix Plot on Genes")
  ggplot_object <-
    cummeRbund::csVolcanoMatrix(object = cummeRbund::genes(object = cuff_set))
  ggplot2::ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggplot2::ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

# Volcano Plot on Genes for each sample pair ------------------------------


for (i in seq_along(along.with = sample_pairs[1L,])) {
  plot_path_pdf <-
    file.path(
      output_directory,
      paste(
        prefix,
        sample_pairs[1L, i],
        sample_pairs[2L, i],
        "genes_volcano.pdf",
        sep = "_"
      )
    )
  plot_path_png <-
    file.path(
      output_directory,
      paste(
        prefix,
        sample_pairs[1L, i],
        sample_pairs[2L, i],
        "genes_volcano.png",
        sep = "_"
      )
    )
  if (file.exists(plot_path_pdf) &&
      (file.info(plot_path_pdf)$size > 0L) &&
      file.exists(plot_path_png) &&
      (file.info(plot_path_png)$size > 0L)) {
    message("Skipping a Volcano Plot on Genes for ",
            sample_pairs[1L, i],
            " versus ",
            sample_pairs[2L, i])
  } else {
    message("Creating a Volcano Plot on Genes for ",
            sample_pairs[1L, i],
            " versus ",
            sample_pairs[2L, i])
    # FIXME: The definition of the generic function "csVolcano" does not
    # include the "alpha" and "showSignificant" options, although the function
    # definition contains them. It does not seem that the option defaults are
    # used.
    # In R/methods-CuffData.R:
    # .volcano <-
    #   function(object, x, y, alpha = 0.05, showSignificant = TRUE,
    #            features = FALSE, xlimits = c(-20, 20), ...)
    # setMethod("csVolcano", signature(object = "CuffData"), .volcano)
    # In R/AllGenerics.R:
    # setGeneric("csVolcano",
    #            function(object, x, y, features=F, ...) standardGeneric("csVolcano"))
    ggplot_object <-
      cummeRbund::csVolcano(
        object = cummeRbund::genes(object = cuff_set),
        x = sample_pairs[1L, i],
        y = sample_pairs[2L, i],
        alpha = 0.05,
        showSignificant = TRUE
      )
    ggplot2::ggsave(
      filename = plot_path_pdf,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ggplot2::ggsave(
      filename = plot_path_png,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    rm(ggplot_object)
  }
  rm(plot_path_pdf, plot_path_png)
}

# Multidimensional Scaling (MDS) Plot on Genes ----------------------------


# Plot only, if the CuffData object contains more than two replicates.
if (replicate_number > 2L) {
  plot_path_pdf <-
    file.path(output_directory, paste0(prefix, "_genes_mds.pdf"))
  plot_path_png <-
    file.path(output_directory, paste0(prefix, "_genes_mds.png"))
  if (file.exists(plot_path_pdf) &&
      (file.info(plot_path_pdf)$size > 0L) &&
      file.exists(plot_path_png) &&
      (file.info(plot_path_png)$size > 0L)) {
    message("Skipping a Multidimensional Scaling Plot on Genes")
  } else {
    # if (have_replicates) {
    message("Creating a Multidimensional Scaling Plot on Genes")
    # Nothing ever is simple. If the set has too many replicates, the
    # standard cummeRbund MDSplot() falls down.
    if (replicate_number <= 24L) {
      ggplot_object <-
        cummeRbund::MDSplot(object = cummeRbund::genes(object = cuff_set),
                            replicates = TRUE)
      ggplot2::ggsave(
        filename = plot_path_pdf,
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height
      )
      ggplot2::ggsave(
        filename = plot_path_png,
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height
      )
      rm(ggplot_object)
    } else {
      # The standard MDSplot has too many replicates.
      gene_rep_fit <-
        stats::cmdscale(
          d = cummeRbund::JSdist(mat = cummeRbund::makeprobs(
            a = cummeRbund::repFpkmMatrix(object = cummeRbund::genes(object = cuff_set))
          )),
          eig = TRUE,
          k = 2
        )
      gene_rep_res <-
        data.frame(
          "names" = rownames(gene_rep_fit$points),
          "M1" = gene_rep_fit$points[, 1L],
          "M2" = gene_rep_fit$points[, 2L],
          stringsAsFactors = FALSE
        )
      ggplot_object <- ggplot2::ggplot(data = gene_rep_res)
      ggplot_object <-
        ggplot_object + ggplot2::theme_bw()  # Add theme black and white.
      ggplot_object <-
        ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(
          x = .data$M1,
          y = .data$M2,
          colour = .data$names
        ))  # Draw points in any case.
      if (replicate_number <= 24L) {
        # Only add text for a sensible number of replicates i.e. less than or
        # equal to 24.
        ggplot_object <-
          ggplot_object + ggplot2::geom_text(
            mapping = ggplot2::aes(
              x = .data$M1,
              y = .data$M2,
              label = .data$names,
              colour = .data$names
            ),
            size = 4
          )
      }
      # Arrange a maximum of 24 replicates in each guide column.
      ggplot_object <-
        ggplot_object + ggplot2::guides(colour = ggplot2::guide_legend(nrow = 24))
      # Use the base plot_witdh and add a quarter of the width for each
      # additional guide legend column.
      plot_width <-
        argument_list$plot_width + (ceiling(x = replicate_number / 24L) - 1L) * argument_list$plot_width * 0.25
      ggplot2::ggsave(
        filename = plot_path_pdf,
        plot = ggplot_object,
        width = plot_width,
        height = argument_list$plot_height
      )
      ggplot2::ggsave(
        filename = plot_path_png,
        plot = ggplot_object,
        width = plot_width,
        height = argument_list$plot_height
      )
      rm(ggplot_object, plot_width, gene_rep_res, gene_rep_fit)
    }
  }
  rm(plot_path_pdf, plot_path_png)
} else {
  message("Skipping Multidimensional Scaling Plot on genes in lack of sufficient replicates")
}

# Principal Component Analysis (PCA) Plot on Genes ------------------------


# TODO: Add also other principal components or even better,
# use plots of the PCA package?
plot_path_pdf <-
  file.path(output_directory, paste0(prefix, "_genes_pca.pdf"))
plot_path_png <-
  file.path(output_directory, paste0(prefix, "_genes_pca.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a Principal Component Analysis Plot (PCA) on Genes")
} else {
  message("Creating a Principal Component Analysis Plot (PCA) on Genes")
  ggplot_object <-
    cummeRbund::PCAplot(
      object = cummeRbund::genes(object = cuff_set),
      x = "PC1",
      y = "PC2",
      replicates = TRUE
    )
  ggplot2::ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggplot2::ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

# Finishing QC plotting ---------------------------------------------------


message("Finishing QC plotting")

# Starting data table splitting -------------------------------------------


# TODO: Process gene and isoform sets and lists. Feature-level information can
# be accessed directly from a CuffData object using the fpkm, repFpkm, count,
# diffData, or annotation methods.

# Split the large and unwieldy differential data tables into
# pairwise comparisons, but read or assemble gene and isoform
# annotation first.

# Gene and Isoform Annotation ---------------------------------------------


gene_annotation_frame <- NULL
isoform_annotation_frame <- NULL
frame_path_genes <-
  file.path(output_directory, paste0(prefix, "_genes_annotation.tsv"))
frame_path_isoforms <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_annotation.tsv"))
if (file.exists(frame_path_genes) &&
    file.info(frame_path_genes)$size > 0L) {
  message("Reading gene annotation from file")
  gene_annotation_frame <-
    read.table(file = frame_path_genes,
               header = TRUE,
               sep = "\t")
  isoform_annotation_frame <-
    read.table(file = frame_path_isoforms,
               header = TRUE,
               sep = "\t")
} else {
  # Unfortunately, cummeRbund does not seem to offer a convenient way to
  # correlate Cufflinks (XLOC) gene identifiers with the original
  # Ensembl (ENSG) gene identifiers.
  # Thus, the following steps are neccessary:
  #
  # 1. Read the reference GTF to correlate Ensembl (ENSG) gene and
  # (ENST) transcript identifiers, as well as official gene symbols.
  # Read only features of type "transcript".
  message("Reading reference transcriptome annotation")
  reference_granges <- rtracklayer::import(
    con = argument_list$gtf_reference,
    format = "gtf",
    genome = argument_list$genome_version,
    feature.type = "transcript"
  )
  # Selecting via S4Vectors::mcols()[] returns a S4Vectors::DataFrame object.
  reference_frame <- unique.data.frame(
    x = data.frame(
      "ensembl_gene_id" = S4Vectors::mcols(x = reference_granges)$gene_id,
      "ensembl_transcript_id" = S4Vectors::mcols(x = reference_granges)$transcript_id,
      stringsAsFactors = TRUE
    )
  )
  rm(reference_granges)

  # 2. Read the assembly GTF, which specifies Cufflinks (XLOC) gene and
  # Cufflinks (TCONS) transcript identifiers, but does no longer provide
  # information about Ensembl (ENSG) gene identifiers. Reference transcriptome
  # loci may have been merged or split by Cuffmerge.
  # The Cuffmerge GTF file has only features of type "exon".
  message("Reading assembled transcriptome annotation")
  assembly_granges <- rtracklayer::import(
    con = argument_list$gtf_assembly,
    format = "gtf",
    genome = argument_list$genome_version,
    feature.type = "exon"
  )
  if ("nearest_ref" %in% names(x = S4Vectors::mcols(x = assembly_granges))) {
    # If a "nearest_ref" variable is defined, the GTF is a Cuffmerge assembly.
    #
    # Example: gene_id "XLOC_000001"; transcript_id "TCONS_00000001"; exon_number "1";
    #          gene_name "DDX11L1"; oId "CUFF.1.2"; nearest_ref "ENST00000450305";
    #          class_code "="; tss_id "TSS1";
    #
    # Selecting via S4Vectors::mcols()[] returns a S4Vectors::DataFrame object.
    assembly_frame <- unique.data.frame(
      x = data.frame(
        "gene_id" = S4Vectors::mcols(x = assembly_granges)$gene_id,
        "transcript_id" = S4Vectors::mcols(x = assembly_granges)$transcript_id,
        "gene_name" = S4Vectors::mcols(x = assembly_granges)$gene_name,
        "ensembl_transcript_id" = S4Vectors::mcols(x = assembly_granges)$nearest_ref,
        stringsAsFactors = TRUE
      )
    )

    # 3. Join the reference and assembly data frames to correlate
    # Cufflinks (XLOC) gene identifiers with Ensembl (ENSG) gene identifiers
    # via Ensembl (ENST) transcript identifiers.
    message("Merging reference and assembled transcriptome annotation")
    ensembl_frame <-
      merge(
        x = assembly_frame,
        y = reference_frame,
        by = "ensembl_transcript_id",
        # Keep all XLOC entries.
        all.x = TRUE
      )
    rm(reference_frame, assembly_frame)

    # 4. Create a new, normalised Ensembl annotation data frame that correlates
    # each Cufflinks (XLOC) gene identifier with comma-separated lists of one or
    # more Ensembl (ENSG) gene and Ensembl (ENST) transcript identifiers.
    message("Aggregating Ensembl annotation by gene_id")
    aggregate_frame <-
      aggregate.data.frame(x = ensembl_frame,
                           by = list(ensembl_frame$gene_id),
                           FUN = paste)
    rm(ensembl_frame)
    # The aggregate frame consists of list objects of character vectors
    # with one or more elements aggregated in each group.
    # For each character vector, the character_to_csv() function finds unique
    # elements and sorts them, before collapsing them into a single,
    # comma-separated value.
    message("Collapsing aggregated Ensembl annotation")
    ensembl_annotation <-
      data.frame(
        "gene_id" = unlist(x = lapply(
          X = aggregate_frame$gene_id, FUN = character_to_csv
        )),
        "transcript_ids" = unlist(
          x = lapply(X = aggregate_frame$transcript_id, FUN = character_to_csv)
        ),
        # "gene_names" = unlist(x = lapply(X = aggregate_frame$gene_name, FUN = character_to_csv)),
        "ensembl_gene_ids" = unlist(
          x = lapply(X = aggregate_frame$ensembl_gene_id, FUN = character_to_csv)
        ),
        "ensembl_transcript_ids" = unlist(
          x = lapply(X = aggregate_frame$ensembl_transcript_id, FUN = character_to_csv)
        ),
        stringsAsFactors = FALSE
      )
    rm(aggregate_frame)
  } else {
    # If a "nearest_ref" variable is missing,
    # the GTF is the original reference without de-novo transcript assembly.
    #
    # Example: gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2";
    #          exon_number "1"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
    #          havana_gene "OTTHUMG00000000961"; havana_gene_version "2";
    #          transcript_name "DDX11L1-002"; transcript_source "havana"; transcript_biotype "processed_transcript";
    #          havana_transcript "OTTHUMT00000362751"; havana_transcript_version "1";
    #          exon_id "ENSE00002234944"; exon_version "1";
    #          tag "basic"; transcript_support_level "1";
    #
    # 3. Joining reference and assembly frames is not required.
    # 4. Create a new, normalised Ensembl annotation data frame.
    ensembl_annotation <- data.frame(
      "gene_id" = reference_frame$ensembl_gene_id,
      "transcript_ids" = reference_frame$ensembl_transcript_id,
      "ensembl_gene_ids" = reference_frame$ensembl_gene_id,
      "ensembl_transcript_ids" = reference_frame$ensembl_transcript_id,
      stringsAsFactors = FALSE
    )
    rm(reference_frame)
  }
  rm(assembly_granges)

  # 5. Create a gene annotation frame, ready for enriching
  # cummeRbund gene information, by merging the "gene_id", "gene_short_name"
  # and "locus" variables of the gene annotation frame with the "transcript_ids",
  # "ensembl_gene_ids" and "ensembl_transcript_ids" of the
  # Ensembl annotation frame.
  message("Creating gene annotation information")
  cufflinks_annotation <-
    cummeRbund::annotation(object = cummeRbund::genes(object = cuff_set))
  gene_annotation_frame <- merge(
    # Merge with a simplified cummeRbund gene annotation data frame, since
    # variables "class_code", "nearest_ref_id", "length" and "coverage" seem
    # empty by design. Remove hidden exon information by finding unique
    # rows only.
    x = unique.data.frame(x = cufflinks_annotation[, c("gene_id", "gene_short_name", "locus"), drop = FALSE]),
    # Merge with the Ensembl annotation data frame that correlates gene (XLOC)
    # identifiers with comma-separated lists of Ensembl Gene and Transcript
    # identifiers.
    y = ensembl_annotation,
    # Merge on gene (XLOC) identifiers, but also keep all rows of the
    # cufflinks_annotation frame, or else, locus information would be lost.
    by = "gene_id",
    all.x = TRUE
  )
  rm(ensembl_annotation, cufflinks_annotation)

  # 6. Create an isoform annotation data frame, ready for enriching
  # cummeRbund isoform information.
  # Simplify the cummeRbund isoform annotation data frame,
  # since variables "CDS_id" and "coverage" seem empty by design.
  # Remove hidden exon information by finding unique rows only.

  message("Creating isoform annotation information")
  cufflinks_annotation <-
    cummeRbund::annotation(object = cummeRbund::isoforms(object = cuff_set))
  isoform_annotation_frame <-
    unique.data.frame(x = cufflinks_annotation[, c(
      "isoform_id",
      "gene_id",
      "gene_short_name",
      "TSS_group_id",
      "class_code",
      "nearest_ref_id",
      "locus",
      "length"
    ), drop = FALSE])
  rm(cufflinks_annotation)

  write.table(x = gene_annotation_frame,
              frame_path_genes,
              sep = "\t",
              row.names = FALSE)
  write.table(
    x = isoform_annotation_frame,
    frame_path_isoforms,
    sep = "\t",
    row.names = FALSE
  )
}
rm(frame_path_genes, frame_path_isoforms)

# Update the cummeRbund SQLite database -----------------------------------


# Push the aggregated Ensembl gene annotation back into the SQLite database.
if ("ensembl_gene_ids" %in% names(x = annotation(object = cummeRbund::genes(object = cuff_set)))) {
  message("Skipping Ensembl annotation for the SQLite database")
} else {
  message("Creating Ensembl annotation for the SQLite database")
  addFeatures(
    object = cuff_set,
    features = data.frame(
      gene_id = gene_annotation_frame$gene_id,
      ensembl_gene_ids = gene_annotation_frame$ensembl_gene_ids,
      ensembl_transcript_ids = gene_annotation_frame$ensembl_transcript_ids,
      stringsAsFactors = FALSE
    ),
    level = "genes"
  )
}

# Differential Genes per sample pair --------------------------------------


# Create an annotated differential genes data frame for each sample pair.
for (i in seq_along(along.with = sample_pairs[1L,])) {
  frame_path <-
    file.path(
      output_directory,
      paste(prefix, sample_pairs[1L, i], sample_pairs[2L, i], "genes_diff.tsv", sep = "_")
    )
  if (file.exists(frame_path) && file.info(frame_path)$size > 0L) {
    message(
      "Skipping a differential data frame on Genes for ",
      sample_pairs[1L, i],
      " versus ",
      sample_pairs[2L, i]
    )
  } else {
    message(
      "Creating a differential data frame on Genes for ",
      sample_pairs[1L, i],
      " versus ",
      sample_pairs[2L, i]
    )
    # The diffData function allows automatic merging with feature annotation,
    # but that includes some empty columns. For cleaner result tables, merge
    # with the smaller gene_annotation_frame established above.
    diff_data_genes <-
      cummeRbund::diffData(
        object = cummeRbund::genes(object = cuff_set),
        x = sample_pairs[1L, i],
        y = sample_pairs[2L, i],
        features = FALSE
      )
    # Remove the second column, which is duplicated as a consequence of a
    # SQL table join between the "genes" and "geneExpDiffData" tables.
    diff_data_genes <- diff_data_genes[, -c(2L)]
    # Calculate ranks for the effect size (log2_fold_change), absolute level
    # and statistical significance (q_value).
    diff_data_genes$rank_log2_fold_change <-
      rank(
        x = -abs(x = diff_data_genes$log2_fold_change),
        ties.method = c("min")
      )
    diff_data_genes$rank_value <-
      rank(
        x = -abs(x = diff_data_genes$value_2 - diff_data_genes$value_1),
        ties.method = c("min")
      )
    diff_data_genes$rank_q_value <-
      rank(x = diff_data_genes$q_value, ties.method = c("min"))
    # Calculate the rank sum for the three ranks.
    # diff_data_genes$rank_sum <-
    #   diff_data_genes$rank_log2_fold_change + diff_data_genes$rank_value +
    #   diff_data_genes$rank_q_value
    # Calculate the maximum rank of value and q-value, as
    # rank_log2_fold_change is dominated by +/- infinity resulting from genes
    # measured only once.
    diff_data_genes$max_rank <-
      pmax(diff_data_genes$rank_value,
           diff_data_genes$rank_q_value)
    gene_merge <-
      merge(
        x = gene_annotation_frame,
        y = diff_data_genes,
        by = "gene_id",
        all = TRUE,
        sort = TRUE
      )
    write.table(
      x = gene_merge,
      file = frame_path,
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )
    rm(gene_merge, diff_data_genes)
  }
  rm(frame_path)
}

# Differential Isoforms per sample pair -----------------------------------


# Create an annotated differential isoforms data frame for each sample pair.
for (i in seq_along(along.with = sample_pairs[1L,])) {
  frame_path <-
    file.path(
      output_directory,
      paste(
        prefix,
        sample_pairs[1L, i],
        sample_pairs[2L, i],
        "isoforms_diff.tsv",
        sep = "_"
      )
    )
  if (file.exists(frame_path) && file.info(frame_path)$size > 0L) {
    message(
      "Skipping a differential data frame on Isoforms for ",
      sample_pairs[1L, i],
      " versus ",
      sample_pairs[2L, i]
    )
  } else {
    message(
      "Creating a differential data frame on Isoforms for ",
      sample_pairs[1L, i],
      " versus ",
      sample_pairs[2L, i]
    )
    # The diffData function allows automatic merging with feature annotation,
    # but that includes some empty columns. For cleaner result tables, merge
    # with the smaller isoform_annotation_frame established above.
    diff_data_isoform <-
      cummeRbund::diffData(
        object = cummeRbund::isoforms(object = cuff_set),
        x = sample_pairs[1L, i],
        y = sample_pairs[2L, i],
        features = FALSE
      )
    # Remove the second column, which is duplicated as a consequence of a
    # SQL table join between the "isoforms" and "isoformsExpDiffData" tables.
    diff_data_isoform <- diff_data_isoform[, -c(2L)]
    # Calculate ranks for the effect size (log2_fold_change), absolute level
    # and statistical significance (q_value).
    diff_data_isoform$rank_log2_fold_change <-
      rank(
        x = -abs(x = diff_data_isoform$log2_fold_change),
        ties.method = c("min")
      )
    diff_data_isoform$rank_value <-
      rank(
        x = -abs(x = diff_data_isoform$value_2 - diff_data_isoform$value_1),
        ties.method = c("min")
      )
    diff_data_isoform$rank_q_value <-
      rank(x = diff_data_isoform$q_value, ties.method = c("min"))
    # Calculate the rank sum for the three ranks.
    # diff_data_isoform$rank_sum <-
    #   diff_data_isoform$rank_log2_fold_change + diff_data_isoform$rank_value +
    #   diff_data_isoform$rank_q_value
    # Calculate the maximum rank of value and q-value, as rank_log2_fold_change
    # is dominated by +/- infinity resulting from genes measured only once.
    diff_data_isoform$max_rank <-
      pmax(diff_data_isoform$rank_value,
           diff_data_isoform$rank_q_value)
    isoform_merge <-
      merge(
        x = isoform_annotation_frame,
        y = diff_data_isoform,
        by = "isoform_id",
        all = TRUE,
        sort = TRUE
      )
    write.table(
      x = isoform_merge,
      file = frame_path,
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )
    rm(diff_data_isoform, isoform_merge)
  }
  rm(frame_path)
}

# Aggregate status information --------------------------------------------


# Aggregate the "status" variable of differential data frames to inform about
# the number of NOTEST, HIDATA, LOWDATA and FAIL states.
# It is silly to redo this outside of the loops splitting the results, but if
# files above were partially written, the status frame would not be complete.

# Status frame to aggregate the test status columns in pairwise comparisons of genes.
frame_path <-
  file.path(output_directory,
            paste0(prefix, "_genes_status.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0L) {
  message("Skipping aggregated status information on Genes")
} else {
  message("Creating aggregated status information on Genes")
  status_frame <-
    data.frame(
      # The sample_pairs character matrix contains pairs in columns,
      # so that the apply function needs to process by columns (i.e. MARGIN = 2L).
      "PAIR" = apply(
        X = sample_pairs,
        MARGIN = 2L,
        FUN = function(x) {
          paste(x, collapse = "__")
        }
      ),
      "OK" = integer(length = length(x = sample_pairs[1L,])),
      "NOTEST" = integer(length = length(x = sample_pairs[1L,])),
      "HIDATA" = integer(length = length(x = sample_pairs[1L,])),
      "LOWDATA" = integer(length = length(x = sample_pairs[1L,])),
      "FAIL" = integer(length = length(x = sample_pairs[1L,])),
      "SUM" = integer(length = length(x = sample_pairs[1L,])),
      row.names = apply(
        X = sample_pairs,
        MARGIN = 2L,
        FUN = function(x) {
          paste(x, collapse = "__")
        }
      )
    )

  for (i in seq_along(along.with = sample_pairs[1L,])) {
    diff_data_genes <-
      cummeRbund::diffData(
        object = cummeRbund::genes(object = cuff_set),
        x = sample_pairs[1L, i],
        y = sample_pairs[2L, i],
        features = FALSE
      )
    # Aggregate the test Status column.
    status_integer <- table(diff_data_genes$status)
    for (status in names(x = status_frame)) {
      if (!is.na(x = status_integer[status])) {
        status_frame[i, status] <- status_integer[status]
      }
    }
    rm(status, status_integer, diff_data_genes)
    status_frame[i, "SUM"] <- sum(status_frame[i, 2L:6L])
  }
  write.table(
    x = status_frame,
    file = frame_path,
    row.names = FALSE,
    sep = "\t"
  )
  rm(status_frame)
}
rm(frame_path)

# Status frame to aggregate the test status columns in pairwise comparisons of isoforms.
frame_path <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_status.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0L) {
  message("Skipping aggregated status information on Isoforms")
} else {
  message("Creating aggregated status information on Isoforms")
  status_frame <-
    data.frame(
      # The sample_pairs character matrix contains pairs in columns,
      # so that the apply function needs to process by columns (i.e. MARGIN = 2L).
      "PAIR" = apply(
        X = sample_pairs,
        MARGIN = 2L,
        FUN = function(x) {
          paste(x, collapse = "__")
        }
      ),
      "OK" = integer(length = length(x = sample_pairs[1L,])),
      "NOTEST" = integer(length = length(x = sample_pairs[1L,])),
      "HIDATA" = integer(length = length(x = sample_pairs[1L,])),
      "LOWDATA" = integer(length = length(x = sample_pairs[1L,])),
      "FAIL" = integer(length = length(x = sample_pairs[1L,])),
      "SUM" = integer(length = length(x = sample_pairs[1L,])),
      row.names = apply(
        X = sample_pairs,
        MARGIN = 2L,
        FUN = function(x) {
          paste(x, collapse = "__")
        }
      )
    )

  for (i in seq_along(along.with = sample_pairs[1L,])) {
    diff_data_isoforms <-
      cummeRbund::diffData(
        object = cummeRbund::isoforms(object = cuff_set),
        x = sample_pairs[1L, i],
        y = sample_pairs[2L, i],
        features = FALSE
      )
    # Aggregate the test Status column.
    status_integer <- table(diff_data_isoforms$status)
    for (status in names(x = status_frame)) {
      if (!is.na(x = status_integer[status])) {
        status_frame[i, status] <- status_integer[status]
      }
    }
    rm(status, status_integer, diff_data_isoforms)
    status_frame[i, "SUM"] <- sum(status_frame[i, 2L:6L])
  }
  write.table(
    x = status_frame,
    file = frame_path,
    row.names = FALSE,
    sep = "\t"
  )
  rm(status_frame)
}
rm(frame_path)

# Matrix of FPKM values per replicate on Genes ----------------------------


frame_path <-
  file.path(output_directory,
            paste0(prefix, "_genes_fpkm_replicates.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0L) {
  message("Skipping a matrix of FPKM values per replicates on genes")
} else {
  message("Creating a matrix of FPKM values per replicates on genes")
  gene_rep_fpkm_matrix <-
    cummeRbund::repFpkmMatrix(object = cummeRbund::genes(object = cuff_set))
  gene_merge <-
    merge(
      x = gene_annotation_frame,
      y = gene_rep_fpkm_matrix,
      by.x = "gene_id",
      by.y = "row.names",
      all = TRUE,
      sort = TRUE
    )
  write.table(
    x = gene_merge,
    file = frame_path,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  rm(gene_merge, gene_rep_fpkm_matrix)
}
rm(frame_path)

# Matrix of count values per replicate on genes.

frame_path <-
  file.path(output_directory,
            paste0(prefix, "_genes_counts_replicates.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0L) {
  message("Skipping a matrix of count values per replicates on genes")
} else {
  message("Creating a matrix of count values per replicates on genes")
  gene_rep_count_matrix <-
    cummeRbund::repCountMatrix(object = cummeRbund::genes(object = cuff_set))
  gene_merge <-
    merge(
      x = gene_annotation_frame,
      y = gene_rep_count_matrix,
      by.x = "gene_id",
      by.y = "row.names",
      all = TRUE,
      sort = TRUE
    )
  write.table(
    x = gene_merge,
    file = frame_path,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  rm(gene_merge, gene_rep_count_matrix)
}
rm(frame_path)

# Matrix of FPKM values per replicate on Isoforms -------------------------


frame_path <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_fpkm_replicates.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0L) {
  message("Skipping a matrix of FPKM per replicates on isoforms")
} else {
  message("Creating a matrix of FPKM per replicates on isoforms")
  isoform_rep_fpkm_matrix <-
    cummeRbund::repFpkmMatrix(object = cummeRbund::isoforms(object = cuff_set))
  isoform_merge <-
    merge(
      x = isoform_annotation_frame,
      y = isoform_rep_fpkm_matrix,
      by.x = "isoform_id",
      by.y = "row.names",
      all = TRUE,
      sort = TRUE
    )
  write.table(
    x = isoform_merge,
    file = frame_path,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  rm(isoform_merge, isoform_rep_fpkm_matrix)
}
rm(frame_path)

# Matrix of count values per replicate on Isoforms ------------------------


frame_path <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_counts_replicates.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0L) {
  message("Skipping a matrix of count values per replicates on isoforms")
} else {
  message("Creating a matrix of count values per replicates on isoforms")
  isoform_rep_count_matrix <-
    cummeRbund::repCountMatrix(object = cummeRbund::isoforms(object = cuff_set))
  isoform_merge <-
    merge(
      x = isoform_annotation_frame,
      y = isoform_rep_count_matrix,
      by.x = "isoform_id",
      by.y = "row.names",
      all = TRUE,
      sort = TRUE
    )
  write.table(
    x = isoform_merge,
    file = frame_path,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  rm(isoform_merge, isoform_rep_count_matrix)
}
rm(frame_path)

# Significance Matrix on Genes --------------------------------------------


plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_genes_significance_matrix.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_genes_significance_matrix.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a significance matrix plot on Genes")
} else {
  message("Creating a significance matrix plot on Genes")
  ggplot_object <- sigMatrix(object = cuff_set, level = "genes")
  ggplot2::ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggplot2::ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

# Significance Matrix on Isoforms -----------------------------------------


plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_significance_matrix.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_significance_matrix.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0L) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0L)) {
  message("Skipping a significance matrix plot on Isoforms")
} else {
  message("Creating a significance matrix plot on Isoforms")
  ggplot_object <- sigMatrix(object = cuff_set, level = "isoforms")
  ggplot2::ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggplot2::ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

rm(gene_annotation_frame, isoform_annotation_frame)

# Get a CuffFeatureSet or CuffGeneSet of significant genes.

# TODO: Comment-out for the moment, as the data is currently not used.
# significant_gene_ids <-
#   getSig(object = cuff_set, level = "genes")
# significant_genes <-
#   cummeRbund::getGenes(object = cuff_set, geneIdList = significant_gene_ids)

# The csHeatmap plot does not seem to be a sensible option for larger sets
# of significant genes.
# for (i in seq_along(along.with = sample_pairs[1L,])) {

# significant_genes_diff <-
# }

# Starting symbolic linking -----------------------------------------------


# Finally, create comparison-specific relative symbolic links for cuffdiff
# results in the rnaseq_process_cuffdiff_* directory to avoid identical file
# names between comparisons.

message("Starting symbolic linking to cuffdiff results")

link_path <-
  file.path(output_directory, paste0(prefix, "_cds_exp_diff.tsv"))
if (!file.exists(link_path)) {
  if (!file.symlink(from = file.path("..", cuffdiff_directory, "cds_exp.diff"),
                    to = link_path)) {
    warning("Encountered an error linking the cds_exp.diff file.")
  }
}
rm(link_path)

link_path <-
  file.path(output_directory, paste0(prefix, "_genes_exp_diff.tsv"))
if (!file.exists(link_path)) {
  if (!file.symlink(from = file.path("..", cuffdiff_directory, "gene_exp.diff"),
                    to = link_path)) {
    warning("Encountered an error linking the gene_exp.diff file.")
  }
}
rm(link_path)

link_path <-
  file.path(output_directory, paste0(prefix, "_isoforms_exp_diff.tsv"))
if (!file.exists(link_path)) {
  if (!file.symlink(from = file.path("..", cuffdiff_directory, "isoform_exp.diff"),
                    to = link_path)) {
    warning("Encountered an error linking the isoform_exp.diff file.")
  }
}
rm(link_path)

link_path <-
  file.path(output_directory, paste0(prefix, "_promoters_diff.tsv"))
if (!file.exists(link_path)) {
  if (!file.symlink(from = file.path("..", cuffdiff_directory, "promoters.diff"),
                    to = link_path)) {
    warning("Encountered an error linking the promoters.diff file.")
  }
}
rm(link_path)

link_path <-
  file.path(output_directory, paste0(prefix, "_splicing_diff.tsv"))
if (!file.exists(link_path)) {
  if (!file.symlink(from = file.path("..", cuffdiff_directory, "splicing.diff"),
                    to = link_path)) {
    warning("Encountered an error linking the splicing.diff file.")
  }
}
rm(link_path)

link_path <-
  file.path(output_directory, paste0(prefix, "_tss_group_exp_diff.tsv"))
if (!file.exists(link_path)) {
  if (!file.symlink(from = file.path("..", cuffdiff_directory, "tss_group_exp.diff"),
                    to = link_path)) {
    warning("Encountered an error linking the tss_group_exp.diff file.")
  }
}
rm(link_path)

# Finishing symbolic linking ----------------------------------------------


message("Finishing symbolic linking to cuffdiff results")

rm(
  cuff_set,
  output_directory,
  cuffdiff_directory,
  sample_pairs,
  sample_number,
  replicate_number,
  have_replicates,
  prefix,
  i,
  argument_list,
  character_to_csv
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

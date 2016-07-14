#! /usr/bin/env Rscript
#
# BSF R script to post-processes Cuffdiff output utilising the cummeRbund
# package. This script creates a cummeRbund SQLite database, writes a set of QC
# plots in PDF and PNG format, splits the large and unwieldy differential data
# tables into pairwise comparisons and creates symbolic links to the original
# Cuffdiff differential tables.-
#
#
# Copyright 2013 -2015 Michael K. Schuster
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

suppressPackageStartupMessages(expr = library(package = "cummeRbund"))
suppressPackageStartupMessages(expr = library(package = "optparse"))
suppressPackageStartupMessages(expr = library(package = "rtracklayer"))

# Get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults.

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
      help = "Comparison name"
    ),
    make_option(
      opt_str = c("--gtf-assembly"),
      default = NULL,
      dest = "gtf_assembly",
      help = "GTF file specifying an assembled and merged transcriptome"
    ),
    make_option(
      opt_str = c("--gtf-reference"),
      default = NULL,
      dest = "gtf_reference",
      help = "GTF file specifying a reference transcriptome"
    ),
    make_option(
      opt_str = c("--genome-version"),
      default = NULL,
      dest = "genome_version",
      help = "Genome version"
    ),
    make_option(
      opt_str = c("--plot-width"),
      default = 7,
      dest = "plot_width",
      help = "Plot width in inches",
      type = "numeric"
    ),
    make_option(
      opt_str = c("--plot-height"),
      default = 7,
      dest = "plot_height",
      help = "Plot height in inches",
      type = "numeric"
    )
  )
))

# Define CuffDiff and output directory names relative to the working directory.

cuffdiff_directory <-
  paste("rnaseq", "cuffdiff", argument_list$comparison_name, sep = "_")
output_directory <-
  paste("rnaseq",
        "process",
        "cuffdiff",
        argument_list$comparison_name,
        sep = "_")

# To avoid name clashes when downloading files, use the output directory name
# also as a prefix for all files therein.

prefix <- output_directory

# Read and index Cuffdiff output and create a CuffSet object. The CuffSet object
# has slots genes, isoforms, TSS and CDS that are each instances of the CuffData
# class. By default, CuffData accessor methods applied to a CuffSet class will
# operate on the ’genes’ slot.

message("Create or load the cummeRbund database")
cuff_set <-
  readCufflinks(
    dir = cuffdiff_directory,
    gtfFile = argument_list$gtf_reference,
    genome = argument_list$genome_version
  )

# Create a new sub-directory for plots if it does not exist.

if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

# Process Cuffdiff run information.

frame_path <-
  file.path(output_directory, paste0(prefix, "_run_information.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
  message("Skipping run information table")
} else {
  message("Creating run information table")
  write.table(
    x = runInfo(object = cuff_set),
    file = frame_path,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
}
rm(frame_path)

# Process Cuffdiff sample information.

sample_frame <- samples(object = cuff_set)
# Get the number of samples.
sample_number <- nrow(x = sample_frame)
# Write the sample_frame table.
frame_path <-
  file.path(output_directory, paste0(prefix, "_samples.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
  message("Skipping sample table")
} else {
  message("Creating sample table")
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
sample_pairs <- combn(x = sample_frame$sample_name, m = 2)
# Write a table of sample pair information by transposing the sample pairs
# array. This table allows the Python web code to link in pairwise plots.
frame_path <-
  file.path(output_directory, paste0(prefix, "_sample_pairs.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
  message("Skipping sample pairs table")
} else {
  message("Create sample pairs table")
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

# Process Cuffdiff replicate information.

replicate_frame <- replicates(object = cuff_set)
# Get the number of replicates.
replicate_number <- nrow(x = replicate_frame)
# Some plots require replicates. Check whether at least one row has a replicate
# column value greater than 0.
have_replicates <-
  (nrow(x = replicate_frame[replicate_frame$replicate > 0,]) > 0)
# Write the replicate_frame.
frame_path <-
  file.path(output_directory, paste0(prefix, "_replicates.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
  message("Skipping replicate table")
} else {
  message("Creating replicate table")
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
# Plot the log10(internal_scale) of the replicate_frame to visualise outliers.
plot_path_pdf <-
  file.path(output_directory, paste0(prefix, "_replicate_scale.pdf"))
plot_path_png <-
  file.path(output_directory, paste0(prefix, "_replicate_scale.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a library scale plot on replicates")
} else {
  message("Creating a library scale plot on replicates")
  ggplot_object <- ggplot(data = replicate_frame)
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = rep_name, y = log10(internal_scale)))
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  plot_width <-
    argument_list$plot_width + (ceiling(x = replicate_number / 24) - 1) * argument_list$plot_width * 0.25
  ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object, plot_width)
}
rm(plot_path_pdf, plot_path_png, replicate_frame)

# Create a set of QC plots.

message("Started QC plotting")

# Create a Dispersion Plot on Genes.

plot_path_pdf <-
  file.path(output_directory, paste0(prefix, "_genes_dispersion.pdf"))
plot_path_png <-
  file.path(output_directory, paste0(prefix, "_genes_dispersion.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a Dispersion Plot on Genes")
} else {
  message("Creating a Dispersion Plot on Genes")
  ggplot_object <- dispersionPlot(object = genes(object = cuff_set))
  ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

# Create a Dispersion Plot on Isoforms.

plot_path_pdf <-
  file.path(output_directory, paste0(prefix, "_isoforms_dispersion.pdf"))
plot_path_png <-
  file.path(output_directory, paste0(prefix, "_isoforms_dispersion.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a Dispersion Plot on Isoforms")
} else {
  message("Creating Dispersion Plot on Isoforms")
  tryCatch({
    ggplot_object <-
      dispersionPlot(object = isoforms(object = cuff_set))
    ggsave(
      filename = plot_path_pdf,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ggsave(
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
  })
}
rm(plot_path_pdf, plot_path_png)

# Create a Squared Coefficient of Variation (SCV) Plot on Genes.
# The plot requires replicates.

plot_path_pdf <-
  file.path(output_directory, paste0(prefix, "_genes_scv.pdf"))
plot_path_png <-
  file.path(output_directory, paste0(prefix, "_genes_scv.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a SCV Plot on Genes")
} else {
  if (have_replicates) {
    message("Creating a SCV Plot on Genes")
    ggplot_object <- fpkmSCVPlot(object = genes(object = cuff_set))
    ggsave(
      filename = plot_path_pdf,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ggsave(
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

# Create a Squared Coefficient of Variation (SCV) Plot on Isoforms.
# The plot requires replicates.

plot_path_pdf <-
  file.path(output_directory, paste0(prefix, "_isoforms_scv.pdf"))
plot_path_png <-
  file.path(output_directory, paste0(prefix, "_isoforms_scv.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a SCV Plot on Isoforms")
} else {
  if (have_replicates) {
    message("Creating a SCV Plot on Isoforms")
    ggplot_object <-
      fpkmSCVPlot(object = isoforms(object = cuff_set))
    ggsave(
      filename = plot_path_pdf,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ggsave(
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

# Create a Density Plot on Genes with and without replicates.

plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_genes_density_wo_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_genes_density_wo_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a Density Plot on Genes without replicates")
} else {
  message("Creating a Density Plot on Genes without replicates")
  ggplot_object <-
    csDensity(object = genes(object = cuff_set), replicates = FALSE)
  ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_genes_density_w_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_genes_density_w_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a Density Plot on Genes with replicates")
} else {
  message("Creating a Density Plot on Genes with replicates")
  ggplot_object <-
    csDensity(object = genes(object = cuff_set), replicates = TRUE)
  ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

# Create a Density Plot on Isoforms with and without replicates.

plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_density_wo_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_density_wo_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a Density Plot on Isoforms without replicates")
} else {
  message("Creating a Density Plot on Isoforms without replicates")
  ggplot_object <-
    csDensity(object = isoforms(object = cuff_set), replicates = FALSE)
  ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_density_w_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_density_w_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a Density Plot on Isoforms with replicates")
} else {
  message("Creating a Density Plot on Isoforms with replicates")
  ggplot_object <-
    csDensity(object = isoforms(object = cuff_set), replicates = TRUE)
  ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

# Create a Box Plot on Genes with and without replicates.

plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_genes_box_w_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_genes_box_w_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a Box Plot on Genes with replicates")
} else {
  message("Creating a Box Plot on Genes with replicates")
  # By default, the csBoxplot function adds a pseudocount of 1e-04 for log10
  # transformed fpkms. In case of a large number of missing values in shallowly
  # sequenced samples, this affects the plot. Hence reproduce the plot here.
  # ggplot_object <- csBoxplot(object = genes(object = cuff_set), replicates = TRUE)
  rep_fpkm_genes <- repFpkm(object = genes(object = cuff_set))
  # Rename the "rep_name" column into "condition".
  colnames(x = rep_fpkm_genes)[colnames(x = rep_fpkm_genes) == "rep_name"] <-
    "condition"
  ggplot_object <- ggplot(data = rep_fpkm_genes)
  ggplot_object <-
    ggplot_object + geom_boxplot(
      mapping = aes(
        x = condition,
        y = log10(fpkm),
        fill = condition
      ),
      size = 0.3,
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  ggplot_object <-
    ggplot_object + theme(legend.text = element_text(size = rel(x = 0.8)))  # Reduce the legend text
  ggplot_object <-
    ggplot_object + scale_fill_hue(l = 50, h.start = 200)
  # Arrange a maximum of 24 replicates in each guide column.
  # ggplot_object <-
  #   ggplot_object + guides(fill = guide_legend(nrow = 24))
  # Use the base plot_witdh and add a quarter of the width for each additional guide legend column.
  plot_width <-
    argument_list$plot_width + (ceiling(x = replicate_number / 24) - 1) * argument_list$plot_width * 0.25
  ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object, plot_width, rep_fpkm_genes)
}
rm(plot_path_pdf, plot_path_png)

plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_genes_box_wo_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_genes_box_wo_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a Box Plot on Genes without replicates")
} else {
  message("Creating a Box Plot on Genes without replicates")
  # ggplot_object <- csBoxplot(object = genes(object = cuff_set), replicates = FALSE)
  fpkm_genes <- fpkm(object = genes(object = cuff_set))
  # Rename the "sample_name" column into "condition".
  colnames(fpkm_genes)[colnames(fpkm_genes) == "sample_name"] <-
    "condition"
  ggplot_object <- ggplot(data = fpkm_genes)
  ggplot_object <-
    ggplot_object + geom_boxplot(
      mapping = aes(
        x = condition,
        y = log10(fpkm),
        fill = condition
      ),
      size = 0.3,
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  ggplot_object <-
    ggplot_object + theme(legend.text = element_text(size = rel(x = 0.8)))  # Reduce the legend text
  ggplot_object <-
    ggplot_object + scale_fill_hue(l = 50, h.start = 200)
  # Arrange a maximum of 24 samples in each guide column.
  # ggplot_object <-
  #   ggplot_object + guides(fill = guide_legend(nrow = 24))
  # Use the base plot_witdh and add a quarter of the width for each additional guide legend column.
  plot_width <-
    argument_list$plot_width + (ceiling(x = sample_number / 24) - 1) * argument_list$plot_width * 0.25
  ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object, plot_width, fpkm_genes)
}
rm(plot_path_pdf, plot_path_png)

# Create a Box Plot on Isoforms with and without replicates.

plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_box_w_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_box_w_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a Box Plot on Isoforms with replicates")
} else {
  message("Creating a Box Plot on Isoforms with replicates")
  # By default, the csBoxplot function adds a pseudocount of 1e-04 for log10
  # transformed fpkms. In case of a large number of missing values in shallowly
  # sequenced samples, this affects the plot. Hence reproduce the plot here.
  # ggplot_object <- csBoxplot(object = isoforms(object = cuff_set), replicates = TRUE)
  rep_fpkm_isoforms <- repFpkm(object = isoforms(object = cuff_set))
  # Rename the "rep_name" column into "condition".
  colnames(rep_fpkm_isoforms)[colnames(rep_fpkm_isoforms) == "rep_name"] <-
    "condition"
  ggplot_object <- ggplot(data = rep_fpkm_isoforms)
  ggplot_object <-
    ggplot_object + geom_boxplot(
      mapping = aes(
        x = condition,
        y = log10(fpkm),
        fill = condition
      ),
      size = 0.3,
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  ggplot_object <-
    ggplot_object + theme(legend.text = element_text(size = rel(x = 0.8)))  # Reduce the legend text
  ggplot_object <-
    ggplot_object + scale_fill_hue(l = 50, h.start = 200)
  # Arrange a maximum of 24 replicates in each guide column.
  # ggplot_object <-
  #   ggplot_object + guides(fill = guide_legend(nrow = 24))
  # Use the base plot_witdh and add a quarter of the width for each additional guide legend column.
  plot_width <-
    argument_list$plot_width + (ceiling(x = replicate_number / 24) - 1) * argument_list$plot_width * 0.25
  ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object, plot_width, rep_fpkm_isoforms)
}
rm(plot_path_pdf, plot_path_png)

plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_box_wo_replicates.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_box_wo_replicates.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a Box Plot on Isoforms without replicates")
} else {
  message("Creating a Box Plot on Isoforms without replicates")
  # ggplot_object <- csBoxplot(object = isoforms(object = cuff_set), replicates = FALSE)
  fpkm_isoforms <- fpkm(object = isoforms(object = cuff_set))
  # Rename the "sample_name" column into "condition".
  colnames(fpkm_isoforms)[colnames(fpkm_isoforms) == "sample_name"] <-
    "condition"
  ggplot_object <- ggplot(data = fpkm_isoforms)
  ggplot_object <-
    ggplot_object + geom_boxplot(
      mapping = aes(
        x = condition,
        y = log10(fpkm),
        fill = condition
      ),
      size = 0.3,
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  ggplot_object <-
    ggplot_object + theme(legend.text = element_text(size = rel(x = 0.8)))  # Reduce the legend text
  ggplot_object <-
    ggplot_object + scale_fill_hue(l = 50, h.start = 200)
  # Arrange a maximum of 24 replicates in each guide column.
  # ggplot_object <-
  #   ggplot_object + guides(fill = guide_legend(nrow = 24))
  # Use the base plot_witdh and add a quarter of the with for each additional guide legend column.
  plot_width <-
    argument_list$plot_width + (ceiling(x = sample_number / 24) - 1) * argument_list$plot_width * 0.25
  ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object, plot_width, fpkm_isoforms)
}
rm(plot_path_pdf, plot_path_png)

# Create a Scatter Matrix Plot on Genes and Isoforms for less than or equal to
# 20 samples.

if (sample_number <= 20) {
  plot_path_pdf <-
    file.path(output_directory,
              paste0(prefix, "_genes_scatter_matrix.pdf"))
  plot_path_png <-
    file.path(output_directory,
              paste0(prefix, "_genes_scatter_matrix.png"))
  if (file.exists(plot_path_pdf) &&
      (file.info(plot_path_pdf)$size > 0) &&
      file.exists(plot_path_png) &&
      (file.info(plot_path_png)$size > 0)) {
    message("Skipping a Scatter Matrix Plot on Genes")
  } else {
    message("Creating a Scatter Matrix Plot on Genes")
    ggplot_object <-
      csScatterMatrix(object = genes(object = cuff_set))
    ggsave(
      filename = plot_path_pdf,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ggsave(
      filename = plot_path_png,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    rm(ggplot_object)
  }
  rm(plot_path_pdf, plot_path_png)
}

if (sample_number <= 20) {
  plot_path_pdf <-
    file.path(output_directory,
              paste0(prefix, "_isoforms_scatter_matrix.pdf"))
  plot_path_png <-
    file.path(output_directory,
              paste0(prefix, "_isoforms_scatter_matrix.png"))
  if (file.exists(plot_path_pdf) &&
      (file.info(plot_path_pdf)$size > 0) &&
      file.exists(plot_path_png) &&
      (file.info(plot_path_png)$size > 0)) {
    message("Skipping a Scatter Matrix Plot on Isoforms")
  } else {
    message("Creating a Scatter Matrix Plot on Isoforms")
    ggplot_object <-
      csScatterMatrix(object = isoforms(object = cuff_set))
    ggsave(
      filename = plot_path_pdf,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ggsave(
      filename = plot_path_png,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    rm(ggplot_object)
  }
  rm(plot_path_pdf, plot_path_png)
}

# Create a Scatter Plot on Genes for each sample pair.

for (i in 1:length(sample_pairs[1,])) {
  plot_path_pdf <-
    file.path(
      output_directory,
      paste(
        prefix,
        sample_pairs[1, i],
        sample_pairs[2, i],
        "genes_scatter.pdf",
        sep = "_"
      )
    )
  plot_path_png <-
    file.path(
      output_directory,
      paste(
        prefix,
        sample_pairs[1, i],
        sample_pairs[2, i],
        "genes_scatter.png",
        sep = "_"
      )
    )
  if (file.exists(plot_path_pdf) &&
      (file.info(plot_path_pdf)$size > 0) &&
      file.exists(plot_path_png) &&
      (file.info(plot_path_png)$size > 0)) {
    message(
      paste(
        "Skipping a Scatter Plot on Genes for",
        sample_pairs[1, i],
        "versus",
        sample_pairs[2, i]
      )
    )
  } else {
    message(
      paste(
        "Creating a Scatter Plot on Genes for",
        sample_pairs[1, i],
        "versus",
        sample_pairs[2, i]
      )
    )
    ggplot_object <-
      csScatter(
        object = genes(object = cuff_set),
        x = sample_pairs[1, i],
        y = sample_pairs[2, i],
        colorByStatus = TRUE
      )
    ggsave(
      filename = plot_path_pdf,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ggsave(
      filename = plot_path_png,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    rm(ggplot_object)
  }
  rm(plot_path_pdf, plot_path_png)
}

# Create a Dendrogram Plot on Genes for time-series analyses. The csDendro
# function returns a dendrogram object that cannot be saved with the ggsave
# function.

plot_path_pdf <-
  file.path(output_directory, paste0(prefix, "_genes_dendrogram.pdf"))
if (file.exists(plot_path_pdf) &&
    file.info(plot_path_pdf)$size > 0) {
  message("Skipping a Dendrogram Plot on Genes [PDF]")
} else {
  message("Creating a Dendrogram Plot on Genes [PDF]")
  pdf(file = plot_path_pdf)
  csDendro(object = genes(object = cuff_set))
  active_device <- dev.off()
  rm(active_device)
}
rm(plot_path_pdf)

plot_path_png <-
  file.path(output_directory, paste0(prefix, "_genes_dendrogram.png"))
if (file.exists(plot_path_png) &&
    file.info(plot_path_png)$size > 0) {
  message("Skipping a Dendrogram Plot on Genes [PNG]")
} else {
  message("Creating a Dendrogram Plot on Genes [PNG]")
  png(filename = plot_path_png)
  csDendro(object = genes(object = cuff_set))
  active_device <- dev.off()
  rm(active_device)
}
rm(plot_path_png)

# Create a MA Plot on genes for each sample pair based on FPKM values.

for (i in 1:length(sample_pairs[1,])) {
  plot_path_pdf <-
    file.path(
      output_directory,
      paste(prefix, sample_pairs[1, i], sample_pairs[2, i], "maplot.pdf", sep = "_")
    )
  plot_path_png <-
    file.path(
      output_directory,
      paste(prefix, sample_pairs[1, i], sample_pairs[2, i], "maplot.png", sep = "_")
    )
  if (file.exists(plot_path_pdf) &&
      (file.info(plot_path_pdf)$size > 0) &&
      file.exists(plot_path_png) &&
      (file.info(plot_path_png)$size > 0)) {
    message(paste(
      "Skipping a MAplot on Genes for",
      sample_pairs[1, i],
      "versus",
      sample_pairs[2, i]
    ))
  } else {
    message(paste(
      "Creating a MAplot on Genes for",
      sample_pairs[1, i],
      "versus",
      sample_pairs[2, i]
    ))
    ggplot_object <-
      MAplot(object = genes(object = cuff_set),
             x = sample_pairs[1, i],
             y = sample_pairs[2, i])
    ggsave(
      filename = plot_path_pdf,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ggsave(
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

# Create a Volcano Matrix Plot on Genes.

plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_genes_volcano_matrix.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_genes_volcano_matrix.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a Volcano Matrix Plot on Genes")
} else {
  message("Creating a Volcano Matrix Plot on Genes")
  ggplot_object <-
    csVolcanoMatrix(object = genes(object = cuff_set))
  ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

# Create a Volcano Plot on genes for each sample pair.

for (i in 1:length(sample_pairs[1,])) {
  plot_path_pdf <-
    file.path(
      output_directory,
      paste(
        prefix,
        sample_pairs[1, i],
        sample_pairs[2, i],
        "genes_volcano.pdf",
        sep = "_"
      )
    )
  plot_path_png <-
    file.path(
      output_directory,
      paste(
        prefix,
        sample_pairs[1, i],
        sample_pairs[2, i],
        "genes_volcano.png",
        sep = "_"
      )
    )
  if (file.exists(plot_path_pdf) &&
      (file.info(plot_path_pdf)$size > 0) &&
      file.exists(plot_path_png) &&
      (file.info(plot_path_png)$size > 0)) {
    message(
      paste(
        "Skipping a Volcano Plot on Genes for",
        sample_pairs[1, i],
        "versus",
        sample_pairs[2, i]
      )
    )
  } else {
    message(
      paste(
        "Creating a Volcano Plot on Genes for",
        sample_pairs[1, i],
        "versus",
        sample_pairs[2, i]
      )
    )
    # FIXME: The definition of the generic function "csVolcano" does not include
    # the "alpha" and "showSignificant" options, although the function
    # definition contains them. It doesn not seem that the option defaults are used.
    # In R/methods-CuffData.R:
    # .volcano<-function(object,x,y,alpha=0.05,showSignificant=TRUE,features=FALSE,xlimits=c(-20,20),...)
    # setMethod("csVolcano",signature(object="CuffData"), .volcano)
    # In R/AllGenerics.R:
    # setGeneric("csVolcano",function(object, x, y, features=F, ...) standardGeneric("csVolcano"))
    ggplot_object <-
      csVolcano(
        object = genes(object = cuff_set),
        x = sample_pairs[1, i],
        y = sample_pairs[2, i],
        alpha = 0.05,
        showSignificant = TRUE
      )
    ggsave(
      filename = plot_path_pdf,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    ggsave(
      filename = plot_path_png,
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
    rm(ggplot_object)
  }
  rm(plot_path_pdf, plot_path_png)
}

# Create a Multidimensional Scaling (MDS) Plot on genes,
# but only if the CuffData object contains more than two samples.

if (sample_number > 2) {
  plot_path_pdf <-
    file.path(output_directory, paste0(prefix, "_genes_mds.pdf"))
  plot_path_png <-
    file.path(output_directory, paste0(prefix, "_genes_mds.png"))
  if (file.exists(plot_path_pdf) &&
      (file.info(plot_path_pdf)$size > 0) &&
      file.exists(plot_path_png) &&
      (file.info(plot_path_png)$size > 0)) {
    message("Skipping a Multidimensional Scaling Plot on Genes")
  } else {
    # if (have_replicates) {
    message("Creating a Multidimensional Scaling Plot on Genes")
    # Nothing ever is simple. If the set has too many replicates, the standard
    # cummeRbund MDSplot() falls down.
    if (replicate_number <= 24) {
      ggplot_object <-
        MDSplot(object = genes(object = cuff_set),
                replicates = TRUE)
      ggsave(
        filename = plot_path_pdf,
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height
      )
      ggsave(
        filename = plot_path_png,
        plot = ggplot_object,
        width = argument_list$plot_width,
        height = argument_list$plot_height
      )
      rm(ggplot_object)
    } else {
      # The standard MDSplot has too many replicates.
      gene_rep_fit <-
        cmdscale(d = JSdist(mat = makeprobs(a = repFpkmMatrix(
          object = genes(object = cuff_set)
        ))),
        eig = TRUE,
        k = 2)
      gene_rep_res <-
        data.frame(
          names = rownames(gene_rep_fit$points),
          M1 = gene_rep_fit$points[, 1],
          M2 = gene_rep_fit$points[, 2],
          stringsAsFactors = FALSE
        )
      ggplot_object <- ggplot(data = gene_rep_res)
      ggplot_object <-
        ggplot_object + theme_bw()  # Add theme black and white.
      ggplot_object <-
        ggplot_object + geom_point(mapping = aes(x = M1, y = M2, color = names))  # Draw points in any case.
      if (replicate_number <= 24) {
        # Only add text for a sensible number of replicates i.e. less than or
        # equal to 24.
        ggplot_object <-
          ggplot_object + geom_text(mapping = aes(
            x = M1,
            y = M2,
            label = names,
            color = names
          ),
          size = 4)
      }
      # Arrange a maximum of 24 replicates in each guide column.
      ggplot_object <-
        ggplot_object + guides(color = guide_legend(nrow = 24))
      # Use the base plot_witdh and add a quarter of the width for each additional guide legend column.
      plot_width <-
        argument_list$plot_width + (ceiling(x = replicate_number / 24) - 1) * argument_list$plot_width * 0.25
      ggsave(
        filename = plot_path_pdf,
        plot = ggplot_object,
        width = plot_width,
        height = argument_list$plot_height
      )
      ggsave(
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
  message("Skipping Multidimensional Scaling Plot on genes in lack of sufficient samples")
}

# Create a Principal Component Analysis (PCA) Plot on Genes.
# TODO: Add also other principal components or even better use plots of the PCA package?

plot_path_pdf <-
  file.path(output_directory, paste0(prefix, "_genes_pca.pdf"))
plot_path_png <-
  file.path(output_directory, paste0(prefix, "_genes_pca.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a Principal Component Analysis Plot (PCA) on Genes")
} else {
  message("Creating a Principal Component Analysis Plot (PCA) on Genes")
  ggplot_object <-
    PCAplot(
      object = genes(object = cuff_set),
      x = "PC1",
      y = "PC2",
      replicates = TRUE
    )
  ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

# Finished plotting.

message("Finished QC plotting")

# TODO: Process gene and isoform sets and lists. Feature-level information can
# be accessed directly from a CuffData object using the fpkm, repFpkm, count,
# diffData, or annotation methods.

# Unfortunately, there does not seem to be a simple way to correlate XLOC and
# Ensembl gene identifiers within cummeRbund. Read the reference GTF to
# correlate Ensembl gene and transcript identifiers, as well as official gene
# symbols.
message("Read reference transcriptome information")
reference_granges <- import(con = argument_list$gtf_reference)
reference_frame <- unique(
  x = data.frame(
    ensembl_gene_id = reference_granges$gene_id,
    ensembl_transcript_id = reference_granges$transcript_id,
    stringsAsFactors = FALSE
  )
)
rm(reference_granges)

# Read the assmbly GTF, which does no longer have information about Ensembl gene
# identifiers, as loci may have been merged or split by Cuffmerge.
message("Read assembled transcriptome information")
assembly_granges <- import(con = argument_list$gtf_assembly)
assembly_frame <- unique(
  x = data.frame(
    gene_id = assembly_granges$gene_id,
    transcript_id = assembly_granges$transcript_id,
    gene_name = assembly_granges$gene_name,
    ensembl_transcript_id = assembly_granges$nearest_ref,
    stringsAsFactors = FALSE
  )
)
rm(assembly_granges)

# Join the reference and assembly data frames to complement Ensembl gene
# identifiers.
message("Merge reference and assembled annotation information")
ensembl_frame <-
  merge(x = reference_frame, y = assembly_frame, by = "ensembl_transcript_id")
rm(reference_frame, assembly_frame)

# Annotate Ensembl gene and transcript identifiers as independent,
# comma-separated lists. Find unique gene (XLOC) identifiers and correlate them
# with Ensembl gene and transcript identifiers.
unique_gene_identifiers <- unique(x = ensembl_frame$gene_id)
message(paste(
  "Create Ensembl annotation information for",
  length(x = unique_gene_identifiers),
  "genes"
))
ensembl_gene_identifiers <-
  character(length = length(x = unique_gene_identifiers))
ensembl_transcript_identifiers <-
  character(length = length(x = unique_gene_identifiers))
for (i in 1:length(x = unique_gene_identifiers)) {
  gene_frame <-
    ensembl_frame[ensembl_frame$gene_id == unique_gene_identifiers[i],]
  ensembl_gene_identifiers[i] <-
    paste(unique(x = gene_frame$ensembl_gene_id), collapse = ',')
  ensembl_transcript_identifiers[i] <-
    paste(unique(x = gene_frame$ensembl_transcript_id), collapse = ',')
  rm(gene_frame)
}
rm(ensembl_frame)

message("Create gene annotation information")
gene_annotation <- annotation(object = genes(object = cuff_set))
gene_frame <- merge(
  # Create a new, simpler data frame for gene annotation ready for merging with
  # other data frames. The class_code, nearest_ref_id, length and coverage
  # fields seem to be empty by design.
  x = data.frame(
    gene_id = gene_annotation$gene_id,
    gene_short_name = gene_annotation$gene_short_name,
    locus = gene_annotation$locus,
    stringsAsFactors = FALSE
  ),
  # Create a new Ensembl data frame that correlates gene (XLOC) identifiers with
  # comma-separated lists of Ensembl Gene and Transcript Identifiers.
  y = data.frame(
    gene_id = unique_gene_identifiers,
    ensembl_gene_id = ensembl_gene_identifiers,
    ensembl_transcript_id = ensembl_transcript_identifiers,
    stringsAsFactors = FALSE
  ),
  # Merge on gene (XLOC) identifiers.
  by = "gene_id"
)

rm(
  gene_annotation,
  unique_gene_identifiers,
  ensembl_gene_identifiers,
  ensembl_transcript_identifiers
)

# Create a new, simpler data frame for isoform annotation ready for merging with
# other data frames. The CDS_id and coverage fields seem to be empty by design.

message("Create isoform annotation information")
isoform_annotation <-
  annotation(object = isoforms(object = cuff_set))
isoform_frame <- data.frame(
  isoform_id = isoform_annotation$isoform_id,
  gene_id = isoform_annotation$gene_id,
  gene_short_name = isoform_annotation$gene_short_name,
  TSS_group_id = isoform_annotation$TSS_group_id,
  class_code = isoform_annotation$class_code,
  nearest_ref_id = isoform_annotation$nearest_ref_id,
  locus = isoform_annotation$locus,
  length = isoform_annotation$length,
  stringsAsFactors = FALSE
)
rm(isoform_annotation)

# Split the large and unwieldy differential data tables into pairwise comparisons.-

for (i in 1:length(sample_pairs[1,])) {
  frame_path <-
    file.path(
      output_directory,
      paste(prefix, sample_pairs[1, i], sample_pairs[2, i], "genes_diff.tsv", sep = "_")
    )
  if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
    message(
      paste(
        "Skipping a differential data frame on Genes for",
        sample_pairs[1, i],
        "versus",
        sample_pairs[2, i]
      )
    )
  } else {
    message(
      paste(
        "Creating a differential data frame on Genes for",
        sample_pairs[1, i],
        "versus",
        sample_pairs[2, i]
      )
    )
    # The diffData function allows automatic merging with feature annotation,
    # but that includes some empty columns. For cleaner result tables, merge
    # with the smaller gene_frame established above. Unfortunately, the
    # 'gene_id' column appears twice as a consequence of a SQL table join of the
    # 'genes' table with the 'geneExpDiffData' table. Separate data.frame
    # columns of the same name interfere with the merge() function, so that the
    # first column needs removing.
    gene_diff <-
      diffData(
        object = genes(object = cuff_set),
        x = sample_pairs[1, i],
        y = sample_pairs[2, i],
        features = FALSE
      )
    gene_diff <- gene_diff[,-1]
    # Calculate ranks for the effect size (log2_fold_change), absolute level and
    # statistical significance (q_value).
    gene_diff$rank_log2_fold_change <-
      rank(
        x = -abs(x = gene_diff$log2_fold_change),
        ties.method = c("min")
      )
    gene_diff$rank_value <-
      rank(
        x = -abs(x = gene_diff$value_2 - gene_diff$value_1),
        ties.method = c("min")
      )
    gene_diff$rank_q_value <-
      rank(x = gene_diff$q_value, ties.method = c("min"))
    # Calculate the rank sum for the three ranks.
    gene_diff$rank_sum <-
      gene_diff$rank_log2_fold_change + gene_diff$rank_value + gene_diff$rank_q_value
    gene_merge <-
      merge(
        x = gene_frame,
        y = gene_diff,
        by.x = "gene_id",
        by.y = "gene_id",
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
    rm(gene_diff, gene_merge)
  }
  rm(frame_path)
}

for (i in 1:length(sample_pairs[1,])) {
  frame_path <-
    file.path(
      output_directory,
      paste(
        prefix,
        sample_pairs[1, i],
        sample_pairs[2, i],
        "isoforms_diff.tsv",
        sep = "_"
      )
    )
  if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
    message(
      paste(
        "Skipping a differential data frame on Isoforms for",
        sample_pairs[1, i],
        "versus",
        sample_pairs[2, i]
      )
    )
  } else {
    message(
      paste(
        "Creating a differential data frame on Isoforms for",
        sample_pairs[1, i],
        "versus",
        sample_pairs[2, i]
      )
    )
    # The diffData function allows automatic merging with feature annotation,
    # but that includes some empty columns. For cleaner result tables, merge
    # with the smaller isoform_frame estbalished above. Unfortunately, the
    # 'isoform_id' column appears twice as a consequence of a SQL table join of
    # the 'isoforms' table with the 'isoformsExpDiffData' table. Separate
    # data.frame columns of the same name interfere with the merge() function,
    # so that the first column needs removing.
    isoform_diff <-
      diffData(
        object = isoforms(object = cuff_set),
        x = sample_pairs[1, i],
        y = sample_pairs[2, i],
        features = FALSE
      )
    isoform_diff <- isoform_diff[,-1]
    # Calculate ranks for the effect size (log2_fold_change), absolute level and statistical significance (q_value).
    isoform_diff$rank_log2_fold_change <-
      rank(
        x = -abs(x = isoform_diff$log2_fold_change),
        ties.method = c("min")
      )
    isoform_diff$rank_value <-
      rank(
        x = -abs(x = isoform_diff$value_2 - isoform_diff$value_1),
        ties.method = c("min")
      )
    isoform_diff$rank_q_value <-
      rank(x = isoform_diff$q_value, ties.method = c("min"))
    # Calculate the rank sum for the three ranks.
    isoform_diff$rank_sum <-
      isoform_diff$rank_log2_fold_change + isoform_diff$rank_value + isoform_diff$rank_q_value
    isoform_merge <-
      merge(
        x = isoform_frame,
        y = isoform_diff,
        by.x = "isoform_id",
        by.y = "isoform_id",
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
    rm(isoform_diff, isoform_merge)
  }
  rm(frame_path)
}

# Matrix of FPKM values per replicate on genes.

frame_path <-
  file.path(output_directory,
            paste0(prefix, "_genes_fpkm_replicates.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
  message("Skipping a matrix of FPKM values per replicates on genes")
} else {
  message("Creating a matrix of FPKM values per replicates on genes")
  gene_rep_fpkm_matrix <-
    repFpkmMatrix(object = genes(object = cuff_set))
  gene_merge <-
    merge(
      x = gene_frame,
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
if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
  message("Skipping a matrix of count values per replicates on genes")
} else {
  message("Creating a matrix of count values per replicates on genes")
  gene_rep_count_matrix <-
    repCountMatrix(object = genes(object = cuff_set))
  gene_merge <-
    merge(
      x = gene_frame,
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

# Matrix of FPKM values per replicate on isoforms

frame_path <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_fpkm_replicates.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
  message("Skipping a matrix of FPKM per replicates on isoforms")
} else {
  message("Creating a matrix of FPKM per replicates on isoforms")
  isoform_rep_fpkm_matrix <-
    repFpkmMatrix(object = isoforms(object = cuff_set))
  isoform_merge <-
    merge(
      x = isoform_frame,
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

# Matrix of count values per replicate on isoforms.

frame_path <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_counts_replicates.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
  message("Skipping a matrix of count values per replicates on isoforms")
} else {
  message("Creating a matrix of count values per replicates on isoforms")
  isoform_rep_count_matrix <-
    repCountMatrix(object = isoforms(object = cuff_set))
  isoform_merge <-
    merge(
      x = isoform_frame,
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

# TODO: Deal with sets of significantly regulated genes, including a sigMatrix
# ggplot plot.
plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_genes_significance_matrix.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_genes_significance_matrix.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a significance matrix plot on Genes")
} else {
  message("Creating a a significance matrix plot on Genes")
  ggplot_object <- sigMatrix(object = cuff_set, level = "genes")
  ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

plot_path_pdf <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_significance_matrix.pdf"))
plot_path_png <-
  file.path(output_directory,
            paste0(prefix, "_isoforms_significance_matrix.png"))
if (file.exists(plot_path_pdf) &&
    (file.info(plot_path_pdf)$size > 0) &&
    file.exists(plot_path_png) &&
    (file.info(plot_path_png)$size > 0)) {
  message("Skipping a significance matrix plot on Isoforms")
} else {
  message("Creating a a significance matrix plot on Isoforms")
  ggplot_object <- sigMatrix(object = cuff_set, level = "isoforms")
  ggsave(
    filename = plot_path_pdf,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  ggsave(
    filename = plot_path_png,
    plot = ggplot_object,
    width = argument_list$plot_width,
    height = argument_list$plot_height
  )
  rm(ggplot_object)
}
rm(plot_path_pdf, plot_path_png)

rm(gene_frame, isoform_frame)

# Get a CuffFeatureSet or CuffGeneSet of significant genes.

# TODO: Comment-out for the moment, as the data is currently not used.
# significant_gene_ids <- getSig(object = cuff_set, level = "genes")
# significant_genes <- getGenes(object = cuff_set, geneIdList = significant_gene_ids)

# The csHeatmap plot does not seem to be a sensible option for larger sets of significant genes.
# for (i in 1:length(sample_pairs[1,])) {

# significant_genes_diff <-
# }

# Finally, create comparison-specific relative symbolic links for cuffdiff
# results in the rnaseq_process_cuffdiff_* directory to avoid identical file
# names between comparisons.

message("Started creating symbolic links to cuffdiff results")

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

message("Finished creating symbolic links to cuffdiff results")

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
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

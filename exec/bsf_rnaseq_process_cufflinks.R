#!/usr/bin/env Rscript
#
# BSF R script to post-processes Cufflinks output by enriching gene and
# transcript tables with Ensembl annotation, either downloaded from BioMart or
# imported from a reference GTF file.
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
        opt_str = c("--gtf-reference"),
        default = NULL,
        dest = "gtf_reference",
        help = "GTF file specifying a reference transcriptome",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--genome-version"),
        default = NULL,
        dest = "genome_version",
        help = "Genome version",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--biomart-instance"),
        default = "ENSEMBL_MART_ENSEMBL",
        dest = "biomart_instance",
        help = "BioMart instance [ENSEMBL_MART_ENSEMBL]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--biomart-data-set"),
        dest = "biomart_data_set",
        help = "BioMart data set",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--biomart-host"),
        dest = "biomart_host",
        default = "www.ensembl.org",
        help = "BioMart host [www.ensembl.org]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--biomart-port"),
        default = 80L,
        dest = "biomart_port",
        help = "BioMart port [80]",
        type = "integer"
      ),
      optparse::make_option(
        opt_str = c("--biomart-path"),
        default = "/biomart/martservice",
        dest = "biomart_path",
        help = "BioMart path [/biomart/martservice]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--pattern-path"),
        default = "/rnaseq_cufflinks_.*$",
        dest = "pattern_path",
        help = "Cufflinks directory path pattern [/rnaseq_cufflinks_.*$]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--pattern-sample"),
        default = ".*/rnaseq_cufflinks_(.*)$",
        dest = "pattern_sample",
        help = "Cufflinks sample name pattern [.*/rnaseq_cufflinks_(.*)$]",
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

# Validate the argument_list.

if (is.null(x = argument_list$biomart_data_set)) {
  # If a --biomart-data-set was not specified, ...
  if (is.null(x = argument_list$gtf_reference)) {
    # ... a --gtf-reference option must be specified.
    stop("Missing --gtf-reference or --biomart-data-set option")
  }
}

suppressPackageStartupMessages(expr = library(package = "tidyverse"))

# Save plots in the following formats.

graphics_formats <- c("pdf" = "pdf", "png" = "png")

gene_annotation_tibble <- NULL
isoform_annotation_tibble <- NULL

if (is.null(x = argument_list$biomart_data_set)) {
  # Import GTF file-based annotation --------------------------------------

  if (is.null(x = argument_list$gtf_reference)) {
    stop("Missing --gtf-reference option")
  }
  # If a --genome-version option was not provided, set it to NA.
  if (is.null(x = argument_list$genome_version)) {
    argument_list$genome_version <- NA_character_
  }

  suppressPackageStartupMessages(expr = library(package = "Biostrings"))
  suppressPackageStartupMessages(expr = library(package = "rtracklayer"))

  message("Import GTF annotation file")
  reference_granges <- rtracklayer::import(
    con = argument_list$gtf_reference,
    format = "gtf",
    genome = argument_list$genome_version,
    feature.type = "transcript"
  )
  reference_mcols <- S4Vectors::mcols(x = reference_granges)

  # Standard GTF attributes
  gene_annotation_tibble <- tibble::tibble(
    "ensembl_gene_id" = reference_mcols$gene_id,
    "ensembl_gene_version" = reference_mcols$gene_version,
    "ensembl_gene_name" = reference_mcols$gene_name,
    "ensembl_gene_source" = reference_mcols$gene_source,
    "ensembl_gene_biotype" = reference_mcols$gene_biotype
    # "havana_gene" = if ("havana_gene" %in% names(x = reference_mcols))
    #   reference_mcols$havana_gene
    # else
    #   NA_character_,
    # "havana_gene_version" = if ("havana_gene_version" %in% names(x = reference_mcols))
    #   reference_mcols$havana_gene_version
    # else
    #   NA_character_
  )
  # Optional GTF attributes
  # Havana annotation is only avaliable for earlier Ensembl releases.
  if ("havana_gene" %in% names(x = reference_mcols)) {
    gene_annotation_tibble$havana_gene <-
      reference_mcols$havana_gene
  }
  if ("havana_gene_version" %in% names(x = reference_mcols)) {
    gene_annotation_tibble$havana_gene_version <-
      reference_mcols$havana_gene_version
  }

  # Gene information is available for each transcript.
  gene_annotation_tibble <-
    dplyr::distinct(.data = gene_annotation_tibble)

  # Standard GTF attributes
  isoform_annotation_tibble <- tibble::tibble(
    "ensembl_transcript_id" = reference_mcols$transcript_id,
    "ensembl_transcript_version" = reference_mcols$transcript_version,
    "ensembl_transcript_name" = reference_mcols$transcript_name,
    "ensembl_transcript_source" = reference_mcols$transcript_source,
    "ensembl_transcript_biotype" = reference_mcols$transcript_biotype,
    "ensembl_gene_id" = reference_mcols$gene_id,
    "ensembl_gene_version" = reference_mcols$gene_version,
    "ensembl_gene_name" = reference_mcols$gene_name,
    "ensembl_gene_source" = reference_mcols$gene_source,
    "ensembl_gene_biotype" = reference_mcols$gene_biotype,
    "ccds_id" = if ("ccds_id" %in% names(x = reference_mcols))
      reference_mcols$ccds_id
    else
      NA_character_,
    "tag" = if ("tag" %in% names(x = reference_mcols))
      reference_mcols$tag
    else
      NA_character_
    # "havana_transcript" = if ("havana_transcript" %in% names(x = reference_mcols))
    #   reference_mcols$havana_transcript
    # else
    #   NA_character_,
    # "havana_transcript_version" = if ("havana_transcript_version" %in% names(x = reference_mcols))
    #   reference_mcols$havana_transcript_version
    # else
    #   NA_character_,
    # "havana_transcript_support_level" = if ("havana_transcript_support_level" %in% names(x = reference_mcols))
    #   reference_mcols$havana_transcript_support_level
    # else
    #   NA_character_,
    # "havana_gene" = if ("havana_gene" %in% names(x = reference_mcols))
    #   reference_mcols$havana_gene
    # else
    #   NA_character_,
    # "havana_gene_version" = if ("havana_gene_version" %in% names(x = reference_mcols))
    #   reference_mcols$havana_gene_version
    # else
    #   NA_character_
  )
  # Optional GTF attributes
  # if ("ccds_id" %in% names(x = reference_mcols)) {
  #   isoform_annotation_tibble$ccds_id <-
  #     reference_mcols$ccds_id
  # }
  # if ("tag" %in% names(x = reference_mcols)) {
  #   isoform_annotation_tibble$tag <-
  #     reference_mcols$tag
  # }
  # Havana annotation is only avaliable for earlier Ensembl releases.
  if ("havana_transcript" %in% names(x = reference_mcols)) {
    isoform_annotation_tibble$havana_transcript <-
      reference_mcols$havana_transcript
  }
  if ("havana_transcript_version" %in% names(x = reference_mcols)) {
    isoform_annotation_tibble$havana_transcript_version <-
      reference_mcols$havana_transcript_version
  }
  if ("havana_transcript_support_level" %in% names(x = reference_mcols)) {
    isoform_annotation_tibble$havana_transcript_support_level <-
      reference_mcols$havana_transcript_support_level
  }
  if ("havana_gene" %in% names(x = reference_mcols)) {
    isoform_annotation_tibble$havana_gene <-
      reference_mcols$havana_gene
  }
  if ("havana_gene_version" %in% names(x = reference_mcols)) {
    isoform_annotation_tibble$havana_gene_version <-
      reference_mcols$havana_gene_version
  }

  rm(reference_mcols, reference_granges)
} else {
  # Import BioMart-based annotation ---------------------------------------

  suppressPackageStartupMessages(expr = library(package = "biomaRt"))

  # Connect to the Ensembl BioMart.
  message("Connect to BioMart")

  ensembl_mart <- biomaRt::useMart(
    biomart = argument_list$biomart_instance,
    dataset = argument_list$biomart_data_set,
    host = argument_list$biomart_host,
    path = argument_list$biomart_path,
    port = argument_list$biomart_port
  )

  message("Loading attribute data from BioMart")

  ensembl_attributes <- biomaRt::listAttributes(
    mart = ensembl_mart,
    page = "feature_page",
    what = c("name", "description", "page")
  )

  # Get Ensembl Gene information. From Ensembl version e75, the attributes
  # "external_gene_id" and "external_gene_db" are called "external_gene_name" and
  # "external_gene_source", respectively.

  if ("external_gene_id" %in% ensembl_attributes$name) {
    # Pre e75.
    ensembl_gene_attributes <- c("ensembl_gene_id",
                                 "external_gene_id",
                                 "external_gene_db",
                                 "gene_biotype")
  } else if ("external_gene_name" %in% ensembl_attributes$name) {
    # Post e75.
    ensembl_gene_attributes <- c("ensembl_gene_id",
                                 "external_gene_name",
                                 "external_gene_db",
                                 "gene_biotype")
  } else {
    stop(
      "Neither external_gene_id nor external_gene_name are available as BioMart attributes."
    )
  }

  message("Loading gene data from BioMart")

  gene_annotation_tibble <-
    tibble::as_tibble(x = biomaRt::getBM(attributes = ensembl_gene_attributes, mart = ensembl_mart))
  rm(ensembl_gene_attributes)

  # Get Ensembl Transcript information.

  if ("external_gene_id" %in% ensembl_attributes$name) {
    # Pre e75.
    ensembl_transcript_attributes <- c(
      "ensembl_transcript_id",
      "external_transcript_id",
      "transcript_db_name",
      "transcript_biotype",
      "ensembl_gene_id",
      "external_gene_id",
      "external_gene_db",
      "gene_biotype"
    )
  } else if ("external_gene_name" %in% ensembl_attributes$name) {
    # Post e75.
    ensembl_transcript_attributes <- c(
      "ensembl_transcript_id",
      "external_transcript_name",
      "transcript_db_name",
      "transcript_biotype",
      "ensembl_gene_id",
      "external_gene_name",
      "external_gene_db",
      "gene_biotype"
    )
  } else {
    stop(
      "Neither external_gene_id nor external_gene_name are available as BioMart attributes."
    )
  }

  message("Loading transcript data from BioMart")
  isoform_annotation_tibble <-
    tibble::as_tibble(biomaRt::getBM(attributes = ensembl_transcript_attributes, mart = ensembl_mart))
  rm(ensembl_transcript_attributes)

  # Destroy and thus discconnect the Ensembl BioMart connection already here.

  rm(ensembl_mart, ensembl_attributes)
}

# Process all "rnaseq_cufflinks_*" directories in the current working directory.

message("Processing Cufflinks reports for sample:")
sample_tibble <- tibble::tibble(
  # R character vector of directory paths.
  "directory_path" = grep(
    pattern = argument_list$pattern_path,
    x = list.dirs(
      path = argument_list$genome_directory,
      full.names = TRUE,
      recursive = FALSE
    ),
    value = TRUE
  ),
  # R character vector of sample names.
  "sample_name" = gsub(
    pattern = argument_list$pattern_sample,
    replacement = "\\1",
    x = .data$directory_path
  )
)

# Add further variables for all FPKM_status bio types and states.
for (biotype in c("genes", "isoforms")) {
  for (fpkm_status in c("OK", "LOWDATA", "HIDATA", "FAIL")) {
    sample_tibble[, paste("FPKM_status", biotype, fpkm_status, sep = "_")] <-
      0L
  }
  rm(fpkm_status)
}
rm(biotype)

for (i in seq_len(length.out = nrow(x = sample_tibble))) {
  message("  ", sample_tibble$sample_name[i])

  # Construct sample-specific prefixes for Cufflinks and Tophat directories.
  prefix_cufflinks <-
    paste("rnaseq", "cufflinks", sample_tibble$sample_name[i], sep = "_")

  output_directory <-
    file.path(argument_list$output_directory, prefix_cufflinks)
  if (!file.exists(output_directory)) {
    dir.create(path = output_directory,
               showWarnings = TRUE,
               recursive = FALSE)
  }

  # Read, summarise, merge and write gene (genes.fpkm_tracking) tables.

  tracking_tibble <-
    dplyr::group_by(
      .data = readr::read_tsv(
        file = file.path(
          argument_list$genome_directory,
          prefix_cufflinks,
          "genes.fpkm_tracking"
        ),
        col_names = TRUE,
        col_types = readr::cols(
          "tracking_id" = readr::col_character(),
          # Not populated for genes and isoforms.
          "class_code" = readr::col_factor(),
          # Not populated for genes and isoforms.
          "nearest_ref_id" = readr::col_character(),
          "gene_id" = readr::col_character(),
          # Not populated for genes and isoforms.
          "gene_short_name" = readr::col_character(),
          # Not populated for genes and isoforms.
          "tss_id" = readr::col_character(),
          "locus" = readr::col_character(),
          # Not populated (i.e. "-") for genes, read as factor..
          "length" = readr::col_factor(),
          # Not populated (i.e. "-") for genes, read as factor.
          "coverage" = readr::col_factor(),
          "FPKM" = readr::col_double(),
          "FPKM_conf_lo" = readr::col_double(),
          "FPKM_conf_hi" = readr::col_double(),
          "FPKM_status" = readr::col_factor()
        )
      ),
      .data$FPKM_status
    )

  # Summarise the FPKM_status variable and drop the grouping in the summary tibble.
  fpkm_status_tibble <-
    dplyr::summarise(
      .data = tracking_tibble,
      "FPKM_status_count" = dplyr::n(),
      .groups = "drop"
    )
  fpkm_status_tibble <-
    dplyr::mutate(
      .data = fpkm_status_tibble,
      "FPKM_status_variable" = paste("FPKM_status_genes", .data$FPKM_status, sep = "_")
    )
  for (j in seq_along(along.with = fpkm_status_tibble$FPKM_status_count)) {
    sample_tibble[i, fpkm_status_tibble$FPKM_status_variable[j]] <-
      fpkm_status_tibble$FPKM_status_count[j]
  }
  rm(j, fpkm_status_tibble)

  readr::write_tsv(
    x = dplyr::left_join(
      x = gene_annotation_tibble,
      y = tracking_tibble,
      by = c("ensembl_gene_id" = "tracking_id")
    ),
    path = file.path(
      output_directory,
      paste(prefix_cufflinks, "genes_fpkm_tracking.tsv", sep = "_")
    )
  )
  rm(tracking_tibble)

  # Read, summarize, merge and write transcript (isoforms.fpkm_tracking) tables.

  tracking_tibble <-
    dplyr::group_by(
      .data = readr::read_tsv(
        file = file.path(
          argument_list$genome_directory,
          prefix_cufflinks,
          "isoforms.fpkm_tracking"
        ),
        col_names = TRUE,
        col_types = readr::cols(
          "tracking_id" = readr::col_character(),
          # Not populated for genes and isoforms.
          "class_code" = readr::col_factor(),
          # Not populated for genes and isoforms.
          "nearest_ref_id" = readr::col_character(),
          "gene_id" = readr::col_character(),
          # Not populated for genes and isoforms.
          "gene_short_name" = readr::col_character(),
          # Not populated for genes and isoforms.
          "tss_id" = readr::col_character(),
          "locus" = readr::col_character(),
          # Not populated for genes.
          "length" = readr::col_integer(),
          # Not populated for genes.
          "coverage" = readr::col_number(),
          "FPKM" = readr::col_double(),
          "FPKM_conf_lo" = readr::col_double(),
          "FPKM_conf_hi" = readr::col_double(),
          "FPKM_status" = readr::col_factor()
        )
      ),
      .data$FPKM_status
    )

  # Summarise the FPKM_status variable and drop the grouping in the summary tibble.
  fpkm_status_tibble <-
    dplyr::summarise(
      .data = tracking_tibble,
      "FPKM_status_count" = dplyr::n(),
      .groups = "drop"
    )
  fpkm_status_tibble <-
    dplyr::mutate(
      .data = fpkm_status_tibble,
      "FPKM_status_variable" = paste("FPKM_status_isoforms", .data$FPKM_status, sep = "_")
    )
  for (j in seq_along(along.with = fpkm_status_tibble$FPKM_status_count)) {
    sample_tibble[i, fpkm_status_tibble$FPKM_status_variable[j]] <-
      fpkm_status_tibble$FPKM_status_count[j]
  }
  rm(j, fpkm_status_tibble)

  readr::write_tsv(
    x = dplyr::left_join(
      x = isoform_annotation_tibble,
      y = tracking_tibble,
      by = c("ensembl_transcript_id" = "tracking_id")
    ),
    path = file.path(
      output_directory,
      paste(prefix_cufflinks, "isoforms_fpkm_tracking.tsv", sep = "_")
    )
  )
  rm(tracking_tibble,
     output_directory,
     prefix_cufflinks)
}
rm(i)

readr::write_tsv(
  x = sample_tibble,
  path = file.path(
    argument_list$output_directory,
    "rnaseq_cufflinks_summary.tsv"
  )
)

rm(
  sample_tibble,
  isoform_annotation_tibble,
  gene_annotation_tibble,
  graphics_formats,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

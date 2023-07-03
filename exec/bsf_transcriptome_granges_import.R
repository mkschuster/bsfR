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


# BSF R Script to import a GenomicRanges::GRanges object.
#
# The --drop-variables option allows for specifying metadata variables
# that should be dropped from the GenomicRanges::GRanges object.
#
# Since Ensembl provides multiple "tag" attributes per GTF line, the value of
# only one of them can be included in the meta data DataFrame. Since a single
# "tag" value renders the annotation misleading, the "tag" variable can be
# dropped.
#
# Please note that Ensembl provides a singe tag attribute in GFF3 files with one
# or more comma-separated tag values, which can be parsed into
# IRanges::CharacterList objects.
#
# Output file:
# - <output_directory>/<transcriptome_version>_granges.rds
# - <output_directory>/<transcriptome_version>.bed

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
        opt_str = "--transcriptome-path",
        dest = "transcriptome_path",
        help = "Transcriptome (GTF) path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--transcriptome-version",
        dest = "transcriptome_version",
        help = "Transcriptome version",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--seqinfo-path",
        dest = "seqinfo_path",
        help = "GenomeInfoDb::Seqinfo path (preferred over --genome-version)",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--genome-version",
        dest = "genome_version",
        help = "Genome version (--seqinfo-path is preferred)",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--drop-variables",
        dest = "drop_variables",
        help = "A comma-separated list of meta data variables to drop",
        type = "character"
      ),
      optparse::make_option(
        opt_str = "--output-directory",
        default = ".",
        dest = "output_directory",
        help = "Output directory path [.]",
        type = "character"
      )
    )
  ))

if (is.null(x = argument_list$transcriptome_path)) {
  stop("Missing --transcriptome-path option")
}

if (is.null(x = argument_list$transcriptome_version)) {
  stop("Missing --transcriptome-version option")
}

if (is.null(x = argument_list$genome_version) &&
    is.null(x = argument_list$seqinfo_path)) {
  stop(
    "Either the --genome-version or the --seqinfo-path option ",
    "is required for output file naming."
  )
}

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "readr"))
suppressPackageStartupMessages(expr = library(package = "stringr"))
# Bioconductor
suppressPackageStartupMessages(expr = library(package = "BiocVersion"))
suppressPackageStartupMessages(expr = library(package = "rtracklayer"))

seqinfo_object <-
  if (is.null(x = argument_list$seqinfo_path)) {
    GenomeInfoDb::Seqinfo(genome = argument_list$genome_version)
  } else {
    readr::read_rds(file = argument_list$seqinfo_path)
  }

if (is.null(x = seqinfo_object)) {
  stop("Could not retrieve a GenomeInfoDb::Seqinfo object")
}

granges_object <-
  rtracklayer::import(con = argument_list$transcriptome_path,
                      genome = seqinfo_object)

# Drop metadata variables if requested and present.
if (!is.null(argument_list$drop_variables)) {
  # Split on "," characters.
  drop_variables_character <-
    stringr::str_split_1(string = argument_list$drop_variables,
                         pattern = stringr::fixed(pattern = ","))

  # Trim remaining white space around "," characters.
  drop_variables_character <-
    stringr::str_trim(string = drop_variables_character)

  # Only subset if the variable is not just a single empty string.
  if ((length(x = drop_variables_character) > 1) ||
      (drop_variables_character[1L] != "")) {
    if (argument_list$verbose) {
      message("Removing metadata variables: ",
              paste(drop_variables_character, collapse = " "))
    }

    mcols_dframe <- S4Vectors::mcols(x = granges_object)

    mcols_dframe <- S4Vectors::subset(
      x = mcols_dframe,
      select = base::setdiff(x = S4Vectors::colnames(x = mcols_dframe),
                             y = drop_variables_character),
      drop = FALSE
    )

    S4Vectors::mcols(x = granges_object) <- mcols_dframe
    rm(mcols_dframe)
  }
  rm(drop_variables_character)
}

if (TRUE) {
  mcols_dframe <- S4Vectors::mcols(x = granges_object)
  if ("type" %in% S4Vectors::colnames(x = mcols_dframe)) {
    # Check if the GTF file has only "exon" features. If so, create also
    # "transcript" and "gene" features.
    if (sum(mcols_dframe$type == "exon") == S4Vectors::nrow(x = mcols_dframe)) {
      # Only "exon records are available.
      transcript_granges <-
        granges_object

      transcript_mcols_dframe <-
        S4Vectors::mcols(x = transcript_granges)

      transcript_mcols_dframe$type <- "transcript"

      if ("exon_id" %in% S4Vectors::colnames(x = transcript_mcols_dframe)) {
        transcript_mcols_dframe$exon_id <- NA_character_
      }

      S4Vectors::mcols(x = transcript_granges) <-
        transcript_mcols_dframe
      rm(transcript_mcols_dframe)

      gene_granges <- transcript_granges

      gene_mcols_dframe <- S4Vectors::mcols(x = gene_granges)

      gene_mcols_dframe$type <- "gene"

      if ("transcript_id" %in% S4Vectors::colnames(x = gene_mcols_dframe)) {
        gene_mcols_dframe$transcript_id <- NA_character_
      }

      S4Vectors::mcols(x = gene_granges) <- gene_mcols_dframe
      rm(gene_mcols_dframe)

      # Merge "exon", "transcript" and "gene" GRanges.
      granges_object <-
        c(granges_object, transcript_granges, gene_granges)

      rm(transcript_granges, gene_granges)

      granges_object <- GenomicRanges::sort(x = granges_object)
    }
  }
  rm(mcols_dframe)
}

readr::write_rds(
  x = granges_object,
  file = file.path(argument_list$output_directory,
                   paste(
                     paste(argument_list$transcriptome_version, "granges", sep = "_"),
                     "rds",
                     sep = "."
                   )),
  compress = "gz"
)

# Re-export the merged GenomicRanges::GRanges objects in GTF format.
# The rtracklayer::export() function checks for a single genome.
merged_seqinfo_object <- GenomicRanges::seqinfo(x = granges_object)
GenomeInfoDb::genome(x = merged_seqinfo_object) <-
  paste(unique(x = GenomeInfoDb::genome(x = merged_seqinfo_object)), collapse = "_")
GenomicRanges::seqinfo(x = granges_object) <- merged_seqinfo_object
rm(merged_seqinfo_object)

rtracklayer::export(object = granges_object,
                    con = file.path(
                      argument_list$output_directory,
                      paste(argument_list$transcriptome_version,
                            "gtf",
                            sep = ".")
                    ))

rm(granges_object, seqinfo_object, argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessioninfo::session_info())

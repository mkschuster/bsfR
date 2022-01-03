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


# BSF R script to merge data frames of each replicate after the Cufflinks
# RNA-seq analysis stage (bsf_rnaseq_process_cufflinks.R) to come up with
# consistent tables for genes (rnaseq_cufflinks_isoforms_fpkm_tracking.tsv) and
# isoforms (rnaseq_cufflinks_isoforms_fpkm_tracking.tsv).
#
# Since FPKM values have not been normalised at this stage, a direct comparison
# between replicates is strictly not possible. Hence, this script is for special
# purposes only and not part of the standard pipeline. Please see the data
# frames after the Cuffdiff RNA-seq analysis stage that provide normalised and
# thus perfectly comaprable data.

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
      )
    )
  ))

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))


#' Process a Cufflinks FPKM tracking table for a particular, named replicate.
#'
#' @noRd
#' @param replicate_name: Replicate name
#' @type replicate_name: char
#' @param object_type: Cufflinks object type i.e. genes or isoforms
#' @type object_type: char
#' @param merge_frame: Merge data frame
#' @type merge_frame: data.frame
#' @return: Merged data frame
#' @rtype: data.frame

process_cufflinks_table <-
  function(replicate_name,
           object_type,
           merge_frame = NULL) {
    if (is.null(x = replicate_name)) {
      stop("Missing replicate_name argument")
    }

    if (is.null(x = object_type)) {
      stop("Missing object_type argument")
    }

    message(paste("Processing", object_type, replicate_name, sep = ' '))

    # Construct replicate-specific prefixes for Cufflinks directories.
    prefix_cufflinks <-
      paste("rnaseq", "cufflinks", replicate_name, sep = "_")

    cufflinks_frame <-
      utils::read.table(
        file = file.path(
          prefix_cufflinks,
          paste(prefix_cufflinks, object_type, "fpkm_tracking.tsv", sep = "_")
        ),
        header = TRUE,
        stringsAsFactors = FALSE
      )

    # Subset the data frame by removing unused columns.
    obsolete_columns <- NULL
    if (object_type == "genes") {
      obsolete_columns <- c(
        "class_code",
        "nearest_ref_id",
        "gene_short_name",
        "tss_id",
        "length",
        "coverage"
      )
    } else if (object_type == "isoforms") {
      obsolete_columns <- c("class_code",
                            "nearest_ref_id",
                            "gene_id",
                            "gene_short_name",
                            "tss_id")
    }

    subset_frame <-
      cufflinks_frame[, !(names(x = cufflinks_frame) %in% obsolete_columns), drop = FALSE]
    rm(cufflinks_frame, obsolete_columns)

    # Select only Ensembl objects i.e., those that have a gene_biotype set.
    ensembl_frame <-
      subset_frame[!is.na(x = subset_frame$gene_biotype), , drop = FALSE]
    rm(subset_frame)

    # Rename columns to include the replicate name.
    if (object_type == "isoforms") {
      names(x = ensembl_frame)[grepl('^coverage$', names(x = ensembl_frame))] <-
        paste('coverage', replicate_name, sep = '_')
    }

    names(x = ensembl_frame)[grepl('^FPKM$', names(x = ensembl_frame))] <-
      paste('FPKM', replicate_name, sep = '_')

    names(x = ensembl_frame)[grepl('^FPKM_conf_lo$', names(x = ensembl_frame))] <-
      paste('FPKM_conf_lo', replicate_name, sep = '_')

    names(x = ensembl_frame)[grepl('^FPKM_conf_hi$', names(x = ensembl_frame))] <-
      paste('FPKM_conf_hi', replicate_name, sep = '_')

    names(x = ensembl_frame)[grepl('^FPKM_status$', names(x = ensembl_frame))] <-
      paste('FPKM_status', replicate_name, sep = '_')

    if (is.null(x = merge_frame)) {
      return(ensembl_frame)
    } else {
      # Merge the data frames. The 'by' parameter defaults to an intersection of the column names.
      return(base::merge.data.frame(
        x = merge_frame,
        y = ensembl_frame,
        all = TRUE,
        sort = TRUE
      ))
    }
  }

# Process all "rnaseq_cufflinks_*" directories in the current working directory.
# List all rnaseq_cufflinks directories via their common prefix and
# parse the sample (or replicate) name simply by removing the prefix.

replicate_names <- sub(
  pattern = "^rnaseq_cufflinks_",
  replacement = "",
  x = grep(
    pattern = '^rnaseq_cufflinks_',
    x = list.dirs(full.names = FALSE, recursive = FALSE),
    value = TRUE
  )
)

# Process genes tables and write the merged data frame onto disk.
merge_frame <- NULL
for (replicate_name in replicate_names) {
  merge_frame <- process_cufflinks_table(
    replicate_name = replicate_name,
    object_type = "genes",
    merge_frame = merge_frame
  )
}
rm(replicate_name)

utils::write.table(
  x = merge_frame,
  file = "rnaseq_cufflinks_genes_fpkm_tracking.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)
rm(merge_frame)

# Process isoforms tables and write the merged data frame onto disk.
merge_frame <- NULL
for (replicate_name in replicate_names) {
  merge_frame <- process_cufflinks_table(
    replicate_name = replicate_name,
    object_type = "isoforms",
    merge_frame = merge_frame
  )
}
rm(replicate_name)

utils::write.table(
  x = merge_frame,
  file = "rnaseq_cufflinks_isoforms_fpkm_tracking.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)
rm(merge_frame)

rm(replicate_names, process_cufflinks_table, argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessioninfo::session_info())

#!/usr/bin/env Rscript
#
# BSF R script to refine the GATK CallableLoci analysis.
#
# The coverage assessment loads (Ensembl) exon information, collates (or
# projects) all (overlapping) exons into transcribed regions on the genome,
# overlaps those with the target regions to get transcribed target regions and
# finally overlaps those with non-callable loci to get the minimal set of
# problematic regions. Each problematic region is annotated with target name, as
# well as exon, transcript and gene information.
#
# This script writes three files:
#
# A tab-separated value file of non-callable loci, annotated with Ensembl gene,
# transcript and exon information, as well as target information if available.
# In contrast to the file above, this is grouped by the non-callable regions so
# that gene, transcript and exon information is collated.
#   variant_calling_diagnose_sample_{sample_name}_non_callable_regions.tsv
#
# A tab-separated value file of non-callable loci, annotated with Ensembl gene,
# transcript and exon information, as well as target information if available.
#   variant_calling_diagnose_sample_{sample_name}_non_callable_loci.tsv
#
# A summary data frame of metrics collected along the procedure.
#   variant_calling_diagnose_sample_{sample_name}_non_callable_summary.tsv
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
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))

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
        opt_str = c("--exons"),
        dest = "exon_path",
        help = "File path to the gene, transcript and exon annotation GTF",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--exon-flanks"),
        default = 0L,
        dest = "exon_flanks",
        help = "Exon flanking regions [0]",
        type = "integer"
      ),
      optparse::make_option(
        opt_str = c("--no-filter"),
        action = "store_true",
        default = FALSE,
        dest = "no_filter",
        help = "Do not filter GTF annotation by 'tag \"basic\";' [FALSE]",
        type = "logical"
      ),
      optparse::make_option(
        opt_str = c("--callable-loci"),
        dest = "callable_loci_path",
        help = "File path to the GATK CallableLoci BED",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--targets"),
        dest = "target_path",
        help = "File path to the enrichment targets BED",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--genome-version"),
        default = "hs37d5",
        dest = "genome_version",
        help = "Genome version [hs37d5]",
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

if (file.size(argument_list$callable_loci_path) > 1e+09) {
  message("The file specified by the --callable-loci option is too large to process.")
  print(x = sessioninfo::session_info())
  quit()
}

suppressPackageStartupMessages(expr = library(package = "bsfR"))
suppressPackageStartupMessages(expr = library(package = "Biostrings"))
suppressPackageStartupMessages(expr = library(package = "GenomicRanges"))

prefix <- "variant_calling_diagnose_sample"

if (!file.exists(argument_list$output_directory)) {
  dir.create(
    path = argument_list$output_directory,
    showWarnings = TRUE,
    recursive = FALSE
  )
}

# Import constrained target GRanges ---------------------------------------


summary_list <-
  bsfR::bsfvc_import_constrained_granges(
    exon_path = argument_list$exon_path,
    exon_flanks = argument_list$exon_flanks,
    exon_basic = !argument_list$no_filter,
    target_path = argument_list$target_path,
    target_flanks = argument_list$exon_flanks,
    genome_version = argument_list$genome_version,
    verbose = argument_list$verbose
  )

# Eventually, delete those relatively big GenomicRanges::GRanges objects.
# summary_list$transcribed_granges <- NULL

# Import GATK CallableLoci GRanges ----------------------------------------


# Import the BED file produced by the GATK CallableLoci analysis.
summary_list$callable_loci_path <-
  argument_list$callable_loci_path

# Store the sample name in the summary list.
summary_list$sample_name <-
  base::gsub(
    pattern = paste(paste0(".*", prefix), "(.*?)_callable_loci.bed$", sep = "_"),
    replacement = "\\1",
    x = summary_list$callable_loci_path
  )

if (argument_list$verbose) {
  message("Importing GATK CallableLoci GRanges: ",
          summary_list$callable_loci_path)
}

callable_granges <-
  rtracklayer::import(
    con = summary_list$callable_loci_path,
    format = "bed",
    genome = argument_list$genome_version
  )

# Rename the "name" variable of the CallableLoci GenomicRanges::GRanges to the
# more informative "mapping_status".
S4Vectors::mcols(x = callable_granges) <-
  S4Vectors::rename(x = S4Vectors::mcols(x = callable_granges), c("name" = "mapping_status"))

# Filter out "CALLABLE" GenomicRanges::GRanges.
non_callable_granges <-
  callable_granges[callable_granges$mapping_status != "CALLABLE", ]
rm(callable_granges)

non_callable_granges$mapping_status <-
  as.factor(x = non_callable_granges$mapping_status)

summary_list$non_callable_number_raw <-
  length(x = non_callable_granges)

summary_list$non_callable_width_raw <-
  sum(GenomicRanges::width(x = non_callable_granges))

if (argument_list$verbose) {
  message("Processing sample name: ", summary_list$sample_name)

  message("Number of non-callable raw GRanges: ",
          summary_list$non_callable_number_raw)

  message(
    "Cumulative width of non-callable raw GRanges: ",
    summary_list$non_callable_width_raw
  )
}

# Overlap non-callable GRanges with constrained GRanges -------------------


# To get accurate non-callable statistics with regards to the target regions
# that are transcribed, non-callable GenomicRanges::GRanges need overlapping
# with the constrained GenomicRanges::GRanges.
overlap_granges <-
  IRanges::pintersect(
    x = IRanges::findOverlapPairs(
      query = non_callable_granges,
      subject = summary_list$constrained_granges,
      ignore.strand = TRUE
    ),
    ignore.strand = TRUE
  )

summary_list[["non_callable_number_constrained.TOTAL"]] <-
  length(x = overlap_granges)

summary_list[["non_callable_width_constrained.TOTAL"]] <-
  sum(GenomicRanges::width(x = overlap_granges))

if (argument_list$verbose) {
  message("Number of non-callable constrained GRanges: ",
          summary_list[["non_callable_number_constrained.TOTAL"]])

  message("Cumulative width of non-callable constrained GRanges: ",
          summary_list[["non_callable_width_constrained.TOTAL"]])
}

# Summarise by mapping status ---------------------------------------------


# Summarise also separately by mapping status. Populate the summary frame with
# columns for each mapping status level, regardless of whether it is associated
# with data or not. Use fixed mapping status levels emitted by GATK
# CallableLoci.
for (level in c(
  "REF_N",
  "NO_COVERAGE",
  "LOW_COVERAGE",
  "EXCESSIVE_COVERAGE",
  "POOR_MAPPING_QUALITY"
)) {
  summary_list[[paste("non_callable_number_constrained", level, sep = ".")]] <-
    0L

  summary_list[[paste("non_callable_width_constrained", level, sep = ".")]] <-
    0L
}
rm(level)

if (length(x = overlap_granges) > 0L) {
  # Count the number of entries for each mapping_status level.
  aggregate_frame <-
    as.data.frame(x = base::table(S4Vectors::mcols(x = overlap_granges)$mapping_status))

  # Assign the result levels (rows) as summary frame columns.
  for (j in seq_len(length.out = nrow(x = aggregate_frame))) {
    summary_list[[paste("non_callable_number_constrained",
                        aggregate_frame[j, 1L, drop = TRUE],
                        sep = ".")]] <-
      aggregate_frame[j, 2L, drop = TRUE]
  }
  rm(j, aggregate_frame)

  # Sum the widths of entries for each mapping_status level.
  aggregate_frame <-
    aggregate.data.frame(
      x = data.frame(width = GenomicRanges::width(x = overlap_granges)),
      by = list(mapping_status = S4Vectors::mcols(x = overlap_granges)$mapping_status),
      FUN = sum
    )

  # Assign the result levels (rows) as summary frame columns.
  for (j in seq_len(length.out = nrow(x = aggregate_frame))) {
    summary_list[[paste("non_callable_width_constrained", aggregate_frame[j, 1L, drop = TRUE], sep = ".")]] <-
      aggregate_frame[j, 2L, drop = TRUE]
  }
  rm(j, aggregate_frame)
}
rm(overlap_granges)

# Annotate non-callable GRanges with target names -------------------------


# Annotate the non-callable GenomicRanges::GRanges with target names if
# available.
diagnose_granges <- NULL
if (!is.null(x = summary_list$target_path)) {
  if (argument_list$verbose) {
    message("Annotate the diagnostic GRanges with target GRanges ...")
  }
  # If the target GenomicRanges::GRanges are available, merge by overlap with
  # the non-callable GenomicRanges::GRanges into a new S4Vectors::DataFrame.
  overlap_dframe <-
    IRanges::mergeByOverlaps(query = summary_list$target_granges, subject = non_callable_granges)
  # Extract the "non_callable_granges".
  diagnose_granges <-
    overlap_dframe[, deparse(expr = quote(expr = non_callable_granges)), drop = TRUE]
  # Annotate with the target_name from the target GenomicRanges::GRanges object.
  S4Vectors::mcols(x = diagnose_granges)$target_name <-
    overlap_dframe$name
  rm(overlap_dframe)
} else {
  # Diagnose GenomicRanges::GRanges without target annotation.
  diagnose_granges <- non_callable_granges
}

# Annotate non-callable GRanges with exon names ---------------------------


if (argument_list$verbose) {
  message("Annotate the diagnostic GRanges with exon GRanges ...")
}
# To annotate non-callable regions, merge by overlap with the exon
# GenomicRanges::GRanges.
overlap_dframe <-
  IRanges::mergeByOverlaps(query = diagnose_granges, subject = summary_list$exon_granges)

# The IRanges::mergeByOverlaps() function returns a S4Vectors::DataFrame of
# query (diagnose GenomicRanges::GRanges) and subject (exon
# GenomicRanges::GRanges) variables, as well as all S4Vectors::mcols() variables
# that were present in either GenomicRanges::GRanges object. To remove these
# redundant S4Vectors::mcols() variables, keep only the GenomicRanges::GRanges
# objects themselves.
#
# The IRanges::mergeByOverlaps() function is rather peculiar, as it constructs
# names for query and subject by quoting and de-parsing.

utils::write.table(
  x = overlap_dframe[, c(deparse(expr = quote(expr = diagnose_granges)),
                         deparse(expr = quote(expr = summary_list$exon_granges))),
                     drop = FALSE],
  file = file.path(
    argument_list$output_directory,
    paste(
      prefix,
      summary_list$sample_name,
      "non_callable_loci.tsv",
      sep = "_"
    )
  ),
  quote = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

# Group GRanges by non-callable regions -----------------------------------


# Since the IRanges::mergeByOverlaps() function provides the Cartesian product
# of diagnosis and exon GenomicRanges::GRanges, the table can be rather long and
# unwieldy. Annotate the diagnosis GenomicRanges::GRanges with exon
# GenomicRanges::GRanges meta information before grouping by diagnosis
# GenomicRanges::GRanges so that there is a single observation for each
# problematic region. Gene and transcript identifiers and names are turned into
# comma-separated lists.

if (argument_list$verbose) {
  message("Group the diagnostic GRanges ...")
}

# Simplify access to the rather complicated merged frame, via quoting and de-parsing.
overlap_exon_dframe <-
  S4Vectors::mcols(x = overlap_dframe[, deparse(expr = quote(expr = summary_list$exon_granges)), drop  = TRUE])

overlap_diagnose_granges <-
  overlap_dframe[, deparse(expr = quote(expr = diagnose_granges)), drop = TRUE]

overlap_diagnose_dframe <-
  S4Vectors::mcols(x = overlap_diagnose_granges)

overlap_diagnose_dframe$gene_id <-
  overlap_exon_dframe$gene_id

overlap_diagnose_dframe$gene_name <-
  overlap_exon_dframe$gene_name

overlap_diagnose_dframe$transcript_id <-
  overlap_exon_dframe$transcript_id

overlap_diagnose_dframe$transcript_name <-
  overlap_exon_dframe$transcript_name

overlap_diagnose_dframe$exon_id <-
  overlap_exon_dframe$exon_id

S4Vectors::mcols(x = overlap_diagnose_granges) <-
  overlap_diagnose_dframe

rm(overlap_diagnose_dframe, overlap_exon_dframe)

if (argument_list$verbose) {
  message("Group annotated diagnostic GRanges by region ...")
}

collapse_character <- function(x) {
  return(base::paste(base::unique(x = base::sort(x = x)), collapse = ","))
}

if (TRUE) {
  # Summarise the "gene_id", "transcript_id" and "exon_id" annotation by
  # non-callable region with Tidyverse code.

  diagnostic_tibble <-
    tibble::tibble(
      # The coercion of the GenomicRanges::GRanges class to the character class
      # provides a compact string (seqnames:start-end), while saving
      # GenomicRanges::GRanges as a S4Vectors::DFrame provides columns
      # "seqnames", "start", "end" and "width".
      "seqnames" = methods::as(
        object = GenomicRanges::seqnames(x = .env$overlap_diagnose_granges),
        Class = "character"
      ),
      "start" = GenomicRanges::start(x = .env$overlap_diagnose_granges),
      "end" = GenomicRanges::end(x = .env$overlap_diagnose_granges),
      "width" = GenomicRanges::width(x = .env$overlap_diagnose_granges),
      "strand" = methods::as(
        object = GenomicRanges::strand(x = .env$overlap_diagnose_granges),
        Class = "character"
      ),
      as.data.frame(x = S4Vectors::mcols(x = .env$overlap_diagnose_granges))
    )

  diagnostic_tibble <-
    dplyr::group_by(.data = diagnostic_tibble,
                    .data$seqnames,
                    .data$start,
                    .data$end,
                    .data$width,
                    .data$strand)

  diagnostic_tibble <- dplyr::summarise(
    .data = diagnostic_tibble,
    "maping_status" = dplyr::first(x = .data$mapping_status),
    "gene_ids" = collapse_character(x = .data$gene_id),
    "gene_names" = collapse_character(x = .data$gene_name),
    "transcript_ids" = collapse_character(x = .data$transcript_id),
    "transcript_names" = collapse_character(x = .data$transcript_name),
    "exon_ids" = collapse_character(x = .data$exon_id),
    .groups = "drop_last"
  )

  readr::write_tsv(x = diagnostic_tibble,
                   file = file.path(
                     argument_list$output_directory,
                     paste(
                       prefix,
                       summary_list$sample_name,
                       "non_callable_regions.tsv",
                       sep = "_"
                     )
                   ))

  rm(diagnostic_tibble)
} else {
  # Coerce the annotated diagnostic GenomicRanges::GRanges object into an
  # IRanges::CompressedManyToOneGrouping object, which specifies which groups
  # contain which indices to the original object.

  collapse_granges <- function(x) {
    group_granges <- overlap_diagnose_granges[x]
    group_mcols <- S4Vectors::mcols(x = group_granges)
    first_grange <- group_granges[1L]

    S4Vectors::mcols(x = first_grange) <- S4Vectors::DataFrame(
      # The "mapping_status" variable should have the same value for all
      # GenomicRanges::GRanges in the group. However, all other variables need
      # collapsing.
      "mapping_status" = group_mcols$mapping_status[1L],
      "gene_ids" = collapse_character(x = group_mcols$gene_id),
      "gene_names" = collapse_character(x = group_mcols$gene_name),
      "transcript_ids" = collapse_character(x = group_mcols$transcript_id),
      "transcript_names" = collapse_character(x = group_mcols$transcript_name),
      "exon_ids" = collapse_character(x = group_mcols$exon_id)
    )
    return(first_grange)
  }

  overlap_diagnose_grouping <-
    methods::as(object = overlap_diagnose_granges, "Grouping")

  # Collapse the IRanges::CompressedManyToOneGrouping object.
  grouped_granges <-
    unlist(x = GenomicRanges::GRangesList(lapply(X = overlap_diagnose_grouping, FUN = collapse_granges)))

  utils::write.table(
    x = grouped_granges,
    file = file.path(
      argument_list$output_directory,
      paste(
        prefix,
        summary_list$sample_name,
        "non_callable_regions.tsv",
        sep = "_"
      )
    ),
    quote = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  rm(grouped_granges,
     overlap_diagnose_grouping,
     collapse_granges)
}

rm(
  collapse_character,
  overlap_diagnose_granges,
  overlap_dframe,
  diagnose_granges,
  non_callable_granges
)

# Write the summary list as tibble -----------------------------------------


# Delete all GRanges objects from the summary list.
summary_list$exon_granges <- NULL
summary_list$transcribed_granges <- NULL
summary_list$target_granges <- NULL
summary_list$constrained_granges <- NULL

readr::write_tsv(
  x = tibble::as_tibble(x = summary_list),
  file = file.path(
    argument_list$output_directory,
    paste(
      prefix,
      summary_list$sample_name,
      "non_callable_summary.tsv",
      sep = "_"
    )
  )
)

rm(summary_list,
   prefix,
   argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessioninfo::session_info())

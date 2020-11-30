#!/usr/bin/env Rscript
#
# BSF R script to refine the GATK Callable Loci analyis.
#
# The coverage assessment loads (Ensembl) exon information, collates (or
# projects) all (overlapping) exons into transcribed regions on the genome,
# overlaps those with the target regions to get transcribed target regions and
# finally overlaps those with non-callable loci to get the minimal set of
# problematic regions. Each problematic region is annotated with target name, as
# well as exon, transcript and gene information. A summary data frame of metrics
# collected along the procedure is also written to disk.
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
        help = "File path to the GATK Callable Loci BED",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--targets"),
        dest = "target_path",
        help = "File path to the enrichment targets BED",
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
  message("The file specified by the --callable-loci option is too large to process")
  print(x = sessionInfo())
  quit()
}

suppressPackageStartupMessages(expr = library(package = "tidyverse"))
suppressPackageStartupMessages(expr = library(package = "rtracklayer"))

# Keep overall statistics in a summary list.
summary_list <- list()

# Import Ensembl annotation -----------------------------------------------


# Import Ensembl gene, transcript and exon annotation as GenomicRanges::GRanges
# vector object, where only exon components are relevant for this analysis.
summary_list$exon_path <- argument_list$exon_path

message("Importing exon GRanges: ", summary_list$exon_path)
exon_granges <-
  rtracklayer::import(con = summary_list$exon_path,
                      genome = "hs37d5",
                      feature.type = "exon")

summary_list$exon_number_raw <- length(x = exon_granges)
message("Number of raw exon GRanges: ", summary_list$exon_number_raw)

summary_list$exon_width_raw <-
  sum(GenomicRanges::width(x = exon_granges))
message("Cumulative width of raw exon GRanges: ",
        summary_list$exon_width_raw)

# Ensembl now annotates a "tag" in GFF files with value "basic" indicating
# standard (basic) transcript models. Can the exon GenomicRanges::GRanges be
# subset by such a tag?
if ("tag" %in% names(x = S4Vectors::mcols(x = exon_granges)) &&
    !argument_list$no_filter) {
  message("Filtering by GTF 'tag = \"basic\"' annotation")
  # Use the %in% operator for character matching as it sets NA values to FALSE,
  # automatically.
  exon_granges <-
    exon_granges[S4Vectors::mcols(x = exon_granges)$tag %in% "basic",]
}

summary_list$exon_number <- length(x = exon_granges)
message("Number of exon GRanges: ", summary_list$exon_number)

summary_list$exon_width <-
  sum(GenomicRanges::width(x = exon_granges))
message("Cumulative width of exon GRanges: ", summary_list$exon_width)

# Apply flanking regions --------------------------------------------------


# Apply flanking regions, by default 0L, to the exon GenomicRanges::GRanges.
exon_granges <-
  GenomicRanges::resize(
    x = exon_granges,
    width = GenomicRanges::width(x = exon_granges) + argument_list$exon_flanks,
    fix = "end"
  )
exon_granges <-
  GenomicRanges::resize(
    x = exon_granges,
    width = GenomicRanges::width(x = exon_granges) + argument_list$exon_flanks,
    fix = "start"
  )

summary_list$exon_flank_width <-
  sum(GenomicRanges::width(x = exon_granges))
message("Cumulative width of exon GRanges with flanks: ",
        summary_list$exon_flank_width)

# Reduce non-redundant Ensembl exons --------------------------------------


# Reduce the non-redundant Ensembl exons to their footprint on the genome to get
# transcribed regions. The revmap column contains the mapping to original
# exon_granges components.
message("Reducing exon GRanges to transcribed GRanges.")
transcribed_granges <-
  GenomicRanges::reduce(
    x = exon_granges,
    drop.empty.ranges = TRUE,
    with.revmap = TRUE,
    ignore.strand = TRUE
  )

summary_list$transcribed_number <-
  length(x = transcribed_granges)
message("Number of transcribed GRanges: ",
        summary_list$transcribed_number)

summary_list$transcribed_width <-
  sum(GenomicRanges::width(x = transcribed_granges))
message("Cumulative width of transcribed GRanges: ",
        summary_list$transcribed_width)

# Read target regions -----------------------------------------------------


# Read the file of targeted regions, if available. Although the GATK Callable
# Loci analysis is generally only run on these target regions, this file
# provides the target (probe) names of the enrichment design.
target_granges <- NULL
constrained_granges <- NULL
if (!is.null(x = argument_list$target_path)) {
  summary_list$target_path <- argument_list$target_path

  message("Importing target range annotation: ",
          summary_list$target_path)
  # The rtrackayer::import() function reads the genome version from the "db"
  # attribute of the BED "track" line.
  target_granges <-
    rtracklayer::import(con = summary_list$target_path)

  summary_list$target_number_raw <-
    length(x = target_granges)
  message("Number of target GRanges: ", summary_list$target_number_raw)

  summary_list$target_width_raw <-
    sum(GenomicRanges::width(x = target_granges))
  message("Cumulative width of target GRanges: ",
          summary_list$target_width_raw)

  message("Overlapping target and transcribed GRanges.")
  overlap_frame <-
    mergeByOverlaps(query = target_granges, subject = transcribed_granges)
  overlap_granges <- overlap_frame$target_granges
  # Adjust start and end to the minimally overlapping regions.
  constrained_granges <-
    GenomicRanges::GRanges(
      seqnames = seqnames(x = overlap_frame$target_granges),
      ranges = IRanges(
        start = pmax(
          start(x = overlap_frame$target_granges),
          start(x = overlap_frame$transcribed_granges)
        ),
        end = pmin(
          end(x = overlap_frame$target_granges),
          end(x = overlap_frame$transcribed_granges)
        )
      )
    )
  rm(overlap_frame)
} else {
  # If target regions are not available, all transcribed GenomicRanges::GRanges
  # count.
  summary_list$target_path <- NA
  message("Not importing target range annotation.")
  target_granges <- transcribed_granges
  summary_list$target_width_raw <-
    summary_list$transcribed_width
  summary_list$target_number_raw <- 0L
  constrained_granges <- transcribed_granges
}

summary_list$target_number_constrained <-
  length(x = constrained_granges)
message("Number of transcribed target GRanges: ",
        summary_list$target_number_constrained)

summary_list$target_width_constrained <-
  sum(GenomicRanges::width(x = constrained_granges))
message(
  "Cumulative width of transcribed target GRanges: ",
  summary_list$target_width_constrained
)

# Read callable loci ------------------------------------------------------


# Import the callable loci BED file produced by the GATK CallableLoci analysis.
summary_list$callable_loci_path <-
  argument_list$callable_loci_path
# Store the sample name in the summary frame.
summary_list$sample_name <-
  gsub(
    pattern = "variant_calling_diagnose_sample_(.*?)_callable_loci.bed",
    replacement = "\\1",
    x = summary_list$callable_loci_path
  )
message("Processing sample name: ", summary_list$sample_name)
callable_granges <-
  rtracklayer::import(con = summary_list$callable_loci_path)
non_callable_granges <-
  callable_granges[callable_granges$name != "CALLABLE",]
non_callable_granges$name <-
  as.factor(x = non_callable_granges$name)
summary_list$non_callable_number_raw <-
  length(x = non_callable_granges)
message("Number of non-callable raw GRanges: ",
        summary_list$non_callable_number_raw)
summary_list$non_callable_width_raw <-
  sum(GenomicRanges::width(x = non_callable_granges))
message("Cumulative width of non-callable raw GRanges: ",
        summary_list$non_callable_width_raw)

# Overlap constrained GRanges ---------------------------------------------


# To get accurate non-callable statistics with regards to the target regions
# that are transcribed, non-callable GenomicRanges::GRanges need overlapping
# with the constrained GenomicRanges::GRanges.
overlap_frame <-
  mergeByOverlaps(query = non_callable_granges, constrained_granges)
# Constrain the GenomicRanges::GRanges to the minimum overlap.
overlap_granges <- GenomicRanges::GRanges(
  seqnames = seqnames(x = overlap_frame$non_callable_granges),
  ranges = IRanges(start = pmax(
    start(x = overlap_frame$non_callable_granges),
    start(x = overlap_frame$constrained_granges)
  ),
  end = pmin(
    end(x = overlap_frame$non_callable_granges),
    end(x = overlap_frame$constrained_granges)
  )),
  mapping_status = overlap_frame$name
)
summary_list[["non_callable_number_constrained.TOTAL"]] <-
  length(x = overlap_granges)
message("Number of non-callable constrained GRanges: ",
        summary_list[["non_callable_number_constrained.TOTAL"]])
summary_list[["non_callable_width_constrained.TOTAL"]] <-
  sum(GenomicRanges::width(x = overlap_granges))
message("Cumulative width of non-callable constrained GRanges: ",
        summary_list[["non_callable_width_constrained.TOTAL"]])
# summary_list[["non_callable_constrained_fraction.TOTAL"]] <-
#   summary_list[["non_callable_width_constrained.TOTAL"]] / summary_list$target_width_constrained

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
  # Count the number of entries for each mapping status level.
  aggregate_frame <-
    as.data.frame(x = table(S4Vectors::mcols(x = overlap_granges)$mapping_status))
  # Assign the result levels (rows) as summary frame columns.
  for (j in seq_len(length.out = nrow(x = aggregate_frame))) {
    summary_list[[paste("non_callable_number_constrained",
                        aggregate_frame[j, 1L, drop = TRUE],
                        sep = ".")]] <-
      aggregate_frame[j, 2L, drop = TRUE]
  }
  # Sum the widths of entries for each mapping_status level.
  aggregate_frame <-
    aggregate.data.frame(
      x = data.frame(width = GenomicRanges::width(x = overlap_granges)),
      by = list(mapping_status = S4Vectors::mcols(x = overlap_granges)$mapping_status),
      FUN = "sum"
    )
  # Assign the result levels (rows) as summary frame columns.
  for (j in seq_len(length.out = nrow(x = aggregate_frame))) {
    summary_list[[paste("non_callable_width_constrained", aggregate_frame[j, 1L, drop = TRUE], sep = ".")]] <-
      aggregate_frame[j, 2L, drop = TRUE]
  }
  rm(j, aggregate_frame)
}
rm(overlap_granges, overlap_frame)

# Annotate non-callable GRanges with target region names ------------------


# Annotate the table of non-callable GenomicRanges::GRanges with target region
# names if available.
diagnose_granges <- NULL
if (!is.null(x = summary_list$target_path)) {
  # If the target GenomicRanges::GRanges are available, merge by overlap with
  # the non-callable GenomicRanges::GRanges into a new S4Vectors::DataFrame.
  overlap_frame <-
    mergeByOverlaps(query = target_granges, subject = non_callable_granges)
  # Extract the "non_callable_granges".
  diagnose_granges <- overlap_frame[, c("non_callable_granges")]
  # Annotate with the target_name from the target GenomicRanges::GRanges object.
  S4Vectors::mcols(x = diagnose_granges)$target_name <-
    overlap_frame$name
  # Rename the "name" column of the S4Vectors::mcols() DataFrame into "mapping_status".
  colnames(x = S4Vectors::mcols(x = diagnose_granges))[colnames(x = S4Vectors::mcols(x = diagnose_granges)) == "name"] <-
    "mapping_status"
  rm(overlap_frame)
} else {
  # Diagnose GenomicRanges::GRanges without target annotation.
  diagnose_granges <- non_callable_granges
}
# To annotate non-callable regions, merge by overlap with the exon
# GenomicRanges::GRanges.
overlap_frame <-
  mergeByOverlaps(query = diagnose_granges, subject = exon_granges)
# The mergeByOverlaps() function returns a S4Vectors::DataFrame of query
# (diagnose GenomicRanges::GRanges) and subject (exon GenomicRanges::GRanges)
# variables, as well as all S4Vectors::mcols() variables that were present in
# either GenomicRanges::GRanges object. To remove these redundant
# S4Vectors::mcols() variables, keep only the GenomicRanges::GRanges objects
# themselves.
utils::write.table(
  x = overlap_frame[, c("diagnose_granges", "exon_granges"), drop = FALSE],
  file = paste(
    "variant_calling_diagnose_sample",
    summary_list$sample_name,
    "non_callable_loci.tsv",
    sep = "_"
  ),
  quote = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

# Group GRanges by non-callable regions -----------------------------------


# Since the mergeByOverlaps() function provides the cartesian product of
# diagnosis and exon GenomicRanges::GRanges, the table can be rather long and
# unwieldy. Annotate the diagnosis GenomicRanges::GRanges with exon
# GenomicRanges::GRanges meta information before grouping by diagnosis
# GenomicRanges::GRanges so that there is a single observation for each
# problematic region. Gene and transcript identifiers and names are turned into
# comma-separated lists.
message("Annotate the diagnostic GRanges.")
overlap_diagnose_granges <- overlap_frame$diagnose_granges
overlap_diagnose_frame <-
  S4Vectors::mcols(x = overlap_diagnose_granges)
overlap_diagnose_frame$gene_id <-
  S4Vectors::mcols(x = overlap_frame$exon_granges)$gene_id
overlap_diagnose_frame$gene_name <-
  S4Vectors::mcols(x = overlap_frame$exon_granges)$gene_name
overlap_diagnose_frame$transcript_id <-
  S4Vectors::mcols(x = overlap_frame$exon_granges)$transcript_id
overlap_diagnose_frame$transcript_name <-
  S4Vectors::mcols(x = overlap_frame$exon_granges)$transcript_name
overlap_diagnose_frame$exon_id <-
  S4Vectors::mcols(x = overlap_frame$exon_granges)$exon_id
S4Vectors::mcols(x = overlap_diagnose_granges) <-
  overlap_diagnose_frame
# This returns a Grouping object (CompressedManyToOneGrouping) of the IRanges
# package, specifying which groups contain which indices to the original object.
message("Group annotated diagnose GRanges by region.")
overlap_diagnose_grouping <-
  methods::as(object = overlap_diagnose_granges, "Grouping")

grouped_granges <- unlist(x = GenomicRanges::GRangesList(lapply(
  X = overlap_diagnose_grouping,
  FUN = function(x) {
    sub_granges <- overlap_diagnose_granges[x]
    sub_mcols <- S4Vectors::mcols(x = sub_granges)
    selected_range <- sub_granges[1L]

    S4Vectors::mcols(x = selected_range) <- S4Vectors::DataFrame(
      "mapping_status" = sub_mcols$mapping_status[1L],

      "gene_ids" = paste(unique(x = sort(x = sub_mcols$gene_id)), collapse = ","),

      "gene_names" = paste(unique(x = sort(x = sub_mcols$gene_name)), collapse = ","),

      "transcript_ids" = paste(unique(x = sort(x = sub_mcols$transcript_id)), collapse = ","),

      "transcript_names" = paste(unique(x = sort(
        x = sub_mcols$transcript_name
      )), collapse = ","),

      "exon_ids" = paste(unique(x = sort(x = sub_mcols$exon_id)), collapse = ",")
    )
    return(selected_range)
  }
)))

utils::write.table(
  x = grouped_granges,
  file = paste(
    "variant_calling_diagnose_sample",
    summary_list$sample_name,
    "non_callable_regions.tsv",
    sep = "_"
  ),
  quote = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

rm(
  grouped_granges,
  overlap_diagnose_grouping,
  overlap_diagnose_frame,
  overlap_diagnose_granges,
  overlap_frame,
  diagnose_granges,
  callable_granges,
  non_callable_granges
)

readr::write_tsv(
  x = tibble::as_tibble(x = summary_list),
  path = paste(
    "variant_calling_diagnose_sample",
    summary_list$sample_name,
    "non_callable_summary.tsv",
    sep = "_"
  )
)

rm(
  i,
  exon_granges,
  transcribed_granges,
  target_granges,
  constrained_granges,
  summary_list,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

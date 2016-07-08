#! /usr/bin/env Rscript
#
# BSF R script to refine the GATK Callable Loci analyis.
# The coverage assessment loads (Ensembl) exon information, collates (or projects)
# all (overlapping) exons into transcribed regions on the genome, overlaps those with
# the target regions to get transcribed target regions and finally overlaps those with
# non-callable loci to get the minimal set of problematic regions.
# Each problematic region is annotated with teh target name and the exon, transcript and
# gene information. A summary data frame of metrics collected along the procedure is also
# written to disk.
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
      opt_str = c("--exons"),
      dest = "exon_path",
      help = "File path to the gene transcript an exon annotation GTF",
      type = "character"
    ),
    make_option(
      opt_str = c("--callable-loci"),
      dest = "callable_loci_path",
      help = "File path to the GATK Callable Loci BED",
      type = "character"
    ),
    make_option(
      opt_str = c("--targets"),
      dest = "target_path",
      help = "File path to the enrichment targets BED",
      type = "character"
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

suppressPackageStartupMessages(expr = library(package = "rtracklayer"))

# Keep overall statistics in a summary data frame.
summary_frame <- data.frame(stringsAsFactors = FALSE)
i <- 1

# Import Ensembl gene, transcript and exon annotation as GRanges vector object,
# where only exon components are relevant for this analysis.
summary_frame[i, "exon_path"] <- argument_list$exon_path
message(paste0("Importing exon ranges: ", summary_frame[i, "exon_path"]))
exon_ranges <-
  import(con = summary_frame[i, "exon_path"],
         genome = "hs37d5",
         feature.type = "exon")
summary_frame[i, "exon_number"] <- length(x = exon_ranges)
message(paste0("Number of exon ranges: ", summary_frame[i, "exon_number"]))
summary_frame[i, "exon_width"] <- sum(width(x = exon_ranges))
message(paste0("Cumulative width of exon ranges: ", summary_frame[i, "exon_width"]))

# Reduce the non-redundant Ensembl exons to their footprint on the genome to get
# transcribed regions. The revmap column contains the mapping to original
# exon_ranges components.
message("Reducing exon ranges to transcribed ranges.")
transcribed_ranges <-
  reduce(
    x = exon_ranges,
    drop.empty.ranges = TRUE,
    with.revmap = TRUE,
    ignore.strand = TRUE
  )
summary_frame[i, "transcribed_number"] <-
  length(x = transcribed_ranges)
message(paste0("Number of transcribed ranges: ", summary_frame[i, "transcribed_number"]))
summary_frame[i, "transcribed_width"] <-
  sum(width(x = transcribed_ranges))
message(paste0("Cumulative width of transcribed ranges: ", summary_frame[i, "transcribed_width"]))

# Read the file of targeted regions, if available. Although the GATK Callable
# Loci analysis is generally only run on these target regions, this file
# provides the target (probe) names of the enrichment design.
target_ranges <- NULL
constrained_ranges <- NULL
if (!is.null(x = argument_list$target_path)) {
  summary_frame[i, "target_path"] <- argument_list$target_path
  message(paste0("Importing target range annotation: ", summary_frame[i, "target_path"]))
  target_ranges <- import(con = summary_frame[i, "target_path"])
  # TODO: Could the genome version be specified to get the sequence lengths from the
  # sequence dictionary (*.dict) file.
  summary_frame[i, "target_number_raw"] <-
    length(x = target_ranges)
  message(paste0("Number of target ranges: ", summary_frame[i, "target_number_raw"]))
  summary_frame[i, "target_width_raw"] <-
    sum(width(x = target_ranges))
  message(paste0("Cumulative width of target ranges: ", summary_frame[i, "target_width_raw"]))
  
  message("Overlapping target and transcribed ranges.")
  overlap_frame <-
    mergeByOverlaps(query = target_ranges, subject = transcribed_ranges)
  overlap_ranges <- overlap_frame[, "target_ranges"]
  # Adjust start and end to the minimally overlapping regions.
  constrained_ranges <-
    GRanges(
      seqnames = seqnames(x = overlap_frame$target_ranges),
      ranges = IRanges(
        start = pmax(
          start(x = overlap_frame$target_ranges),
          start(x = overlap_frame$transcribed_ranges)
        ),
        end = pmin(
          end(x = overlap_frame$target_ranges),
          end(x = overlap_frame$transcribed_ranges)
        )
      )
    )
  rm(overlap_frame)
} else {
  # If target regions are not avaiable, all transcribed GRanges count.
  summary_frame[i, "target_path"] <- NA
  message("Not importing target range annotation.")
  target_ranges <- transcribed_ranges
  summary_frame[i, "target_width_raw"] <-
    summary_frame[i, "transcribed_width"]
  summary_frame[i, "target_number_raw"] <- 0
  constrained_ranges <- transcribed_ranges
}
summary_frame[i, "target_number_constrained"] <-
  length(x = constrained_ranges)
message(paste0("Number of transcribed target ranges: ", summary_frame[i, "target_number_constrained"]))
summary_frame[i, "target_width_constrained"] <-
  sum(width(x = constrained_ranges))
message(paste0("Cumulative width of transcribed target ranges: ", summary_frame[i, "target_width_constrained"]))

# Import the callable loci BED file produced by the GATK CallableLoci analysis.
summary_frame[i, "callable_loci_path"] <-
  argument_list$callable_loci_path
# Store the sample name in the summary frame.
summary_frame[i, "sample_name"] <-
  gsub(pattern = "variant_calling_diagnose_sample_(.*?)_callable_loci.bed",
       replacement = "\\1",
       x = summary_frame[i, "callable_loci_path"])
message(paste0("Processing sample name: ", summary_frame[i, "sample_name"]))
callable_ranges <-
  import(con = summary_frame[i, "callable_loci_path"])
non_callable_ranges <-
  callable_ranges[callable_ranges$name != "CALLABLE",]
non_callable_ranges$name <-
  as.factor(x = non_callable_ranges$name)
summary_frame[i, "non_callable_number_raw"] <-
  length(x = non_callable_ranges)
message(paste0("Number of non-callable raw ranges: ", summary_frame[i, "non_callable_number_raw"]))
summary_frame[i, "non_callable_width_raw"] <-
  sum(width(x = non_callable_ranges))
message(paste0("Cumulative width of non-callable raw ranges: ", summary_frame[i, "non_callable_width_raw"]))

# To get accurate non-callable statistics with regards to the target regions
# that are transcribed, non-callable GRanges need overlapping with the
# constrained GRanges.
overlap_frame <-
  mergeByOverlaps(query = non_callable_ranges, constrained_ranges)
# Constrain the GRanges to the minimum overlap.
overlap_ranges <- GRanges(
  seqnames = seqnames(x = overlap_frame$non_callable_ranges),
  ranges = IRanges(start = pmax(
    start(x = overlap_frame$non_callable_ranges),
    start(x = overlap_frame$constrained_ranges)
  ),
  end = pmin(
    end(x = overlap_frame$non_callable_ranges),
    end(x = overlap_frame$constrained_ranges)
  )),
  mapping_status = overlap_frame$name
)
summary_frame[i, "non_callable_number_constrained.TOTAL"] <-
  length(x = overlap_ranges)
message(paste0(
  "Number of non-callable constrained ranges: ",
  summary_frame[i, "non_callable_number_constrained.TOTAL"]))
summary_frame[i, "non_callable_width_constrained.TOTAL"] <-
  sum(width(x = overlap_ranges))
message(paste0(
  "Cumulative width of non-callable constrained ranges: ",
  summary_frame[i, "non_callable_width_constrained.TOTAL"]
))
# summary_frame[i, "non_callable_constrained_fraction.TOTAL"] <-
#   summary_frame[i, "non_callable_width_constrained.TOTAL"] / summary_frame[i, "target_width_constrained"]

# Summarise also separately by mapping status.
# Populate the summary frame with columns for each mapping status level,
# regardless of whether it is associated with data or not.
for (level in levels(x = mcols(x = overlap_ranges)$mapping_status)) {
  summary_frame[i, paste("non_callable_number_constrained", level, sep = ".")] <- 0
  summary_frame[i, paste("non_callable_width_constrained", level, sep = ".")] <- 0
}
rm(level)
# Count the number of entries for each mapping status level.
aggregate_frame <-
  as.data.frame(x = table(mcols(x = overlap_ranges)$mapping_status))
# Assign the result levels (rows) as summary frame columns.
for (j in 1:nrow(x = aggregate_frame)) {
  summary_frame[i, paste("non_callable_number_constrained", aggregate_frame[j, 1], sep = ".")] <-
    aggregate_frame[j, 2]
}
# Sum the widths of entries for each mapping_status level.
aggregate_frame <-
  aggregate.data.frame(
    x = data.frame(width = width(x = overlap_ranges)),
    by = list(mapping_status = mcols(x = overlap_ranges)$mapping_status),
    FUN = "sum"
  )
# Assign the result levels (rows) as summary frame columns.
for (j in 1:nrow(x = aggregate_frame)) {
  summary_frame[i, paste("non_callable_width_constrained", aggregate_frame[j, 1], sep = ".")] <-
    aggregate_frame[j, 2]
}
rm(j, aggregate_frame, overlap_ranges, overlap_frame)

# Annotate the table of non-callable GRanges with target region names if available.
diagnose_ranges <- NULL
if (!is.null(x = summary_frame[i, "target_path"])) {
  # If the target GRanges are available, merge by overlap with the non-callable GRanges into a new DataFrame.
  overlap_frame <-
    mergeByOverlaps(query = target_ranges, subject = non_callable_ranges)
  # Extract the "non_callable_ranges".
  diagnose_ranges <- overlap_frame[, c("non_callable_ranges")]
  # Annotate with the target_name from the target GRanges object.
  mcols(x = diagnose_ranges)$target_name <-
    overlap_frame[, "name"]
  # Rename the "name" column of the mcols() DataFrame into "mapping_status".
  colnames(x = mcols(x = diagnose_ranges))[colnames(x = mcols(x = diagnose_ranges)) == "name"] <-
    "mapping_status"
  rm(overlap_frame)
} else {
  # Diagnose GRanges without target annotation.
  diagnose_ranges <- non_callable_ranges
}
# To annotate non-callable regions, merge by overlap with the exon GRanges.
overlap_frame <-
  mergeByOverlaps(query = diagnose_ranges, subject = exon_ranges)
# Remove redundant columns, as the GRanges objects contain some of the columns internally.
overlap_frame <-
  overlap_frame[, c("diagnose_ranges", "exon_ranges")]
write.table(
  x = overlap_frame,
  file = paste(
    "variant_calling_diagnose_sample",
    summary_frame[i, "sample_name"],
    "non_callable_loci.tsv",
    sep = "_"
  ),
  quote = TRUE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)
rm(overlap_frame,
   diagnose_ranges,
   callable_ranges,
   non_callable_ranges)

write.table(
  x = summary_frame,
  file = paste(
    "variant_calling_diagnose_sample",
    summary_frame[i, "sample_name"],
    "non_callable_summary.tsv",
    sep = "_"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

rm(
  i,
  exon_ranges,
  transcribed_ranges,
  target_ranges,
  constrained_ranges,
  summary_frame,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

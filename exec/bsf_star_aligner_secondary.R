#!/usr/bin/env Rscript
#
# BSF R script to count STAR aligner secondary alignments per genome tile.
#
# The size of the tiles is configurable, resuts are returned as
# RangedSummarizedExperiment objects saved to a
# star_aligner_secondary_ranged_summarized_experiment.R file.
#
#
# Copyright 2013 - 2019 Michael K. Schuster
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
      opt_str = c("--bsgenome"),
      default = "BSgenome.Hsapiens.UCSC.hg38",
      dest = "bsgenome",
      help = "Bioconductor BSgenome package [BSgenome.Hsapiens.UCSC.hg38]",
      type = "character"
    ),
    make_option(
      opt_str = c("--directory"),
      default = ".",
      dest = "directory",
      help = "Directory of star_aligner_merged_*.bam files [.]",
      type = "character"
    ),
    make_option(
      opt_str = c("--tile-width"),
      default = 1000000L,
      dest = "tile_width",
      help = "Tile width [1,000,000]",
      type = "integer"
    ),
    make_option(
      opt_str = c("--threads"),
      default = 1L,
      dest = "threads",
      help = "Number of parallel processing threads [1]",
      type = "integer"
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

# Check the input.
if (is.null(x = argument_list$directory)) {
  stop("Missing --directory option")
}

suppressPackageStartupMessages(expr = library(package = "BiocParallel"))  # for BiocParallel::register()
suppressPackageStartupMessages(expr = library(package = "GenomicAlignments"))  # for GenomicAlignments::summarizeOverlaps()
suppressPackageStartupMessages(expr = library(package = "GenomicRanges"))  # for GenomicRanges::tileGenome()
suppressPackageStartupMessages(expr = library(package = "Rsamtools"))  # for Rsamtools::BamFileList()
suppressPackageStartupMessages(expr = library(package = argument_list$bsgenome, character.only = TRUE))  # BSgenome package

# Set the number of parallel threads in the MulticoreParam instance.

BiocParallel::register(BPPARAM = MulticoreParam(workers = argument_list$threads))

# Create genome tiles -----------------------------------------------------


message("Creating genome tiles")
granges_list <-
  tileGenome(seqlengths = seqinfo(x = get(x = argument_list$bsgenome)),
             tilewidth = argument_list$tile_width)

# Create a BamFileList ----------------------------------------------------


message("Creating BamFileList object")
sample_frame <-
  DataFrame(
    bam_path = base::list.files(
      path = argument_list$directory,
      pattern = "^star_aligner_merge_.*\\.bam$",
      full.names = TRUE
    )
  )
if (nrow(x = sample_frame) == 0L) {
  stop("Could not find any files in the supplied directory.")
}
sample_frame$bai_path <-
  gsub(pattern = "\\.bam$",
       replacement = ".bai",
       x = sample_frame$bam_path)
sample_frame$sample <-
  gsub(pattern = "^star_aligner_merge_(.*)\\.bam$",
       replacement = "\\1",
       x = sample_frame$bam_path)
bam_file_list <- BamFileList(
  file = as.character(x = sample_frame$bam_path),
  index = as.character(x = sample_frame$bai_path),
  yieldSize = 2000000L
)
names(x = bam_file_list) <-
  as.character(x = sample_frame$sample)

# Create a RangedSummarizedExperiment -------------------------------------


message("Creating a RangedSummarizedExperiment object")
ranged_summarized_experiment <-
  summarizeOverlaps(
    features = granges_list,
    reads = bam_file_list,
    mode = "Union",
    # For the purpose of assigning secondary reads to tiles the strand is not essential.
    ignore.strand = TRUE,
    # Count reads overlapping features.
    inter.feature = FALSE,
    # include only reads that represent secondary alignments and pass the vendor quality filter.
    param = ScanBamParam(
      flag = scanBamFlag(
        isSecondaryAlignment = TRUE,
        isNotPassingQualityControls = FALSE
      )
    )
  )
colData(x = ranged_summarized_experiment) <-
  sample_frame

# Save the RangedSummarizedExperiment -------------------------------------


message("Saving the RangedSummarizedExperiment object")
file_path <- "star_aligner_secondary_ranged_summarized_experiment.R"
save(ranged_summarized_experiment, file = file_path)

rm(
  ranged_summarized_experiment,
  file_path,
  bam_file_list,
  sample_frame,
  granges_list,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

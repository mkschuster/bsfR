#! /usr/bin/env Rscript
#
# BSF R script to aggregate and plot Picard tools IlluminaSamDemux and
# Illumina2bam tools BamIndexDecoder metrics files.
#
#
# Copyright 2013 - 2018 Michael K. Schuster
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
  parse_args(object = OptionParser(
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
        opt_str = c("--directory-path"),
        dest = "directory_path",
        help = "Directory path of EXPERIMENT_LANE_metrics.tsv files",
        type = "character"
      ),
      make_option(
        opt_str = c("--file-path"),
        dest = "file_path",
        help = "File path of a EXPERIMENT_LANE_metrics.tsv file",
        type = "character"
      ),
      make_option(
        opt_str = c("--plot-factor"),
        default = 0.5,
        dest = "plot_factor",
        help = "Plot width increase per 24 samples [0.5]",
        type = "numeric"
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

if (is.null(x = argument_list$directory_path) &
    is.null(x = argument_list$file_path)) {
  stop("Missing --directory-path or --file-path option")
}

suppressPackageStartupMessages(expr = library(package = "ggplot2"))
suppressPackageStartupMessages(expr = library(package = "methods"))
suppressPackageStartupMessages(expr = library(package = "reshape2"))

# Save plots in the following formats.

graphics_formats <- c("pdf", "png")

# General Picard metrics class.

setClass(
  Class = "PicardMetrics",
  slots = c(
    program = "character",
    cmd_line = "character",
    start_date = "character",
    metrics = "data.frame"
  )
)

base::invisible(x = setMethod(
  f = "initialize",
  signature = "PicardMetrics",
  definition = function(.Object,
                        program = NULL,
                        cmd_line = NULL,
                        start_date = NULL,
                        metrics = NULL ,
                        ...) {
    .Object <- callNextMethod(
      .Object,
      program = program,
      cmd_line = cmd_line,
      start_date = start_date,
      metrics = metrics,
      ...
    )
  }
))

# Illumina2bam tools BamIndexDecoder BarcodeMetrics class.

setClass(Class = "PicardBarcodeMetric",
         contains = "PicardMetrics")

base::invisible(x = setGeneric(
  name = "process_metrics",
  def = function(.Object, name, ...) {
    # Melt the data frame to have series for total and pass_filter matches only.
    data_frame <- .Object@metrics
    data_frame[data_frame$BARCODE_NAME == "", "BARCODE_NAME"] <-
      "Unassigned"
    # Increase the plot width per 24 samples.
    plot_width <-
      argument_list$plot_width + (ceiling(x = nrow(x = data_frame) / 24L) - 1L) * argument_list$plot_width * argument_list$plot_factor
    
    # Fractions of barcode matches ----------------------------------------
    molten_frame <-
      reshape2::melt(
        data = data_frame,
        id.vars = "BARCODE_NAME",
        measure.vars = c("PCT_MATCHES", "PF_PCT_MATCHES"),
        variable.name = "Type",
        value.name = "Fraction"
      )
    # Plot the molten data frame.
    ggplot_object <- ggplot2::ggplot(data = molten_frame)
    ggplot_object <-
      ggplot_object + ggplot2::geom_point(mapping = aes(x = BARCODE_NAME, y = Fraction, colour = Type))
    ggplot_object <-
      ggplot_object + ggplot2::ggtitle(label = paste("De-Multiplexing Statistics", name))
    ggplot_object <-
      ggplot_object + ggplot2::theme(axis.text.x = element_text(
        size = 8.0,
        hjust = 1.0,
        vjust = 0.5,
        angle = 90.0
      ))
    ggplot_object <-
      ggplot_object + ggplot2::xlab(label = "Sample Name")
    ggplot_object <- ggplot_object + ggplot2::ylim(0.0, NA)
    for (graphics_format in graphics_formats) {
      ggplot2::ggsave(
        filename = paste(
          paste(name, "metrics_fraction", sep = "_"),
          graphics_format,
          sep = "."
        ),
        plot = ggplot_object,
        width = plot_width,
        height = argument_list$plot_height,
        limitsize = FALSE
      )
    }
    rm(graphics_format, ggplot_object, molten_frame)
    
    # Numbers of barcode matches ------------------------------------------
    molten_frame <- reshape2::melt(
      data = data_frame,
      id.vars = "BARCODE_NAME",
      measure.vars = c(
        "READS",
        "PF_READS",
        "PERFECT_MATCHES",
        "PF_PERFECT_MATCHES",
        "ONE_MISMATCH_MATCHES",
        "PF_ONE_MISMATCH_MATCHES"
      ),
      variable.name = "Type",
      value.name = "Number"
    )
    # Plot the molten data frame.
    ggplot_object <- ggplot2::ggplot(data = molten_frame)
    ggplot_object <-
      ggplot_object + ggplot2::geom_point(mapping = aes(x = BARCODE_NAME, y = Number, colour = Type))
    ggplot_object <-
      ggplot_object + ggplot2::ggtitle(label = paste("De-Multiplexing Statistics", name))
    ggplot_object <-
      ggplot_object + ggplot2::theme(axis.text.x = element_text(
        size = 8.0,
        hjust = 1.0,
        vjust = 0.5,
        angle = 90.0
      ))
    ggplot_object <-
      ggplot_object + ggplot2::xlab(label = "Sample Name")
    ggplot_object <- ggplot_object + ggplot2::ylim(0, NA)
    for (graphics_format in graphics_formats) {
      ggplot2::ggsave(
        filename = paste(
          paste(name, "metrics_number", sep = "_"),
          graphics_format,
          sep = "."
        ),
        plot = ggplot_object,
        width = plot_width,
        height = argument_list$plot_height,
        limitsize = FALSE
      )
    }
    rm(graphics_format, ggplot_object, molten_frame)
    rm(plot_width, data_frame)
  }
))

# Picard AlignmentsSummary class.

setClass(Class = "PicardAlignmentSummary",
         contains = "PicardMetrics")

base::invisible(x = setMethod(
  f = "process_metrics",
  signature = c(".Object" = "PicardAlignmentSummary", "name" = "character"),
  definition = function(.Object, name, ...) {
    print(x = paste("The plot method is not implemented."))
  }
))

# Correlate Picard metrics classes with R classes.

metrics_class <- list(
  "picard.analysis.AlignmentSummaryMetrics" = "PicardAlignmentSummary",
  "picard.illumina.ExtractIlluminaBarcodes$BarcodeMetric" = "PicardBarcodeMetric",
  "uk.ac.sanger.npg.picard.IndexDecoder$BarcodeMetric" = "PicardBarcodeMetric"
)


#' Process a Picard-style metrics file.
#'
#' @param file_path
#'
#' @return PicardMetrics sub-class
#' @export
#'
#' @examples
process_metrics_file <- function(file_path) {
  if (is.null(x = file_path)) {
    stop("Missing file_path argument.")
  }
  
  # picard.illumina.ExtractIlluminaBarcodes$BarcodeMetric -----------------
  
  ## htsjdk.samtools.metrics.StringHeader
  # IlluminaSamDemux INPUT=/data/groups/lab_bsf/sequences/BSF_0491_HWM3JBBXX/BSF_0491_HWM3JBBXX_1.bam OUTPUT_DIR=BSF_0491_HWM3JBBXX_1_samples OUTPUT_PREFIX=BSF_0491_HWM3JBBXX_1 LIBRARY_PARAMS=BSF_0491_HWM3JBBXX_1_library.tsv METRICS_FILE=BSF_0491_HWM3JBBXX_1_metrics.tsv READ_STRUCTURE=51T8B TMP_DIR=[illumina_demultiplex_sam_lane_BSF_0491_HWM3JBBXX_1_temporary] COMPRESSION_LEVEL=9 CREATE_MD5_FILE=true    OUTPUT_FORMAT=bam BARCODE_TAG_NAME=BC BARCODE_QUALITY_TAG_NAME=QT MAX_MISMATCHES=1 MIN_MISMATCH_DELTA=1 MAX_NO_CALLS=2 MINIMUM_BASE_QUALITY=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
  ## htsjdk.samtools.metrics.StringHeader
  # Started on: Sat Jul 14 15:59:38 CEST 2018
  
  ## METRICS CLASS        picard.illumina.ExtractIlluminaBarcodes$BarcodeMetric
  # BARCODE
  # BARCODE_WITHOUT_DELIMITER
  # BARCODE_NAME
  # LIBRARY_NAME
  # READS
  # PF_READS
  # PERFECT_MATCHES
  # PF_PERFECT_MATCHES
  # ONE_MISMATCH_MATCHES
  # PF_ONE_MISMATCH_MATCHES
  # PCT_MATCHES
  # RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT
  # PF_PCT_MATCHES
  # PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT
  # PF_NORMALIZED_MATCHES
  
  # BamIndexDecoder metrics -----------------------------------------------
  
  ## net.sf.picard.metrics.StringHeader
  # uk.ac.sanger.npg.picard.BamIndexDecoder INPUT=/data/groups/lab_bsf/sequences/BSF_0490_HWMMFBBXX/BSF_0490_HWMMFBBXX_1.bam OUTPUT_DIR=BSF_0490_HWMMFBBXX_1_samples OUTPUT_PREFIX=BSF_0490_HWMMFBBXX_1 OUTPUT_FORMAT=bam BARCODE_FILE=BSF_0490_HWMMFBBXX_1_barcode.tsv METRICS_FILE=BSF_0490_HWMMFBBXX_1_metrics.tsv TMP_DIR=[bam_index_decoder_lane_BSF_0490_HWMMFBBXX_1_temporary] VERBOSITY=WARNING COMPRESSION_LEVEL=9 CREATE_MD5_FILE=true    BARCODE_TAG_NAME=BC BARCODE_QUALITY_TAG_NAME=QT MAX_MISMATCHES=1 MIN_MISMATCH_DELTA=1 MAX_NO_CALLS=2 CONVERT_LOW_QUALITY_TO_NO_CALL=false MAX_LOW_QUALITY_TO_CONVERT=15 USE_HASH_DEMULTIPLEXING=false BARCODE_STARTBASE_INDEX=0 QUIET=false VALIDATION_STRINGENCY=STRICT MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false
  ## net.sf.picard.metrics.StringHeader
  # Started on: Fri Jul 13 16:44:46 CEST 2018
  
  ## METRICS CLASS        uk.ac.sanger.npg.picard.IndexDecoder$BarcodeMetric
  # BARCODE
  # BARCODE_NAME
  # LIBRARY_NAME
  # SAMPLE_NAME
  # DESCRIPTION
  # READS
  # PF_READS
  # PERFECT_MATCHES
  # PF_PERFECT_MATCHES
  # ONE_MISMATCH_MATCHES
  # PF_ONE_MISMATCH_MATCHES
  # PCT_MATCHES
  # RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT
  # PF_PCT_MATCHES
  # PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT
  # PF_NORMALIZED_MATCHES
  
  # CollectAlignmentSummaryMetrics ----------------------------------------
  
  ## htsjdk.samtools.metrics.StringHeader
  # picard.analysis.CollectAlignmentSummaryMetrics METRIC_ACCUMULATION_LEVEL=[READ_GROUP, ALL_READS] INPUT=/data/groups/lab_bsf/sequences/BSF_0149_C6H48ANXX/BSF_0149_C6H48ANXX_8.bam OUTPUT=BSF_0149_C6H48ANXX_8_metrics.tsv TMP_DIR=[bam_index_decoder_BSF_0149_C6H48ANXX_8_temporary] VERBOSITY=WARNING    MAX_INSERT_SIZE=100000 ADAPTER_SEQUENCE=[AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG] IS_BISULFITE_SEQUENCED=false ASSUME_SORTED=true STOP_AFTER=0 QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false
  ## htsjdk.samtools.metrics.StringHeader
  # Started on: Mon Jun 22 15:03:59 CEST 2015
  
  ## METRICS CLASS        picard.analysis.AlignmentSummaryMetrics
  # CATEGORY
  # TOTAL_READS
  # PF_READS
  # PCT_PF_READS
  # PF_NOISE_READS
  # PF_READS_ALIGNED
  # PCT_PF_READS_ALIGNED
  # PF_ALIGNED_BASES
  # PF_HQ_ALIGNED_READS
  # PF_HQ_ALIGNED_BASES
  # PF_HQ_ALIGNED_Q20_BASES
  # PF_HQ_MEDIAN_MISMATCHES
  # PF_MISMATCH_RATE
  # PF_HQ_ERROR_RATE
  # PF_INDEL_RATE
  # MEAN_READ_LENGTH
  # READS_ALIGNED_IN_PAIRS
  # PCT_READS_ALIGNED_IN_PAIRS
  # BAD_CYCLES
  # STRAND_BALANCE
  # PCT_CHIMERAS
  # PCT_ADAPTER
  # SAMPLE
  # LIBRARY
  # READ_GROUP
  
  # TODO: It would be more scalable to parse line by line, rather than fixed lines (i.e. 2, 4, 6, ...).
  
  metrics_lines <- readLines(con = file_path)
  
  # The second line contains the program class and the command line.
  
  matching_text_2 <- regmatches(
    x = metrics_lines[2],
    m = regexec(pattern = "#[[:space:]]+([^[:space:]]+)[[:space:]]+(.*)$",
                text = metrics_lines[2])
  )
  
  # The fourth line contains the start date
  matching_text_4 <- regmatches(
    x = metrics_lines[4],
    m = regexec(pattern = "^# Started on:[[:space:]]+(.*)$",
                text = metrics_lines[4])
  )
  
  # The sixth line constains the metrics class.
  matching_text_6 <- regmatches(
    x = metrics_lines[6],
    m = regexec(pattern = "^## METRICS CLASS[[:space:]]+(.*)$",
                text = metrics_lines[6])
  )
  
  metrics_frame <- read.table(
    file = file_path,
    header = TRUE,
    sep = "\t",
    skip = 6,
    stringsAsFactors = FALSE
  )
  
  # Create a new Picard Report object depending on the metrics class.
  
  picard_report <- new(
    Class = metrics_class[[matching_text_6[[1]][2]]],
    program = matching_text_2[[1]][2],
    cmd_line = matching_text_2[[1]][3],
    start_date = matching_text_4[[1]][2],
    metrics = metrics_frame
  )
  
  rm(metrics_lines,
     matching_text_2,
     matching_text_4,
     matching_text_6)
  
  return(picard_report)
}

metrics_files <- NULL
if (is.null(x = argument_list$directory_path)) {
  argument_list$directory_path = dirname(path = argument_list$file_path)
  metrics_files <- list(basename(path = argument_list$file_path))
} else {
  metrics_files <-
    dir(path = argument_list$directory_path, pattern = "_metrics.tsv")
}
metrics_reports <-
  vector(mode = "list", length = length(x = metrics_files))
names(x = metrics_reports) <-
  sub(pattern = "^(.*)_metrics.tsv$",
      replacement = "\\1",
      x = metrics_files)
for (i in seq_along(along.with = metrics_files)) {
  metrics_reports[[i]] <-
    process_metrics_file(file_path = file.path(argument_list$directory, metrics_files[i]))
  # Create a plot.
  process_metrics(metrics_reports[[i]], names(metrics_reports)[i])
}
rm(i, metrics_files)

rm(
  metrics_reports,
  process_metrics_file,
  metrics_class,
  graphics_formats,
  argument_list
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

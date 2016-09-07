#! /usr/bin/env Rscript
#
# BSF R script to aggregate Trimmomatic output files.
#
#
# Copyright 2013 - 2016 Michael K. Schuster
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

suppressPackageStartupMessages(expr = library(package = "ggplot2"))
suppressPackageStartupMessages(expr = library(package = "optparse"))
suppressPackageStartupMessages(expr = library(package = "reshape2"))

#' Process trimmomatic STDOUT files and return a named character vector.
#' Names are: "line", "input", "both", "first", "second", "dropped"
#'
#' @param file_path: File path
#' @return: Named character vector

process_stdout <- function(file_path) {
  if (is.null(x = file_path)) {
    stop("Missing file_path argument.")
  }
  trimmomatic_lines <- readLines(con = file_path)
  # Match the relevant line, which looks like:
  # Input Read Pairs: 7999402 Both Surviving: 6629652 (82.88%) Forward Only Surviving: 1271736 (15.90%) Reverse Only Surviving: 18573 (0.23%) Dropped: 79441 (0.99%)
  trimmomatic_matches <-
    regexec(pattern = "^Input Read Pairs: ([[:digit:]]+) Both Surviving: ([[:digit:]]+) .*Forward Only Surviving: ([[:digit:]]+).*Reverse Only Surviving: ([[:digit:]]+).*Dropped: ([[:digit:]]+) ",
            text = trimmomatic_lines)
  # Get the sub-strings corresponding to the matches.
  trimmomatic_strings <-
    regmatches(x = trimmomatic_lines, m = trimmomatic_matches)
  # Select only lines with matches from the list via a logical vector.
  trimmomatic_filter <-
    as.logical(x = lapply(
      X = trimmomatic_matches,
      FUN = function(x) {
        x[1] > -1
      }
    ))
  trimmomatic_filtered <- trimmomatic_strings[trimmomatic_filter]
  # trimmomatic_frame = NULL
  if (length(x = trimmomatic_filtered) > 0) {
    # Annotate the list of strings with names.
    names(x = trimmomatic_filtered[[1]]) <-
      c("line", "input", "both", "first", "second", "dropped")
    # trimmomatic_frame <- data.frame(trimmomatic_filtered[[1]], stringsAsFactors = FALSE)
  }
  rm(trimmomatic_lines,
     trimmomatic_matches,
     trimmomatic_strings,
     trimmomatic_filter)
  return(trimmomatic_filtered[[1]])
}

#' Process trimlog files.
#'
#' Trimmomatic trimlog files are tab-separated value (TSV) files with the
#' following fields:
#' 1: read name
#' 2: surviving sequence length
#' 3: location of the first surviving base, or the amount trimmed from the start
#' 4: location of the last surviving base in the original read
#' 5: amount trimmed from the end
#' @param file_path: File path
#' @param number: Maximum number of (paired) reads to read
#' @return:

process_trimlog <- function(file_path, number = -1L) {
  file_prefix <- sub(
    pattern = "_trim_log\\.tsv(\\.gz)?",
    replacement = "",
    x = basename(path = file_path)
  )
  
  initialise_read_frame <- function(length, read_number = 1L) {
    return(
      data.frame(
        "position" = seq.int(from = 0L, to = length - 1L),
        "read" = as.factor(x = paste0('R', read_number)),
        "surviving" = integer(length = length),
        "frequency_5" = integer(length = length),
        "frequency_3" = integer(length = length),
        stringsAsFactors = FALSE
      )
    )
  }
  
  update_read_frame <- function(read_frame, trimlog_list) {
    # Increment the surviving position.
    read_frame$surviving[as.integer(x = trimlog_list[[1]][2]) + 1L] <-
      read_frame$surviving[as.integer(x = trimlog_list[[1]][2]) + 1L] + 1L
    # Increment the frequency_5 position.
    read_frame$frequency_5[as.integer(x = trimlog_list[[1]][3]) + 1L] <-
      read_frame$frequency_5[as.integer(x = trimlog_list[[1]][3]) + 1L] + 1L
    # Increment the frequency_3 position.
    read_frame$frequency_3[as.integer(x = trimlog_list[[1]][4]) + 1L] <-
      read_frame$frequency_3[as.integer(x = trimlog_list[[1]][4]) + 1L] + 1L
    return(read_frame)
  }
  
  # The data frames need to be intialised with the correct length, which is only
  # available from reading the trim log file.
  counter_1 <- 0L
  counter_2 <- 0L
  read_frame_1 <- NULL
  read_frame_2 <- NULL
  
  trimlog_connection <- file(description = file_path, open = "rt")
  # Since trim log files may start with completely trimmed reads that allow no conclusion
  # about read lengths, the file needs searching for meaningful coordinates first.
  while (TRUE) {
    trimlog_line <- readLines(con = trimlog_connection, n = 1L)
    if (length(x = trimlog_line) == 0L) {
      break()
    }
    trimlog_list <-
      strsplit(x = trimlog_line, split = "[[:space:]]")
    # The read length is the sum of the end trimming position and the amount trimmed from the end.
    # The data frame length is one longer to store the end position.
    frame_length <-
      as.integer(x = trimlog_list[[1]][4]) + as.integer(x = trimlog_list[[1]][5]) + 1L
    
    # Check for read 1.
    # TODO: The endsWith() function would be available in R base with version 3.3.0.
    if (substr(
      x = trimlog_list[[1]][1],
      start = nchar(x = trimlog_list[[1]][1]) - 1,
      stop = nchar(x = trimlog_list[[1]][1])
    ) == "/1") {
      counter_1 <- counter_1 + 1L  # Increment the counter for read 1.
      # Initialise only with meaningful coordinates.
      if (is.null(x = read_frame_1) && frame_length > 1) {
        read_frame_1 <-
          initialise_read_frame(length = frame_length, read_number = 1L)
      }
    }
    
    # Check for read 2.
    # TODO: The endsWith() function would be available in R base with version 3.3.0.
    if (substr(
      x = trimlog_list[[1]][1],
      start = nchar(x = trimlog_list[[1]][1]) - 1,
      stop = nchar(x = trimlog_list[[1]][1])
    ) == "/2") {
      counter_2 <- counter_2 + 1L  # Increment the counter for read 2.
      # Initialise only with meaningful coordinates.
      if (is.null(x = read_frame_2) && frame_length > 1) {
        read_frame_2 <-
          initialise_read_frame(length = frame_length, read_number = 2L)
      }
    }
    
    # No decision before the first read has been seen twice ...
    if (counter_1 > 2L) {
      # Break out, if ...
      if (!is.null(x = read_frame_1)) {
        # ... the data frame for read 1 is initialised ...
        if (counter_2 > 0L) {
          # ... and after two frist reads a second read seems to exist ...
          if (!is.null(x = read_frame_2)) {
            # ... and the data frame for read 2 is initialised.
            break()
          }
        } else {
          # ... and after two first reads no second read seems to exist.
          break()
        }
      }
    }
  }
  
  # Re-postion the connection to the start, reset the read counters and re-read the entire file.
  seek(con = trimlog_connection, where = 0L)
  counter_1 <- 0L
  counter_2 <- 0L
  while (TRUE) {
    trimlog_line <- readLines(con = trimlog_connection, n = 1L)
    if (length(x = trimlog_line) == 0L) {
      break()
    }
    trimlog_list <-
      strsplit(x = trimlog_line, split = "[[:space:]]")
    
    # Check if first or second read.
    # TODO: The endsWith() function would be available in R base with version 3.3.0.
    # if (endsWith(x = trimlog_list[[1]][1], "/1")) {
    if (substr(
      x = trimlog_list[[1]][1],
      start = nchar(x = trimlog_list[[1]][1]) - 1,
      stop = nchar(x = trimlog_list[[1]][1])
    ) == "/1") {
      counter_1 <- counter_1 + 1L  # Increment the counter for read 1.
      read_frame_1 <-
        update_read_frame(read_frame = read_frame_1, trimlog_list = trimlog_list)
    }
    
    # TODO: The endsWith() function would be available in R base with version 3.3.0.
    # if (endsWith(x = trimlog_list[[1]][1], "/2")) {
    if (substr(
      x = trimlog_list[[1]][1],
      start = nchar(x = trimlog_list[[1]][1]) - 1,
      stop = nchar(x = trimlog_list[[1]][1])
    ) == "/2") {
      counter_2 <- counter_2 + 1L  # Increment the counter for read 2.
      read_frame_2 <-
        update_read_frame(read_frame = read_frame_2, trimlog_list = trimlog_list)
    }
    rm(trimlog_list, trimlog_line)
    
    # if ((counter_1 %% 10000L) == 0L) {
    #   print(x = counter_1)
    # }
    
    # Break after reaching the maximum nuber of reads to process.
    if (number > 0L) {
      if (counter_2 > 0L) {
        if (counter_2 >= number) {
          # If read 2 exists, break after reaching its counter and thus completing the pair.
          break()
        }
      } else {
        if (counter_1 >= number) {
          # Otherwise, simply break after reaching the first read's counter.
          break()
        }
      }
    }
  }
  close(con = trimlog_connection)
  rm(trimlog_connection,
     update_read_frame,
     initialise_read_frame)
  
  complete_read_frame <- function(read_frame, counter) {
    # Calculate the 5-prime coverage, which is the
    # cumulative sum of 5-prime frequencies.
    read_frame$coverage_5 <-
      cumsum(x = read_frame$frequency_5)
    # Calculate the 3-prime coverage, which is the
    # total count minus the
    # cumulative sum of 3-prime frequencies.
    read_frame$coverage_3 <-
      counter -
      cumsum(x = read_frame$frequency_3)
    # Calculate the total coverage, which is the
    # 5-prime cumulative sum minus the
    # 3-prime cumulative sum.
    read_frame$coverage <-
      cumsum(x = read_frame$frequency_5) -
      cumsum(x = read_frame$frequency_3)
    # Also set the counter, which is a bit of a waste,
    # as it applies to all positions equally.
    read_frame$counter <-
      counter
    return(read_frame)
  }
  
  trimmomatic_frame <-
    complete_read_frame(read_frame = read_frame_1, counter = counter_1)
  if (!is.null(x = read_frame_2)) {
    trimmomatic_frame <-
      rbind.data.frame(
        trimmomatic_frame,
        complete_read_frame(read_frame = read_frame_2, counter = counter_2)
      )
  }
  rm(read_frame_1,
     read_frame_2,
     complete_read_frame,
     counter_1,
     counter_2)
  
  # Write the summary frame to disk.
  write.table(
    x = trimmomatic_frame,
    file = paste(file_prefix, "summary.tsv", sep = "_"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  
  # Melt the data frame on measure variables "coverage_5", "coverage_3" and "coverage"
  # and plot faceted by column on the "read" factor.
  molten_frame <- melt(
    data = trimmomatic_frame,
    id.vars = c("position", "read"),
    measure.vars = c("coverage_5", "coverage_3", "coverage"),
    variable.name = "type",
    value.name = "reads"
  )
  ggplot_object <- ggplot(data = molten_frame)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Coverage Summary")
  ggplot_object <- ggplot_object + facet_grid(facets = read ~ .)
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = position, y = reads, colour = type))
  ggsave(filename = paste(file_prefix, "coverage.png", sep = "_"),
         plot = ggplot_object)
  
  # Melt the data frame on measure variables "frequency_5" and "frequency_3"
  # and plot faceted by column on the "read" factor.
  molten_frame <- melt(
    data = trimmomatic_frame,
    id.vars = c("position", "read"),
    measure.vars = c("frequency_5", "frequency_3"),
    variable.name = "type",
    value.name = "reads"
  )
  ggplot_object <- ggplot(data = molten_frame)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Frequency Summary")
  ggplot_object <- ggplot_object + facet_grid(facets = read ~ .)
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = position, y = reads, colour = type))
  ggsave(filename = paste(file_prefix, "frequency.png", sep = "_"),
         plot = ggplot_object)
  
  # Melt the data frame on measure variables "surviving"
  # and plot faceted by column on the "read" factor.
  # molten_frame <- melt(
  #   data = trimmomatic_frame,
  #   id.vars = c("position", "read"),
  #   measure.vars = c("surviving"),
  #   variable.name = "type",
  #   value.name = "reads"
  # )
  # Rather than melting with a single measure variable, just select from the data frame.
  molten_frame <-
    trimmomatic_frame[, c("position", "read", "surviving")]
  ggplot_object <- ggplot(data = molten_frame)
  ggplot_object <-
    ggplot_object + ggtitle(label = "Surviving Sequence")
  ggplot_object <- ggplot_object + facet_grid(facets = read ~ .)
  # ggplot_object <-
  #   ggplot_object + geom_point(mapping = aes(x = position, y = reads, colour = type))
  ggplot_object <-
    ggplot_object + geom_point(mapping = aes(x = position, y = surviving)) + ylab(label = "reads")
  ggsave(filename = paste(file_prefix, "surviving.png", sep = "_"),
         plot = ggplot_object)
  
  rm(ggplot_object, molten_frame, trimmomatic_frame, file_prefix)
  return()
}

#' Process Trimmomatic summary data frame files, produced by this script by
#' reading Trimmomatic trimlog files.
#'
#' @param directory_path: Directory path
#' @return:

process_summary <- function(directory_path) {
  file_list <- list.files(
    path = directory_path,
    pattern = '^trimmomatic_.*_summary.tsv',
    full.names = FALSE,
    recursive = FALSE
  )
  summary_frame <- data.frame(
    "position" = integer(),
    "surviving" = integer(),
    "frequency_5" = integer(),
    "frequency_3" = integer(),
    "coverage_5" = integer(),
    "coverage_3" = integer(),
    "coverage" = integer(),
    "sample" = character(),
    stringsAsFactors = TRUE
  )
  for (file_path in file_list) {
    sample_frame <-
      read.table(
        file = file_path,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = TRUE
      )
    sample_frame$sample <-
      sub(pattern = '^trimmomatic_(.*)_summary.tsv',
          replacement = "\\1",
          x = file_path)
    summary_frame <- rbind.data.frame(summary_frame, sample_frame)
    rm(sample_frame)
  }
  write.table(
    x = summary_frame,
    file = "trimmomatic_summary.tsv",
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  rm(file_list, summary_frame)
  return()
}

# Specify the desired options in a list.
# By default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action = "store_true", default = FALSE,
# help = "Show this help message and exit").

# Get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults.

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
        opt_str = c("--directory"),
        help = "Trimmomatic directory",
        type = "character"
      ),
      make_option(
        opt_str = c("--file_path"),
        help = "Trimmomatic trimlog file path",
        type = "character"
      ),
      make_option(
        opt_str = c("--number"),
        help = "Maximum number of Trimmomatic trimlog lines, or -1 (default) for unlimited",
        default = -1L,
        type = "integer"
      ),
      make_option(
        opt_str = c("--summary_path"),
        help = "Trimmomatic summary data frame directory path",
        type = "character"
      )
    )
  ))

if (!is.null(x = argument_list$file_path)) {
  return_value <-
    process_trimlog(file_path = argument_list$file_path, number = argument_list$number)
  rm(return_value)
}

if (!is.null(x = argument_list$summary_path)) {
  return_value <-
    process_summary(directory_path = argument_list$summary_path)
  rm(return_value)
}

if (!is.null(x = argument_list$directory)) {
  trimmomatic_paths <-
    list.files(
      path = argument_list$directory,
      pattern = "\\.bsf_run_trimmomatic_.*\\.err",
      all.files = TRUE,
      full.names = TRUE
    )
  # trimmomatic_reports <- vector(mode = "list", length = length(x = trimmomatic_paths))
  
  trimmomatic_list <-
    lapply(X = trimmomatic_paths, FUN = process_stdout)
  # print(x = paste("trimmomatic_list class:", class(x = trimmomatic_list)))
  # print(x = paste("trimmomatic_list length:", length(x = trimmomatic_list)))
  # print(x = data.frame(matrix(data = unlist(x = trimmomatic_list)), byrow = TRUE), stringsAsFactors=FALSE)
  # print(x = lapply(X = trimmomatic_list, FUN = class))
  # print(x = trimmomatic_list)
  
  # trimmomatic_frame <- rbind.data.frame(trimmomatic_list)
  # trimmomatic_frame <- do.call(rbind.data.frame, trimmomatic_list)
  # trimmomatic_frame <- data.frame(
  #   line = trimmomatic_list$line,
  #   input = trimmomatic_list$input,
  #   both = trimmomatic_list$both,
  #   first = trimmomatic_list$first,
  #   second = trimmomatic_list$second,
  #   dropped = trimmomatic_list$dropped)
  # trimmomatic_frame <- rbind.data.frame(trimmomatic_list)
  # print(x = trimmomatic_frame)
  # print(x = paste("trimmomatic_frame names:", names(x = trimmomatic_frame)))
  # print(head(x = trimmomatic_frame))
  # print(row.names(x = trimmomatic_frame))
  
  # for (trimmomatic_file in trimmomatic_files) {
  # trimmomatic_lines <- readLines(con = trimmomatic_file)
  
  # The relevant line is:
  # Input Read Pairs: 7999402 Both Surviving: 6629652 (82.88%) Forward Only Surviving: 1271736 (15.90%) Reverse Only Surviving: 18573 (0.23%) Dropped: 79441 (0.99%)
  # trimmomatic_matches <- regexec(pattern = "^Input Read Pairs: ([[:digit:]]+) Both Surviving: ([[:digit:]]+) .*Forward Only Surviving: ([[:digit:]]+).*Reverse Only Surviving: ([[:digit:]]+).*Dropped: ([[:digit:]]+) ",
  #         text = trimmomatic_lines)
  # trimmomatic_values <- regmatches(x = trimmomatic_lines, m = matches)
  # The second line contains the program class and the command line.
  # print(x = trimmomatic_values)
  # }
  rm(trimmomatic_paths)
}

rm(argument_list,
   process_summary,
   process_trimlog,
   process_stdout)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

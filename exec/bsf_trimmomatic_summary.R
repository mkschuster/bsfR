#!/usr/bin/env Rscript
#
# BSF R script to aggregate Trimmomatic output files.
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
        opt_str = c("--directory-path"),
        dest = "directory_path",
        help = "Trimmomatic STDERR directory path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--file-path"),
        dest = "file_path",
        help = "Trimmomatic trimlog file path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--number"),
        help = "Maximum number of Trimmomatic trimlog lines, or -1 for unlimited [-1]",
        default = -1L,
        type = "integer"
      ),
      optparse::make_option(
        opt_str = c("--chunk-size"),
        default = 10000L,
        dest = "chunk_size",
        help = "Number of lines to process per chunk [10000]",
        type = "integer"
      ),
      optparse::make_option(
        opt_str = c("--summary-path"),
        dest = "summary_path",
        help = "Trimmomatic summary data frame directory path",
        type = "character"
      )
    )
  ))

suppressPackageStartupMessages(expr = library(package = "tidyverse"))

#' Process Trimmomatic STDOUT files and return a named character vector.
#' Names are: "line", "input", "both", "first", "second", "dropped"
#'
#' @noRd
#' @param file_path A \code{character} scalar with the file path.
#' @return A named \code{character} vector of parsed components.
process_stdout <- function(file_path) {
  if (is.null(x = file_path)) {
    stop("Missing file_path argument.")
  }
  trimmomatic_lines <- base::readLines(con = file_path)
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
        return(x[1] > -1)
      }
    ))
  trimmomatic_filtered <- trimmomatic_strings[trimmomatic_filter]
  if (length(x = trimmomatic_filtered) > 0) {
    # Annotate the list of strings with names.
    names(x = trimmomatic_filtered[[1]]) <-
      c("line", "input", "both", "first", "second", "dropped")
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
#' following variables:
#'
#' 1: read name
#' 2: surviving sequence length
#' 3: location of the first surviving base, or the amount trimmed from the start
#' 4: location of the last surviving base in the original read
#' 5: amount trimmed from the end
#'
#' @noRd
#' @param file_path A \code{character} scalar with the file path.
#' @param number An \code{integer} scalar indicating the maximum number of
#'   (paired) reads to read.
#' @return
process_trimlog <- function(file_path, number = -1L) {
  file_prefix <- sub(
    pattern = "_trim_log\\.tsv(\\.gz)?",
    replacement = "",
    x = base::basename(path = file_path)
  )

  initialise_summary_tibble <- function(length, read_number = 1L) {
    return(
      tibble::tibble(
        "position" = seq.int(from = 0L, to = length - 1L),
        "read" = paste0('R', read_number),
        "surviving" = integer(length = length),
        "frequency_5" = integer(length = length),
        "frequency_3" = integer(length = length)
      )
    )
  }

  update_summary_tibble <-
    function(summary_tibble, trimlog_character) {
      # Subset the character matrix into columns 2 (surviving), 3 (frequency_5)
      # and 4 (frequency_3). Increment all matrix elements by 1L to convert
      # 0-based sequence indices to 1-based R vector indices.
      trimlog_integer <-
        matrix(data = as.integer(x = trimlog_character[, c(2L, 3L, 4L)]),
               nrow = nrow(x = trimlog_character)) + 1L

      for (i in seq_len(length.out = nrow(x = trimlog_integer))) {
        # Increment the surviving position.
        summary_tibble$surviving[trimlog_integer[i, 1L]] <-
          summary_tibble$surviving[trimlog_integer[i, 1L]] + 1L
        # Increment the frequency_5 position.
        summary_tibble$frequency_5[trimlog_integer[i, 2L]] <-
          summary_tibble$frequency_5[trimlog_integer[i, 2L]] + 1L
        # Increment the frequency_3 position.
        summary_tibble$frequency_3[trimlog_integer[i, 3L]] <-
          summary_tibble$frequency_3[trimlog_integer[i, 3L]] + 1L
      }
      rm(i, trimlog_integer)

      return(summary_tibble)
    }

  # The data frames need to be initialised with the correct length, which is only
  # available from reading the trim log file.
  counter_1 <- 0L
  counter_2 <- 0L
  summary_tibble_1 <- NULL
  summary_tibble_2 <- NULL

  trimlog_connection <-
    base::file(description = file_path, open = "rt")
  # Since trim log files may start with completely trimmed reads that allow no
  # conclusion about read lengths, the file needs searching for meaningful
  # coordinates first.
  while (TRUE) {
    trimlog_line <- base::readLines(con = trimlog_connection, n = 1L)
    if (length(x = trimlog_line) == 0L) {
      break()
    }
    trimlog_character <-
      stringr::str_split(string = trimlog_line,
                         pattern = stringr::fixed(pattern = " "))[[1L]]
    # The read length is the sum of the end trimming position and the amount
    # trimmed from the end. The tibble length is one longer to store the end
    # position.
    tibble_length <-
      as.integer(x = trimlog_character[4L]) + as.integer(x = trimlog_character[5L]) + 1L

    # Check for read 1.
    if (endsWith(x = trimlog_character[1L], "/1")) {
      counter_1 <- counter_1 + 1L  # Increment the counter for read 1.
      # Initialise only with meaningful coordinates.
      if (is.null(x = summary_tibble_1) && tibble_length > 1) {
        summary_tibble_1 <-
          initialise_summary_tibble(length = tibble_length, read_number = 1L)
      }
    } else {
      # Check for read 2.
      if (endsWith(x = trimlog_character[1L], "/2")) {
        counter_2 <- counter_2 + 1L  # Increment the counter for read 2.
        # Initialise only with meaningful coordinates.
        if (is.null(x = summary_tibble_2) && tibble_length > 1) {
          summary_tibble_2 <-
            initialise_summary_tibble(length = tibble_length, read_number = 2L)
        }
      } else {
        # This must be Trimmomatic data in SE mode, which lacks /1 and /2 suffices.
        counter_1 <-
          counter_1 + 1L  # Increment the counter for read 1.
        # Initialise only with meaningful coordinates.
        if (is.null(x = summary_tibble_1) && tibble_length > 1) {
          summary_tibble_1 <-
            initialise_summary_tibble(length = tibble_length, read_number = 1L)
        }
      }
    }

    # No decision before the first read has been seen twice ...
    if (counter_1 > 2L) {
      # Break out, if ...
      if (!is.null(x = summary_tibble_1)) {
        # ... the data frame for read 1 is initialised ...
        if (counter_2 > 0L) {
          # ... and after two frist reads a second read seems to exist ...
          if (!is.null(x = summary_tibble_2)) {
            # ... and the data frame for read 2 is initialised.
            break()
          }
        } else {
          # ... and after two first reads no second read seems to exist.
          break()
        }
      }
    }
    rm(tibble_length, trimlog_character, trimlog_line)
  }

  # Re-postion the connection to the start, reset the read counters and re-read the entire file.
  base::seek(con = trimlog_connection, where = 0L)
  counter_1 <- 0L
  counter_2 <- 0L
  while (TRUE) {
    trimlog_lines <-
      base::readLines(con = trimlog_connection, n = argument_list$chunk_size)
    if (length(x = trimlog_lines) == 0L) {
      break()
    }
    trimlog_character <-
      stringr::str_split_fixed(
        string = trimlog_lines,
        pattern = stringr::fixed(pattern = " "),
        n = 5L
      )

    reads_1 <- base::endsWith(x = trimlog_character[, 1L], "/1")
    reads_2 <- base::endsWith(x = trimlog_character[, 1L], "/2")
    # Check for read 1.
    if (any(reads_1)) {
      counter_1 <- counter_1 + length(x = which(x = reads_1))
      summary_tibble_1 <-
        update_summary_tibble(summary_tibble = summary_tibble_1,
                              trimlog_character = trimlog_character[reads_1, ])
    }
    # Check for read 2.
    if (any(reads_2)) {
      counter_2 <- counter_2 + length(x = which(x = reads_2))
      summary_tibble_2 <-
        update_summary_tibble(summary_tibble = summary_tibble_2,
                              trimlog_character = trimlog_character[reads_2, ])
    }
    # This must be Trimmomatic data in SE mode, which lacks /1 and /2 suffices.
    if (!any(reads_1, reads_2)) {
      counter_1 <- counter_1 + nrow(x = trimlog_character)
      summary_tibble_1 <-
        update_summary_tibble(summary_tibble = summary_tibble_1,
                              trimlog_character = trimlog_character)
    }
    rm(reads_2, reads_1, trimlog_character, trimlog_lines)

    # FIXME: For debugging only!
    print(x = paste("Counter 1:", counter_1))

    # Break after reaching the maximum number of reads to process.
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
  base::close(con = trimlog_connection)
  rm(trimlog_connection,
     update_summary_tibble,
     initialise_summary_tibble)

  finalise_summary_tibble <- function(summary_tibble, counter) {
    return(
      dplyr::mutate(
        .data = summary_tibble,
        # Calculate the 5-prime coverage, which is the cumulative sum of 5-prime
        # frequencies.
        "coverage_5" = cumsum(x = .data$frequency_5),
        # Calculate the 3-prime coverage, which is the total count minus the
        # cumulative sum of 3-prime frequencies.
        "coverage_3" = counter - cumsum(x = .data$frequency_3),
        # Calculate the total coverage, which is the 5-prime cumulative sum
        # minus the 3-prime cumulative sum.
        "coverage" = cumsum(x = .data$frequency_5) - cumsum(x = .data$frequency_3),
        # Also set the counter, which is a bit of a waste, as it applies to all
        # positions equally.
        "counter" = counter
      )
    )
  }

  summary_tibble <-
    finalise_summary_tibble(summary_tibble = summary_tibble_1, counter = counter_1)

  if (!is.null(x = summary_tibble_2)) {
    summary_tibble <-
      dplyr::bind_rows(
        summary_tibble,
        finalise_summary_tibble(summary_tibble = summary_tibble_2, counter = counter_2)
      )
  }
  rm(
    summary_tibble_1,
    summary_tibble_2,
    finalise_summary_tibble,
    counter_1,
    counter_2
  )

  # Write the summary frame to disk.
  readr::write_tsv(x = summary_tibble,
                   file = paste(paste(file_prefix, "summary", sep = "_"), "tsv", sep = "."))

  # Coverage Summary Plot -------------------------------------------------


  # Pivot the data frame on measure variables "coverage_5", "coverage_3" and
  # "coverage" and plot faceted by column on the "read" factor.
  ggplot_object <- ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = summary_tibble,
      cols = c(.data$coverage_5, .data$coverage_3, .data$coverage),
      names_to = "type",
      values_to = "reads"
    )
  )
  ggplot_object <-
    ggplot_object + ggplot2::facet_grid(cols = ggplot2::vars(read))
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$position,
        y = .data$reads,
        colour = .data$type
      ),
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Position",
      y = "Reads",
      colour = "Type",
      title = "Coverage Summary"
    )
  ggplot2::ggsave(filename = paste(paste(file_prefix, "coverage", sep = "_"), "png", sep = "."),
                  plot = ggplot_object)

  # Frequency Summary Plot ------------------------------------------------


  # Pivot the data frame on measure variables "frequency_5" and "frequency_3"
  # and plot faceted by column on the "read" factor.
  ggplot_object <- ggplot2::ggplot(
    data = tidyr::pivot_longer(
      data = summary_tibble,
      cols = c(.data$frequency_5, .data$frequency_3),
      names_to = "type",
      values_to = "reads"
    )
  )
  ggplot_object <-
    ggplot_object + ggplot2::facet_grid(cols = ggplot2::vars(read))
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$position,
        y = .data$reads,
        colour = .data$type
      ),
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Position",
      y = "Reads",
      colour = "Type",
      title = "Frequency Summary"
    )
  ggplot2::ggsave(filename = paste(paste(file_prefix, "frequency", sep = "_"), "png", sep = "."),
                  plot = ggplot_object)

  # Surviving Sequence Plot -----------------------------------------------


  # Select measure variable "surviving" and variables "position" and "read" and
  # plot faceted by column on the "read" factor.
  ggplot_object <-
    ggplot2::ggplot(
      data = dplyr::select(.data = summary_tibble, .data$position, .data$read, .data$surviving)
    )
  ggplot_object <-
    ggplot_object + ggplot2::facet_grid(cols = ggplot2::vars(read))
  ggplot_object <-
    ggplot_object + ggplot2::geom_point(
      mapping = ggplot2::aes(x = .data$position, y = .data$surviving),
      alpha = I(1 / 3)
    )
  ggplot_object <-
    ggplot_object + ggplot2::labs(x = "Position", y = "Reads", title = "Surviving Sequence")
  ggplot2::ggsave(filename = paste(paste(file_prefix, "surviving", sep = "_"), "png", sep = "."),
                  plot = ggplot_object)

  rm(ggplot_object, summary_tibble, file_prefix)
}

#' Process Trimmomatic summary data frame files, produced by this script by
#' reading Trimmomatic trimlog files.
#'
#' @param directory_path A \code{character} scalar with a directory path.
#' @return
#' @noRd
process_summary <- function(directory_path) {
  file_list <- base::list.files(
    path = directory_path,
    pattern = '^trimmomatic_.*_summary.tsv',
    full.names = FALSE,
    recursive = FALSE
  )
  summary_tibble <- NULL
  for (file_path in file_list) {
    sample_tibble <- readr::read_tsv(
      file = file_path,
      col_types = readr::cols(
        "position" = readr::col_integer(),
        "surviving" = readr::col_integer(),
        "frequency_5" = readr::col_integer(),
        "frequency_3" = readr::col_integer(),
        "coverage_5" = readr::col_integer(),
        "coverage_3" = readr::col_integer(),
        "coverage" = readr::col_integer(),
        "counter" = readr::col_integer()
        # "sample" = readr::col_character()
      )
    )
    sample_tibble$sample <-
      sub(pattern = '^trimmomatic_(.*)_summary.tsv',
          replacement = "\\1",
          x = file_path)
    summary_tibble <-
      dplyr::bind_rows(summary_tibble, sample_tibble)
    rm(sample_tibble)
  }
  readr::write_tsv(x = summary_tibble, file = "trimmomatic_summary.tsv")
  rm(file_path, file_list, summary_tibble)
}

if (!is.null(x = argument_list$file_path)) {
  process_trimlog(file_path = argument_list$file_path, number = argument_list$number)
}

if (!is.null(x = argument_list$summary_path)) {
  process_summary(directory_path = argument_list$summary_path)
}

if (!is.null(x = argument_list$directory_path)) {
  trimmomatic_paths <-
    base::list.files(
      path = argument_list$directory_path,
      pattern = "\\.bsf_run_trimmomatic_.*\\.err",
      all.files = TRUE,
      full.names = TRUE
    )
  # trimmomatic_reports <- vector(mode = "list", length = length(x = trimmomatic_paths))

  trimmomatic_list <-
    lapply(X = trimmomatic_paths, FUN = process_stdout)
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
  rm(trimmomatic_list, trimmomatic_paths)
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

print(x = sessionInfo())

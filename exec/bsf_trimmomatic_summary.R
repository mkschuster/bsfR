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
        opt_str = c("--file-path"),
        dest = "file_path",
        help = "Trimmomatic trim log file path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--number"),
        help = "Maximum number of Trimmomatic trim log lines, or -1 for unlimited [-1]",
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
        opt_str = c("--stderr-path"),
        dest = "stderr_path",
        help = "Trimmomatic STDERR directory path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--stderr-pattern-file"),
        default = "^trimmomatic_read_group_.*\\.err$",
        dest = "stderr_pattern_file",
        help = "Trimmomatic STDERR file pattern [trimmomatic_read_group_*_[0-9]+.err]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--stderr-pattern-read-group"),
        default = "^trimmomatic_read_group_(.*)_[0-9]+\\.err$",
        dest = "stderr_pattern_read_group",
        help = "Trimmomatic STDERR read group pattern [^trimmomatic_read_group_(.*)_[0-9]+\\.err$]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--summary-path"),
        dest = "summary_path",
        help = "Trimmomatic summary data frame directory path",
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

suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
suppressPackageStartupMessages(expr = library(package = "tidyverse"))

# Save plots in the following formats.
graphics_formats <- c("pdf" = "pdf", "png" = "png")
# Maximum size for the PNG device in inches.
graphics_maximum_size_png <- 100.0

#' Process Trimmomatic STDERR files and return a tibble.
#'
#' @noRd
#' @param file_path A \code{character} scalar with the file path.
#' @return A \code{tibble} of parsed components.
#' \describe{
#' \item{input}{Number of input reads.}
#' \item{both}{Number of first and second reads.}
#' \item{first}{Number of first reads.}
#' \item{second}{Number of second reads.}
#' \item{dropped}{Number of dropped reads.}
#' \item{file_path}{STDERR file path.}
#' }
process_stderr <- function(file_path) {
  if (is.null(x = file_path)) {
    stop("Missing file_path argument.")
  }
  trimmomatic_tibble <- NULL
  trimmomatic_lines <- base::readLines(con = file_path)

  # Match the relevant line, which looks like:
  # Input Read Pairs: 7999402 Both Surviving: 6629652 (82.88%) Forward Only Surviving: 1271736 (15.90%) Reverse Only Surviving: 18573 (0.23%) Dropped: 79441 (0.99%)
  trimmomatic_matches <-
    base::regexec(pattern = "^Input Read Pairs: ([[:digit:]]+) Both Surviving: ([[:digit:]]+) .*Forward Only Surviving: ([[:digit:]]+).*Reverse Only Surviving: ([[:digit:]]+).*Dropped: ([[:digit:]]+) ",
                  text = trimmomatic_lines)

  # Get a list of sub-strings corresponding to the matches.
  trimmomatic_strings <-
    base::regmatches(x = trimmomatic_lines, m = trimmomatic_matches)

  # Select only lines with matches (> -1L) from the list via a logical vector.
  trimmomatic_filtered <-
    trimmomatic_strings[purrr::map_lgl(.x = trimmomatic_matches, .f = ~ .[1L] > -1L)]

  if (length(x = trimmomatic_filtered) > 0L) {
    # Annotate the character vector with names, but skip the first element, which is the original matched line.
    trimmomatic_list <-
      as.list(x = trimmomatic_filtered[[1L]][2L:6L])
    names(x = trimmomatic_list) <-
      c("input", "both", "first", "second", "dropped")
    trimmomatic_list$file_path <- file_path
    trimmomatic_tibble <- tibble::as_tibble(x = trimmomatic_list)
    rm(trimmomatic_list)
  }

  rm(
    trimmomatic_filtered,
    trimmomatic_strings,
    trimmomatic_matches,
    trimmomatic_lines
  )

  return(trimmomatic_tibble)
}

#' Process trim log files.
#'
#' Trimmomatic trim log files are tab-separated value (TSV) files with the
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
process_trim_log <- function(file_path, number = -1L) {
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
    function(summary_tibble, trim_log_character) {
      # Subset the character matrix into columns 2 (surviving), 3 (frequency_5)
      # and 4 (frequency_3). Increment all matrix elements by 1L to convert
      # 0-based sequence indices to 1-based R vector indices.
      trim_log_integer <-
        matrix(
          data = as.integer(x = trim_log_character[, c(2L, 3L, 4L)]),
          nrow = nrow(x = trim_log_character)
        ) + 1L

      for (i in seq_len(length.out = nrow(x = trim_log_integer))) {
        # Increment the surviving position.
        summary_tibble$surviving[trim_log_integer[i, 1L]] <-
          summary_tibble$surviving[trim_log_integer[i, 1L]] + 1L
        # Increment the frequency_5 position.
        summary_tibble$frequency_5[trim_log_integer[i, 2L]] <-
          summary_tibble$frequency_5[trim_log_integer[i, 2L]] + 1L
        # Increment the frequency_3 position.
        summary_tibble$frequency_3[trim_log_integer[i, 3L]] <-
          summary_tibble$frequency_3[trim_log_integer[i, 3L]] + 1L
      }
      rm(i, trim_log_integer)

      return(summary_tibble)
    }

  # The data frames need to be initialised with the correct length, which is only
  # available from reading the trim log file.
  counter_1 <- 0L
  counter_2 <- 0L
  summary_tibble_1 <- NULL
  summary_tibble_2 <- NULL

  trim_log_connection <-
    base::file(description = file_path, open = "rt")
  # Since trim log files may start with completely trimmed reads that allow no
  # conclusion about read lengths, the file needs searching for meaningful
  # coordinates first.
  while (TRUE) {
    trim_log_line <- base::readLines(con = trim_log_connection, n = 1L)
    if (length(x = trim_log_line) == 0L) {
      break()
    }
    trim_log_character <-
      stringr::str_split(string = trim_log_line,
                         pattern = stringr::fixed(pattern = " "))[[1L]]
    # The read length is the sum of the end trimming position and the amount
    # trimmed from the end. The tibble length is one longer to store the end
    # position.
    tibble_length <-
      as.integer(x = trim_log_character[4L]) + as.integer(x = trim_log_character[5L]) + 1L

    # Check for read 1.
    if (endsWith(x = trim_log_character[1L], "/1")) {
      counter_1 <- counter_1 + 1L  # Increment the counter for read 1.
      # Initialise only with meaningful coordinates.
      if (is.null(x = summary_tibble_1) && tibble_length > 1) {
        summary_tibble_1 <-
          initialise_summary_tibble(length = tibble_length, read_number = 1L)
      }
    } else {
      # Check for read 2.
      if (endsWith(x = trim_log_character[1L], "/2")) {
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
    rm(tibble_length, trim_log_character, trim_log_line)
  }

  # Re-position the connection to the start, reset the read counters and re-read
  # the entire file.
  base::seek(con = trim_log_connection, where = 0L)
  counter_1 <- 0L
  counter_2 <- 0L
  while (TRUE) {
    trim_log_lines <-
      base::readLines(con = trim_log_connection, n = argument_list$chunk_size)
    if (length(x = trim_log_lines) == 0L) {
      break()
    }
    trim_log_character <-
      stringr::str_split_fixed(
        string = trim_log_lines,
        pattern = stringr::fixed(pattern = " "),
        n = 5L
      )

    reads_1 <- base::endsWith(x = trim_log_character[, 1L], "/1")
    reads_2 <- base::endsWith(x = trim_log_character[, 1L], "/2")
    # Check for read 1.
    if (any(reads_1)) {
      counter_1 <- counter_1 + length(x = which(x = reads_1))
      summary_tibble_1 <-
        update_summary_tibble(summary_tibble = summary_tibble_1,
                              trim_log_character = trim_log_character[reads_1, ])
    }
    # Check for read 2.
    if (any(reads_2)) {
      counter_2 <- counter_2 + length(x = which(x = reads_2))
      summary_tibble_2 <-
        update_summary_tibble(summary_tibble = summary_tibble_2,
                              trim_log_character = trim_log_character[reads_2, ])
    }
    # This must be Trimmomatic data in SE mode, which lacks /1 and /2 suffices.
    if (!any(reads_1, reads_2)) {
      counter_1 <- counter_1 + nrow(x = trim_log_character)
      summary_tibble_1 <-
        update_summary_tibble(summary_tibble = summary_tibble_1,
                              trim_log_character = trim_log_character)
    }
    rm(reads_2, reads_1, trim_log_character, trim_log_lines)

    # Break after reaching the maximum number of reads to process.
    if (number > 0L) {
      if (counter_2 > 0L) {
        if (counter_2 >= number) {
          # If read 2 exists, break after reaching its counter and thus
          # completing the pair.
          break()
        }
      } else {
        if (counter_1 >= number) {
          # Otherwise, simply break after reaching the read 1 counter.
          break()
        }
      }
    }
  }
  base::close(con = trim_log_connection)
  rm(trim_log_connection,
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

  for (graphics_format in graphics_formats) {
    ggplot2::ggsave(
      filename = paste(
        paste(file_prefix, "coverage", sep = "_"),
        graphics_format,
        sep = "."
      ),
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
  }
  rm(graphics_format, ggplot_object)

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

  for (graphics_format in graphics_formats) {
    ggplot2::ggsave(
      filename = paste(
        paste(file_prefix, "frequency", sep = "_"),
        graphics_format,
        sep = "."
      ),
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
  }
  rm(graphics_format, ggplot_object)

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

  for (graphics_format in graphics_formats) {
    ggplot2::ggsave(
      filename = paste(
        paste(file_prefix, "surviving", sep = "_"),
        graphics_format,
        sep = "."
      ),
      plot = ggplot_object,
      width = argument_list$plot_width,
      height = argument_list$plot_height
    )
  }
  rm(graphics_format, ggplot_object)

  rm(summary_tibble, file_prefix)
}

#' Process Trimmomatic summary data frame files, produced by this script by
#' reading Trimmomatic trim log files.
#'
#' @param directory_path A \code{character} scalar with a directory path.
#' @return
#' @noRd
process_summary <- function(directory_path) {
  process_read_group_tibble <- function(file_path) {
    read_group_tibble <- readr::read_tsv(
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
      )
    )

    read_group_tibble$read_group <-
      base::sub(pattern = '^trimmomatic_read_group_(.*)_summary.tsv',
                replacement = "\\1",
                x = file_path)

    return(read_group_tibble)
  }

  file_paths <- base::list.files(
    path = directory_path,
    pattern = '^trimmomatic_read_group_.*_summary.tsv',
    full.names = FALSE,
    recursive = FALSE
  )

  summary_tibble <-
    purrr::map_dfr(.x = file_paths, .f = process_read_group_tibble)

  readr::write_tsv(x = summary_tibble, file = "trimmomatic_summary.tsv")

  rm(summary_tibble, file_paths, process_read_group_tibble)
}

# Process Trimmomatic trim log files.
if (!is.null(x = argument_list$file_path)) {
  process_trim_log(file_path = argument_list$file_path, number = argument_list$number)
}

# Process Trimmomatic summary files that were created from trim log files with
# this script via the --file-path option before.
if (!is.null(x = argument_list$summary_path)) {
  process_summary(directory_path = argument_list$summary_path)
}

# Process STDERR files with Trimmomatic statistics ------------------------


if (!is.null(x = argument_list$stderr_path)) {
  summary_tibble <-
    purrr::map_dfr(
      .x = base::list.files(
        path = argument_list$stderr_path,
        pattern = argument_list$stderr_pattern_file,
        full.names = TRUE,
        recursive = TRUE
      ),
      .f = process_stderr
    )

  summary_tibble <-
    readr::type_convert(
      df = summary_tibble,
      col_types = readr::cols(
        "input" = readr::col_integer(),
        "both" = readr::col_integer(),
        "first" = readr::col_integer(),
        "second" = readr::col_integer(),
        "dropped" = readr::col_integer()
      )
    )

  summary_tibble <-
    dplyr::mutate(
      .data = summary_tibble,
      "read_group" = base::gsub(
        pattern = argument_list$stderr_pattern_read_group,
        replacement = "\\1",
        x = base::basename(path = .data$file_path)
      ),
      .before = .data$input
    )

  readr::write_tsv(x = summary_tibble, file = "trimmomatic_statistics.tsv")

  # Scale the plot width with the number of read groups, by adding a quarter of
  # the original width for each 24 read groups.
  # Because read group names are quite long, extend already for the first column.
  plot_width <-
    argument_list$plot_width + (ceiling(x = nrow(x = summary_tibble) / 24L) - 1L) * argument_list$plot_width * 0.5

  # Plot the absolute numbers of reads ------------------------------------

  ggplot_object <-
    ggplot2::ggplot(
      data = tidyr::pivot_longer(
        data = summary_tibble,
        cols = c(
          .data$input,
          .data$both,
          .data$first,
          .data$second,
          .data$dropped
        ),
        names_to = "state",
        values_to = "numbers"
      )
    )

  ggplot_object <-
    ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(
      x = .data$read_group,
      y = .data$numbers,
      colour = .data$state
    ))

  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Read Group",
      y = "Numbers",
      colour = "Status",
      title = "Trimmomatic Read Numbers by Read Group"
    )

  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = ggplot2::rel(x = 0.7),
        hjust = 0.0,
        vjust = 0.5,
        angle = 90.0
      ),
      legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7))
    )

  for (graphics_format in graphics_formats) {
    ggplot2::ggsave(
      filename = paste("trimmomatic_statistics_number", graphics_format, sep = "."),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object)

  # Plot the fractions of reads -------------------------------------------

  # Calculate fractions on the basis of "input" reads.
  summary_tibble <-
    dplyr::mutate(.data = summary_tibble, dplyr::across(
      .cols = c(.data$both, .data$first, .data$second, .data$dropped),
      .fns = ~ .x / .data$input
    ))

  summary_tibble <-
    dplyr::select(.data = summary_tibble,
                  c(
                    .data$read_group,
                    .data$both,
                    .data$first,
                    .data$second,
                    .data$dropped
                  ))

  summary_tibble <-
    tidyr::pivot_longer(
      data = summary_tibble,
      cols = c(.data$both, .data$first, .data$second, .data$dropped),
      names_to = "state",
      values_to = "fractions"
    )

  ggplot_object <- ggplot2::ggplot(data = summary_tibble)

  ggplot_object <-
    ggplot_object + ggplot2::geom_point(mapping = ggplot2::aes(
      x = .data$read_group,
      y = .data$fractions,
      colour = .data$state
    ))

  ggplot_object <-
    ggplot_object + ggplot2::labs(
      x = "Read Group",
      y = "Fractions",
      colour = "Status",
      title = "Trimmomatic Read Fractions by Read Group"
    )

  ggplot_object <-
    ggplot_object + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = ggplot2::rel(x = 0.7),
        hjust = 0.0,
        vjust = 0.5,
        angle = 90.0
      ),
      legend.text = ggplot2::element_text(size = ggplot2::rel(x = 0.7))
    )

  for (graphics_format in graphics_formats) {
    ggplot2::ggsave(
      filename = paste(
        "trimmomatic_statistics_fractions",
        graphics_format,
        sep = "."
      ),
      plot = ggplot_object,
      width = if (graphics_format == "png" &&
                  plot_width > graphics_maximum_size_png)
        graphics_maximum_size_png
      else
        plot_width,
      height = argument_list$plot_height,
      limitsize = FALSE
    )
  }
  rm(graphics_format, ggplot_object, plot_width, summary_tibble)
}

rm(
  argument_list,
  process_summary,
  process_trim_log,
  process_stderr,
  graphics_maximum_size_png,
  graphics_formats
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessioninfo::session_info())

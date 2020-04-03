#!/usr/bin/env Rscript
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --distribution=block
#SBATCH --mem=2048
#SBATCH --time=12:00:00
#SBATCH --partition=shortq
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --error=.bsf_rnaseq_deseq_summary_%j.err
#SBATCH --output=.bsf_rnaseq_deseq_summary_%j.out
#
# BSF R script to collate DESeq2 analysis summary tables.
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

suppressPackageStartupMessages(expr = library(package = "bsfR"))
suppressPackageStartupMessages(expr = library(package = "tidyverse"))

# For the moment, no other prefix-based directory gets created under the output
# directory.
output_directory <-
  file.path(argument_list$output_directory)
if (!file.exists(output_directory)) {
  dir.create(path = output_directory,
             showWarnings = TRUE,
             recursive = FALSE)
}

# Parse rnaseq_deseq_[design]_contrasts_summary.tsv files -----------------


summary_tibble <- NULL
# Get a list of all rnaseq_deseq_.*_contrasts_summary.tsv files recursively,
# extract directory names, base names and finally match design names.
design_names <-
  base::gsub(
    pattern = "^rnaseq_deseq_(.*?)$",
    replacement = "\\1",
    x = base::basename(path = base::dirname(
      path = base::list.files(
        path = argument_list$genome_directory,
        pattern = "rnaseq_deseq_.*_contrasts_summary.tsv$",
        full.names = TRUE,
        recursive = TRUE
      )
    ))
  )

for (design_name in design_names) {
  summary_tibble <-
    dplyr::bind_rows(
      summary_tibble,
      bsfR::bsfrd_read_contrast_tibble(
        genome_directory = argument_list$genome_directory,
        design_name = design_name,
        summary = TRUE,
        verbose = argument_list$verbose
      )
    )
}
rm(design_name, design_names)

# Drop the "Exclude" variable.
summary_tibble <- dplyr::select(.data = summary_tibble, -Exclude)

# Replace NA and "" values in the Numerator and Denominator with character "1".
summary_tibble <-
  dplyr::mutate_at(
    .tbl = summary_tibble,
    .vars = c("Numerator", "Denominator"),
    .funs = list( ~ replace(., is.na(x = .) | . == "", "1"))
  )

# Summarise by Numerator and Denominator ----------------------------------


# Check for unique keys.
key_tibble <-
  dplyr::select(.data = summary_tibble, Design, Numerator, Denominator)
duplicated_tibble <-
  summary_tibble[duplicated(x = key_tibble) |
                   duplicated(x = key_tibble, fromLast = TRUE),]
if (nrow(x = duplicated_tibble)) {
  print(x = "Duplicated Design, Numerator and Denominator rows:")
  print(x = duplicated_tibble)
}
rm(duplicated_tibble, key_tibble)

# Remove the "Label" variable,
# then spread "Significant" values on the "Design" key.
readr::write_tsv(
  x = tidyr::spread(
    data = dplyr::select(.data = summary_tibble, -Label),
    key = Design,
    value = Significant
  ),
  path = file.path(output_directory, "rnaseq_deseq_summary_design.tsv"),
  col_names = TRUE
)

# Summarise by Label ------------------------------------------------------


# Check for unique keys.
key_tibble <- dplyr::select(.data = summary_tibble, Design, Label)
duplicated_tibble <-
  summary_tibble[duplicated(x = key_tibble) |
                   duplicated(x = key_tibble, fromLast = TRUE),]
if (nrow(x = duplicated_tibble)) {
  print(x = "Duplicated Design and Label rows:")
  print(x = duplicated_tibble)
}
rm(duplicated_tibble, key_tibble)

# Remove the "Numerator" and "Denominator" variables,
# then spread "Significant" values on the "Design" key.
readr::write_tsv(
  x = tidyr::spread(
    data = dplyr::select(.data = summary_tibble, -Numerator, -Denominator),
    key = Design,
    value = Significant
  ),
  path = file.path(output_directory, "rnaseq_deseq_summary_labels.tsv"),
  col_names = TRUE
)

rm(summary_tibble, output_directory, argument_list)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

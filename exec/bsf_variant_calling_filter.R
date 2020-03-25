#!/usr/bin/env Rscript
#
# BSF R script to filter snpEff annotated, multi-sample variant calling format
# (VCF) files. While SNPEFF_IMPACT values "MODIFIER" and "LOW" are always
# filtered out, the minor allele frequency (MAF) threshold on the basis of the
# 1000 Genomes calculated allele frequencies (CAF) variable can be configured.
# Individual samples can be selected from the set, as can be a recurrence
# threshold for each given variant.
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
        opt_str = c("--vcf"),
        dest = "vcf_file_path",
        help = "Multi-sample VCF file path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--genome-assembly"),
        dest = "genome_assembly",
        help = "Genome assembly version (e.g. 'b37', 'hg38', 'hg19', ...)",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--set-name"),
        dest = "set_name",
        help = "Set name to automatically name output files",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--maf"),
        default = 0.01,
        dest = "maf",
        help = "Minor allele frequency (MAF) threshold [0.01]",
        type = "numeric"
      ),
      optparse::make_option(
        opt_str = c("--recurrence"),
        default = 1L,
        dest = "recurrence",
        help = "Variant recurrence threshold [1]",
        type = "integer"
      ),
      optparse::make_option(
        opt_str = c("--samples"),
        dest = "samples",
        help = "Comma-separated sample name filter [NULL]",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--chunk-size"),
        default = 10000L,
        dest = "chunk_size",
        help = "Chunk size, i.e. number of VCF lines to process at once [10000]",
        type = "integer"
      )
    )
  ))

# Check for missing options.

if (is.null(x = argument_list$vcf_file_path)) {
  stop("Missing --vcf option")
}

if (is.null(x = argument_list$genome_assembly)) {
  stop("Missing --genome-assembly option")
}

if (is.null(x = argument_list$set_name)) {
  stop("Missing --set-name option")
}

# Check that the --set-name option only contains word characters.
if (grepl(pattern = "\\W",
          x = argument_list$set_name,
          perl = TRUE)) {
  stop("Only word characters [0-9A-Z_a-z] are allowed for the --set-name option.")
}

suppressPackageStartupMessages(expr = library(package = "tidyverse"))

# Pre-process the sample names to get a character vector.
selected_sample_names <- character()
if (length(x = argument_list$samples) > 0L) {
  selected_sample_names <-
    unlist(x = stringr::str_split(
      string = argument_list$samples,
      pattern = stringr::fixed(pattern = ",")
    ))
}

# Check the recurrence option.
selected_sample_names_length <- length(x = selected_sample_names)
if (selected_sample_names_length > 0L &&
    selected_sample_names_length < argument_list$recurrence) {
  stop(
    sprintf(
      fmt = "The --recurrence option (%d) is set higher than the number of selected samples (%d).",
      selected_sample_names_length,
      argument_list$recurrence
    )
  )
}
rm(selected_sample_names_length)

# Begin with the VariantAnnotation package.
suppressPackageStartupMessages(expr = library(package = "VariantAnnotation"))

# Define filter functions for the FilterRules object required by the filterVcf() function.

#' Filter function object for FilterRules to filter out variants with
#' INFO variable SNPEFF_IMPACT equal to "MODIFIER" or "LOW".
#'
#' @param x: Fully parsed VCF object
#' @type x: VCF
#' @return: Logical vector of lines to keep
#' @rtype: logical

filter_info_SNPEFF_IMPACT <- function(x) {
  character_vector <- info(x = x)$SNPEFF_IMPACT
  logical_vector <-
    !(character_vector == "MODIFIER" | character_vector == "LOW")
  # Keep all records with NA.
  logical_vector[is.na(x = logical_vector)] <- TRUE
  return(logical_vector)
}

#' Helper function object for the filter_info_CAF() function to process
#' INFO CAF character vectors, checking that the minimum of the
#' double-converted CAF values is less than the minor allele frequency (MAF)
#' specified by the "--maf" command line option.
#'
#' @param x: Fully parsed VCF object
#' @type x: VCF
#' @return: Logical vector indicating minimum CAF less than MAF threshold
#' @rtype: logical

process_CAF_element <- function(x) {
  # If the character vector is empty, retain the record in any case.
  if (length(x = x) == 0L) {
    return(TRUE)
  }
  # TODO: Report a BioCoductor VariantAnnotation bug?
  #   Square brackets in the INFO variable are lost, if the last value of a list is missing (i.e. '.').
  # if (length(x = which(x = grepl(pattern = "^\\.$", x = x, perl = TRUE)))) { print(x = x) }
  # Replace square brackets seemingly specific to the GATK bundle CAF variable.
  x <-
    gsub(
      pattern = "[\\[\\]]",
      replacement = "",
      x = x,
      perl = TRUE
    )
  # Replace missing values ('.') with '0', which prevents the line from being filtered.
  x[x == "."] <- "0"
  return(min(as.double(x = x)) < argument_list$maf)
}

#' Filter function object for the FilterRules object to filter out variants with
#' a calcualted allele frequency (CAF) greater than the minor allele frequency (MAF)
#' specified by the "--maf" option.
#'
#' @param x: Fully parsed VCF object
#' @type x: VCF
#' @return: Logical vector of lines to keep
#' @rtype: logical

filter_info_CAF <- function(x) {
  character_list <- info(x = x)$db_snp.CAF
  logical_vector <-
    unlist(x = lapply(X = character_list, FUN = process_CAF_element))
  # Keep all records with NA.
  logical_vector[is.na(x = logical_vector)] <- TRUE
  return(logical_vector)
}

# Filer out variants, which are COMMON.
#' Filter function object for the FilterRules object to filter out variants with
#' an INFO variable COMMON equal to 1.
#'
#' @param x: Fully parsed VCF object
#' @type x: VCF
#' @return: Logical vector of lines to keep
#' @rtype: logical

filter_info_COMMON <- function(x) {
  integer_vector <- info(x = x)$db_snp.COMMON
  logical_vector <- integer_vector != 1L
  # Keep all records with NA.
  logical_vector[is.na(x = logical_vector)] <- TRUE
  return(logical_vector)
}

# Filter the VCF file by SNPEFF_IMPACT and CAF ----------------------------


# Write the filtered VCF file into the current working directory.
filtered_vcf_path <-
  sub(
    pattern = "\\.gz$",
    replacement = "",
    x = sub(
      pattern = "_annotated",
      replacement = "_filtered",
      x = base::basename(path = argument_list$vcf_file_path)
    )
  )
filtered_tbi_path <-
  paste(filtered_vcf_path, "gz", "tbi", sep = ".")

# Filter the VCF file by the SNPEFF_IMPACT and CAF variables.
if (!(file.exists(filtered_tbi_path) &&
      (file.info(filtered_tbi_path)$size > 0L))) {
  message("Filtering the annotated variants file.")
  filterVcf(
    file = TabixFile(
      file = argument_list$vcf_file_path,
      yieldSize = argument_list$chunk_size
    ),
    genome = argument_list$genome_assembly,
    destination = filtered_vcf_path,
    index = TRUE,
    filters = FilterRules(exprs = list(
      filter_info_CAF,
      # filter_info_COMMON,
      filter_info_SNPEFF_IMPACT
    ))
  )
} else {
  message("Skipping filtering the annotated variants file.")
}

# The filterVcf() method reindexes the filtered VCF file so that "gz" gets reappended.
filtered_vcf_path <- paste(filtered_vcf_path, "gz", sep = ".")

# Iterate over the filtered VCF file --------------------------------------


# Now, read the filtered VCF file in chunks and optionally select for samples.
message(
  "Iterating over filtered VCF file for samples: ",
  paste(selected_sample_names, collapse = ", ")
)

# Create the TSV file path.
selected_tsv_path <-
  sub(
    pattern = "\\.vcf\\.gz$",
    replacement = ".tsv",
    x = sub(
      pattern = "_annotated",
      # The underscore-separated list of sample names does not scale.
      # replacement = paste(c("_selected", selected_sample_names), collapse = "_"),
      replacement = paste(c("_selected", argument_list$set_name), collapse = "_"),
      x = base::basename(path = argument_list$vcf_file_path)
    )
  )
selected_tsv_first_chunk <- TRUE
filtered_vcf_file <- TabixFile(file = filtered_vcf_path,
                               yieldSize = argument_list$chunk_size)
open(con = filtered_vcf_file)
sum_records_read <- 0L
sum_records_written <- 0L
while (nrow(
  x = vcf_object <- readVcf(
    file = filtered_vcf_file,
    genome = argument_list$genome_assembly,
    param = ScanVcfParam(samples = selected_sample_names)
  )
)) {
  sum_records_read <- sum_records_read + nrow(x = vcf_object)
  message(sprintf(fmt = "Number of VCF records read: %i (%i total)", nrow(x = vcf_object), sum_records_read))
  # TODO: Get a list of samples, create a data frame and write it to disk.

  if (selected_tsv_first_chunk) {
    # Get the VCF header.
    vcf_header_object <- header(x = vcf_object)
    # Write a DataFrame containing VCF header information.
    filtered_vcf_annotation_path <-
      sub(
        pattern = "\\.vcf\\.gz$",
        replacement = ".tsv",
        x = sub(
          pattern = "_annotated",
          replacement = paste(
            c("_information_variables", argument_list$set_name),
            collapse = "_"
          ),
          x = base::basename(path = argument_list$vcf_file_path)
        )
      )
    message("Writing VCF information file: ",
            filtered_vcf_annotation_path)
    info_frame <- info(x = vcf_header_object)
    # Add another "Variable" column to properly list row names and reorder the data frame.
    info_frame[["Variable"]] <- row.names(x = info_frame)
    info_frame <-
      info_frame[, c("Variable", "Number", "Type", "Description")]
    write.table(
      x = info_frame,
      file = filtered_vcf_annotation_path,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )
    rm(filtered_vcf_annotation_path, info_frame)
    # Write a DataFrame containing VCF header genotype information.
    filtered_vcf_annotation_path <-
      sub(
        pattern = "\\.vcf\\.gz$",
        replacement = ".tsv",
        x = sub(
          pattern = "_annotated",
          replacement = paste(
            c("_genotype_variables", argument_list$set_name),
            collapse = "_"
          ),
          x = base::basename(path = argument_list$vcf_file_path)
        )
      )
    message("Writing VCF genotype file: ",
            filtered_vcf_annotation_path)
    geno_frame <- geno(x = vcf_header_object)
    # Add another "Variable" column to properly list row names and reorder the data frame.
    geno_frame[["Variable"]] <- row.names(x = geno_frame)
    geno_frame <-
      geno_frame[, c("Variable", "Number", "Type", "Description")]
    write.table(
      x = geno_frame,
      file = filtered_vcf_annotation_path,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )
    rm(filtered_vcf_annotation_path, geno_frame)
    rm(vcf_header_object)
  }

  # Firstly, retrieve the GRanges object via rowRanges() and convert into a
  # data frame, whereby several columns need adjusting.
  row_ranges_frame <-
    as.data.frame(x = rowRanges(x = vcf_object, fixed = TRUE))
  # The original VCF ID variable is only available in form of row names.
  row_ranges_frame$identifier <- row.names(x = row_ranges_frame)
  # The VCF REF variable is a DNAStringSet object
  # that needs converting into a character vector.
  row_ranges_frame$REF <- as.character(x = row_ranges_frame$REF)
  # The ALT column contains a DNAStringSetList object,
  # which individual components need converting into character vectors,
  # in which empty strings need to be replaced with '*' characters (VCF 4.2),
  # before the character vector can be collapsed into a comma-separated string.
  row_ranges_frame$ALT <-
    unlist(x = lapply(
      X = row_ranges_frame$ALT,
      FUN = function(x) {
        paste(sub(
          pattern = "^$",
          replacement = "*",
          x = as.character(x),
          perl = TRUE
        ),
        collapse = ",")
      }
    ))
  # Reorder columns, but drop "strand" and "paramRangeID".
  row_ranges_frame <-
    row_ranges_frame[c("seqnames",
                       "start",
                       "end",
                       "width",
                       "identifier",
                       "REF",
                       "ALT",
                       "QUAL",
                       "FILTER")]

  # Secondly, get the VCF INFO frame.
  info_frame <- info(x = vcf_object)
  # Convert List objects into character vector objects of comma-separated strings.
  info_frame$AC <-
    unlist(x = lapply(
      X = info_frame$AC,
      FUN = function(x) {
        paste(x, collapse = ",")
      }
    ))
  info_frame$AF <-
    unlist(x = lapply(
      X = info_frame$AF,
      FUN = function(x) {
        paste(x, collapse = ",")
      }
    ))
  info_frame$MLEAC <-
    unlist(x = lapply(
      X = info_frame$MLEAC,
      FUN = function(x) {
        paste(x, collapse = ",")
      }
    ))
  info_frame$MLEAF <-
    unlist(x = lapply(
      X = info_frame$MLEAF,
      FUN = function(x) {
        paste(x, collapse = ",")
      }
    ))
  # FIXME: Because of a bug in VariantAnnotation, the closing square bracket in the
  # CAF variable is missing, if the last value is unknown (i.e. '.'). Bring it back.
  info_frame$db_snp.CAF <-
    unlist(x = lapply(
      X = info_frame$db_snp.CAF,
      FUN = function(x) {
        paste(sub(
          pattern = "^\\.$",
          replacement = ".]",
          x = x
        ), collapse = ",")
      }
    ))
  info_frame$db_snp.CLNDBN <-
    unlist(x = lapply(
      X = info_frame$db_snp.CLNDBN,
      FUN = function(x) {
        paste(x, collapse = ",")
      }
    ))
  info_frame$db_snp.CLNDSDBID <-
    unlist(x = lapply(
      X = info_frame$db_snp.CLNDSDBID,
      FUN = function(x) {
        paste(x, collapse = ",")
      }
    ))
  info_frame$db_snp.CLNHGVS <-
    unlist(x = lapply(
      X = info_frame$db_snp.CLNHGVS,
      FUN = function(x) {
        paste(x, collapse = ",")
      }
    ))

  # Thirdly, process sample-specific genotype information and build up a new sample frame.
  # Get sample meta information as a DataFrame.
  column_data_frame <- colData(x = vcf_object)
  # Get sample-specific information as a SimpleList of matrix objects.
  genotype_list <- geno(x = vcf_object)
  # Create a new data frame for sample-specific information by expanding the
  # column information for each sample name..
  sample_frame <- NULL
  for (sample_name in row.names(x = column_data_frame)) {
    for (variable_name in names(x = genotype_list)) {
      column_name <- paste(sample_name, variable_name, sep = '.')
      # Exclude the column "SB" (Fisher strand bias),
      # which is an array of dimensions [<variant number>, <sample number>, 4L].
      if (variable_name == "SB") {
        next()
      }
      if (is.null(x = sample_frame)) {
        # If the sample_frame does not exist, create it.
        sample_frame <-
          data.frame(column_name = genotype_list[[variable_name]][, sample_name])
        # Change the "column_name" header into its real value.
        names(x = sample_frame) <- column_name
      } else {
        # If the sample frame exists, add a column to it.
        sample_frame[[column_name]] <-
          genotype_list[[variable_name]][, sample_name]
      }
      # The AD and PL variables are lists that need collapsing.
      if (variable_name %in% c("AD", "PL")) {
        sample_frame[[column_name]] <- unlist(x = lapply(
          X = sample_frame[[column_name]],
          FUN = function(x) {
            paste(x, collapse = ",")
          }
        ))
      }
    }
  }
  rm(column_name, variable_name, sample_name)

  # Calculate the recurrence of variant sites and store it in an integer vector.
  # Check the "GT" matrix line for line, if a genotype has a number from 1 to 9,
  # i.e. not ".", unphased "./.", "0/0", any other ploidy such as "0/0/0/0" or
  # phased ".:." "0:0", "0:0:0:0", etc.
  message("Identifying variant sites and calculating their recurrence.")
  info_frame[["Recurrence"]] <-
    integer(length = nrow(x = info_frame))
  for (i in seq_len(length.out = dim(genotype_list[["GT"]])[1L])) {
    info_frame[i, "Recurrence"] <-
      sum(grepl(pattern = "[1-9]", x = genotype_list[["GT"]][i, ]))
  }
  rm(i, genotype_list, column_data_frame)

  # Simply combine the row ranges, info and column data frames via cbind, as they are all in the same order.
  # message("Merging filtered VCF file.")
  combined_frame <-
    cbind(row_ranges_frame, info_frame, sample_frame)
  # Remove the three sub-frames at this stage to save memory.
  rm(row_ranges_frame, info_frame, sample_frame)
  # Select only rows which Recurrence variable is equal to or more than the recurrence threshold option.
  combined_frame <-
    combined_frame[combined_frame$Recurrence >= argument_list$recurrence,]

  sum_records_written <-
    sum_records_written + nrow(x = combined_frame)
  message(
    sprintf(fmt = "Number of VCF records written: %i (%i total)", nrow(x = combined_frame), sum_records_written)
  )
  write.table(
    x = combined_frame,
    file = selected_tsv_path,
    append = !selected_tsv_first_chunk,
    quote = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = selected_tsv_first_chunk
  )
  # For subsequent chunks, the colum names need not writing.
  selected_tsv_first_chunk <- FALSE
}
if (sum_records_read > 0L) {
  message(
    sprintf(
      fmt = "Read %i VCF records, wrote %i (%.2f %%).",
      sum_records_read,
      sum_records_written,
      sum_records_written * 100.0 / sum_records_read
    )
  )
}
rm(sum_records_read, sum_records_written)

# Close the TabixFile for the filtered VCF file.
close(con = filtered_vcf_file)

rm(
  filtered_tbi_path,
  filtered_vcf_path,
  filtered_vcf_file,
  selected_sample_names,
  vcf_object,
  selected_tsv_path,
  selected_tsv_first_chunk,
  combined_frame,
  argument_list,
  filter_info_SNPEFF_IMPACT,
  filter_info_CAF,
  filter_info_COMMON,
  process_CAF_element
)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

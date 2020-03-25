#!/usr/bin/env Rscript
#
# BSF R script to link genetic variants in a VCF file to SNPedia annotation
# in a GFF3 file on the basis of variant identifiers
# (i.e. NCBI dbSNP reference SNP rs identifiers).
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
        opt_str = c("--snpedia-gff"),
        dest = "snpedia_gff",
        help = "SNPedia GFF file path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--vcf-path"),
        dest = "vcf_path",
        help = "Input VCF file path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--tsv-path"),
        dest = "tsv_path",
        help = "Output TSV file path",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--genome-assembly"),
        dest = "genome_assembly",
        help = "Genome assembly version (e.g. 'b37', 'hg38', 'hg19', ...)",
        type = "character"
      ),
      optparse::make_option(
        opt_str = c("--chunk-size"),
        default = 10000L,
        dest = "chunk_size",
        help = "Chunk size, i.e. number of VCF lines to process at once [10000]",
        type = "integer"
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

if (is.null(x = argument_list$snpedia_gff)) {
  stop("Missing --snpedia-gff option")
}

if (is.null(x = argument_list$vcf_path)) {
  stop("Missing --vcf-path option")
}

if (is.null(x = argument_list$tsv_path)) {
  stop("Missing --tsv-path option")
}

if (is.null(x = argument_list$genome_assembly)) {
  stop("Missing --genome-assembly option")
}

suppressPackageStartupMessages(expr = library(package = "tidyverse"))
suppressPackageStartupMessages(expr = library(package = "rtracklayer"))
suppressPackageStartupMessages(expr = library(package = "VariantAnnotation"))

# Rewrite SNPedia GFF3 file -----------------------------------------------


# Unfortunately, the SNPedia file is *not* GFF3 compliant, since it
# contains unescaped tabs and semicolons in the "attributes" column 9.
#
# See the GFF3 format specification for details.
# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
#
# The following characters needRFC 3986 Percent-Encoding in column 9.
# https://tools.ietf.org/html/rfc3986#section-2.1
#
# ; semicolon (%3B)
# = equals (%3D)
# & ampersand (%26)
# , comma (%2C)

message("Rewriting the SNPedia GFF3 file")
temporary_path <-
  tempfile(pattern = "variant_calling_snpedia_", fileext = ".gff")

snpedia_connection <-
  file(description = argument_list$snpedia_gff, open = "rt")
temporary_connection <-
  file(description = temporary_path, open = "wt")

while (TRUE) {
  snpedia_character <- readLines(con = snpedia_connection, n = 100L)
  if (length(x = snpedia_character) == 0L) {
    rm(snpedia_character)
    break()
  }
  snpedia_list <- lapply(
    X = stringr::str_split(
      string = snpedia_character,
      pattern = stringr::fixed(pattern = "\t"),
      n = 9L,
      omit_empty = FALSE
    ),
    FUN = function(character_1) {
      if (length(x = character_1) == 0L) {
        # Fix empty lines that come out as NA.
        return("")
      }
      # Replace tab characters with their percent encoded value %09.
      character_1[9L] <-
        gsub(pattern = "\t",
             replacement = "%09",
             x = character_1[9L])
      # Replace semicolon characters within quoted values only
      # with their precent encoded value %3B.
      #
      # Anchor the regular expression with a look ahead on the closing quote,
      # avoiding an equal character just before the closing quote.
      # Regex ';(?=[^=]+?")' did not work with constructs that have an equal sign within quotes.
      # 'Name=rs781362878;ID=rs781362878;Note="CLNSIG=255;CLNDBN=Familial hypercholesterolemia"'
      # The following regular expression over-converts empty "Note" attributes (i.e. %3BNote="").
      character_1[9L] <-
        gsub(
          pattern = ';(?=[^"]+?[^=]")',
          replacement = "%3B",
          x = character_1[9L],
          perl = TRUE
        )
      # Revert the semicolon in over-converted empty "Note" attributes.
      character_1[9L] <-
        gsub(pattern = '%3BNote=""',
             replacement = ';Note=""',
             x = character_1[9L])
      # Replace any NA values resulting from shorter comment lines.
      character_1 <- character_1[!is.na(x = character_1)]
      return(character_1)
    }
  )
  snpedia_character <-
    unlist(x = lapply(
      X = snpedia_list,
      FUN = function(character_1) {
        return(paste(character_1, collapse = "\t"))
      }
    ))
  writeLines(text = snpedia_character, con = temporary_connection)
  rm(snpedia_character, snpedia_list)
}
close(con = snpedia_connection)
close(con = temporary_connection)
rm(snpedia_connection, temporary_connection)

message("Loading rewritten SNPedia GFF3 file")
snpedia_granges <- import.gff3(con = temporary_path)
snpedia_frame <- S4Vectors::mcols(x = snpedia_granges)

file.remove(temporary_path)
rm(temporary_path)

# Iterate over the VCF file -----------------------------------------------


message("Iterating over the VCF file")
vcf_file <- TabixFile(file = argument_list$vcf_path,
                      yieldSize = argument_list$chunk_size)
open(con = vcf_file)
output_frame <- NULL
sum_records_read <- 0L
sum_records_written <- 0L
while (nrow(x = vcf_object <- readVcf(file = vcf_file,
                                      genome = argument_list$genome_assembly))) {
  sum_records_read <- sum_records_read + nrow(x = vcf_object)
  # message(sprintf(fmt = "Number of VCF records read: %i (%i total)", nrow(x = vcf_object), sum_records_read))

  # Since the INFO variable "VEP_Existing_variation" contains ampersand-separated values,
  # they need further splitting. Convert into a character vector of variation identifiers.
  # info(x = vcf_object)$VEP_Existing_variation returns a CharacterList so that
  # unlist() needs calling before stringr::str_split().

  existing_variation <-
    unlist(
      x = stringr::str_split(
        string = unlist(x = info(x = vcf_object)$VEP_Existing_variation),
        pattern = stringr::fixed(pattern = "&")
      ),
      recursive = TRUE
    )
  # Find those SNPedia identifiers that match an existing variation in the VCF file
  # and extract the corresponding SNPedia records.
  sub_output_frame <-
    snpedia_frame[snpedia_frame$ID %in% existing_variation, ]
  sum_records_written <-
    sum_records_written + nrow(x = sub_output_frame)
  if (is.null(x = output_frame)) {
    output_frame <- sub_output_frame
  } else {
    output_frame <- rbind(output_frame, sub_output_frame)
  }
  rm(sub_output_frame, existing_variation)
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
close(con = vcf_file)

rm(snpedia_frame, snpedia_granges)

write.table(
  x = output_frame,
  file = argument_list$tsv_path,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

rm(output_frame)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = "Remaining objects:")
  print(x = ls())
}

print(x = sessionInfo())

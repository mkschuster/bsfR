#!/usr/bin/env Rscript
#
# Copyright 2013 - 2022 Michael K. Schuster
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

# Description -------------------------------------------------------------


# BSF R script to create internal data sets stored in the R/sysdata.rda file.

# Library Import ----------------------------------------------------------


# CRAN r-lib
suppressPackageStartupMessages(expr = library(package = "devtools"))
suppressPackageStartupMessages(expr = library(package = "sessioninfo"))
suppressPackageStartupMessages(expr = library(package = "usethis"))
# CRAN Tidyverse
suppressPackageStartupMessages(expr = library(package = "tibble"))

# Genome tibble template --------------------------------------------------


bsfg_genome_tibble <- tibble::tibble(
  "scientific_name" = character(),
  "genus" = character(),
  "species" = character(),
  "taxon" = integer(),
  "assembly_version_ncbi" = character(),
  "assembly_version_ucsc" = character(),
  "assembly_report_url" = character(),
  "bsgenome_ncbi" = character(),
  "bsgenome_ucsc" = character()
)

# Homo sapiens b37 GRCh37 -------------------------------------------------


bsfg_genome_tibble <- tibble::add_row(
  .data = bsfg_genome_tibble,
  "scientific_name" = "Homo sapiens",
  "genus" = "Homo",
  "species" = "sapiens",
  "taxon" = 9606L,
  "assembly_version_ncbi" = "b37",
  "assembly_version_ucsc" = "hg19",
  # "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.1_GRCh37/GCA_000001405.1_GRCh37_assembly_report.txt",
  "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_assembly_report.txt",
  "bsgenome_ncbi" = "BSgenome.Hsapiens.1000genomes.hs37d5",
  "bsgenome_ucsc" = "BSgenome.Hsapiens.UCSC.hg19"
)

# Homo sapiens hg19 GRCh37 ------------------------------------------------


bsfg_genome_tibble <- tibble::add_row(
  .data = bsfg_genome_tibble,
  "scientific_name" = "Homo sapiens",
  "genus" = "Homo",
  "species" = "sapiens",
  "taxon" = 9606L,
  "assembly_version_ncbi" = "GRCh37",
  "assembly_version_ucsc" = "hg19",
  # "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.1_GRCh37/GCA_000001405.1_GRCh37_assembly_report.txt",
  "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_assembly_report.txt",
  "bsgenome_ncbi" = NA_character_,
  "bsgenome_ucsc" = "BSgenome.Hsapiens.UCSC.hg19"
)

# Homo sapiens hg38 GRCh38 ------------------------------------------------


bsfg_genome_tibble <- tibble::add_row(
  .data = bsfg_genome_tibble,
  "scientific_name" = "Homo sapiens",
  "genus" = "Homo",
  "species" = "sapiens",
  "taxon" = 9606L,
  "assembly_version_ncbi" = "GRCh38",
  "assembly_version_ucsc" = "hg38",
  # "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_report+ucsc_names.txt",
  # "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_report.txt",
  "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_assembly_report.txt",
  "bsgenome_ncbi" = "BSgenome.Hsapiens.NCBI.GRCh38",
  "bsgenome_ucsc" = "BSgenome.Hsapiens.UCSC.hg38"
)

# Mus musculus mm10 GRCm38 ------------------------------------------------


bsfg_genome_tibble <- tibble::add_row(
  .data = bsfg_genome_tibble,
  "scientific_name" = "Mus musculus",
  "genus" = "Mus",
  "species" = "musculus",
  "taxon" = 10090L,
  "assembly_version_ncbi" = "GRCm38",
  "assembly_version_ucsc" = "mm10",
  # "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.2_GRCm38/GCA_000001635.2_GRCm38_assembly_report.txt",
  "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.20_GRCm38/GCF_000001635.20_GRCm38_assembly_report.txt",
  "bsgenome_ncbi" = NA_character_,
  "bsgenome_ucsc" = "BSgenome.Mmusculus.UCSC.mm10"
)

# Mus musculus mm39 GRCm39 ------------------------------------------------


bsfg_genome_tibble <- tibble::add_row(
  .data = bsfg_genome_tibble,
  "scientific_name" = "Mus musculus",
  "genus" = "Mus",
  "species" = "musculus",
  "taxon" = 10090L,
  "assembly_version_ncbi" = "GRCm39",
  "assembly_version_ucsc" = "mm39",
  # "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/GCA_000001635.9_GRCm39_assembly_report.txt",
  "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_assembly_report.txt",
  "bsgenome_ncbi" = NA_character_,
  "bsgenome_ucsc" = NA_character_,
  # The Bioconductor BSgenome is not available, yet.
  # "bsgenome_ucsc" = "BSgenome.Mmusculus.UCSC.mm39"
)

# Gallus gallus galGal5 ---------------------------------------------------


bsfg_genome_tibble <- tibble::add_row(
  .data = bsfg_genome_tibble,
  "scientific_name" = "Gallus gallus",
  "genus" = "Gallus",
  "species" = "gallus",
  "taxon" = 9031L,
  "assembly_version_ncbi" = "Gallus_gallus-5.0",
  "assembly_version_ucsc" = "galGal5",
  # "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002/315/GCA_000002315.3_Gallus_gallus-5.0/GCA_000002315.3_Gallus_gallus-5.0_assembly_report.txt",
  "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.4_Gallus_gallus-5.0/GCF_000002315.4_Gallus_gallus-5.0_assembly_report.txt",
  "bsgenome_ncbi" = NA_character_,
  "bsgenome_ucsc" = "BSgenome.Ggallus.UCSC.galGal5"
)

# Rattus norvegicus rn6 Rnor_6.0 ----------------------------------------


bsfg_genome_tibble <- tibble::add_row(
  .data = bsfg_genome_tibble,
  "scientific_name" = "Rattus norvegicus",
  "genus" = "Rattus",
  "species" = "norvegicus",
  "taxon" = 10116L,
  "assembly_version_ncbi" = "Rnor_6.0",
  "assembly_version_ucsc" = "rn6",
  # "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/895/GCA_000001895.4_Rnor_6.0/GCA_000001895.4_Rnor_6.0_assembly_report.txt",
  "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/895/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_assembly_report.txt",
  "bsgenome_ncbi" = NA_character_,
  "bsgenome_ucsc" = NA_character_
)

# Rattus norvegicus rn7 mRatBN7.2 ---------------------------------------


bsfg_genome_tibble <- tibble::add_row(
  .data = bsfg_genome_tibble,
  "scientific_name" = "Rattus norvegicus",
  "genus" = "Rattus",
  "species" = "norvegicus",
  "taxon" = 10116L,
  "assembly_version_ncbi" = "mRatBN7.2",
  "assembly_version_ucsc" = "rn7",
  # "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/015/227/675/GCA_015227675.2_mRatBN7.2/GCA_015227675.2_mRatBN7.2_assembly_report.txt",
  "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/675/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_assembly_report.txt",
  "bsgenome_ncbi" = NA_character_,
  "bsgenome_ucsc" = NA_character_
)

# Sarcophilus harrisii sarHar1 --------------------------------------------


bsfg_genome_tibble <- tibble::add_row(
  .data = bsfg_genome_tibble,
  "scientific_name" = "Sarcophilus harrisii",
  "genus" = "Sarcophilus",
  "species" = "harrisii",
  "taxon" = 9305L,
  "assembly_version_ncbi" = "DEVIL7.0",
  "assembly_version_ucsc" = "sarHar1",
  # "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/189/315/GCA_000189315.1_Devil_ref_v7.0/GCA_000189315.1_Devil_ref_v7.0_assembly_report.txt",
  "assembly_report_url" = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/189/315/GCF_000189315.1_Devil_ref_v7.0/GCF_000189315.1_Devil_ref_v7.0_assembly_report.txt",
  "bsgenome_ncbi" = NA_character_,
  "bsgenome_ucsc" = NA_character_
)

# Install the bsfg_genome_tibble as package-internal data set.

usethis::use_data(bsfg_genome_tibble,
                  internal = TRUE,
                  overwrite = TRUE)

rm(bsfg_genome_tibble)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessioninfo::session_info())

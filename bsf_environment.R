#! /usr/bin/env Rscript
#
# BSF R script to quickly setup the BSF environment after an R upgrade.
#
#
# Copyright 2013 - 2017 Michael K. Schuster
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

base::source(file = "http://bioconductor.org/biocLite.R")

# Let the user choose mirror options or ...
# utils::chooseCRANmirror()
# utils::chooseBioCmirror()
# ... set preferences via the options() function in e.g. .Rprofile.
# options(
#  repos = c("Austria (Vienna) [https]" = "https://cran.wu.ac.at/"),
#  BioC_mirror = c("United Kingdom (Hinxton) [https]" = "http://mirrors.ebi.ac.uk/bioconductor")
# )


# Update CRAN packages ----------------------------------------------------

message("Updating CRAN packages ...")
utils::update.packages(lib.loc = .Library, ask = FALSE)

# Update Bioconductor packages --------------------------------------------

message("Updating Bioconductor packages ...")
BiocInstaller::biocUpdatePackages(pkgs = c(),
                                  lib.loc = .Library,
                                  ask = FALSE)

# Check and install CRAN packages -----------------------------------------
#
# CRAN packages for GATK
# https://github.com/broadgsa/gatk-protected/blob/master/public/gatk-engine/src/main/resources/org/broadinstitute/gatk/engine/recalibration/BQSR.R
# gplots:
# reshape:
# grid: R system library, no need to install explicitly
# tools: R system library, no need to install explicitly

message("Checking CRAN packages ...")
cran_packages <- c(
  "devtools",
  # Development tools, including installation from GitHub repositories, etc.
  "ggplot2",
  # for almost eveything plot related, really ...
  "ggrepel",
  # for ggplot2 extension functions geom_text_repel() and geom_label_repel() to repel labes form data points.
  "gplots",
  # for GATK
  "gsalib",
  # for GATK
  "hexbin",
  # for ggplot2 functions geom_binhex() and stat_bin_hex()
  "optparse",
  # for option parsing
  "reshape",
  # for GATK
  "VGAM",  # for Bioconductor monocle
  "simpleCache", # FIXME: LOLA seems to depend on it, yet only suggest it. Why is it not declared?
  "wordcloud" # FIXME: RnBeads depends on it.
)
for (i in seq_along(along.with = cran_packages)) {
  if (suppressPackageStartupMessages(expr = base::require(
    package = cran_packages[i],
    lib.loc = .Library,
    quietly = TRUE,
    character.only = TRUE
  ))) {
    message(paste0("Skipping package ", cran_packages[i]))
    # base::detach(
    #   name = paste("package", cran_packages[i], sep = ":"),
    #   unload = TRUE,
    #   character.only = TRUE
    # )
  } else {
    message(paste0("Installing package ", cran_packages[i]))
    # utils::install.packages(pkgs = cran_packages[i],
    #                  lib = .Library)
  }
}
rm(i, cran_packages)

# Check and install Bioconductor packages ---------------------------------

message("Checking Bioconductor packages ...")
bioconductor_packages <- c(
  "BiocInstaller",
  "BiocParallel",
  # Biostrings genomes
  "BSgenome.Hsapiens.1000genomes.hs37d5",
  # For bsf_variant_calling_coverage.R
  "BSgenome.Hsapiens.UCSC.hg19",
  "BSgenome.Hsapiens.UCSC.hg19.masked",
  "BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Hsapiens.UCSC.hg38.masked",
  "BSgenome.Mmusculus.UCSC.mm10",
  "BSgenome.Mmusculus.UCSC.mm10.masked",
  #
  "ChIPpeakAnno",
  "ComplexHeatmap",
  "cummeRbund",
  # for Cufflinks mRNA-seq data processing
  "DESeq2",
  # for RNA-seq analysis
  "DiffBind",
  # for differential ChIP-seq analysis
  "goseq",
  # for Gene Ontology annotation
  "monocle",
  # For single cell RNA-seq
  "RnBeads",
  "RnBeads.hg38",
  "doParallel",  # FIXME: RnBeads seems to depend on it, yet only suggest it. Why is it not declared?
  "LOLA", # FIXME: Same as above. Sigh.
  # For Illumina Sequence Anaylsis Viewer information
  "savR",
  # For Meth-seq anaylsis (FDb.InfiniumMethylation.hg19)
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  # UCSC gene set
  "TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts",
  # UCSC gene set
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  # UCSC gene set
  "TxDb.Mmusculus.UCSC.mm10.ensGene",
  # Ensembl gene set
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  # UCSC gene set
  "VariantAnnotation"
)
for (i in seq_along(along.with = bioconductor_packages)) {
  if (suppressPackageStartupMessages(expr = base::require(
    package = bioconductor_packages[i],
    lib.loc = .Library,
    quietly = TRUE,
    character.only = TRUE
  ))) {
    message(paste0("Skipping package ", bioconductor_packages[i]))
    if (!bioconductor_packages[i] %in% c("BiocInstaller")) {
      # base::detach(
      #   name = paste("package", bioconductor_packages[i], sep = ":"),
      #   unload = TRUE,
      #   character.only = TRUE
      # )
    }
  } else {
    message(paste0("Installing package ", bioconductor_packages[i]))
    # BiocInstaller::biocLite(pkgs = bioconductor_packages[i],
    #          lib = .Library)
  }
}
rm(bioconductor_packages, i)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

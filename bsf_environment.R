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

source(file = "http://bioconductor.org/biocLite.R")

# Let the user choose mirror options or ...
# chooseCRANmirror()
# chooseBioCmirror()
# ... set preferences via the options() function in e.g. .Rprofile.
# options(
#  repos = c("Austria (Vienna) [https]" = "https://cran.wu.ac.at/"),
#  BioC_mirror = c("United Kingdom (Hinxton) [https]" = "http://mirrors.ebi.ac.uk/bioconductor")
# )


# Update all CRAN and Bioconductor packages.

update.packages()
biocUpdatePackages(pkgs = c())

# CRAN packages for BSF R scripts.
#
# CRAN packages for GATK
# https://github.com/broadgsa/gatk-protected/blob/master/public/gatk-engine/src/main/resources/org/broadinstitute/gatk/engine/recalibration/BQSR.R
# gplots:
# reshape:
# grid: R system library, no need to install explicitly
# tools: R system library, no need to install explicitly

install.packages(pkgs = c(
  "devtools",  # Development tools, including installation from GitHub repositories, etc.
  "ggplot2",  # for almost eveything plot related, really ...
  "ggrepel",  # for ggplot2 extension functions geom_text_repel() and geom_label_repel() to repel labes form data points.
  "gplots",  # for GATK
  "gsalib",  # for GATK
  "hexbin",  # for ggplot2 functions geom_binhex() and stat_bin_hex()
  "optparse",  # for option parsing
  "reshape"  # for GATK
))

# For the BSF RNA-Seq pipeline
# cummeRbund
# DESeq
#
# For the BSF ChIP-Seq pipeline
# DiffBind

biocLite(pkgs = c(
  "BiocParallel",
  # Biostrings genomes
  "BSgenome.Hsapiens.1000genomes.hs37d5", # For bsf_variant_calling_coverage.R
  "BSgenome.Hsapiens.UCSC.hg19",
  "BSgenome.Hsapiens.UCSC.hg19.masked",
  "BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Hsapiens.UCSC.hg38.masked",
  "BSgenome.Mmusculus.UCSC.mm10",
  "BSgenome.Mmusculus.UCSC.mm10.masked",
  "ChIPpeakAnno",
  "ComplexHeatmap",
  "cummeRbund",  # for Cufflinks data processing
  "DESeq2",  # for RNA-seq analysis
  "DiffBind",  # for differential ChIP-seq analysis
  "goseq",  # for Gene Ontology annotation
  "RnBeads",  # For Meth-seq anaylsis
  "TxDb.Hsapiens.UCSC.hg19.knownGene", # UCSC gene set
  "TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts",  # UCSC gene set
  "TxDb.Hsapiens.UCSC.hg38.knownGene",  # UCSC gene set
  "TxDb.Mmusculus.UCSC.mm10.ensGene",  # Ensembl gene set
  "TxDb.Mmusculus.UCSC.mm10.knownGene",  # UCSC gene set
  "VariantAnnotation"
  ))

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

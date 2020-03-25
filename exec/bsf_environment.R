#!/usr/bin/env Rscript
#
# BSF R script to quickly setup the BSF environment after an R upgrade.
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

# base::source(file = "http://bioconductor.org/biocLite.R")

# Let the user choose mirror options or ...
# utils::chooseCRANmirror()
# utils::chooseBioCmirror()
# ... set preferences via the options() function in e.g. .Rprofile.
# options(
#  repos = c("Austria (Vienna) [https]" = "https://cran.wu.ac.at/"),
#  BioC_mirror = c("United Kingdom (Hinxton) [https]" = "http://mirrors.ebi.ac.uk/bioconductor")
#  BioC_mirror = c("Department of Statistics, TU Dortmund [https]" = "https://bioconductor.statistik.tu-dortmund.de")
# )

# R 3.6.0 implemented a staged install that is not compatible with the BeeGFS (formerly FhGFS) file system.
# A move of a directory on top of an existing one leads to an operating system error.
# https://stat.ethz.ch/pipermail/r-devel/2019-May/077737.html
# As a workaround, set the following Bash environment variable.
# declare -x R_INSTALL_STAGED='false';

# Global install
# library_location <- .Library
# Private install
# library_location <- .libPaths()[1]
library_location <- NULL

# Update CRAN packages ----------------------------------------------------


message("Updating CRAN packages ...")
utils::update.packages(lib.loc = library_location, ask = FALSE)

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
  # Hadley Wickham's Tidyverse
  #    ggplot2  Create Elegant Data Visualisations Using the Grammar of Graphics
  #    purrr Functional Programming Tools
  #    tibble Simple Data Frames
  #    dplyr A Grammar of Data Manipulation
  #    tidyr Easily Tidy Data with 'spread()' and 'gather()' Functions
  #    stringr Simple, Consistent Wrappers for Common String Operations
  #    readr Read Rectangular Text Data
  #    forcats Tools for Working with Categorical Variables (Factors)
  "tidyverse", # Easily Install and Load the 'Tidyverse'

  # For package development
  "devtools", # Tools to Make Developing R Packages Easier
  "roxygen2", # In-Line Documentation for R
  "knitr", # A General-Purpose Package for Dynamic Report Generation in R
  "testthat", # Unit Testing for R

  # For bsf_rnaseq_deseq_analysis.R (DESeq2)
  "ashr", # Methods for Adaptive Shrinkage, using Empirical Bayes
  "caret", # Classification and Regression Training
  "enrichR", # Provides an R Interface to 'Enrichr'

  "argparser", # Command-Line Argument Parser

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
  "reshape", # Flexibly Reshape Data

  "spp", # ChIP-Seq Processing Pipeline
  # For attac-seq pipleine
  # for GATK
  "VGAM",
  "Seurat",
  # for Bioconductor monocle
  "simpleCache",
  # FIXME: LOLA seems to depend on it, yet only suggest it. Why is it not declared?
  "wordcloud",
  # FIXME: RnBeads depends on it.
  "RMariaDB",  # For GenomicFeatures::makeTxDbFromEnsembl(). FIXME: This may require a new installation of MySQ or MariaDB libraries.
  "BiocManager" # Access the Bioconductor Project Package Repository
)
for (i in seq_along(along.with = cran_packages)) {
  if (base::requireNamespace(package = cran_packages[i],
                             lib.loc = library_location,
                             quietly = TRUE)) {
    message("Skipping package ", cran_packages[i])
  } else {
    message("Installing package ", cran_packages[i])
    utils::install.packages(pkgs = cran_packages[i],
                            lib = library_location)
  }
}
rm(i, cran_packages)

# Update Bioconductor packages --------------------------------------------

message("Updating Bioconductor packages ...")
BiocManager::install(lib.loc = library_location,
                     update = TRUE,
                     ask = FALSE)

# Check and install Bioconductor packages ---------------------------------

message("Checking Bioconductor packages ...")
bioconductor_packages <- c(
  # Data packages

  # Biostrings genomes
  # "BSgenome.Ggallus.UCSC.galGal5",
  # "BSgenome.Hsapiens.1000genomes.hs37d5", # For bsf_variant_calling_coverage.R
  # "BSgenome.Hsapiens.UCSC.hg19",
  # "BSgenome.Hsapiens.UCSC.hg19.masked",
  # "BSgenome.Hsapiens.UCSC.hg38",
  # "BSgenome.Hsapiens.UCSC.hg38.masked",
  # "BSgenome.Mmusculus.UCSC.mm10",
  # "BSgenome.Mmusculus.UCSC.mm10.masked",

  # "org.Gg.eg.db", # Genome wide annotation for Chicken
  # "org.Hs.eg.db", # Genome wide annotation for Human
  # "org.Mm.eg.db", # Genome wide annotation for Mouse

  # "Homo.sapiens", # Annotation package for the Homo.sapiens object
  # "Mus.musculus", # Annotation package for the Mus.musculus object

  # "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  # "IlluminaHumanMethylation450kmanifest",
  # "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
  # "IlluminaHumanMethylationEPICanno.ilm10b3.hg19",
  # "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  # "IlluminaHumanMethylationEPICmanifest",

  # TxDb objects
  # "TxDb.Hsapiens.UCSC.hg19.knownGene",
  # "TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts",
  # "TxDb.Hsapiens.UCSC.hg38.knownGene",
  # "TxDb.Mmusculus.UCSC.mm10.ensGene",
  # "TxDb.Mmusculus.UCSC.mm10.knownGene",

  "BiocParallel",

  "AnnotationHub", # Client to access AnnotationHub resources
  "ChIPpeakAnno",
  # for structural variant calling
  "CODEX", # A Normalization and Copy Number Variation Detection Method for Whole Exome Sequencing
  "ComplexHeatmap", # Make Complex Heatmaps
  "EnhancedVolcano", # Publication-ready volcano plots with enhanced colouring and labeling
  # for Cufflinks mRNA-seq data processing
  "cummeRbund", # Analysis, exploration, manipulation, and visualization of Cufflinks high-throughput sequencing data.
  # for RNA-seq analysis
  "DESeq2", # Differential gene expression analysis based on the negative binomial distribution
  "apeglm", # Approximate posterior estimation for GLM coefficients
  "vsn", # Variance stabilization and calibration for microarray data
  # for differential ChIP-seq analysis
  "ChIPQC", # Quality metrics for ChIPseq data
  "DiffBind", # Differential Binding Analysis of ChIP-Seq Peak Data
  "BCRANK",  # Predicting binding site consensus from ranked DNA sequences
  "rGADEM",  # de novo motif discovery
  "DNAcopy",  # DNA copy number data analysis
  "goseq", # Gene Ontology analyser for RNA-seq and other length biased data
  "heatmaps", # Flexible Heatmaps for Functional Genomics and Sequence Features
  # For PWM to JASPAR conversion.
  "universalmotif",
  # For RnBeads
  "impute",
  # for Gene Ontology annotation
  "monocle",
  "PSCBS",
  "PureCN",
  # For single cell RNA-seq
  # "RnBeads",
  # "RnBeads.hg19",
  # "RnBeads.hg38",
  # "RnBeads.mm10",
  "doParallel",
  # FIXME: RnBeads seems to depend on it, yet only suggest it. Why is it not declared?
  "LOLA",
  # FIXME: Same as above. Sigh.
  # For Illumina Sequence Analysis Viewer information
  "savR",
  # For single cell analysis.
  "SingleR", # Reference-Based Single-Cell RNA-Seq Annotation
  "scater", # Single-Cell Analysis Toolkit for Gene Expression Data in R
  # For Meth-seq analysis (FDb.InfiniumMethylation.hg19)

  "topGO",  # Enrichment analysis for Gene Ontology

  "VariantAnnotation" # Annotation of Genetic Variants
)
for (i in seq_along(along.with = bioconductor_packages)) {
  if (base::requireNamespace(package = bioconductor_packages[i],
                             lib.loc = library_location,
                             quietly = FALSE)) {
    message("Skipping package ", bioconductor_packages[i])
  } else {
    message("Installing package ", bioconductor_packages[i])
    BiocManager::install(
      pkgs = bioconductor_packages[i],
      lib = library_location,
      update = FALSE,
      ask = FALSE
    )
  }
}
rm(i, bioconductor_packages, library_location)

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

print(x = sessionInfo())

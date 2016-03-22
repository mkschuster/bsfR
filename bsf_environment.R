#! /usr/bin/env Rscript
#
# BSF R script to quickly setup the BSF environment after an R upgrade.
#
#
# Copyright 2013 -2015 Michael K. Schuster
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
# optparse: option parsing
# ggplot2: for eveything graphics, really
#
# CRAN packages for GATK
# https://github.com/broadgsa/gatk-protected/blob/master/public/gatk-engine/src/main/resources/org/broadinstitute/gatk/engine/recalibration/BQSR.R
# gplots:
# reshape:
# grid: R system library, no need to install explicitly
# tools: R system library, no need to install explicitly
# gsalib:

install.packages(pkgs = c(
  "optparse",
  "ggplot2",
  "gplots",
  "reshape",
  "gsalib"
))

# For the BSF RNA-Seq pipeline
# cummeRbund
# DESeq
#
# For the BSF ChIP-Seq pipeline
# DiffBind

biocLite(pkgs = c("DESeq", "cummeRbund", "DiffBind", "VariantAnnotation"))

# For RnBeads
# The IlluminaHumanMethylation450k.db is > 60 MB
# TODO: This should be taken care of by an eventual RnBeads package.
# biocLite(pkgs = c("foreach", "mclust", "RPMM", "fields",
#                   "matrixStats", "IlluminaHumanMethylation450k.db", "methylumi",
#                   "ggbio", "GEOquery", "GOstats", "wordcloud"))

message("All done")

# Finally, print all objects that have not been removed from the environment.
if (length(x = ls())) {
  print(x = ls())
}

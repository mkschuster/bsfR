#! /usr/bin/env Rscript
#
# Rscript to quickly setup the BSF environment after an R upgrade.
#
# Copyright 2013 -2015 Michael K. Schuster
#
# Biomedical Sequencing Facility (BSF), part of the genomics core facility
# of the Research Center for Molecular Medicine (CeMM) of the
# Austrian Academy of Sciences and the Medical University of Vienna (MUW).
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

# CRAN package for BSF R scripts.
install.packages("optparse")

# For the BSF RNA-Seq pipeline
biocLite(pkgs = c("DESeq", "cummeRbund"))

# For the BSF ChIP-Seq pipeline
biocLite(pkgs = c("DiffBind"))

# For RnBeads
# The IlluminaHumanMethylation450k.db is > 60 MB
# TODO: This should be taken care of by an eventual RnBeads package.
biocLite(pkgs = c("foreach", "mclust", "RPMM", "fields",
                  "matrixStats", "IlluminaHumanMethylation450k.db", "methylumi",
                  "ggbio", "GEOquery", "GOstats", "wordcloud"))

# CRAN package for GATK
install.packages("ggplot2")
install.packages("gsalib")

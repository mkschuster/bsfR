#! /usr/bin/env Rscript

# Rscript to quickly setup the BSF environment after an R upgrade.

source(file="http://bioconductor.org/biocLite.R")

# For the BSF RNA-Seq pipeline
biocLite(pkgs=c("DESeq", "cummeRbund"))

# For the BSF ChIP-Seq pipeline
biocLite(pkgs=c("DiffBind"))

# For RnBeads
# The IlluminaHumanMethylation450k.db is > 60 MB
# TODO: This should be taken care of by an eventual RnBeads package.
biocLite(pkgs=c("foreach", "mclust", "RPMM", "fields",
           "matrixStats", "IlluminaHumanMethylation450k.db", "methylumi",
           "ggbio", "GEOquery", "GOstats", "wordcloud"))

#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(package="optparse"))
suppressPackageStartupMessages(library(package="cummeRbund"))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")

option_list = list(
  make_option(opt_str=c("-v", "--verbose"), action="store_true",
              default=TRUE,
              help="Print extra output [default]"),
  make_option(opt_str=c("-q", "--quietly"), action="store_false",
              default=FALSE,
              dest="verbose",
              help="Print little output"),
  make_option(opt_str=c("-c", "--comparison"),
              dest="comparison",
              help="Comparison name"),
  make_option(opt_str=c("-g", "--genome_directory"),
              dest="genome_directory",
              help="Genome-specific analysis directory")
)

# Get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,

opt = parse_args(OptionParser(option_list=option_list))

prefix = paste("rnaseq", "cuffdiff", opt$comparison, sep="_")

# Read and index Cuffdiff output and create a CuffSet object..

message("Create or load the cummeRbund database")
cuff_set = readCufflinks(dir=file.path(opt$genome_directory, prefix,
                                       fsep=.Platform$file.sep))

# Create a new sub-directory for plots.

prefix = paste("rnaseq", "process", "cuffdiff", opt$comparison, sep="_")

output_directory = file.path(opt$genome_directory, prefix, fsep=.Platform$file.sep)

if (!file.exists(output_directory)) {
  dir.create(output_directory, showWarnings=TRUE, recursive=FALSE)
}

# Some plots require replicates.

have_replicates = (nrow(replicates(cuff_set)) > 2)

# Create Dispersion Plots for Genes and Isoforms.

plot_name = paste(prefix, "_genes_dispersion.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
message("Creating dispersion plot for genes")
dispersionPlot(object=genes(object=cuff_set))
active_device = dev.off()

plot_name = paste(prefix, "_isoforms_dispersion.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
message("Creating dispersion plot for isoforms")
dispersionPlot(object=isoforms(object=cuff_set))
active_device = dev.off()

# Create Squared Coefficient of Variation (SCV) Plots for Genes and Isoforms.
# The plot requires replicates

plot_name = paste(prefix, "_genes_scv.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
if (have_replicates) {
  message("Creating SCV plot for genes")
  fpkmSCVPlot(object=genes(object=cuff_set))
} else {
  message("Skipping SCV plot for genes")  
}
active_device = dev.off()

plot_name = paste(prefix, "_isoforms_scv.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
if (have_replicates) {
  message("Creating SCV plot for isoforms")
  fpkmSCVPlot(object=isoforms(object=cuff_set))
} else {
  message("Skipping SCV plot for isoforms")
}
active_device = dev.off()

# Create Density Plots for Genes with and without Replicates.

plot_name = paste(prefix, "_genes_density_wo_replicates.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
message("Creating Density plot for genes without replicates")
csDensity(object=genes(object=cuff_set), replicates=FALSE)
active_device = dev.off()

plot_name = paste(prefix, "_genes_density_w_replicates.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
message("Creating Density plot for genes with replicates")
csDensity(object=genes(object=cuff_set), replicates=TRUE)
active_device = dev.off()

# Create Density Plots for Isoforms with and without Replicates.

plot_name = paste(prefix, "_isoforms_density_wo_replicates.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
message("Creating Density plot for isoforms without replicates")
csDensity(object=isoforms(object=cuff_set), replicates=FALSE)
active_device = dev.off()

plot_name = paste(prefix, "_isoforms_density_w_replicates.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
message("Creating Density plot for isoforms with replicates")
csDensity(object=isoforms(object=cuff_set), replicates=TRUE)
active_device = dev.off()

# Create Box Plots for Genes with and without Replicates.

plot_name = paste(prefix, "_genes_box_w_replicates.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
message("Creating Box plot for genes with replicates")
csBoxplot(object=genes(object=cuff_set), replicates=TRUE)
active_device = dev.off()

plot_name = paste(prefix, "_genes_box_wo_replicates.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
message("Creating Box plot for genes without replicates")
csBoxplot(object=genes(object=cuff_set), replicates=FALSE)
active_device = dev.off()

# Create a Scatter Matrix Plot of pairwise comparisons.

plot_name = paste(prefix, "_genes_scatter_matrix.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
message("Creating Scatter Matrix plot for genes")
csScatterMatrix(object=genes(object=cuff_set))
active_device = dev.off()

plot_name = paste(prefix, "_isoforms_scatter_matrix.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
message("Creating Scatter Matrix plot for isoforms")
csScatterMatrix(object=isoforms(object=cuff_set))
active_device = dev.off()

# TODO: Scatter Plot by sample pair.
# TODO: Dendrograms are only meaningful for time-series analyses.
# TODO: MAplot for genes with and without count data.

# Create a MAplot for a gene sample pair
sample_list = samples(object=genes(object=cuff_set))

plot_name = paste(prefix, "_maplot.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
message("Creating MAplot")
MAplot(object=genes(object=cuff_set), x=sample_list[1], y=sample_list[2])
active_device = dev.off()

# Create a Volcano Matrix Plot of pairwise comparisons.

plot_name = paste(prefix, "_genes_volcano_matrix.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
message("Creating Volcano Matrix plot for genes")
csVolcanoMatrix(object=genes(object=cuff_set))
active_device = dev.off()

# TODO: Volcano plots by sample pair ...

# Create a Multidimensional Scaling (MDS) Plot

plot_name = paste(prefix, "_genes_mds.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
if (have_replicates) {
  message("Creating Multidimensional Scaling Plot for genes")
  MDSplot(object=genes(object=cuff_set), replicates=TRUE)
} else {
  message("Skipping Multidimensional Scaling Plot for genes")
}
active_device = dev.off()

# Create a Principal Component Analysis (PCA) Plot

plot_name = paste(prefix, "_genes_pca.png", sep="")
png(filename=file.path(output_directory, plot_name, fsep=.Platform$file.sep))
message("Creating Principal Component Analysis Plot (PCA) for genes")
PCAplot(object=genes(object=cuff_set), x="PC1", y="PC2", replicates=TRUE)
active_device = dev.off()

# TODO: Process gene and isoform sets and lists.
# Show Cuffdiff run options
# runInfo(object=cuff_set)
# Show replicate names, total map masses and normalised map masses.
# replicates(object=cuff_set)

# Create a new, simpler data frame for gene annotation and merge both data frames

message("Creating FPKM per replicates for genes")
# gene_replicates = replicates(object=genes(object=cuff_set))
# gene_rep_fpkms = repFpkm(object=genes(object=cuff_set))
# gene_fpkm_matrix = fpkmMatrix(object=genes(object=cuff_set))
gene_annotation = annotation(object=genes(object=cuff_set))
gene_rep_fpkm_matrix = repFpkmMatrix(object=genes(object=cuff_set))
gene_frame = data.frame(gene_annotation$gene_short_name,
                        gene_annotation$locus,
                        row.names=gene_annotation$gene_id)
colnames(gene_frame)[1] = "gene_short_name"
colnames(gene_frame)[2] = "locus"
gene_merge = merge(x=gene_frame, y=gene_rep_fpkm_matrix, by="row.names",
                   all=TRUE, sort=TRUE)
file_name = paste(prefix, "_genes_fpkm_replicates.txt", sep="")
write.table(x=gene_merge,
            file=file.path(output_directory, file_name, fsep=.Platform$file.sep),
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

message("Creating FPKM per replicates for isoforms")
isoform_annotation = annotation(object=isoforms(object=cuff_set))
isoform_frame = data.frame(isoform_annotation$gene_short_name,
                           isoform_annotation$nearest_ref_id,
                           isoform_annotation$locus,
                           isoform_annotation$length,
                           row.names=isoform_annotation$isoform_id)
colnames(isoform_frame)[1] = "gene_short_name"
colnames(isoform_frame)[2] = "nearest_ref_id"
colnames(isoform_frame)[3] = "locus"
colnames(isoform_frame)[4] = "length"
isoform_rep_fpkm_matrix = repFpkmMatrix(object=isoforms(object=cuff_set))
isoform_merge = merge(x=isoform_frame, y=isoform_rep_fpkm_matrix, by="row.names",
                      all=TRUE, sort=TRUE)
file_name = paste(prefix, "_isoforms_fpkm_replicates.txt", sep="")
write.table(x=isoform_merge,
            file=file.path(output_directory, file_name, fsep=.Platform$file.sep),
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Finally, create comparison-specific symbolic links for cuffdiff results in the
# rnaseq_process_cuffdiff_* directory to avoid identical file names between comparisons.

result = file.symlink(from=file.path(opt$genome_directory,
                                     paste("rnaseq", "cuffdiff", opt$comparison, sep="_"),
                                     "cds_exp.diff",
                                     fsep=.Platform$file.sep),
                      to=file.path(opt$genome_directory,
                                   paste("rnaseq", "process", "cuffdiff", opt$comparison, sep="_"),
                                   paste(prefix, "_cds_exp_diff.txt", sep=""),
                                   fsep=.Platform$file.sep))

result = file.symlink(from=file.path(opt$genome_directory,
                                     paste("rnaseq", "cuffdiff", opt$comparison, sep="_"),
                                     "gene_exp.diff",
                                     fsep=.Platform$file.sep),
                      to=file.path(opt$genome_directory,
                                   paste("rnaseq", "process", "cuffdiff", opt$comparison, sep="_"),
                                   paste(prefix, "_gene_exp_diff.txt", sep=""),
                                   fsep=.Platform$file.sep))

result = file.symlink(from=file.path(opt$genome_directory,
                                     paste("rnaseq", "cuffdiff", opt$comparison, sep="_"),
                                     "isoform_exp.diff",
                                     fsep=.Platform$file.sep),
                      to=file.path(opt$genome_directory,
                                   paste("rnaseq", "process", "cuffdiff", opt$comparison, sep="_"),
                                   paste(prefix, "_isoform_exp_diff.txt", sep=""),
                                   fsep=.Platform$file.sep))

result = file.symlink(from=file.path(opt$genome_directory,
                                     paste("rnaseq", "cuffdiff", opt$comparison, sep="_"),
                                     "promoters.diff",
                                     fsep=.Platform$file.sep),
                      to=file.path(opt$genome_directory,
                                   paste("rnaseq", "process", "cuffdiff", opt$comparison, sep="_"),
                                   paste(prefix, "_promoters_diff.txt", sep=""),
                                   fsep=.Platform$file.sep))

result = file.symlink(from=file.path(opt$genome_directory,
                                     paste("rnaseq", "cuffdiff", opt$comparison, sep="_"),
                                     "splicing.diff",
                                     fsep=.Platform$file.sep),
                      to=file.path(opt$genome_directory,
                                   paste("rnaseq", "process", "cuffdiff", opt$comparison, sep="_"),
                                   paste(prefix, "_splicing_diff.txt", sep=""),
                                   fsep=.Platform$file.sep))

result = file.symlink(from=file.path(opt$genome_directory,
                                     paste("rnaseq", "cuffdiff", opt$comparison, sep="_"),
                                     "tss_group_exp.diff",
                                     fsep=.Platform$file.sep),
                      to=file.path(opt$genome_directory,
                                   paste("rnaseq", "process", "cuffdiff", opt$comparison, sep="_"),
                                   paste(prefix, "_tss_group_exp_diff.txt", sep=""),
                                   fsep=.Platform$file.sep))

message("All done")

#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(package="optparse"))
suppressPackageStartupMessages(library(package="cummeRbund"))

# Specify the desired options in a list.
# By default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit").

option_list = list(
  make_option(opt_str=c("--verbose", "-v"),
              action="store_true",
              default=TRUE,
              help="Print extra output [default]"),
  make_option(opt_str=c("--quiet", "-q"),
              action="store_false",
              default=FALSE,
              dest="verbose",
              help="Print little output"),
  make_option(opt_str=c("--comparison-name"),
              dest="comparison_name",
              help="Comparison name"),
  make_option(opt_str=c("--genome-directory"),
              default=NULL,
              dest="genome_directory",
              help="Genome-specific analysis directory"),
  make_option(opt_str=c("--gtf-file"),
              default=NULL,
              dest="gtf_file",
              help="GTF file specifying a reference transcriptome"),
  make_option(opt_str=c("--genome-version"),
              default=NULL,
              dest="genome_version",
              help="Genome version")
)

# Get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults.

opt = parse_args(object=OptionParser(option_list=option_list))

# Define CuffDiff and output directory names relative to the working directory.

cuffdiff_directory = paste("rnaseq", "cuffdiff", opt$comparison_name, sep="_")
output_directory = paste("rnaseq", "process", "cuffdiff", opt$comparison_name, sep="_")

# To avoid name clashes when downloading files, use the output directory name also as a prefix for all files therein.

prefix = output_directory

# Read and index Cuffdiff output and create a CuffSet object.
# The CuffSet object has slots genes, isoforms, TSS and CDS that are each instances of teh CuffData class.
# By default, CuffData accessor methods applied to a CuffSet class will operate on the ’genes’ slot.

message("Create or load the cummeRbund database")
cuff_set = readCufflinks(dir=cuffdiff_directory, gtfFile=opt$gtf_file, genome=opt$genome_version)

# Create a new sub-directory for plots if it does not exist.

if (!file.exists(output_directory)) {
  dir.create(output_directory, showWarnings=TRUE, recursive=FALSE)
}

# Store a table with Cuffdiff run information.

frame_path = file.path(output_directory, paste0(prefix, "_run_information.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
  message("Skipping run information table")
} else {
  message("Creating run information table")
  write.table(x=runInfo(object=cuff_set), file=frame_path, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}

# Store a table with sample information and define all possible pairwise sample comparisons or sample pairs.

sample_frame = samples(object=cuff_set)
frame_path = file.path(output_directory, paste0(prefix, "_samples.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
  message("Skipping sample table")
} else {
  message("Creating sample table")
  write.table(x=sample_frame, file=frame_path, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}
sample_number = nrow(x=sample_frame)
sample_pairs = combn(x=sample_frame$sample_name, m=2)
rm(sample_frame)

# Store a table of sample pair information by transposing the sample pairs array.
# This table allows the Python web code to link in pairwise plots.

frame_path = file.path(output_directory, paste0(prefix, "_sample_pairs.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
  message("Skipping sample pairs table")
} else {
  message("Create sample pairs table")
  write.table(x=aperm(a=sample_pairs), file=frame_path, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}

# Store a table of replicate information and find out, whether the analysis has replicates.

replicate_frame = replicates(object=cuff_set)
frame_path = file.path(output_directory, paste0(prefix, "_replicates.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
  message("Skipping replicate table")
} else {
  message("Creating replicate table")
  write.table(x=replicate_frame, file=frame_path, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}
# Some plots require replicates. Check that there is at least one row with a replicate column value greater than 0.
replicate_number = nrow(x=replicate_frame)
have_replicates = (nrow(x=replicate_frame[replicate_frame$replicate > 0, ]) > 0)
rm(replicate_frame)

rm(frame_path)

# Create a set of QC plots.

message("Started plotting.")

# Create a Dispersion Plot on Genes.

plot_path = file.path(output_directory, paste0(prefix, "_genes_dispersion.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a Dispersion Plot on Genes")
} else {
  message("Creating a Dispersion Plot on Genes")
  dispersionPlot(object=genes(object=cuff_set))
  ggsave(filename=plot_path)
}

# Create a Dispersion Plot on Isoforms.

plot_path = file.path(output_directory, paste0(prefix, "_isoforms_dispersion.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a Dispersion Plot on Isoforms")
} else {
  message("Creating Dispersion Plot on Isoforms")
  dispersionPlot(object=isoforms(object=cuff_set))
  ggsave(filename=plot_path)
}

# Create a Squared Coefficient of Variation (SCV) Plot on Genes.
# The plot requires replicates.

plot_path = file.path(output_directory, paste0(prefix, "_genes_scv.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a SCV Plot on Genes")
} else {
  if (have_replicates) {
    message("Creating a SCV Plot on Genes")
    fpkmSCVPlot(object=genes(object=cuff_set))
    ggsave(filename=plot_path)
  } else {
    message("Skipping a SCV Plot on Genes in lack of replicates")  
  }
}

# Create a Squared Coefficient of Variation (SCV) Plot on Isoforms.
# The plot requires replicates.

plot_path = file.path(output_directory, paste0(prefix, "_isoforms_scv.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a SCV Plot on Isoforms")
} else {
  if (have_replicates) {
    message("Creating a SCV Plot on Isoforms")
    fpkmSCVPlot(object=isoforms(object=cuff_set))
    ggsave(filename=plot_path)
  } else {
    message("Skipping a SCV Plot on Isoforms in lack of replicates")
  }
}

# Create a Density Plot on Genes with and without replicates.

plot_path = file.path(output_directory, paste0(prefix, "_genes_density_wo_replicates.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a Density Plot on Genes without replicates")
} else {
  message("Creating a Density Plot on Genes without replicates")
  csDensity(object=genes(object=cuff_set), replicates=FALSE)
  ggsave(filename=plot_path)
}

plot_path = file.path(output_directory, paste0(prefix, "_genes_density_w_replicates.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a Density Plot on Genes with replicates")
} else {
  message("Creating a Density Plot on Genes with replicates")
  csDensity(object=genes(object=cuff_set), replicates=TRUE)
  ggsave(filename=plot_path)
}

# Create a Density Plot on Isoforms with and without replicates.

plot_path = file.path(output_directory, paste0(prefix, "_isoforms_density_wo_replicates.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a Density Plot on Isoforms without replicates")
} else {
  message("Creating a Density Plot on Isoforms without replicates")
  csDensity(object=isoforms(object=cuff_set), replicates=FALSE)
  ggsave(filename=plot_path)
}

plot_path = file.path(output_directory, paste0(prefix, "_isoforms_density_w_replicates.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a Density Plot on Isoforms with replicates")
} else {
  message("Creating a Density Plot on Isoforms with replicates")
  csDensity(object=isoforms(object=cuff_set), replicates=TRUE)
  ggsave(filename=plot_path)
}

# Create a Box Plot on Genes with and without replicates.

plot_path = file.path(output_directory, paste0(prefix, "_genes_box_w_replicates.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a Box Plot on Genes with replicates")
} else {
  message("Creating a Box Plot on Genes with replicates")
  csBoxplot(object=genes(object=cuff_set), replicates=TRUE)
  ggsave(filename=plot_path)
}

plot_path = file.path(output_directory, paste0(prefix, "_genes_box_wo_replicates.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a Box Plot on Genes without replicates")
} else {
  message("Creating a Box Plot on Genes without replicates")
  csBoxplot(object=genes(object=cuff_set), replicates=FALSE)
  ggsave(filename=plot_path)
}

# Create a Scatter Matrix Plot on Genes and Isoforms.

plot_path = file.path(output_directory, paste0(prefix, "_genes_scatter_matrix.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a Scatter Matrix Plot on Genes")
} else {
  message("Creating a Scatter Matrix Plot on Genes")
  csScatterMatrix(object=genes(object=cuff_set))
  ggsave(filename=plot_path)
}

plot_path = file.path(output_directory, paste0(prefix, "_isoforms_scatter_matrix.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a Scatter Matrix Plot on Isoforms")
} else {
  message("Creating a Scatter Matrix Plot on Isoforms")
  csScatterMatrix(object=isoforms(object=cuff_set))
  ggsave(filename=plot_path)
}

# Create a Scatter Plot on Genes for each sample pair.

for (i in 1:length(sample_pairs[1,])) {
  plot_path = file.path(output_directory, paste(prefix, sample_pairs[1, i], sample_pairs[2, i], "genes_scatter.png", sep="_"))
  if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
    message(paste("Skipping a Scatter Plot on Genes for", sample_pairs[1, i], "versus", sample_pairs[2, i]))
  } else {
    message(paste("Creating a Scatter Plot on Genes for", sample_pairs[1, i], "versus", sample_pairs[2, i]))
    csScatter(object=genes(object=cuff_set), x=sample_pairs[1, i], y=sample_pairs[2, i])
    ggsave(filename=plot_path)
  }
}

# Create a Dendrogram Plot on Genes for time-series analyses.
# The csDendro function returns an dendrogram object that cannot be saved with the ggsave.

plot_path = file.path(output_directory, paste0(prefix, "_genes_dendrogram.png"))
png(filename=plot_path)
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a Dendrogram Plot on Genes")
} else {
  message("Creating a Dendrogram Plot on Genes")
  csDendro(object=genes(object=cuff_set))
  # ggsave(filename=plot_path)
}
active_device = dev.off()

# Create a MA Plot on genes for each sample pair based on FPKM values.

for (i in 1:length(sample_pairs[1,])) {
  plot_path = file.path(output_directory, paste(prefix, sample_pairs[1, i], sample_pairs[2, i], "maplot.png", sep="_"))
  if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
    message(paste("Skipping a MAplot on Genes for", sample_pairs[1, i], "versus", sample_pairs[2, i]))
  } else {
    message(paste("Creating a MAplot on Genes for", sample_pairs[1, i], "versus", sample_pairs[2, i]))
    MAplot(object=genes(object=cuff_set), x=sample_pairs[1, i], y=sample_pairs[2, i])
    ggsave(filename=plot_path)
  }
}

# TODO: Create a MAplot on genes for each sample pair based on count data.

# Create a Volcano Matrix Plot on Genes.

plot_path = file.path(output_directory, paste0(prefix, "_genes_volcano_matrix.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a Volcano Matrix Plot on Genes")
} else {
  message("Creating a Volcano Matrix Plot on Genes")
  csVolcanoMatrix(object=genes(object=cuff_set))
  ggsave(filename=plot_path)
}

# Create a Volcano Plot on genes for each sample pair.

for (i in 1:length(sample_pairs[1,])) {
  plot_path = file.path(output_directory, paste(prefix, sample_pairs[1, i], sample_pairs[2, i], "genes_volcano.png", sep="_"))
  if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
    message(paste("Skipping a Volcano Plot on Genes for", sample_pairs[1, i], "versus", sample_pairs[2, i]))
  } else {
    message(paste("Creating a Volcano Plot on Genes for", sample_pairs[1, i], "versus", sample_pairs[2, i]))
    csVolcano(object=genes(object=cuff_set), x=sample_pairs[1, i], y=sample_pairs[2, i])
    ggsave(filename=plot_path)
  }
}

# Create a Multidimensional Scaling (MDS) Plot on genes.

plot_path = file.path(output_directory, paste0(prefix, "_genes_mds.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a Multidimensional Scaling Plot on Genes")
} else {
# if (have_replicates) {
    message("Creating a Multidimensional Scaling Plot on Genes")
    # Nothing ever is simple. If the set has too many replicates, the standard cummeRbund MDSplot() falls down.
    if (replicate_number <= 24) {
      MDSplot(object=genes(object=cuff_set), replicates=TRUE)
      ggsave(filename=plot_path)
      plot_path = file.path(output_directory, paste0(prefix, "_genes_mds.pdf"))
      ggsave(filename=plot_path)
    } else {
      # The standard MDSplot has too many replicates.
      gene_rep_fit = cmdscale(d=JSdist(mat=makeprobs(a=repFpkmMatrix(object=genes(object=cuff_set)))), eig=TRUE, k=2)
      gene_rep_res = data.frame(names=rownames(gene_rep_fit$points), M1=gene_rep_fit$points[,1], M2=gene_rep_fit$points[,2])
      gene_rep_mdsplot = ggplot(data=gene_rep_res)
      gene_rep_mdsplot = gene_rep_mdsplot + theme_bw()  # Add theme black and white.
      gene_rep_mdsplot = gene_rep_mdsplot + geom_point(aes(x = M1, y = M2, color = names))  # Draw points in any case.
      if (replicate_number <= 24) {
        # Only add text for a sensible number of replicates i.e. less than or equal to 24.
        gene_rep_mdsplot = gene_rep_mdsplot + geom_text(aes(x = M1, y = M2, label = names, color = names, size=4))
      }
      # Arrange a maximum of 24 replicates in each guide column.
      gene_rep_mdsplot = gene_rep_mdsplot + guides(col=guide_legend(ncol=ceiling(x=replicate_number / 24)))
      ggsave(filename=plot_path, plot=gene_rep_mdsplot)
      plot_path = file.path(output_directory, paste0(prefix, "_genes_mds.pdf"))
      ggsave(filename=plot_path, plot=gene_rep_mdsplot)      
      rm(gene_rep_fit, gene_rep_res, gene_rep_mdsplot)
    }
#  } else {
#    message("Skipping Multidimensional Scaling Plot on genes in lack of replicates")
#  }
}

# Create a Principal Component Analysis (PCA) Plot on Genes.
# TODO: Add also other principal components.

plot_path = file.path(output_directory, paste0(prefix, "_genes_pca.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a Principal Component Analysis Plot (PCA) on Genes")
} else {
  message("Creating a Principal Component Analysis Plot (PCA) on Genes")
  PCAplot(object=genes(object=cuff_set), x="PC1", y="PC2", replicates=TRUE)
  ggsave(filename=plot_path)
}

# Finished plotting.

message("Finished plotting.")
rm(plot_path, active_device)

# TODO: Process gene and isoform sets and lists.
# Feature-level information can be accessed directly from a CuffData object using the
# fpkm, repFpkm, count, diffData, or annotation methods

# Create a new, simpler data frame for gene annotation ready for merging with other data frames.
# The class_code, nearest_ref_id, length and coverage fields seem to be empty by design.

message("Create gene annotation information")
gene_annotation = annotation(object=genes(object=cuff_set))
gene_frame = data.frame(gene_annotation$gene_id,
                        gene_annotation$gene_short_name,
                        gene_annotation$locus)
colnames(gene_frame)[1] = "gene_id"
colnames(gene_frame)[2] = "gene_short_name"
colnames(gene_frame)[3] = "locus"
rm(gene_annotation)

message("Create isoform annotation information")
isoform_annotation = annotation(object=isoforms(object=cuff_set))
isoform_frame = data.frame(isoform_annotation$isoform_id,
                           isoform_annotation$gene_short_name,
                           isoform_annotation$nearest_ref_id,
                           isoform_annotation$locus,
                           isoform_annotation$length)
colnames(isoform_frame)[1] = "isoform_id"
colnames(isoform_frame)[2] = "gene_short_name"
colnames(isoform_frame)[3] = "nearest_ref_id"
colnames(isoform_frame)[4] = "locus"
colnames(isoform_frame)[5] = "length"
rm(isoform_annotation)

# Split the large and unwieldy differential data tables into pairwise comparisons.-

for (i in 1:length(sample_pairs[1,])) {
  frame_path = file.path(output_directory, paste(prefix, sample_pairs[1, i], sample_pairs[2, i], "genes_diff.tsv", sep="_"))
  if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
    message(paste("Skipping a differential data frame on Genes for", sample_pairs[1, i], "versus", sample_pairs[2, i]))
  } else {
    message(paste("Creating a differential data frame on Genes for", sample_pairs[1, i], "versus", sample_pairs[2, i]))
    # The diffData function allows merging in feature annotation, but that includes some empty columns.
    # Therefore, merge with the gene_frame.
    # Unfortunately, the gene_id column seems to apper twice, which interferes with merge, so that the first column needs removing.
    gene_diff = diffData(object=genes(object=cuff_set), x=sample_pairs[1, i], y=sample_pairs[2, i], features=FALSE)
    gene_diff = gene_diff[, -1]
    gene_merge = merge(x=gene_frame, y=gene_diff, by.x="gene_id", by.y="gene_id", all=TRUE, sort=TRUE)
    write.table(x=gene_merge, file=frame_path, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    rm(gene_diff, gene_merge)
  }
}

for (i in 1:length(sample_pairs[1,])) {
  frame_path = file.path(output_directory, paste(prefix, sample_pairs[1, i], sample_pairs[2, i], "isoforms_diff.tsv", sep="_"))
  if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
    message(paste("Skipping a differential data frame on Isoforms for", sample_pairs[1, i], "versus", sample_pairs[2, i]))
  } else {
    message(paste("Creating a differential data frame on Isoforms for", sample_pairs[1, i], "versus", sample_pairs[2, i]))
    # The diffData function allows merging in feature annotation, but that includes some empty columns.
    # Therefore, merge with the isoform_frame.
    # Unfortunately, the isoform_id column seems to apper twice, which interferes with merge, so that the first column needs removing.
    isoform_diff = diffData(object=isoforms(object=cuff_set), x=sample_pairs[1, i], y=sample_pairs[2, i], features=FALSE)
    isoform_diff = isoform_diff[, -1]
    isoform_merge = merge(x=isoform_frame, y=isoform_diff, by.x="isoform_id", by.y="isoform_id", all=TRUE, sort=TRUE)
    write.table(x=isoform_merge, file=frame_path, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    rm(isoform_diff, isoform_merge)
  }
}

frame_path = file.path(output_directory, paste0(prefix, "_genes_fpkm_replicates.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
  message("Skipping a matrix of FPKM per replicates for genes")
} else {
  message("Creating a matrix of FPKM per replicates for genes")
  gene_rep_fpkm_matrix = repFpkmMatrix(object=genes(object=cuff_set))
  gene_merge = merge(x=gene_frame, y=gene_rep_fpkm_matrix, by.x="gene_id", by.y="row.names", all=TRUE, sort=TRUE)
  write.table(x=gene_merge, file=frame_path, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
  rm(gene_merge)
}
# repCountMatrix would be a similar function for normalised count rather than FPKM values.

frame_path = file.path(output_directory, paste0(prefix, "_isoforms_fpkm_replicates.tsv"))
if (file.exists(frame_path) && file.info(frame_path)$size > 0) {
  message("Skipping a matrix of FPKM per replicates for isoforms")
} else {
  message("Creating a matrix of FPKM per replicates for isoforms")
  isoform_rep_fpkm_matrix = repFpkmMatrix(object=isoforms(object=cuff_set))
  isoform_merge = merge(x=isoform_frame, y=isoform_rep_fpkm_matrix, by.x="isoform_id", by.y="row.names", all=TRUE, sort=TRUE)
  write.table(x=isoform_merge, file=frame_path, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
  rm(isoform_merge)
}

# TODO: Deal with sets of significantly regulated genes, including a sigMatrix ggplot plot.
plot_path = file.path(output_directory, paste0(prefix, "_genes_significance_matrix.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a significance matrix plot on Genes")
} else {
  message("Creating a a significance matrix plot on Genes")
  sigMatrix(object=cuff_set, level="genes")
  ggsave(filename=plot_path)
}

plot_path = file.path(output_directory, paste0(prefix, "_isoforms_significance_matrix.png"))
if (file.exists(plot_path) && file.info(plot_path)$size > 0) {
  message("Skipping a significance matrix plot on Isoforms")
} else {
  message("Creating a a significance matrix plot on Isoforms")
  sigMatrix(object=cuff_set, level="isoforms")
  ggsave(filename=plot_path)
}

rm(frame_path, plot_path, gene_frame, isoform_frame)

# Get a CuffFeatureSet or CuffGeneSet of significant genes.

# TODO: Comment-out for the moment, as the data is currently not used.
# significant_gene_ids = getSig(object=cuff_set, level='genes')
# significant_genes = getGenes(object=cuff_set, geneIdList=significant_gene_ids)

# The csHeatmap plot does not seem to be a sensible option for larger sets of significant genes.
# for (i in 1:length(sample_pairs[1,])) {
  
  # significant_genes_diff = 
# }

# Finally, create comparison-specific relative symbolic links for cuffdiff results in the
# rnaseq_process_cuffdiff_* directory to avoid identical file names between comparisons.

message("Started creating symbolic links to cuffdiff results")

link_path = file.path(output_directory, paste0(prefix, "_cds_exp_diff.tsv"))
if (! file.exists(link_path)) {
  if (! file.symlink(from=file.path("..", cuffdiff_directory, "cds_exp.diff"), to=link_path)) {
    message("Encountered an error linking the cds_exp.diff file.")
  }
}

link_path = file.path(output_directory, paste0(prefix, "_genes_exp_diff.tsv"))
if (! file.exists(link_path)) {
  if (! file.symlink(from=file.path("..", cuffdiff_directory, "gene_exp.diff"), to=link_path)) {
    message("Encountered an error linking the gene_exp.diff file.")
  }
}

link_path = file.path(output_directory, paste0(prefix, "_isoforms_exp_diff.tsv"))
if (! file.exists(link_path)) {
  if (! file.symlink(from=file.path("..", cuffdiff_directory, "isoform_exp.diff"), to=link_path)) {
    message("Encountered an error linking the isoform_exp.diff file.")
  }
}

link_path = file.path(output_directory, paste0(prefix, "_promoters_diff.tsv"))
if (! file.exists(link_path)) {
  if (! file.symlink(from=file.path("..", cuffdiff_directory, "promoters.diff"), to=link_path)) {
    message("Encountered an error linking the promoters.diff file.")
  }
}

link_path = file.path(output_directory, paste0(prefix, "_splicing_diff.tsv"))
if (! file.exists(link_path)) {
  if (! file.symlink(from=file.path("..", cuffdiff_directory, "splicing.diff"), to=link_path)) {
    message("Encountered an error linking the splicing.diff file.")
  }
}

link_path = file.path(output_directory, paste0(prefix, "_tss_group_exp_diff.tsv"))
if (! file.exists(link_path)) {
  if (! file.symlink(from=file.path("..", cuffdiff_directory, "tss_group_exp.diff"), to=link_path)) {
    message("Encountered an error linking the tss_group_exp.diff file.")
  }
}

message("Finished creating symbolic links to cuffdiff results")
rm(link_path)

rm(cuff_set, output_directory, sample_pairs, prefix, i)
message("All done")
ls()

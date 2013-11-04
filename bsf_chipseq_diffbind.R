#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(package="optparse"))
suppressPackageStartupMessages(library(package="DiffBind"))

option_list = list(
  make_option(opt_str=c("-v", "--verbose"), action="store_true",
              default=TRUE,
              help="Print extra output [default]"),
  make_option(opt_str=c("-q", "--quietly"), action="store_false",
              default=FALSE,
              dest="verbose",
              help="Print little output"),
  make_option(opt_str=c("-f", "--factor"),
              dest="factor",
              help="ChIP factor"),
  make_option(opt_str=c("-s", "--sample_annotation"),
              dest="sample_annotation",
              help="Sample annotation sheet"),
  make_option(opt_str=c("-g", "--genome_directory"),
              dest="genome_directory",
              help="Genome-specific analysis directory")
)

# Get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,

opt = parse_args(OptionParser(option_list=option_list))

prefix = paste("chipseq", "diffbind", opt$factor, sep="_")

output_directory = file.path(opt$genome_directory, prefix, fsep=.Platform$file.sep)

if (!file.exists(output_directory)) {
  dir.create(output_directory, showWarnings=TRUE, recursive=FALSE)
}

setwd(output_directory)

message("Create a DiffBind dba object.")
diffbind_dba = dba(sampleSheet=opt$sample_annotation, bCorPlot=FALSE)

message("Plotting correlation heatmap based on peak caller score data.")
file_name = paste(prefix, "_correlation_peak_caller_score.png", sep="")
png(filename=file.path(output_directory, file_name, fsep=.Platform$file.sep))
return_value = dba.plotHeatmap(DBA=diffbind_dba, margin=25)
active_device = dev.off()

message("Counting reads.")
diffbind_dba = dba.count(DBA=diffbind_dba, bCorPlot=FALSE)

message("Plotting correlation heatmap based on read counts.")
file_name = paste(prefix, "_correlation_read_counts.png", sep="")
png(filename=file.path(output_directory, file_name, fsep=.Platform$file.sep))
return_value = dba.plotHeatmap(DBA=diffbind_dba, margin=25)
active_device = dev.off()

message("Establishing contrast by tissue.")
diffbind_dba = dba.contrast(DBA=diffbind_dba, categories=DBA_TISSUE, minMembers=2)

message("Running differential binding affinity analysis.")
diffbind_dba = dba.analyze(DBA=diffbind_dba, bCorPlot=FALSE)

message("Plotting correlation heatmap based on differential binding affinity analysis.")
file_name = paste(prefix, "_correlation_analysis.png", sep="")
png(filename=file.path(output_directory, file_name, fsep=.Platform$file.sep))
return_value = dba.plotHeatmap(DBA=diffbind_dba, margin=25)
active_device = dev.off()

process_per_contrast = function(contrast, group1, group2) {
  
  # Process per row of a contrasts data frame obtained via dba.show()
  # contrast the row.names() string of the data frame indicating the contrast number
  # group1 Group1 value
  # group2 Group2 value
  
  message(sprintf("Write differentially bound sites for factor %s and contrast %s versus %s to disk.",
                  opt$factor, group1, group2))
  # TODO: The report function is quite peculiar in that it insists on a prefix DBA_
  # for the file name.
  diffbind_report = dba.report(DBA=diffbind_dba,
                               contrast=as.integer(contrast),
                               bNormalized=TRUE,
                               bCalled=TRUE,
                               bCounts=TRUE,
                               bCalledDetail=TRUE,
                               file=sprintf("%s_report_%s__%s", opt$factor, group1, group2),
                               DataType=DBA_DATA_FRAME)
  
  message(sprintf("Creating MA plot for factor %s and contrast %s versus %s.",
                  opt$factor, group1, group2))
  file_name = sprintf("%s_ma_plot_%s__%s.png", prefix, group1, group2)
  png(filename=file.path(output_directory, file_name, fsep=.Platform$file.sep))
  dba.plotMA(DBA=diffbind_dba, bNormalized=TRUE, bXY=FALSE, contrast=as.integer(contrast))
  active_device = dev.off()
  
  message(sprintf("Creating scatter plot for factor %s and contrast %s versus %s.",
                  opt$factor, group1, group2))
  file_name = sprintf("%s_scatter_plot_%s__%s.png", prefix, group1, group2)
  png(filename=file.path(output_directory, file_name, fsep=.Platform$file.sep))
  dba.plotMA(DBA=diffbind_dba, bNormalized=TRUE, bXY=TRUE, contrast=as.integer(contrast))
  active_device = dev.off()
  
  message(sprintf("Creating PCA plot for factor %s and contrast %s versus %s.",
                  opt$factor, group1, group2))
  file_name = sprintf("%s_pca_plot_%s__%s.png", prefix, group1, group2)
  png(filename=file.path(output_directory, file_name, fsep=.Platform$file.sep))
  dba.plotPCA(DBA=diffbind_dba, attributes=DBA_CONDITION)
  active_device = dev.off()
  
  message(sprintf("Creating Box plot for factor %s and contrast %s versus %s.",
                  opt$factor, group1, group2))
  file_name = sprintf("%s_box_plot_%s__%s.png", prefix, group1, group2)
  png(filename=file.path(output_directory, file_name, fsep=.Platform$file.sep))
  diffbind_pvals = dba.plotBox(DBA=diffbind_dba, bNormalized=TRUE)
  active_device = dev.off()
  
  return(TRUE)
}

# Get a data frame with all contrasts to apply the above function to each row.
contrasts = dba.show(DBA=diffbind_dba, bContrasts=TRUE)

# Replace '!' characters with 'not_'.
contrasts$Group1 = gsub(pattern="!", replacement="not_", x=contrasts$Group1)
contrasts$Group2 = gsub(pattern="!", replacement="not_", x=contrasts$Group2)

# Write the contrasts to disk.
# write.table(x=contrasts,
#             file=paste(prefix, "contrasts.txt", sep="_"),
#             col.names=TRUE, row.names=TRUE, sep="\t")

write.csv(x=contrasts,
          file=paste(prefix, "contrasts.csv", sep="_"),
          col.names=TRUE, row.names=FALSE)

# Apply the function to each row and discard the return value.
return_value = mapply(FUN=process_per_contrast, row.names(contrasts), contrasts$Group1, contrasts$Group2)
rm(return_value)
rm(file_name)

# dba.overlap(DBA=diffbind_dba, mode=DBA_OLAP_RATE)

# message("Creating a Venn diagram")
# file_name = paste(prefix, "_box_plot.png", sep="")
# png(filename=file.path(output_directory, file_name, fsep=.Platform$file.sep))
# file_name = dba.plotVenn(DBA=diffbind_dba)
# active_device = dev.off()

message("Save workspace image for manual post-processing.")
save.image()

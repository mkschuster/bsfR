#! /usr/bin/env Rscript

suppressPackageStartupMessages(expr = library(package = "biomaRt"))
suppressPackageStartupMessages(expr = library(package = "ggplot2"))
suppressPackageStartupMessages(expr = library(package = "optparse"))

# Process Tophat and Cufflinks directories for a particular, named replicate.
#
# @param replicate_name: Replicate name
# @type replicate_name: char

process_replicate <- function(replicate_name) {
  
  if (is.null(x = replicate_name)) {
    stop("Missing replicate_name argument")    
  }
  
  message(paste0("Processing replicate ", replicate_name))
  
  # Construct replicate-specific prefixes for Cufflinks and Tophat directories.
  prefix_cufflinks <- paste("rnaseq", "cufflinks", replicate_name, sep = "_")
  prefix_tophat <- paste("rnaseq", "tophat", replicate_name, sep = "_")
  
  # Read, merge, write and delete gene tables.
  
  file_path <- file.path(prefix_cufflinks, paste(prefix_cufflinks, "genes_fpkm_tracking.tsv", sep = "_"))
  if (! (file.exists(file_path) && (file.info(file_path)$size > 0))) {
    cufflinks_genes <- read.table(file = file.path(prefix_cufflinks, "genes.fpkm_tracking"), header = TRUE)
    cufflinks_ensembl <- merge(x = ensembl_genes, y = cufflinks_genes,
                               by.x = "ensembl_gene_id", by.y = "tracking_id",
                               all.y = TRUE, sort = TRUE)
    write.table(x = cufflinks_ensembl, file = file_path, col.names = TRUE, row.names = FALSE, sep = "\t")
    rm(cufflinks_genes, cufflinks_ensembl)
  }
  rm(file_path)
  
  # Read, merge, write and delete transcript tables.
  
  file_path <- file.path(prefix_cufflinks, paste(prefix_cufflinks, "isoforms_fpkm_tracking.tsv", sep = "_"))
  if (! (file.exists(file_path) && (file.info(file_path)$size > 0))) {
    cufflinks_transcripts <- read.table(file = file.path(prefix_cufflinks, "isoforms.fpkm_tracking"), header = TRUE)
    cufflinks_ensembl <- merge(x = ensembl_transcripts, y = cufflinks_transcripts,
                               by.x = "ensembl_transcript_id", by.y = "tracking_id",
                               all.y = TRUE, sort = TRUE)
    write.table(x = cufflinks_ensembl, file = file_path, col.names = TRUE, row.names = FALSE, sep = "\t")
    rm(cufflinks_transcripts, cufflinks_ensembl)
  }
  rm(file_path)
  
  # Finally, create a replicate-specific symbolic link to the Tophat aligned BAM file and its index.
  
  file_path <- file.path("..", prefix_tophat, "accepted_hits.bam")
  link_path <- file.path(prefix_cufflinks, paste(prefix_tophat, "accepted_hits.bam", sep = "_"))
  if (! file.exists(link_path)) {
    if (! file.symlink(from = file_path, to = link_path)) {
      warning("Encountered an error linking the accepted_hits.bam file.")
    }
  }
  rm(file_path, link_path)
  
  file_path <- file.path("..", prefix_tophat, "accepted_hits.bam.bai")
  link_path <- file.path(prefix_cufflinks, paste(prefix_tophat, "accepted_hits.bam.bai", sep = "_"))
  if (! file.exists(link_path)) {
    if (! file.symlink(from = file_path, to = link_path)) {
      warning("Encountered an error linking the accepted_hits.bam.bai file.")
    }
  }
  rm(file_path, link_path)
  
  file_path <- file.path("..", prefix_tophat, "unmapped.bam")
  link_path <- file.path(prefix_cufflinks, paste(prefix_tophat, "unmapped.bam", sep = "_"))
  if (! file.exists(link_path)) {
    if (! file.symlink(from = file_path, to = link_path)) {
      warning("Encountered an error linking the unmapped.bam file.")
    }
  }
  rm(file_path, link_path)
  
  file_path <- file.path("..", prefix_tophat, "align_summary.txt")
  link_path <- file.path(prefix_cufflinks, paste(prefix_tophat, "align_summary.txt", sep = "_"))
  if (! file.exists(link_path)) {
    if (! file.symlink(from = file_path, to = link_path)) {
      warning("Encountered an error linking the align_summary.txt file.")
    }
  }
  rm(file_path, link_path)
  
  file_path <- "transcripts.gtf"
  link_path <- file.path(prefix_cufflinks, paste(prefix_cufflinks, "transcripts.gtf", sep = "_"))
  if (! file.exists(link_path)) {
    if (! file.symlink(from = file_path, to = link_path)) {
      warning("Encountered an error linking the transcripts.gtf file.")
    }
  }
  rm(file_path, link_path)
  
  file_path <- "skipped.gtf"
  link_path <- file.path(prefix_cufflinks, paste(prefix_cufflinks, "skipped.gtf", sep = "_"))
  if (! file.exists(link_path)) {
    if (! file.symlink(from = file_path, to = link_path)) {
      warning("Encountered an error linking the skipped.gtf file.")
    }
  }
  rm(file_path, link_path)
  
  rm(prefix_cufflinks, prefix_tophat)
  
  return()
}

# Parse Tophat align_summary.txt files and return a data frame.
#
# @param replicate_names: list of replicate names
# @type replicate_names: list
# @return: Data frame with alignment summary statistics
# @rtype: data.frame

process_align_summary <- function(replicate_names) {
  
  if (is.null(x = replicate_names)) {
    stop("Missing replicate_names argument")    
  }
  
  list_length = length(x = replicate_names)
  summary_frame = data.frame(
    replicate = replicate_names,
    input = integer(length = list_length),
    mapped = integer(length = list_length),
    multiple = integer(length = list_length),
    above = integer(length = list_length),
    threshold = integer(length = list_length),
    row.names = replicate_names,
    stringsAsFactors = FALSE)
  rm(list_length)
  
  for (i in 1:nrow(x = summary_frame)) {
    
    prefix_tophat <- paste("rnaseq", "tophat", row.names(x = summary_frame)[i], sep = "_")
    file_path <- file.path(prefix_tophat, "align_summary.txt")
    
    if (! file.exists(file_path)) {
      warning(paste0("Missing Tophat alignment summary file ", file_path))
      rm(prefix_tophat, file_path)
      next
    }
    
    align_summary <- readLines(con = file_path)
    
    # This is the layout of a Tophat align_summary.txt file.
    #
    #   [1] "Reads:"                                                                          
    #   [2] "          Input     :  21791622"                                                 
    #   [3] "           Mapped   :  21518402 (98.7% of input)"                                
    #   [4] "            of these:   2010462 ( 9.3%) have multiple alignments (8356 have >20)"
    #   [5] "98.7% overall read mapping rate."                                                
    
    # Parse the second line of "input" reads.
    summary_frame[i, "input"] <- as.integer(
      x = sub(
        pattern = "[[:space:]]+Input[[:space:]]+:[[:space:]]+([[:digit:]]+)",
        replacement = "\\1",
        x = align_summary[2]))
    
    # Parse the third line of "mapped" reads.
    summary_frame[i, "mapped"] <- as.integer(
      x = sub(
        pattern = "[[:space:]]+Mapped[[:space:]]+:[[:space:]]+([[:digit:]]+) .*",
        replacement = "\\1",
        x = align_summary[3]))
    
    # Get the number of "multiply" aligned reads from the fourth line.
    summary_frame[i, "multiple"] <- as.integer(
      x = sub(
        pattern = ".+:[[:space:]]+([[:digit:]]+) .*",
        replacement = "\\1",
        x = align_summary[4]))
    
    # Get the number of multiply aligned reads "above" the multiple alignment threshold.
    summary_frame[i, "above"] <- as.integer(
      x = sub(
        pattern = ".+alignments \\(([[:digit:]]+) have.+",
        replacement = "\\1",
        x = align_summary[4]))
    
    # Get the multiple alignment "threshold".
    summary_frame[i, "threshold"] <- as.integer(
      x = sub(
        pattern = ".+ >([[:digit:]]+)\\)",
        replacement = "\\1",
        x = align_summary[4]))
    
    rm(align_summary, prefix_tophat)
  }
  rm(i)
  
  return(summary_frame)
}

option_list <- list(
  make_option(opt_str = c("--verbose", "-v"),
              action = "store_true",
              default = TRUE,
              help = "Print extra output [default]"),
  make_option(opt_str = c("--quiet", "-q"),
              action = "store_false",
              default = FALSE,
              dest = "verbose",
              help = "Print little output"),
  make_option(opt_str = c("--biomart-instance"),
              default = "ENSEMBL_MART_ENSEMBL",
              dest = "biomart_instance",
              help = "BioMart instance"),
  make_option(opt_str = c("--biomart-data-set"),
              dest = "biomart_data_set",
              help = "BioMart data set"),
  make_option(opt_str = c("--biomart-host"),
              dest = "biomart_host",
              default = "www.ensembl.org",
              help = "BioMart host"),
  make_option(opt_str = c("--biomart-port"),
              default = 80,
              dest = "biomart_port",
              help = "BioMart port",
              type = "integer"),
  make_option(opt_str = c("--biomart-path"),
              dest = "biomart_path",
              help = "BioMart path"),
  make_option(opt_str = c("--sample"),
              dest = "sample",
              help = "Sample name"),
  make_option(opt_str = c("--plot-width"),
              default = 7,
              dest = "plot_width",
              help = "Plot width in inches"),
  make_option(opt_str = c("--plot-height"),
              default = 7,
              dest = "plot_height",
              help = "Plot height in inches")
)

# Get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,

opt <- parse_args(object = OptionParser(option_list = option_list))

if (is.null(x = opt$biomart_data_set)) {
  stop("Missing --data_set option")
}

# Connect to the Ensembl BioMart.

ensembl_mart <- useMart(
  biomart = opt$biomart_instance,
  dataset = opt$biomart_data_set,
  host = opt$biomart_host,
  port = opt$biomart_port)

message("Loading attribute data from BioMart")

ensembl_attributes <- listAttributes(
  mart = ensembl_mart,
  page = "feature_page",
  what = c("name", "description", "page"))

# Get Ensembl Gene information.
# From Ensembl version e75, the attributes "external_gene_id" and "external_gene_db" are called
# "external_gene_name" and "external_gene_source", respectively.

if (any(ensembl_attributes$name == "external_gene_id")) {
  # Pre e75.
  ensembl_gene_attributes <- c(
    "ensembl_gene_id",
    "external_gene_id",
    "external_gene_db",
    "gene_biotype")
} else if (any(ensembl_attributes$name == "external_gene_name")) {
  # Post e75.
  ensembl_gene_attributes <- c(
    "ensembl_gene_id",
    "external_gene_name",
    "external_gene_db",
    "gene_biotype")
} else {
  stop("Neither external_gene_id nor external_gene_name are available as BioMart attributes.")
}

message("Loading gene data from BioMart")

ensembl_genes <- getBM(attributes = ensembl_gene_attributes, mart = ensembl_mart)
rm(ensembl_gene_attributes)

# Get Ensembl Transcript information.

if (any(ensembl_attributes$name == "external_gene_id")) {
  # Pre e75.
  ensembl_transcript_attributes <- c(
    "ensembl_transcript_id",
    "external_transcript_id",
    "transcript_db_name",
    "transcript_biotype",
    "ensembl_gene_id",
    "external_gene_id",
    "external_gene_db",
    "gene_biotype")
} else if (any(ensembl_attributes$name == "external_gene_name")) {
  # Post e75.
  ensembl_transcript_attributes <- c(
    "ensembl_transcript_id",
    "external_transcript_name",
    "transcript_db_name",
    "transcript_biotype",
    "ensembl_gene_id",
    "external_gene_name",
    "external_gene_db",
    "gene_biotype")
} else {
  stop("Neither external_gene_id nor external_gene_name are available as BioMart attributes.")
}

message("Loading transcript data from BioMart")
ensembl_transcripts <- getBM(attributes = ensembl_transcript_attributes, mart = ensembl_mart)
rm(ensembl_transcript_attributes)

# Destroy and discconnect the Ensembl BioMart connection already here.

rm(ensembl_mart)

if (is.null(x = opt$sample)) {
  
  # If a --sample option was not provided, process all "rnaseq_cufflinks_*" directories in the
  # current working directory. List all rnaseq_cufflinks directories via their common prefix and
  # parse the sample (or replicate) name simply by removing the prefix.
  
  replicate_names <- sub(
    pattern = "^rnaseq_cufflinks_",
    replacement = "",
    x = grep(pattern = '^rnaseq_cufflinks_', x = list.dirs(full.names = FALSE, recursive = FALSE), value = TRUE))
  
  for (replicate_name in replicate_names) {
    process_replicate(replicate_name = replicate_name)
  }
  rm(replicate_name)
  
  summary_frame <- process_align_summary(replicate_names)
  
  # Write the alignment summary frame to disk and create a chart.
  
  file_path <- "rnaseq_tophat_alignment_summary.tsv"
  # if (! (file.exists(file_path) && (file.info(file_path)$size > 0))) {
  write.table(x = summary_frame, file = file_path, col.names = TRUE, row.names = FALSE, sep = "\t")
  # }
  rm(file_path)
  
  ggplot_object <- ggplot(data = summary_frame)
  ggplot_object <- ggplot_object + geom_point(mapping = aes(x = replicate, y = mapped / input, color = replicate))
  ggplot_object <- ggplot_object + guides(col = guide_legend(ncol = ceiling(x = (nrow(x = summary_frame) / 24))))
  ggplot_object <- ggplot_object + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  ggplot_object <- ggplot_object + theme(legend.text = element_text(size = rel(x = 0.7)))  # Reduce the legend text
  
  ggsave(filename = "rnaseq_tophat_alignment_summary.png",
         plot = ggplot_object, width = opt$plot_width, height = opt$plot_height)
  ggsave(filename = "rnaseq_tophat_alignment_summary.pdf",
         plot = ggplot_object, width = opt$plot_width, height = opt$plot_height)
  
  rm(ggplot_object, summary_frame, replicate_names)
} else {
  # Process the single sample (replicate).
  process_replicate(replicate_name = opt$sample)
}

rm(ensembl_attributes, ensembl_transcripts, ensembl_genes,
   option_list, opt, process_replicate)

message("All done")
# TODO: For debugging and development.
ls()
save.image()

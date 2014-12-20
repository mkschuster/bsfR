#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(package="optparse"))
suppressPackageStartupMessages(library(package="biomaRt"))

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
  make_option(opt_str=c("--biomart-instance"),
              default="ENSEMBL_MART_ENSEMBL",
              dest="biomart_instance",
              help="BioMart instance"),
  make_option(opt_str=c("--biomart-data-set"),
              dest="biomart_data_set",
              help="BioMart data set"),
  make_option(opt_str=c("--biomart-host"),
              dest="biomart_host",
              default="www.ensembl.org",
              help="BioMart host"),
  make_option(opt_str=c("--biomart-port"),
              default=80,
              dest="biomart_port",
              help="BioMart port",
              type="integer"),
  make_option(opt_str=c("--biomart-path"),
              dest="biomart_path",
              help="BioMart path"),
  make_option(opt_str=c("--sample"),
              dest="sample",
              help="Sample name"),
  make_option(opt_str=c("--genome-directory"),
              dest="genome_directory",
              help="Genome-specific analysis directory")
)

# Get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,

opt = parse_args(object=OptionParser(option_list=option_list))

if (is.null(x=opt$biomart_data_set)) {
  stop("Missing --data_set option")
}

if (is.null(x=opt$sample)) {
  stop("Missing --sample option")
}

# Construct a sample-specific prefix.

prefix = paste("rnaseq", "cufflinks", opt$sample, sep="_")

# Connect to the Ensembl BioMart.

ensembl_mart = useMart(biomart=opt$biomart_instance,
                       dataset=opt$biomart_data_set,
                       host=opt$biomart_host,
                       port=opt$biomart_port)

# Read the Cufflinks genes.fpkm_tracking file into a data_frame.

cufflinks_genes = read.table(file=file.path(opt$genome_directory, prefix,
                                            "genes.fpkm_tracking"),
                             header=TRUE)

message("Load gene data from BioMart")

# Get Ensembl Gene stable identifier, gene symbol and the associated database.
# The following attributes are not needed for the moment.
# "chromosome_name", "start_position", "end_position", "strand"

# From Ensembl version e75, the attributes "external_gene_id" and "external_gene_db" are called
# "external_gene_name" and "external_gene_source", respectively.

ensembl_genes = getBM(attributes=c("ensembl_gene_id", "external_gene_name", "external_gene_db"),
                      mart=ensembl_mart)

# Merge the data frames ...

cufflinks_ensembl = merge(x=ensembl_genes, y=cufflinks_genes,
                          by.x="ensembl_gene_id", by.y="tracking_id",
                          all.y=TRUE, sort=TRUE)

# Write back to disk.

write.table(x=cufflinks_ensembl,
            file=file.path(opt$genome_directory, prefix,
                           paste(prefix, "genes_fpkm_tracking.tsv", sep="_")),
            col.names=TRUE, row.names=FALSE, sep="\t")

rm(cufflinks_genes, ensembl_genes, cufflinks_ensembl)

# Read the Cufflinks isoforms.fpkm_tracking file into a data_frame.

cufflinks_transcripts = read.table(file=file.path(opt$genome_directory, prefix,
                                                  "isoforms.fpkm_tracking"),
                                   header=TRUE)

message("Load transcript data from BioMart")

# Get Ensembl Transcript stable identifier, gene symbol and the associated database.
# The following attributes are not needed
# "chromosome_name", "start_position", "end_position", "strand"

ensembl_transcripts = getBM(attributes=c("ensembl_transcript_id", "external_gene_name", "external_gene_db"), 
                            mart=ensembl_mart)

# Merge the data frames ...

cufflinks_ensembl = merge(x=ensembl_transcripts, y=cufflinks_transcripts,
                          by.x="ensembl_transcript_id", by.y="tracking_id",
                          all.y=TRUE, sort=TRUE)

# Write back to disk.

write.table(x=cufflinks_ensembl,
            file=file.path(opt$genome_directory, prefix,
                           paste(prefix, "isoforms_fpkm_tracking.tsv", sep="_")),
            col.names=TRUE, row.names=FALSE, sep="\t")

rm(cufflinks_transcripts, ensembl_transcripts, cufflinks_ensembl)

# Finally, create a replicate-specific symbolic link to the Tophat aligned BAM file and its index.

file_path = file.path(opt$genome_directory,
                      paste("rnaseq", "tophat", opt$sample, sep="_"),
                      "accepted_hits.bam")
link_path = file.path(opt$genome_directory,
                      prefix,
                      paste("rnaseq", "tophat", opt$sample, "accepted_hits.bam", sep="_"))
if (! file.exists(link_path)) {
  if (! file.symlink(from=file_path, to=link_path)) {
    warning("Encountered an error linking the accepted_hits.bam file.")
  }
}
rm(file_path, link_path)

file_path = file.path(opt$genome_directory,
          paste("rnaseq", "tophat", opt$sample, sep="_"),
          "accepted_hits.bam.bai")
link_path = file.path(opt$genome_directory,
                      prefix,
                      paste("rnaseq", "tophat", opt$sample, "accepted_hits.bam.bai", sep="_"))
if (! file.exists(link_path)) {
  if (! file.symlink(from=file_path, to=link_path)) {
    warning("Encountered an error linking the accepted_hits.bam.bai file.")
  }
}
rm(file_path, link_path)

file_path = file.path(opt$genome_directory,
                      paste("rnaseq", "tophat", opt$sample, sep="_"),
                      "unmapped.bam")
link_path = file.path(opt$genome_directory,
                      prefix,
                      paste("rnaseq", "tophat", opt$sample, "unmapped.bam", sep="_"))
if (! file.exists(link_path)) {
  if (! file.symlink(from=file_path, to=link_path)) {
    warning("Encountered an error linking the unmapped.bam file.")
  }
}
rm(file_path, link_path)

rm(opt, option_list, prefix, ensembl_mart)
message("All done")
ls()

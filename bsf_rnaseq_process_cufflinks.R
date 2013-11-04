#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(package="optparse"))
suppressPackageStartupMessages(library(package="biomaRt"))

option_list = list(
  make_option(opt_str=c("-v", "--verbose"), action="store_true",
              default=TRUE,
              help="Print extra output [default]"),
  make_option(opt_str=c("-q", "--quietly"), action="store_false",
              default=FALSE,
              dest="verbose",
              help="Print little output"),
  make_option(opt_str=c("-b", "--biomart"),
              dest="biomart_instance",
              default="ENSEMBL_MART_ENSEMBL",
              help="BioMart instance"),
  make_option(opt_str=c("-w", "--data_set"),
              dest="biomart_data_set",
              help="BioMart data set"),
  make_option(opt_str=c("-a", "--host"),
              dest="biomart_host",
              default="www.ensembl.org",
              help="BioMart host"),
  make_option(opt_str=c("-p", "--port"),
              dest="biomart_port",
              help="BioMart port"),
  make_option(opt_str=c("-f", "--path"),
              dest="biomart_path",
              help="BioMart path"),
  make_option(opt_str=c("-s", "--sample"),
              dest="sample",
              help="Sample name"),
  make_option(opt_str=c("-g", "--genome_directory"),
              dest="genome_directory",
              help="Genome-specific analysis directory")
)

# Get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,

opt = parse_args(OptionParser(option_list=option_list))

# TODO: Check that there is a value for the --data_set option.
# TODO: Check that there is a value for the --sample option.

# Construct a sample-specific prefix.

prefix = paste("rnaseq", "cufflinks", opt$sample, sep="_")

# Connect to the Ensembl BioMart.

# port=opt$biomart_port
# path=opt$biomart_path

ensembl_mart = useMart(biomart=opt$biomart_instance,
                       dataset=opt$biomart_data_set,
                       host=opt$biomart_host)

# Read the Cufflinks file into a data_frame.

cufflinks_genes = read.table(file=file.path(opt$genome_directory, prefix,
                                            "genes.fpkm_tracking",
                                            fsep=.Platform$file.sep),
                             header=TRUE)

message("Load gene data from BioMart")

# Get Ensembl Gene stable identifier, gene symbol and the associated database.
# The following attributes are not needed for the moment.
# "chromosome_name", "start_position", "end_position", "strand"

ensembl_genes = getBM(attributes=c("ensembl_gene_id", "external_gene_id", "external_gene_db"),
                      mart=ensembl_mart)

# Merge the data frames ...

cufflinks_ensembl = merge(x=ensembl_genes, y=cufflinks_genes,
                          by.x="ensembl_gene_id", by.y="tracking_id",
                          all.y=TRUE, sort=TRUE)

# Write back to disk.

write.table(x=cufflinks_ensembl,
            file=file.path(opt$genome_directory, prefix,
                           paste(prefix, "genes_fpkm_tracking.txt", sep="_"),
                           fsep=.Platform$file.sep),
            col.names=TRUE, row.names=FALSE, sep="\t")

rm(ensembl_genes, cufflinks_genes, cufflinks_ensembl)

# Read the Cufflinks file into a data_frame.

cufflinks_transcripts = read.table(file=file.path(opt$genome_directory, prefix,
                                                  "isoforms.fpkm_tracking",
                                                  fsep=.Platform$file.sep),
                                   header=TRUE)

message("Load transcript data from BioMart")

# Get Ensembl Transcript stable identifier, gene symbol and the associated database.
# The following attributes are not needed
# "chromosome_name", "start_position", "end_position", "strand"

ensembl_transcripts = getBM(attributes=c("ensembl_transcript_id", "external_gene_id", "external_gene_db"), 
                      mart=ensembl_mart)

# Merge the data frames ...

cufflinks_ensembl = merge(x=ensembl_transcripts, y=cufflinks_transcripts,
                          by.x="ensembl_transcript_id", by.y="tracking_id",
                          all.y=TRUE, sort=TRUE)

# Write back to disk.

write.table(x=cufflinks_ensembl,
            file=file.path(opt$genome_directory, prefix,
                           paste(prefix, "isoforms_fpkm_tracking.txt", sep="_"),
                           fsep=.Platform$file.sep),
            col.names=TRUE, row.names=FALSE, sep="\t")

rm(ensembl_transcripts, cufflinks_transcripts, cufflinks_ensembl)

#! /usr/bin/env Rscript

# Sort, index and idxstats BAM files
library(optparse)
library(TOP)

# ================ Process the command-line arguments ================

option_list <- list(
  make_option("--bam_file", action="store", default=NULL, type='character',
              help="BAM file name [default: %default]"),
  make_option("--sorted_bam_file", action="store", default=NULL, type='double',
              help="Sorted BAM file name if sort=TRUE [default: %default]"),
  make_option("--sort", action="store_true", default=FALSE,
              help="If TRUE, sort the BAM file [default: %default]"),
  make_option("--index", action="store_true", default=TRUE,
              help="If TRUE, index the BAM file [default: %default]"),
  make_option("--idxstats", action="store_true", default=TRUE,
              help="If TRUE, get the idxstats of the BAM file [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))
bam_file            <- opt$bam_file
sorted_bam_file     <- opt$sorted_bam_file
sort                <- opt$sort
bfile               <- opt$index
idxstats            <- opt$idxstats

# ================ Sort, index and idxstats BAM files ================
sort_index_idxstats_bam(bam_file,
                        sorted_bam_file = sorted_bam_file,
                        sort=FALSE, index=TRUE, idxstats=TRUE)


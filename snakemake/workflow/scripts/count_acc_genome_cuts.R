#! /usr/bin/env Rscript

# Count DNase- or ATAC-seq genome-wide cleavage
library(optparse)
library(TOP)
options(scipen=999) # suppress scientific notations

# ================ Process the command-line arguments ================
option_list <- list(
  make_option("--bam", action="store", default=NULL, type="character",
              help="DNase-seq or ATAC-seq alignment BAM filename."),
  make_option("--chrom_size", action="store", default=NULL, type="character",
              help="Filename listing the size for each chromosome."),
  make_option("--shift_ATAC", action="store_true", default=FALSE,
              help="Shift ATAC reads."),
  make_option("--shift_ATAC_bases", action="store", default="4,-4", type="character",
              help="Shift ATAC bases"),
  make_option("--outdir", action="store", default=NULL, type="character",
              help="Output directory."),
  make_option("--outname", action="store", default=NULL, type="character",
              help="Output filename prefix."),
  make_option("--bedtools", action="store", default="bedtools", type="character",
              help="Path to bedtools executable."),
  make_option("--bedGraphToBigWig", action="store", default="bedGraphToBigWig", type="character",
              help="Path to UCSC bedGraphToBigWig executable.")
)

opt <- parse_args(OptionParser(option_list=option_list))
bam_file               <- opt$bam
chrom_size_file        <- opt$chrom_size
shift_ATAC             <- opt$shift_ATAC
shift_ATAC_bases       <- opt$shift_ATAC_bases
outdir                 <- opt$outdir
outname                <- opt$outname
bedtools_path          <- opt$bedtools
bedGraphToBigWig_path  <- opt$bedGraphToBigWig

# ================ Count cuts along the genome ================

if(!is.null(shift_ATAC_bases)){
  shift_ATAC_bases <- as.integer(unlist(strsplit(shift_ATAC_bases, ",|[.]")))
}else{
  shift_ATAC_bases <- c(4,-4)
}

count_genome_cuts(bam_file,
                  chrom_size_file,
                  shift_ATAC=shift_ATAC,
                  shift_ATAC_bases=shift_ATAC_bases,
                  outdir=outdir,
                  outname=outname,
                  bedtools_path=bedtools_path,
                  bedGraphToBigWig_path=bedGraphToBigWig_path)


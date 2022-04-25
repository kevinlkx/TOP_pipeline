#! /usr/bin/env Rscript

# Get candidate TF binding sites using FIMO motif match result
library(optparse)
library(TOP)
library(GenomicRanges)
library(data.table)
options(scipen=999) # suppress scientific notations

# ================ Process the command-line arguments ================

option_list <- list(
  make_option("--fimo", action="store", default=NULL, type='character',
              help="Filename of FIMO result"),
  make_option("--out", action="store", default=NULL, type='character',
              help="Output filename"),
  make_option("--flank", action="store", default=100, type='integer',
              help="Size (bp) of flanking region on each side of motif [default: %default]"),
  make_option("--thresh_pValue", action="store", default=1e-5, type='double',
              help="FIMO p-value threshold [default: %default]"),
  make_option("--thresh_pwmscore", action="store", default=0, type='double',
              help="FIMO PMW score threshold [default: %default]"),
  make_option("--blacklist", action="store", default=NULL, type='character',
              help="Filename of the blacklist regions."),
  make_option("--mapability", action="store", default=NULL, type='character',
              help="Filename of the mapability reference file in bigWig format."),
  make_option("--bigWigAverageOverBed", action="store", default="bigWigAverageOverBed", type='character',
              help="Path to bigWigAverageOverBed executable. Only needed for computing mapability.")
  )

opt <- parse_args(OptionParser(option_list=option_list))
fimo_file                  <- opt$fimo
out_file                   <- opt$out
flank                      <- opt$flank
thresh_pValue              <- opt$thresh_pValue
thresh_pwmscore            <- opt$thresh_pwmscore
blacklist_file             <- opt$blacklist
mapability_file            <- opt$mapability
bigWigAverageOverBed_path  <- opt$bigWigAverageOverBed

# ================ Get candidate sites ================
if( !file.exists(fimo_file) || file.size(fimo_file) == 0 ){
  cat(paste(fimo_file, 'file does not exist or is empty.
            Probably no motif matches were found.\n'))
  res <- file.create(out_file) # create an empty file for snakemake downstream pipeline.

}else{

  sites.df <- process_candidate_sites(fimo_file,
                                      flank = flank,
                                      thresh_pValue = thresh_pValue,
                                      thresh_pwmscore = thresh_pwmscore,
                                      blacklist_file = blacklist_file,
                                      mapability_file = mapability_file,
                                      thresh_mapability=0.8,
                                      bigWigAverageOverBed_path=bigWigAverageOverBed_path)

  if(!dir.exists(dirname(out_file))){
    dir.create(dirname(out_file), showWarnings = F, recursive = T)
  }

  cat(nrow(sites.df), 'candidate sites. \n')
  if(nrow(sites.df) > 0){
    fwrite(sites.df, out_file, sep = '\t')
    cat(nrow(sites.df), 'candidate sites written in', out_file, '\n')
  }else{
    res <- file.create(out_file)
  }
}


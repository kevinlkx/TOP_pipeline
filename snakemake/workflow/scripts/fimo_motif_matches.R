#! /usr/bin/env Rscript

# Scan for motif matches using FIMO
library(optparse)
library(TOP)
options(scipen=999) # suppress scientific notations

# ================ Process the command-line arguments ================

option_list <- list(
  make_option("--motif_file", action="store", default=NULL, type='character',
              help="Filename of motif file"),
  make_option("--sequence_file", action="store", default=NULL, type='character',
              help="Filename of FASTA sequence file"),
  make_option("--outname", action="store", default="fimo.txt", type='character',
              help="Output file name [default: %default]"),
  make_option("--outdir", action="store", default="./", type='character',
              help="Output directory [default: %default]"),
  make_option("--thresh_pValue", action="store", default=1e-5, type='double',
              help="FIMO p-value threshold [default: %default]"),
  make_option("--background", action="store", default='default', type='character',
              help="background model [default: %default]"),
  make_option("--skip_matched_sequence", action="store_true", default=TRUE,
              help="If TRUE, skip matched sequence[default: %default]"),
  make_option("--max_strand", action="store_true", default=FALSE,
              help="If TRUE, only keep the strand with the better motif match [default: %default]"),
  make_option("--options", action="store", default="", type='character',
              help="Other options for FIMO"),
  make_option("--fimo_path", action="store", default="fimo", type='character',
              help="Path to FIMO executable.")
  )

opt <- parse_args(OptionParser(option_list=option_list))
motif_file                 <- opt$motif_file
sequence_file              <- opt$sequence_file
outname                    <- opt$outname
outdir                     <- opt$outdir
thresh_pValue              <- opt$thresh_pValue
background                 <- opt$background
skip_matched_sequence      <- opt$skip_matched_sequence
max_strand                 <- opt$max_strand
options                    <- opt$options
fimo_path                  <- opt$fimo_path

# ================ Scan for motif matches using FIMO ================

if(dir.exists(outdir))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

fimo_motif_matches(motif_file,
                   sequence_file,
                   outname = outname,
                   outdir = outdir,
                   thresh_pValue = thresh_pValue,
                   background = background,
                   skip_matched_sequence = skip_matched_sequence,
                   max_strand = max_strand,
                   options = options,
                   verbosity = 2,
                   fimo_path = fimo_path)


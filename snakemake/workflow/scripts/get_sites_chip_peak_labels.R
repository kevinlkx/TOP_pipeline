#! /usr/bin/env Rscript

# get ChIP peak labels for candidate sites
library(optparse)
library(TOP)
library(data.table)

# ================ Process the command-line arguments ================
option_list <- list(
  make_option('--tf_name', action='store', default=NULL, type='character',
              help='TF name.'),
  make_option('--cell_type', action='store', default=NULL, type='character',
              help='Cell type.'),
  make_option('--sites', action='store', default=NULL, type='character',
              help='Filename for candidate sites (BED format).'),
  make_option('--metadata', action='store', default=NULL, type='character',
              help='Filename for training metadata table.'),
  make_option('--chip_peakfile_colname', action='store', default='chip_peak_file', type='character',
              help='Column name for ChIP-seq peak bed narrowPeak files.'),
  make_option('--chip_peak_dir', action='store', default=NULL, type='character',
              help='Directory for ChIP-seq peak bed narrowPeak files directory.'),
  make_option('--outdir', action='store', default=NULL, type='character',
              help='Output directory.'),
  make_option('--outname', action='store', default=NULL, type='character',
              help='Output filename prefix.')
)

opt <- parse_args(OptionParser(option_list=option_list))
tf_name                <- opt$tf_name
cell_type              <- opt$cell_type
sites_file             <- opt$sites
metadata_file          <- opt$metadata
chip_peakfile_colname  <- opt$chip_peakfile_colname
chip_peak_dir          <- opt$chip_peak_dir
outdir                 <- opt$outdir
outname                <- opt$outname

# ================ get ChIP peak labels for candidate sites =================

chip_labels_file <- file.path(outdir, paste0(outname, '.chip.peak.labels.rds'))

if( !file.exists(sites_file) || file.size(sites_file) == 0 ){
  cat(paste(sites_file, 'file does not exist or is empty. Skipped. \n'))
  res <- file.create(chip_labels_file) # create an empty file for snakemake downstream pipeline.
}else{
  sites.df <- as.data.frame(fread(sites_file))
  metadata <- as.data.frame(fread(metadata_file, sep = '\t'))
  chip_peak_sampleID <- metadata[metadata$tf_name == tf_name & metadata$cell_type == cell_type, chip_peakfile_colname]

  if( length(chip_peak_sampleID) == 0 || is.na(chip_peak_sampleID) || chip_peak_sampleID == '' ){
    cat(paste('No ChIP-seq peaks for ', tf_name, 'in', cell_type, '\n'))
    res <- file.create(chip_labels_file)
  }else{
    cat('Labeling candidate sites :', tf_name, 'in', cell_type, '\nChIP-seq peak file:',chip_peak_sampleID, '\n')

    sites_chip_labels.df <- add_chip_peak_labels_to_sites(sites.df,
                                                          chip_peak_sampleID = chip_peak_sampleID,
                                                          chip_peak_dir = chip_peak_dir)
    saveRDS(sites_chip_labels.df, chip_labels_file)
    cat('Done. Results saved at:', chip_labels_file, '\n')

  }

}



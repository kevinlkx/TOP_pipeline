#! /usr/bin/env Rscript

# Count ChIP-seq read coverage for candidate sites
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
  make_option('--chrom_size', action='store', default=NULL, type='character',
              help='Filename listing the size for each chromosome.'),
  make_option('--chip_bam_dir', action='store', default=NULL, type='character',
              help='Directory for ChIP-seq bam files'),
  make_option('--chip_file_colname', action='store', default='chip_file', type='character',
              help='Column name for ChIP-seq bam files.'),
  make_option('--ref_size', action='store', default=2e7, type='integer',
              help='Reference library size for ChIP-seq.'),
  make_option('--outdir', action='store', default=NULL, type='character',
              help='Output directory.'),
  make_option('--outname', action='store', default=NULL, type='character',
              help='Output filename prefix.'),
  make_option('--bedtools', action='store', default='bedtools', type='character',
              help='Path to bedtools executable.')
)

opt <- parse_args(OptionParser(option_list=option_list))
tf_name                <- opt$tf_name
cell_type              <- opt$cell_type
sites_file             <- opt$sites
metadata_file          <- opt$metadata
chrom_size_file        <- opt$chrom_size
chip_bam_dir           <- opt$chip_bam_dir
chip_file_colname      <- opt$chip_file_colname
ref_size               <- opt$ref_size
outdir                 <- opt$outdir
outname                <- opt$outname
bedtools_path          <- opt$bedtools

# ============ Counting ChIP-seq read coverage for candidate sites =============

chip_counts_file <- file.path(outdir, paste0(outname, '.normalized_chipcounts.rds'))

if( !file.exists(sites_file) || file.size(sites_file) == 0 ){
  cat(paste(sites_file, 'file does not exist or is empty. Probably no candidate sites were found.\n'))
  res <- file.create(chip_counts_file) # create an empty file for snakemake downstream pipeline.
}else{
  sites.df <- as.data.frame(fread(sites_file))

  # Find ChIP-seq bam files of the TF in the cell type of interest using the metadata table
  metadata <- as.data.frame(fread(metadata_file, sep = '\t'))
  chip_files <- metadata[metadata$tf_name == tf_name & metadata$cell_type == cell_type, chip_file_colname]
  # Split the ChIP-seq bam files, which were separated by ';' in the metadata table.
  chip_samples <- unlist(strsplit(chip_files, ';'))

  if ( length(chip_samples) == 0 || is.na(chip_samples) || chip_samples == '' ) {
    cat(paste('No ChIP-seq files for ', tf_name, 'in', cell_type, '\n'))
    res <- file.create(chip_counts_file)
  }else{
    cat('Counting ChIP-seq read coverage for:', tf_name, 'in', cell_type, '\nChIP-seq samples:',chip_samples, '\n')
    chip_bam_files <- file.path(chip_bam_dir, paste0(chip_samples, '.bam'))
    chip_idxstats_files <- file.path(chip_bam_dir, paste0(chip_samples, '.bam.idxstats.txt'))

    sites_chip_counts.df <- count_normalize_chip(sites.df,
                                          chip_bam_files = chip_bam_files,
                                          chip_idxstats_files = chip_idxstats_files,
                                          chrom_size_file = chrom_size_file,
                                          ref_size = ref_size,
                                          bedtools_path = bedtools_path)

    saveRDS(sites_chip_counts.df, chip_counts_file)
    cat('Done. Results saved at:', chip_counts_file, '\n')
  }

}

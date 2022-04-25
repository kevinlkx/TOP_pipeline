#! /usr/bin/env Rscript

# Combine candidate sites, DNase (or ATAC) bins and ChIP data for a TF x cell type combo
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
  make_option("--thresh_pValue", action="store", default='1e-5', type='character',
              help="FIMO p-value threshold [default: %default]"),
  make_option('--acc_bam_dir', action='store', default=NULL, type='character',
              help='Directory of DNase-seq or ATAC-seq bam file.'),
  make_option('--acc_counts_dir', action='store', default=NULL, type='character',
              help='Directory of DNase-seq or ATAC-seq counts.'),
  make_option('--chip_counts_dir', action='store', default=NULL, type='character',
              help='Directory of ChIP counts.'),
  make_option('--chip_peak_dir', action='store', default=NULL, type='character',
              help='Directory of ChIP peaks file.'),
  make_option('--ref_size', action='store', default=1e8, type='integer',
              help='Reference library size for DNase-seq or ATAC-seq.'),
  make_option("--bin", action="store", default='M5', type='character',
              help="MILLIPEDE binning scheme."),
  make_option("--transform", action="store", default='asinh', type='character',
              help="Transformation of binned DNase counts. (asinh or log2)"),
  make_option('--outdir', action='store', default=NULL, type='character',
              help='Output directory.')
)

opt <- parse_args(OptionParser(option_list=option_list))
tf_name                <- opt$tf_name
cell_type              <- opt$cell_type
sites_file             <- opt$sites
metadata_file          <- opt$metadata
thresh_pValue          <- opt$thresh_pValue
acc_bam_dir            <- opt$acc_bam_dir
acc_counts_dir         <- opt$acc_counts_dir
chip_counts_dir        <- opt$chip_counts_dir
chip_peak_dir         <- opt$chip_peak_dir
ref_size               <- opt$ref_size
bin_method             <- opt$bin
transform              <- opt$transform
outdir                 <- opt$outdir


# ================ Combine sites, bins and ChIP data ================

if( !file.exists(sites_file) || file.size(sites_file) == 0 ){
  cat(paste(sites_file, 'file does not exist or is empty.\n'))
} else {

  # Select metadata for the TF x cell type combo
  metadata <- as.data.frame(fread(metadata_file, sep = '\t'))
  selected_metadata <- metadata[which(metadata$tf_name == tf_name & metadata$cell_type == cell_type), ]
  pwm_id <- selected_metadata$pwm_id
  acc_file <- grep("acc_file|dnase_file|atac_file", colnames(selected_metadata), value = T, ignore.case = T)
  if(length(acc_file) > 1){
    stop('More than 1 column of DNase (or ATAC) file names! Check the metadata!')
  }
  acc_samples <- unlist(strsplit(selected_metadata[, acc_file],';'))

  if(!dir.exists(outdir)){dir.create(outdir, showWarnings = F, recursive = T)}

  # Load candidate sites
  cat('Load candidate sites...\n')
  sites.df <- as.data.frame(fread(sites_file))

  # Load ChIP counts data
  chip_counts_file <- file.path(chip_counts_dir, paste0(tf_name, '_', pwm_id, '_', thresh_pValue, '.', cell_type, '.normalized_chipcounts.rds'))
  if (file.exists(chip_counts_file) & file.size(chip_counts_file) > 0){
    cat('Load ChIP counts...\n')
    chip_counts.df <- readRDS(chip_counts_file)
  }else{
    cat(tf_name, 'in',cell_type, 'ChIP counts not available!\n')
    chip_labels.df <- cbind(sites.df, chip = NA)
  }

  # Load ChIP peak labels data
  chip_labels_file <- file.path(chip_counts_dir, paste0(tf_name, '_', pwm_id, '_', thresh_pValue, '.', cell_type, '.chip.peak.labels.rds'))
  if (file.exists(chip_labels_file) & file.size(chip_labels_file) > 0){
    cat('Load ChIP peak labels...\n')
    chip_labels.df <- readRDS(chip_labels_file)
  }else{
    cat(tf_name, 'in',cell_type, 'ChIP peak labels not available!\n')
    chip_labels.df <- cbind(sites.df, chip_label = NA)
  }

  # DNase or ATAC samples
  if ( length(acc_samples) == 1 ) {

    cat('Combining sites, DNase or ATAC bins and ChIP data for', tf_name, 'in', cell_type, '...\n')
    sites_bins_file <- paste0(acc_counts_dir, '/', acc_samples, '/', pwm_id, '_', thresh_pValue, '.', bin_method, '_bins.rds')
    if(!file.exists(sites_bins_file)){stop('Bins data not available!')}
    sites_bins.df <- readRDS(sites_bins_file)
    if(!all.equal(sites_bins.df$name, chip_counts.df$name)){stop('Sites do not match!')}
    if(!all.equal(sites_bins.df$name, chip_labels.df$name)){stop('Sites do not match!')}
    combined_data.df <- data.frame(sites_bins.df, chip = chip_counts.df$chip, chip_label = chip_labels.df$chip_label)
    combined_data_file <- file.path(outdir, paste0(tf_name, '_', pwm_id, '_', thresh_pValue, '.', cell_type, '.', bin_method, '.combined.data.rds'))
    saveRDS(combined_data.df, combined_data_file)

  } else if ( length(acc_samples) > 1 ) {
    # Combine data for each replicate
    cat('Combining sites, DNase or ATAC bins and ChIP data of', tf_name, 'in', cell_type, 'with', length(acc_samples), 'replicates...\n')
    for (i in 1:length(acc_samples)) {
      cat('Combining data for Rep', i, ': ', acc_samples[i], '...\n')
      sites_bins_file <- paste0(acc_counts_dir, '/', acc_samples[1], '/', pwm_id, '_', thresh_pValue, '.', bin_method, '_bins.rds')
      if(!file.exists(sites_bins_file)){stop('Bins data not available!')}
      sites_bins.df <- readRDS(sites_bins_file)
      if(!all.equal(sites_bins.df$name, chip_counts.df$name)){stop('Sites do not match!')}
      if(!all.equal(sites_bins.df$name, chip_labels.df$name)){stop('Sites do not match!')}
      combined_data.df <- data.frame(sites_bins.df, chip = chip_counts.df$chip, chip_label = chip_labels.df$chip_label)
      combined_data_file <- file.path(outdir, paste0(tf_name, '_', pwm_id, '_', thresh_pValue, '.', cell_type, '.', bin_method, '.rep', i, '.combined.data.rds'))
      saveRDS(combined_data.df, combined_data_file)
    }

    # Merge replicates then combine data
    cat('Merge samples:', acc_samples, '...\n')
    acc_counts_files <- paste0(acc_counts_dir, '/', acc_samples, '/', pwm_id, '_', thresh_pValue, '.countmatrix.rds')
    acc_idxstats_files <- paste0(acc_bam_dir, '/', acc_samples, '.bam.idxstats.txt')
    bins.df <- merge_normalize_bin_transform_counts(acc_counts_files,
                                                    acc_idxstats_files,
                                                    ref_size = ref_size,
                                                    bin_method,
                                                    transform)
    sites_bins.df <- data.frame(sites.df, bins.df)
    colnames(sites_bins.df) <- c(colnames(sites.df), paste0('bin', 1:ncol(bins.df)))
    if(!all.equal(sites_bins.df$name, chip_counts.df$name)){stop('Sites do not match!')}
    if(!all.equal(sites_bins.df$name, chip_labels.df$name)){stop('Sites do not match!')}
    combined_data.df <- data.frame(sites_bins.df, chip = chip_counts.df$chip, chip_label = chip_labels.df$chip_label)
    combined_data_file <- file.path(outdir, paste0(tf_name, '_', pwm_id, '_', thresh_pValue, '.', cell_type, '.', bin_method, '.combined.data.rds'))
    saveRDS(combined_data.df, combined_data_file)

  }else{
    cat('No DNase or ATAC samples!\n')
  }

}


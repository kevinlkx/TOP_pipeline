#! /usr/bin/env Rscript

# Predict TF occupancy
library(optparse)
library(TOP)
library(data.table)
# ================ Process the command-line arguments ================
option_list <- list(
  make_option('--metadata', action='store', default=NULL, type='character',
              help='Filename for training metadata table.'),
  make_option('--data_dir', action='store', default=NULL, type='character',
              help='Directory of combined data for each TF x cell type combo.'),
  make_option("--thresh_pValue", action="store", default='1e-5', type='character',
              help="FIMO p-value threshold [default: %default]"),
  make_option("--bin", action="store", default='M5', type='character',
              help="MILLIPEDE binning scheme."),
  make_option('--outdir', action='store', default=NULL, type='character',
              help='Output directory.'),
  make_option('--outname', action='store', default=NULL, type='character',
              help='Output filename prefix.')
)

opt <- parse_args(OptionParser(option_list=option_list))
metadata_file          <- opt$metadata
model_coef_file        <- opt$model
data_dir               <- opt$data_dir
thresh_pValue          <- opt$thresh_pValue
bin_method             <- opt$bin
outdir                 <- opt$outdir

metadata_file          <- '/datacommons/harteminklab/kl124/TOP/data/ENCODE/metadata/processed_JASPARver1/hg38/ATAC_JASPARver1_training_motifs_predict_data_table.tsv'
model_coef_file        <- '/datacommons/harteminklab/kl124/TOP/output/TOP_postfit/hg38/ATAC/JASPARver1/priorVar1/M5/TOP_posterior_mean_coef.rds'
data_dir               <- '/hpc/home/kl124/work/TOP/processed_data_202112/atac_sites_counts/hg38/JASPARver1/'
thresh_pValue          <- '1e-5'
bin_method             <- 'M5'
outdir                 <- '/hpc/home/kl124/work/TOP/TOP_predictions_202112/'


# ================ Load data ================

if( !dir.exists(outdir) ){
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
}

metadata <- read.table(metadata_file, header = T, sep = "\t", stringsAsFactors = FALSE)
tf_list <- c('ATF1', 'NRF1')
celltype_list <- c('K562', 'A549', 'GM12878', 'HepG2')

metadata <- metadata[which( (metadata$tf_name %in% tf_list) & (metadata$cell_type %in% celltype_list) ), ]

TOP_coef <- readRDS(model_coef_file)
summary(TOP_coef)


# ================ Load model ================

for( i in 1:nrow(metadata)) {
  tf_name <- metadata[i,'tf_name']
  pwm_id <- metadata[i,'pwm_id']
  cell_type <- metadata[i,'cell_type']

  cat(tf_name, 'in', cell_type, '...\n')

  acc_file_col <- grep('acc_file|atac_file|dnase_file', colnames(metadata), ignore.case = TRUE, value = TRUE)
  sampleIDs <- unlist(strsplit(metadata[i,acc_file_col], split = ';'))

  if(length(sampleIDs) ==0){
    message(paste('No data for', tf_name, 'in', cell_type, '!'))
  }else if(length(sampleIDs) == 1){
    data_file <- paste0(data_dir, '/',sampleIDs, '/', pwm_id, '_', thresh_pValue, '.', bin_method, '_bins.rds')
    data <- readRDS(data_file)
    res <- predict_TOP(data,
                       tf_name = 'CTCF',
                       cell_type = 'K562',
                       TOP_coef = TOP_coef,
                       level = 'best',
                       logistic_model = FALSE,
                       transform = 'asinh')
    model_level <- res$level
    data$predicted.occupancy <- res$predicted

    fwrite(data, file.path(outdir, paste0(tf_name, '_', pwm_id, '_', thresh_pValue, '.', cell_type, '.', model_level, '.level.predictions.txt.gz')))

  }else{
    for (j in 1:length(sampleIDs)){
      data_file <- paste0(data_dir, '/',sampleIDs[j], '/', pwm_id, '_', thresh_pValue, '.', bin_method, '_bins.rds')
      data <- readRDS(data_file)
      res <- predict_TOP(data,
                         tf_name = 'CTCF',
                         cell_type = 'K562',
                         TOP_coef = TOP_coef,
                         level = 'best',
                         logistic_model = FALSE,
                         transform = 'asinh')
      model_level <- res$level
      data$predicted.occupancy <- res$predicted

      fwrite(data, file.path(outdir, paste0(tf_name, '_', pwm_id, '_', thresh_pValue, '.', cell_type, '.rep', j,'.', model_level, '.level.predictions.txt.gz')))
    }
  }

}


# Choose trained regression coefficients (bottom level)




#! /usr/bin/env Rscript

# Predict binding probability for multiple TF motifs across multiple cell types
library(optparse)
library(TOP)
library(data.table)

# ================ Process the command-line arguments ================
option_list <- list(
  make_option('--metadata', action='store', default=NULL, type='character',
              help='Filename for training metadata table.'),
  make_option('--tf_list', action='store', default=NULL, type='character',
              help='TF name list.'),
  make_option('--celltype_list', action='store', default=NULL, type='character',
              help='Cell type list.'),
  make_option('--model', action='store', default=NULL, type='character',
              help='TOP model file.'),
  make_option('--use_model', action='store', default=NULL, type='character',
              help='Use model: "ATAC", "DukeDNase", "UwDNase"'),
  make_option('--data_dir', action='store', default=NULL, type='character',
              help='Directory of combined data for each TF x cell type combo.'),
  make_option('--sites_dir', action='store', default=NULL, type='character',
              help='Directory for candidate sites.'),
  make_option('--genomecount_dir', action='store', default=NULL, type='character',
              help='Directory for genome counts'),
  make_option('--idxstats_dir', action='store', default=NULL, type='character',
              help='Directory of idxstats files.'),
  make_option("--thresh_pValue", action="store", default='1e-5', type='character',
              help="FIMO p-value threshold [default: %default]"),
  make_option("--bin", action="store", default='M5', type='character',
              help="MILLIPEDE binning scheme."),
  make_option('--outdir', action='store', default=NULL, type='character',
              help='Output directory.')
)

opt <- parse_args(OptionParser(option_list=option_list))
metadata_file          <- opt$metadata
tf_list                <- opt$tf_list
celltype_list          <- opt$celltype_list
model_coef_file        <- opt$model
use_model              <- opt$use_model
data_dir               <- opt$data_dir
sites_dir              <- opt$sites_dir
genomecount_dir        <- opt$genomecount_dir
idxstats_dir           <- opt$idxstats_dir
thresh_pValue          <- opt$thresh_pValue
bin_method             <- opt$bin
outdir                 <- opt$outdir

# ================ Load data ================

metadata <- read.table(metadata_file, header = T, sep = "\t", stringsAsFactors = FALSE)

if( !is.null(tf_list) ) {
  tf_list <- gsub('\\s', '', tf_list)
  tf_list <- unlist(strsplit(tf_list, split = ','))
}else{
  tf_list <- sort(unique(metadata$tf_name))
}

if( !is.null(celltype_list) ) {
  celltype_list <- gsub('\\s', '', celltype_list)
  celltype_list <- unlist(strsplit(celltype_list, split = ','))
}else{
  celltype_list <- sort(unique(metadata$cell_type))
}

cat('TFs:', tf_list, '\n')
cat('Cell types:', celltype_list, '\n')

metadata <- metadata[which( (metadata$tf_name %in% tf_list) & (metadata$cell_type %in% celltype_list) ), ]

# ================ Check and prepare input data ================
acc_file_col <- grep('acc_file|atac_file|dnase_file', colnames(metadata), ignore.case = TRUE, value = TRUE)
sample_ids <- unique(unlist(strsplit(metadata[,acc_file_col], split = ';')))
cat(length(sample_ids), 'samples. \n')
cat(length(unique(metadata$pwm_id)), 'pwm ids. \n')
pwm_id_list <- setdiff(unique(metadata$pwm_id), c('MA1125.1', 'MA1596.1', 'MA1587.1'))
for(pwm_id in pwm_id_list){
  for(sample_id in sample_ids){
    data_file <- paste0(data_dir, '/',sample_id, '/', pwm_id, '_', thresh_pValue, '.', bin_method, '_bins.rds')
    if(!file.exists(data_file)){
      cat('Input data for', pwm_id, 'in', sample_id, 'is not available. \n')
      # if(use_model  == 'ATAC'){
      #   ref_size <- 5e7
      # }else if(use_model %in% c('DukeDNase', 'UwDNase')){
      #   ref_size < 1e8
      # }else{
      #   stop('use_model needs to be "ATAC", "DukeDNase", "UwDNase"!')
      # }
      # cmd <- paste("Rscript get_acc_counts_bins.R",
      #              "--sites", paste0(sites_dir, "/",pwm_id, "_", thresh_pValue, ".candidate_sites.txt.gz"),
      #              "--genomecount_dir", genomecount_dir,
      #              "--genomecount_name", sample_id,
      #              "--idxstats", paste0(idxstats_dir, "/",sample_id, ".bam.idxstats.txt"),
      #              "--bin", bin_method,
      #              "--ref_size", ref_size,
      #              "--transform asinh",
      #              "--outdir", paste0(data_dir, "/", sample_id),
      #              "--outname", paste0(pwm_id, "_", thresh_pValue))
      # system(cmd)
    }
  }
}

cat('Finished checking input data.\n')

# ================ Predict occupancy ================

if( !dir.exists(outdir) ){
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
}

for( i in 1:nrow(metadata)) {

  tf_name <- metadata[i,'tf_name']
  pwm_id <- metadata[i,'pwm_id']
  cell_type <- metadata[i,'cell_type']
  sample_ids <- unique(unlist(strsplit(metadata[i,acc_file_col], split = ';')))

  sites_file <- paste0(sites_dir, "/",pwm_id, "_", thresh_pValue, ".candidate_sites.txt.gz")
  if( !file.exists(sites_file) || file.size(sites_file)==0){
    cat('No sites in', sites_file, '\nSkipped. \n')
    next;
  }

  cat('Predict', tf_name, 'in', cell_type, '...\n')
  cat(length(sample_ids), 'samples. \n')

  if(length(sample_ids) == 1){
    sample_id <- sample_ids
    data_file <- paste0(data_dir, '/',sample_id, '/', pwm_id, '_', thresh_pValue, '.', bin_method, '_bins.rds')

    if(!file.exists(data_file) || file.size(data_file) == 0){
      cat('No data in', data_file, '\nSkipped. \n')
    }else{
      data <- readRDS(data_file)
      top_res <- predict_TOP(data,
                             tf_name = toupper(tf_name),
                             cell_type = cell_type,
                             use_model = use_model,
                             level = 'best',
                             logistic_model = TRUE)

      top_predictions <- top_res$predictions
      top_predictions$predicted[top_predictions$predicted < 0] <- 0

      predicted <- top_predictions[, c('chr', 'start', 'end', 'name', 'pwm.score', 'strand', 'p.value', 'predicted')]
      predicted$predicted <- round(predicted$predicted,3)

      model_level <- top_res$level

      dir.create(paste0(outdir, '/', cell_type), showWarnings = FALSE, recursive = TRUE)
      outfile <- paste0(tf_name, '_', pwm_id, '_', thresh_pValue, '.', cell_type, '.', model_level, '.level.predicted.binding.probability.txt.gz')
      fwrite(predicted, file.path(paste0(outdir, '/', cell_type), outfile), sep = '\t')
    }

  }else if(length(sample_ids) > 1){
    predicted_reps <- data.frame()

    for (j in 1:length(sample_ids)){
      data_file <- paste0(data_dir, '/',sample_ids[j], '/', pwm_id, '_', thresh_pValue, '.', bin_method, '_bins.rds')

      if(!file.exists(data_file) || file.size(data_file) == 0){
        cat('No data in', data_file, '\nSkipped. \n')
      }else{
        data <- readRDS(data_file)
        top_res <- predict_TOP(data,
                               tf_name = toupper(tf_name),
                               cell_type = cell_type,
                               use_model = use_model,
                               level = 'best',
                               logistic_model = TRUE)

        top_predictions <- top_res$predictions
        top_predictions$predicted[top_predictions$predicted < 0] <- 0

        if(j == 1){
          predicted_reps <- top_predictions[, c('chr', 'start', 'end', 'name', 'pwm.score', 'strand', 'p.value', 'predicted')]
        }else{
          predicted_reps <- cbind(predicted_reps, top_predictions$predicted)
        }
      }
    }

    if(nrow(predicted_reps) > 0){
      colnames(predicted_reps) <- c(c('chr', 'start', 'end', 'name', 'pwm.score', 'strand', 'p.value'), paste0('predicted_rep', 1:length(sample_ids)))
      predicted_reps$predicted_mean <- rowMeans(predicted_reps[,paste0('predicted_rep', 1:length(sample_ids))])
      predicted_reps[, grep('predicted', colnames(predicted_reps))] <- round(predicted_reps[, grep('predicted', colnames(predicted_reps))], 3)
      model_level <- top_res$level
      dir.create(paste0(outdir, '/', cell_type), showWarnings = FALSE, recursive = TRUE)
      outfile <- paste0(tf_name, '_', pwm_id, '_', thresh_pValue, '.', cell_type, '.', model_level, '.level.predicted.binding.probability.txt.gz')
      fwrite(predicted_reps, file.path(paste0(outdir, '/', cell_type), outfile), sep = '\t')
    }
  }

}

cat('Finished predicting TF binding probability. \n')




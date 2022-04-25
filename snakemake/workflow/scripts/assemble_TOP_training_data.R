#! /usr/bin/env Rscript

# Assemble training data for selected TF x cell type combos.
library(optparse)
library(doParallel)
library(TOP)

# ================ Process the command-line arguments ================
option_list <- list(
  make_option('--metadata', action='store', default=NULL, type='character',
              help='Filename for training metadata table.'),
  make_option('--tf_list', action='store', default=NULL, type='character',
              help='List of selected training TFs.'),
  make_option('--celltype_list', action='store', default=NULL, type='character',
              help='List of selected training cell types.'),
  make_option('--combined_data_dir', action='store', default=NULL, type='character',
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
tf_list                <- opt$tf_list
celltype_list          <- opt$celltype_list
combined_data_dir      <- opt$combined_data_dir
thresh_pValue          <- opt$thresh_pValue
bin_method             <- opt$bin
outdir                 <- opt$outdir
outname                <- opt$outname

# ================ Assemble TOP training data with ChIP counts ================

metadata <- read.table(metadata_file, header = T, sep = "\t", stringsAsFactors = FALSE)

if( !is.null(tf_list) ) {
  tf_list <- unlist(strsplit(tf_list, split = ','))
}else{
  tf_list <- sort(unique(metadata$tf_name))
}

if( !is.null(celltype_list) ) {
  celltype_list <- unlist(strsplit(celltype_list, split = ','))
}else{
  celltype_list <- sort(unique(metadata$cell_type))
}

metadata <- metadata[which( (metadata$tf_name %in% tf_list) & (metadata$cell_type %in% celltype_list) ), ]
metadata$pwm_name <- paste(metadata$tf_name, metadata$pwm_id, thresh_pValue, sep = '_')
metadata$data_file <- file.path(combined_data_dir,
                                paste0(metadata$pwm_name, '.', metadata$cell_type, '.', bin_method, '.combined.data.rds'))

tf_cell_table <- as.data.frame(metadata[, c('tf_name', 'cell_type', 'data_file')])

if( !dir.exists(outdir) ){
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
}

cat("Assemble TOP training data with ChIP counts ... \n")
# Load data, split into partitions, and combine all tf x cell type combos
assembled_training_data <- assemble_training_data(tf_cell_table,
                                                  logistic_model = FALSE,
                                                  chip_col = 'chip',
                                                  training_chrs = paste0('chr', seq(1,21,2)),
                                                  n_partitions = 10)
# summary(assembled_training_data)
saveRDS(assembled_training_data, file.path(outdir, paste0(outname, '.chip.training.data.all.partitions.rds')))

# Save a table listing all TF and cell type combinations.
tf_cell_combos <- extract_tf_cell_combos(assembled_training_data)
write.table(tf_cell_combos, file.path(outdir, paste0(outname, '.chip.tf_cell_combos.txt')),
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
cat("Assemble training data saved at:", outdir, "\n")

# ================ Assemble TOP training data with ChIP labels ================

cat("Assemble TOP training data with ChIP labels ... \n")
assembled_training_data <- assemble_training_data(tf_cell_table,
                                                  logistic_model = TRUE,
                                                  chip_col = 'chip_label',
                                                  training_chrs = paste0('chr', seq(1,21,2)),
                                                  n_cores = 10,
                                                  n_partitions = 10)
# summary(assembled_training_data)
saveRDS(assembled_training_data, file.path(outdir, paste0(outname, '.chip_label.training.data.all.partitions.rds')))

# Save a table listing all TF and cell type combinations.
tf_cell_combos <- extract_tf_cell_combos(assembled_training_data)
write.table(tf_cell_combos, file.path(outdir, paste0(outname, '.chip_label.tf_cell_combos.txt')),
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

cat("Assemble training data saved at:", outdir, "\n")


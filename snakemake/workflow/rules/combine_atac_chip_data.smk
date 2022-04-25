
# List ATAC result filenames, which will be part of the input of the combine_atac_chip_data rule,
# to check the availability of those files.
def get_atac_files(wildcards):
    metadata_table = pd.read_csv(METADATA_TABLE_FILE, sep='\t')
    metadata_table = metadata_table[(metadata_table["tf_name"] == wildcards.tf_name) &
                                    (metadata_table["cell_type"] == wildcards.cell_type)]
    sample_names = ';'.join(metadata_table["atac_file"].tolist()).split(';')
    atac_countmatrix_files = expand(ATAC_COUNTMATRIX_FILE, atac_sample = sample_names, allow_missing=True)
    atac_bins_files = expand(ATAC_BINS_FILE, atac_sample = sample_names, allow_missing=True)
    atac_files = atac_countmatrix_files + atac_bins_files
    return atac_files

# combine ATAC and ChIP training data for a TF in a cell type
rule combine_atac_chip_data:
    input:
        sites = SITES_FILE,
        metadata = METADATA_TABLE_FILE,
        atac = get_atac_files,
        chip_counts = CHIP_COUNTS_FILE,
        chip_label = CHIP_PEAKLABELS_FILE
    output:
        touch(COMBINED_ATAC_CHIP_DATA_FILE)
    params:
        bin = BIN,
        transform = 'asinh',
        thresh_pValue = THRESH_PVALUE,
        acc_bam_dir = ATAC_BAM_DIR,
        acc_counts_dir = ATAC_SITES_COUNTS_DIR,
        chip_counts_dir = CHIP_COUNTS_DIR,
        chip_peak_dir = CHIP_PEAK_DIR,
        atac_ref_size = ATAC_REF_SIZE,
        outdir = COMBINED_ATAC_CHIP_DATA_DIR
    log:
        f'{LOG_DIR}/combine_atac_chip_data/{VER_GENOME}/{DATA_TYPE}/{MOTIF_SET}/{{tf_name}}_{{pwm_id}}_{THRESH_PVALUE}.{{cell_type}}.log'
    shell:
        '''
        module load R
        Rscript scripts/combine_acc_chip_data.R \
        --tf_name {wildcards.tf_name} \
        --cell_type {wildcards.cell_type} \
        --sites {input.sites} \
        --metadata {input.metadata} \
        --thresh_pValue {params.thresh_pValue} \
        --acc_bam_dir {params.acc_bam_dir} \
        --acc_counts_dir {params.acc_counts_dir} \
        --chip_counts_dir {params.chip_counts_dir} \
        --chip_peak_dir {params.chip_peak_dir} \
        --ref_size {params.atac_ref_size} \
        --bin {params.bin} \
        --transform {params.transform} \
        --outdir {params.outdir} \
        &> {log}
        '''


rule all_combined_atac_chip_data:
    input:
        ALL_COMBINED_ATAC_CHIP_DATA_FILES


rule clean_combined_atac_chip_data:
    shell:
        '''
        rm -rf {COMBINED_ATAC_CHIP_DATA_DIR}
        '''

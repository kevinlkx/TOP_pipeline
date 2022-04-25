
# get normalized_countmatrix.rds filenames, which will be checked for the input of combine_dnase_chip_data rule
def get_dnase_files(wildcards):
    metadata_table = pd.read_csv(METADATA_TABLE_FILE, sep='\t')
    metadata_table = metadata_table[(metadata_table["tf_name"] == wildcards.tf_name) &
                                    (metadata_table["cell_type"] == wildcards.cell_type)]
    sample_names = ';'.join(metadata_table["dnase_file"].tolist()).split(';')
    dnase_countmatrix_files = expand(DNASE_COUNTMATRIX_FILE, dnase_sample = sample_names, allow_missing=True)
    dnase_bins_files = expand(DNASE_BINS_FILE, dnase_sample = sample_names, allow_missing=True)
    dnase_files = dnase_countmatrix_files + dnase_bins_files

# combine DNase and ChIP training data for a TF in a cell type
rule combine_dnase_chip_data:
    input:
        sites = SITES_FILE,
        metadata = METADATA_TABLE_FILE,
        dnase = get_dnase_files,
        chip_counts = CHIP_COUNTS_FILE,
        chip_label = CHIP_PEAKLABELS_FILE
    output:
        touch(COMBINED_DNASE_CHIP_DATA_FILE)
    params:
        bin = BIN,
        transform = 'asinh',
        thresh_pValue = THRESH_PVALUE,
        acc_bam_dir = DNASE_BAM_DIR,
        acc_counts_dir = DNASE_SITES_COUNTS_DIR,
        chip_counts_dir = CHIP_COUNTS_DIR,
        chip_peak_dir = CHIP_PEAK_DIR,
        dnase_ref_size = DNASE_REF_SIZE,
        outdir = COMBINED_DNASE_CHIP_DATA_DIR
    log:
        f'{LOG_DIR}/combine_dnase_chip_data/{VER_GENOME}/{DATA_TYPE}/{MOTIF_SET}/{{tf_name}}_{{pwm_id}}_{THRESH_PVALUE}.{{cell_type}}.log'
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
        --ref_size {params.dnase_ref_size} \
        --bin {params.bin} \
        --transform {params.transform} \
        --outdir {params.outdir} \
        &> {log}
        '''

rule all_combined_dnase_chip_data:
    input:
        ALL_COMBINED_DNASE_CHIP_DATA_FILES


rule clean_combined_dnase_chip_data:
    shell:
        '''
        rm -rf {COMBINED_DNASE_CHIP_DATA_DIR}
        '''

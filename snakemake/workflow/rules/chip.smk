
# List ChIP-seq bam files, which will be part of the input for count_normalize_chip
# to check the availability of those files.
def get_chip_files(wildcards):
    metadata_table = pd.read_csv(METADATA_TABLE_FILE, sep='\t')
    metadata_table = metadata_table[(metadata_table["tf_name"] == wildcards.tf_name) &
                                    (metadata_table["cell_type"] == wildcards.cell_type)]
    sample_names = ';'.join(metadata_table["chip_file"].tolist()).split(';')
    chip_bam_files = expand(CHIP_BAM_FILE, chip_sample = sample_names, allow_missing=True)
    chip_idxstats_files = expand(CHIP_IDXSTATS_FILE, chip_sample = sample_names, allow_missing=True)
    return chip_idxstats_files

# Count ChIP-seq read coverage for candidate sites
# this requires getting the idxstats for (all) the ChIP-seq bam files first,
# as snakemake won't not know beforehand which ChIP-seq bam files belong to the TF x cell type combo.
rule get_sites_chip_counts:
    input:
        sites = SITES_FILE,
        chrom_size = CHROM_SIZE_FILE,
        metadata = METADATA_TABLE_FILE,
        chip_files = get_chip_files
    output:
        touch(CHIP_COUNTS_FILE)
    params:
        ref_size = CHIP_REF_SIZE,
        thresh_pValue = THRESH_PVALUE,
        chip_bam_dir = CHIP_BAM_DIR,
        chip_counts_dir = CHIP_COUNTS_DIR,
        chip_file_colname = 'chip_file'
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/get_sites_chip_counts/{VER_GENOME}/{{tf_name}}_{{pwm_id}}_{THRESH_PVALUE}.{{cell_type}}.log'
    shell:
        '''
        module load R
        Rscript scripts/get_sites_chip_counts.R \
        --tf_name {wildcards.tf_name} \
        --cell_type {wildcards.cell_type} \
        --sites {input.sites} \
        --metadata {input.metadata} \
        --chrom_size {input.chrom_size} \
        --chip_bam_dir {params.chip_bam_dir} \
        --chip_file_colname {params.chip_file_colname} \
        --ref_size {params.ref_size} \
        --outdir {params.chip_counts_dir} \
        --outname {wildcards.tf_name}_{wildcards.pwm_id}_{params.thresh_pValue}.{wildcards.cell_type} \
        &> {log}
        '''


# get ChIP peaks labels for candidate sites
rule get_sites_chip_peak_labels:
    input:
        sites = SITES_FILE,
        metadata = METADATA_TABLE_FILE
    output:
        touch(CHIP_PEAKLABELS_FILE)
    params:
        thresh_pValue = THRESH_PVALUE,
        chip_peakfile_colname = 'chip_peak_file',
        chip_peak_dir = CHIP_PEAK_DIR,
        chip_counts_dir = CHIP_COUNTS_DIR
    log:
        f'{LOG_DIR}/get_sites_chip_peaks_labels/{VER_GENOME}/{{tf_name}}_{{pwm_id}}_{THRESH_PVALUE}.{{cell_type}}.log'
    shell:
        '''
        module load R
        Rscript scripts/get_sites_chip_peak_labels.R \
        --tf_name {wildcards.tf_name} \
        --cell_type {wildcards.cell_type} \
        --sites {input.sites} \
        --metadata {input.metadata} \
        --chip_peak_dir {params.chip_peak_dir} \
        --chip_peakfile_colname {params.chip_peakfile_colname} \
        --outdir {params.chip_counts_dir} \
        --outname {wildcards.tf_name}_{wildcards.pwm_id}_{params.thresh_pValue}.{wildcards.cell_type} \
        &> {log}
        '''


rule list_chip_samples:
    input: ALL_CHIP_BAM_FILES
    shell:
        ' echo "ChIP samples: {CHIP_SET}" '


rule all_chip_files:
    input:
        ALL_CHIP_BAM_FILES,
        ALL_CHIP_BAI_FILES,
        ALL_CHIP_IDXSTATS_FILES


rule all_chip_counts:
    input:
        ALL_CHIP_COUNTS_FILES


rule all_chip_peak_labels:
    input:
        ALL_CHIP_PEAKLABELS_FILES

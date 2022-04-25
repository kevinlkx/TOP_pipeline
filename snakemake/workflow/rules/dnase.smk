
# Count DNase cuts along the genome
rule count_dnase_genome_cuts:
    input:
        bam = DNASE_BAM_FILE,
        chrom_size = CHROM_SIZE_FILE,
    output:
        fwd = DNASE_GENOMECOUNTS_FWD_FILE,
        rev = DNASE_GENOMECOUNTS_REV_FILE
    params:
        data_type = 'DNase'
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/count_dnase_genome_cuts/{VER_GENOME}/{{dnase_sample}}.log'
    shell:
        '''
        module load R
        Rscript scripts/count_acc_genome_cuts.R \
        --bam {input.bam} \
        --chrom_size {input.chrom_size} \
        --data_type {params.data_type} \
        --outdir {DNASE_GENOMECOUNTS_DIR} \
        --outname {wildcards.dnase_sample} \
        &> {log}
        '''


# Get DNase-seq count matrices for the candidate sites,
# then normalize, bin and transform the counts
rule get_dnase_counts_bins:
    input:
        sites = SITES_FILE,
        fwd_genomecount = DNASE_GENOMECOUNTS_FWD_FILE,
        rev_genomecount = DNASE_GENOMECOUNTS_REV_FILE,
        idxstats = DNASE_IDXSTATS_FILE,
    output:
        touch(DNASE_COUNTMATRIX_FILE),
        touch(DNASE_BINS_FILE)
    params:
        transform = 'asinh',
        ref_size = DNASE_REF_SIZE,
        bin = BIN,
        genomecount_dir = DNASE_GENOMECOUNTS_DIR
    log:
        f'{LOG_DIR}/get_sites_dnase_counts_bins/{VER_GENOME}/{{dnase_sample}}/{{pwm_id}}_{THRESH_PVALUE}.log'
    shell:
        '''
        module load R
        Rscript scripts/get_acc_counts_bins.R \
        --sites {input.sites} \
        --genomecount_dir {params.genomecount_dir} \
        --genomecount_name {wildcards.dnase_sample} \
        --idxstats {input.idxstats} \
        --bin {params.bin} \
        --ref_size {params.ref_size} \
        --transform {params.transform} \
        --outdir {DNASE_SITES_COUNTS_DIR}/{wildcards.dnase_sample} \
        --outname {wildcards.pwm_id}_{THRESH_PVALUE} \
        &> {log}
        '''


rule list_dnase_samples:
    input: ALL_DNASE_BAM_FILES
    shell:
        ' echo "DNase samples: {DNASE_SET}" '


rule all_dnase_files:
    input:
        ALL_DNASE_BAM_FILES,
        ALL_DNASE_BAI_FILES,
        ALL_DNASE_IDXSTATS_FILES


rule all_dnase_genomecounts:
    input:
        ALL_DNASE_BAM_FILES,
        ALL_DNASE_IDXSTATS_FILES,
        ALL_DNASE_GENOMECOUNTS_FWD_FILES,
        ALL_DNASE_GENOMECOUNTS_REV_FILES


rule all_dnase_countmatrices:
    input:
        ALL_DNASE_COUNTMATRIX_FILES


rule all_dnase_bins:
    input:
        ALL_DNASE_IDXSTATS_FILES,
        ALL_DNASE_BINS_FILES

rule clean_dnase:
    shell:
        '''
        rm -rf {OUT_DIR}/atac_genome_counts/{VER_GENOME}
        rm -rf {OUT_DIR}/atac_sites_counts/{VER_GENOME}/{MOTIF_SET}
        '''

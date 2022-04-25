
# Count ATAC cuts along the genome
rule count_atac_genome_cuts:
    input:
        bam = ATAC_BAM_FILE,
        chrom_size = CHROM_SIZE_FILE
    output:
        fwd = ATAC_GENOMECOUNTS_FWD_FILE,
        rev = ATAC_GENOMECOUNTS_REV_FILE
    params:
        shift_ATAC = '--shift_ATAC',
        shift_ATAC_bases = '4,-4'
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/count_atac_genome_cuts/{VER_GENOME}/{{atac_sample}}.log'
    shell:
        '''
        module load R
        Rscript scripts/count_acc_genome_cuts.R \
        --bam {input.bam} \
        --chrom_size {input.chrom_size} \
        {params.shift_ATAC} \
        --shift_ATAC_bases {params.shift_ATAC_bases} \
        --outdir {ATAC_GENOMECOUNTS_DIR} \
        --outname {wildcards.atac_sample} \
        &> {log}
        '''


# Get ATAC-seq count matrices for the candidate sites,
# then normalize, bin and transform the counts
rule get_atac_counts_bins:
    input:
        sites = SITES_FILE,
        fwd_genomecount = ATAC_GENOMECOUNTS_FWD_FILE,
        rev_genomecount = ATAC_GENOMECOUNTS_REV_FILE,
        idxstats = ATAC_IDXSTATS_FILE
    output:
        touch(ATAC_COUNTMATRIX_FILE),
        touch(ATAC_BINS_FILE)
    params:
        transform = 'asinh',
        ref_size = ATAC_REF_SIZE,
        bin = BIN,
        genomecount_dir = ATAC_GENOMECOUNTS_DIR,
    log:
        f'{LOG_DIR}/get_sites_atac_counts_bins/{VER_GENOME}/{{atac_sample}}/{{pwm_id}}_{THRESH_PVALUE}.log'
    shell:
        '''
        module load R
        Rscript scripts/get_acc_counts_bins.R \
        --sites {input.sites} \
        --genomecount_dir {params.genomecount_dir} \
        --genomecount_name {wildcards.atac_sample} \
        --idxstats {input.idxstats} \
        --bin {params.bin} \
        --ref_size {params.ref_size} \
        --transform {params.transform} \
        --outdir {ATAC_SITES_COUNTS_DIR}/{wildcards.atac_sample} \
        --outname {wildcards.pwm_id}_{THRESH_PVALUE} \
        &> {log}
        '''


rule list_atac_samples:
    input: ALL_ATAC_BAM_FILES
    shell:
        ' echo "ATAC samples: {ATAC_SET}" '


rule all_atac_files:
    input:
        ALL_ATAC_BAM_FILES,
        ALL_ATAC_BAI_FILES,
        ALL_ATAC_IDXSTATS_FILES


rule all_atac_genomecounts:
    input:
        ALL_ATAC_BAM_FILES,
        ALL_ATAC_IDXSTATS_FILES,
        ALL_ATAC_GENOMECOUNTS_FWD_FILES,
        ALL_ATAC_GENOMECOUNTS_REV_FILES


rule all_atac_countmatrices:
    input:
        ALL_ATAC_COUNTMATRIX_FILES


rule all_atac_bins:
    input:
        ALL_ATAC_IDXSTATS_FILES,
        ALL_ATAC_BINS_FILES

rule clean_atac:
    shell:
        '''
        rm -rf {ATAC_GENOMECOUNTS_DIR}
        rm -rf {ATAC_SITES_COUNTS_DIR}
        '''

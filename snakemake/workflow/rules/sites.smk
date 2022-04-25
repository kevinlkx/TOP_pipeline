# Scan for motif matches using FIMO
rule fimo_motif_matches:
    input:
        pwm = PWM_FILE,
        ref_genome = REF_GENOME_FILE
    output:
        touch(FIMO_FILE)
    params:
        thresh_pValue = THRESH_PVALUE,
        background = 'default',
        outdir = FIMO_DIR
    log:
        f'{LOG_DIR}/fimo/{VER_GENOME}/{{pwm_id}}_{THRESH_PVALUE}.log'
    shell:
        '''
        module load R
        Rscript scripts/fimo_motif_matches.R \
        --motif_file {input.pwm} \
        --sequence_file {input.ref_genome} \
        --outname {output} \
        --outdir {params.outdir} \
        --thresh_pValue {params.thresh_pValue} \
        --background {params.background} \
        --max_strand \
        --skip_matched_sequence \
        &> {log}
        '''

# Get candidate sites using FIMO motif matches
rule candidate_sites:
    input:
        fimo = FIMO_FILE,
        blacklist = BLACKLIST_FILE
    output:
        touch(SITES_FILE)
    params:
        flank = 100,
        thresh_pValue = THRESH_PVALUE,
        thresh_pwmscore = THRESH_PWMSCORE
    conda:
        '../envs/top.yaml'
    log:
        f'{LOG_DIR}/candidate_sites/{VER_GENOME}/{{pwm_id}}_{THRESH_PVALUE}.log'
    shell:
        '''
        module load R
        Rscript scripts/get_candidate_sites.R \
        --fimo {input.fimo} \
        --flank {params.flank} \
        --thresh_pValue {params.thresh_pValue} \
        --thresh_pwmscore {params.thresh_pwmscore} \
        --blacklist {input.blacklist} \
        --out {output} \
        &> {log}
        '''


# list all PWM IDs
rule list_pwms:
    input: ALL_PWM_FILES
    shell:
        ' echo "PWMs: {PWMID_SET}" '

# get all fimo files
rule all_fimo:
    input: ALL_FIMO_FILES

# get all candidate sites
rule all_sites:
    input: ALL_SITES_FILES

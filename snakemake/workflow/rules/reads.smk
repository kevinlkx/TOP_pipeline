
# # Generate index and idxstats for sorted reads
# rule samtools_index_idxstats:
#     input:
#         '{bam_dir}/{sample}.bam',
#     output:
#         bai = '{bam_dir}/{sample}.bam.bai',
#         idxstats = '{bam_dir}/{sample}.bam.idxstats.txt'
#     conda:
#         '../envs/top.yaml'
#     shell:
#         '''
#         samtools index {input}
#         samtools idxstats {input} > {output.idxstats}
#         '''

# Sort, index and idxstats for BAM alignments
rule index_idxstats_bam:
    input:
        '{bam_dir}/{sample}.bam',
    output:
        bai = '{bam_dir}/{sample}.bam.bai',
        idxstats = '{bam_dir}/{sample}.bam.idxstats.txt'
    conda:
        '../envs/top.yaml'
    shell:
        '''
        module load R
        Rscript scripts/sort_index_idxstats_bam.R \
        --bam_file {input} \
        --index --idxstats
        '''

# Download ENCODE bam file
rule download_encode_bam_file:
    output:
        '{bam_dir}/{sample}.bam'
    shell:
        '''
        cd {wildcards.bam_dir}
        curl -O -L -sS https://www.encodeproject.org/files/{wildcards.sample}/@@download/{wildcards.sample}.bam
        '''

# Download ENCODE bed file
rule download_encode_bed_file:
    output:
        '{bed_dir}/{sample}.bed.gz'
    shell:
        '''
        cd {wildcards.bed_dir}
        curl -O -L -sS https://www.encodeproject.org/files/{wildcards.sample}/@@download/{wildcards.sample}.bed.gz
        '''

# Download ENCODE bigwig file
rule download_encode_bigwig_file:
    output:
        '{bigwig_dir}/{sample}.bigwig.gz'
    shell:
        '''
        cd {wildcards.bigwig_dir}
        curl -O -L -sS https://www.encodeproject.org/files/{wildcards.sample}/@@download/{wildcards.sample}.bigwig.gz
        '''

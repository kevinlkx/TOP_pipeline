#!/bin/bash

#SBATCH --job-name=predict_UWdnase
#SBATCH --output=predict_UWdnase_%A_%a.out
#SBATCH --mem=15G

## Parameters
module load R/4.1.1-rhel8

METADATA='/datacommons/harteminklab/kl124/TOP/data/ENCODE/metadata/predictions/ENCODE2012_hg19/UWDNase_JASPAR2022NR_all_motifs_predict_data_table.tsv'
MODEL='/datacommons/harteminklab/kl124/TOP/output/TOP_postfit/hg19/UwDnase/JASPARver1/priorVar1/M5/TOP_posterior_mean_coef.rds'
DATADIR='/hpc/home/kl124/work/TOP/processed_ENCODE2012_data/dnase_sites_counts/hg19/JASPAR2022'
SITESDIR='/hpc/home/kl124/work/TOP/processed_ENCODE2012_data/candidate_sites/hg19/JASPAR2022'
GENOMECOUNTSDIR='/hpc/home/kl124/work/TOP/processed_ENCODE2012_data/dnase_genome_counts/hg19'
BAMDIR='/hpc/home/kl124/work/data/ENCODE/DNase-seq/hg19'
OUTDIR='/hpc/home/kl124/datacommons/TOP/TOP_predictions_202112/hg19/UWDNase/JASPAR2022/predicted_occupancy'

echo metadata=${METADATA}
echo model=${MODEL}
echo datadir=${DATADIR}
echo sites_dir=${SITESDIR}
echo genomecount_dir=${GENOMECOUNTSDIR}
echo bamdir=${BAMDIR}
echo outdir=${OUTDIR}

tf_list_file='/datacommons/harteminklab/kl124/TOP/data/ENCODE/metadata/predictions/hg38/JASPAR2022_CORE_vertebrates_nr_filtered_TF_motif_table.tsv'

echo "File listing TFs and cell types: ${tf_list_file}"
echo "Run #${SLURM_ARRAY_TASK_ID} on the tf list."
idx_line=${SLURM_ARRAY_TASK_ID}
tf_name=$(awk -v i=${idx_line} -v j=1 'FNR == i {print $j}' ${tf_list_file})
pwm_id=$(awk -v i=${idx_line} -v j=2 'FNR == i {print $j}' ${tf_list_file})
echo "${tf_name}, ${pwm_id}"

# if [[ "${pwm_id}" =~ ^('MA1125.1'|'MA1596.1'|'MA1587.1'|'MA1723.1'|'MA1107.2')$ ]]; then
#    echo "Try large memory for motif ${pwm_id}."
#    exit 1
# fi

cd /hpc/group/harteminklab/kl124/projects/TOP/TOP_workflow/snakemake/workflow/scripts/

start=`date +%s`

Rscript predict_occupancy_motifs_cells.R \
  --metadata ${METADATA} \
  --tf_list ${tf_name} \
  --model ${MODEL} \
  --data_dir ${DATADIR} \
  --sites_dir ${SITESDIR} \
  --genomecount_dir ${GENOMECOUNTSDIR} \
  --idxstats_dir ${BAMDIR} \
  --ref_size '1e8' \
  --thresh_pValue '1e-5' \
  --bin 'M5' \
  --outdir ${OUTDIR}

end=`date +%s`

runtime=$((end-start))
echo "Run time: ${runtime} seconds"
echo "Results saved at: ${OUTDIR}"


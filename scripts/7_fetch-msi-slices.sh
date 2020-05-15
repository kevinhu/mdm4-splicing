#!/bin/sh

# qsub -t 1-329 fetch-msi-slices.sh -l h_rt=8:00:00

source /broad/software/scripts/useuse
use Anaconda3

source activate /home/unix/khu/.conda/envs/telomerehunter_env

export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

BAM_FILES=/broad/hptmp/khu/wgs_paths.txt
ACH_IDS=/broad/hptmp/khu/wgs_ids.txt
BAM_PATH=$(awk "NR==$SGE_TASK_ID" $BAM_FILES)
ACH_ID=$(awk "NR==$SGE_TASK_ID" $ACH_IDS)

echo $BAM_PATH

INTERVALS=/broad/hptmp/khu/MSI_exon_bounds.bed

cd /broad/hptmp/khu/ccle_wgs_bams

samtools view -b -h $BAM_PATH -L $INTERVALS -M > /broad/hptmp/khu/ccle_wgs_bams/$ACH_ID.bam
samtools index /broad/hptmp/khu/ccle_wgs_bams/$ACH_ID.bam

echo "COMPLETED"

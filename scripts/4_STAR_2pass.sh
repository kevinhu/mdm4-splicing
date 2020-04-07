#!/bin/sh

#for i in {1..30} ; do qsub -l h_vmem=48G -l h_rt=24:00:00 ./4_STAR_2pass.sh $i ; done

source /broad/software/scripts/useuse

use STAR

FASTQ_FILES=/broad/hptmp/khu/fastq_samples.txt
sample_id=$(awk "NR==$1" $FASTQ_FILES)

echo $sample_id

fastqDir=/broad/hptmp/khu/zipped_fastqs
outputDir=/broad/hptmp/khu/STAR_2pass
genomeDir=/broad/hptmp/khu/STAR_index

mkdir ${outputDir}/${sample_id}/

echo ${fastqDir}/${sample_id}_R1.fq.gz
echo ${fastqDir}/${sample_id}_R2.fq.gz

STAR --runMode alignReads \
    --genomeDir $genomeDir \
    --readFilesCommand zcat \
    --readFilesIn ${fastqDir}/${sample_id}_R1.fq.gz ${fastqDir}/${sample_id}_R2.fq.gz \
    --outFileNamePrefix ${outputDir}/${sample_id}/ \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic \
    --twopass1readsN 1000000000000 \
    --sjdbOverhang 100 \
    --runThreadN 4

echo "COMPLETED"

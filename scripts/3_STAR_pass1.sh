#!/bin/sh

#for i in {1..21} ; do qsub -l h_vmem=32G -l h_rt=24:00:00 ./3_STAR_pass1.sh $i ; done

source /broad/software/scripts/useuse

use STAR

FASTQ_FILES=/broad/hptmp/khu/fastq_samples.txt
sample_id=$(awk "NR==$1" $FASTQ_FILES)

echo $sample_id

fastqDir=/broad/hptmp/khu/zipped_fastqs
outputDir=/broad/hptmp/khu/STAR_pass1
genomeDir=/broad/hptmp/khu/STAR_b37

mkdir ${outputDir}/${sample_id}/

STAR --runMode alignReads \
    --genomeDir $genomeDir \
    --readFilesIn ${fastqDir}/${sample_id}_R1.fq.gz ${fastqDir}/${sample_id}_R2.fq.gz \
    --outFileNamePrefix ${outputDir}/${sample_id}/ \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN 4 \

echo "COMPLETED"


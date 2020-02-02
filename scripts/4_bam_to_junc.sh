#!/bin/sh

#for i in {1..21} ; do ./4_bam_to_junc.sh $i ; done

FASTQ_FILES=/Users/khu/Desktop/github/MDM4-splicing/data/raw/fastq_samples.txt
sample_id=$(awk "NR==$1" $FASTQ_FILES)

echo $sample_id

/Users/khu/Desktop/github/leafcutter/scripts/bam2junc.sh \
    /Volumes/aulyd_xk/STAR_pass1/${sample_id}/Aligned.sortedByCoord.out.bam \
    /Volumes/aulyd_xk/leafcutter_juncs/${sample_id}.junc

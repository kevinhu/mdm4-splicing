#!/bin/sh

cd /Users/khu/Desktop/github/MDM4-splicing/data/raw/leafcutter_outputs/

/Users/khu/Desktop/github/leafcutter/scripts/leafcutter_ds.R \
    --num_threads 4 \
    --exon_file=/Users/khu/Desktop/github/leafcutter/leafcutter/data/gencode19_exons.txt.gz \
    --min_samples_per_intron 3 \
    --output_prefix RPL22_oe \
    --min_coverage 16 \
    --min_samples_per_intron 3 \
    /Volumes/aulyd_xk/leafcutter_juncs/leafcutter_clusters_perind_numers.counts.gz \
    /Users/khu/Desktop/github/MDM4-splicing/data/raw/sample_groups/RPL22_oe.txt

/Users/khu/Desktop/github/leafcutter/scripts/leafcutter_ds.R \
    --num_threads 4 \
    --exon_file=/Users/khu/Desktop/github/leafcutter/leafcutter/data/gencode19_exons.txt.gz \
    --min_samples_per_intron 3 \
    --output_prefix RPL22L1_oe \
    --min_coverage 16 \
    --min_samples_per_intron 3 \
    /Volumes/aulyd_xk/leafcutter_juncs/leafcutter_clusters_perind_numers.counts.gz \
    /Users/khu/Desktop/github/MDM4-splicing/data/raw/sample_groups/RPL22L1_oe.txt

/Users/khu/Desktop/github/leafcutter/scripts/leafcutter_ds.R \
    --num_threads 4 \
    --exon_file=/Users/khu/Desktop/github/leafcutter/leafcutter/data/gencode19_exons.txt.gz \
    --min_samples_per_intron 3 \
    --output_prefix sh704 \
    --min_coverage 16 \
    --min_samples_per_intron 3 \
    /Volumes/aulyd_xk/leafcutter_juncs/leafcutter_clusters_perind_numers.counts.gz \
    /Users/khu/Desktop/github/MDM4-splicing/data/raw/sample_groups/sh704.txt

/Users/khu/Desktop/github/leafcutter/scripts/leafcutter_ds.R \
    --num_threads 4 \
    --exon_file=/Users/khu/Desktop/github/leafcutter/leafcutter/data/gencode19_exons.txt.gz \
    --min_samples_per_intron 3 \
    --output_prefix sh705 \
    --min_coverage 16 \
    --min_samples_per_intron 3 \
    /Volumes/aulyd_xk/leafcutter_juncs/leafcutter_clusters_perind_numers.counts.gz \
    /Users/khu/Desktop/github/MDM4-splicing/data/raw/sample_groups/sh705.txt

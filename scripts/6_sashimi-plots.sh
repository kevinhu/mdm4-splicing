#!/bin/sh

lncap_gfp="../data/raw/STAR_BAMs/LNCaP_GFP_1.bam,../data/raw/STAR_BAMs/LNCaP_GFP_2.bam,../data/raw/STAR_BAMs/LNCaP_GFP_3.bam"
lncap_rpl22="../data/raw/STAR_BAMs/LNCaP_RPL22_1.bam,../data/raw/STAR_BAMs/LNCaP_RPL22_2.bam,../data/raw/STAR_BAMs/LNCaP_RPL22_3.bam"

cal851_gfp="../data/raw/STAR_BAMs/CAL851_GFP_1.bam,../data/raw/STAR_BAMs/CAL851_GFP_2.bam,../data/raw/STAR_BAMs/CAL851_GFP_3.bam"
cal851_rpl22l1="../data/raw/STAR_BAMs/CAL851_RPL22L1_1.bam,../data/raw/STAR_BAMs/CAL851_RPL22L1_2.bam,../data/raw/STAR_BAMs/CAL851_RPL22L1_3.bam"

lncap_luc="../data/raw/STAR_BAMs/LNCaP_shLuc_1.bam,../data/raw/STAR_BAMs/LNCaP_shLuc_2.bam,../data/raw/STAR_BAMs/LNCaP_shLuc_3.bam"
lncap_rpl22l1_1="../data/raw/STAR_BAMs/LNCaP_sh704_1.bam,../data/raw/STAR_BAMs/LNCaP_sh704_2.bam,../data/raw/STAR_BAMs/LNCaP_sh704_3.bam"
lncap_rpl22l1_2="../data/raw/STAR_BAMs/LNCaP_sh705_1.bam,../data/raw/STAR_BAMs/LNCaP_sh705_2.bam,../data/raw/STAR_BAMs/LNCaP_sh705_3.bam"

ncih2110_gfp="../data/raw/STAR_BAMs/NCIH2110_GFP_1.bam,../data/raw/STAR_BAMs/NCIH2110_GFP_2.bam,../data/raw/STAR_BAMs/NCIH2110_GFP_3.bam"
ncih2110_rpl22_1="../data/raw/STAR_BAMs/NCIH2110_RPL22-1A1_1.bam,../data/raw/STAR_BAMs/NCIH2110_RPL22-1A1_2.bam,../data/raw/STAR_BAMs/NCIH2110_RPL22-1A1_3.bam"
ncih2110_rpl22_2="../data/raw/STAR_BAMs/NCIH2110_RPL22-4A1_1.bam,../data/raw/STAR_BAMs/NCIH2110_RPL22-4A1_2.bam,../data/raw/STAR_BAMs/NCIH2110_RPL22-4A1_3.bam"

zr751_gfp="../data/raw/STAR_BAMs/ZR751_GFP_1.bam,../data/raw/STAR_BAMs/ZR751_GFP_2.bam,../data/raw/STAR_BAMs/ZR751_GFP_3.bam"
zr751_rpl22_1="../data/raw/STAR_BAMs/ZR751_RPL22-1A1_1.bam,../data/raw/STAR_BAMs/ZR751_RPL22-1A1_2.bam,../data/raw/STAR_BAMs/ZR751_RPL22-1A1_3.bam"
zr751_rpl22_2="../data/raw/STAR_BAMs/ZR751_RPL22-4A1_1.bam,../data/raw/STAR_BAMs/ZR751_RPL22-4A1_2.bam,../data/raw/STAR_BAMs/ZR751_RPL22-4A1_3.bam"

sashimi_dir="../plots/sashimi_plots/"

colors="#393e46,#e84545"
exon_scaledown=1
intron_scaledown=10

rm -rf $sashimi_dir*

# A3SS events

rmats2sashimiplot \
	--b1 $lncap_gfp \
	--b2 $lncap_rpl22 \
	-t A3SS \
	-e ../data/raw/sashimi_events/rpl22_oe_a3ss.txt \
	--l1 control --l2 treatment \
	--exon_s $exon_scaledown \
	--intron_s $intron_scaledown \
	-o $sashimi_dir/rpl22_oe_a3ss \
	--group-info ../data/raw/sashimi_groups/rpl22_oe.gf \
	--color $colors

rmats2sashimiplot \
	--b1 $ncih2110_gfp \
	--b2 $ncih2110_rpl22_1 \
	-t A3SS \
	-e ../data/raw/sashimi_events/rpl22_a_ko1_a3ss.txt \
	--l1 control --l2 treatment \
	--exon_s $exon_scaledown \
	--intron_s $intron_scaledown \
	-o $sashimi_dir/rpl22_a_ko1_a3ss \
	--group-info ../data/raw/sashimi_groups/rpl22_a_ko1.gf \
	--color $colors

rmats2sashimiplot \
	--b1 $ncih2110_gfp \
	--b2 $ncih2110_rpl22_2 \
	-t A3SS \
	-e ../data/raw/sashimi_events/rpl22_a_ko2_a3ss.txt \
	--l1 control --l2 treatment \
	--exon_s $exon_scaledown \
	--intron_s $intron_scaledown \
	-o $sashimi_dir/rpl22_a_ko2_a3ss \
	--group-info ../data/raw/sashimi_groups/rpl22_a_ko2.gf \
	--color $colors

rmats2sashimiplot \
	--b1 $zr751_gfp \
	--b2 $zr751_rpl22_1 \
	-t A3SS \
	-e ../data/raw/sashimi_events/rpl22_b_ko1_a3ss.txt \
	--l1 control --l2 treatment \
	--exon_s $exon_scaledown \
	--intron_s $intron_scaledown \
	-o $sashimi_dir/rpl22_b_ko1_a3ss \
	--group-info ../data/raw/sashimi_groups/rpl22_b_ko1.gf \
	--color $colors

rmats2sashimiplot \
	--b1 $zr751_gfp \
	--b2 $zr751_rpl22_2 \
	-t A3SS \
	-e ../data/raw/sashimi_events/rpl22_b_ko2_a3ss.txt \
	--l1 control --l2 treatment \
	--exon_s $exon_scaledown \
	--intron_s $intron_scaledown \
	-o $sashimi_dir/rpl22_b_ko2_a3ss \
	--group-info ../data/raw/sashimi_groups/rpl22_b_ko2.gf \
	--color $colors

# SE events

rmats2sashimiplot \
	--b1 $lncap_gfp \
	--b2 $lncap_rpl22 \
	-t SE \
	-e ../data/raw/sashimi_events/rpl22_oe_se.txt \
	--l1 control --l2 treatment \
	--exon_s $exon_scaledown \
	--intron_s $intron_scaledown \
	-o $sashimi_dir/rpl22_oe_se \
	--group-info ../data/raw/sashimi_groups/rpl22_oe.gf \
	--color $colors

rmats2sashimiplot \
	--b1 $ncih2110_gfp \
	--b2 $ncih2110_rpl22_1 \
	-t SE \
	-e ../data/raw/sashimi_events/rpl22_a_ko1_se.txt \
	--l1 control --l2 treatment \
	--exon_s $exon_scaledown \
	--intron_s $intron_scaledown \
	-o $sashimi_dir/rpl22_a_ko1_se \
	--group-info ../data/raw/sashimi_groups/rpl22_a_ko1.gf \
	--color $colors

rmats2sashimiplot \
	--b1 $ncih2110_gfp \
	--b2 $ncih2110_rpl22_2 \
	-t SE \
	-e ../data/raw/sashimi_events/rpl22_a_ko2_se.txt \
	--l1 control --l2 treatment \
	--exon_s $exon_scaledown \
	--intron_s $intron_scaledown \
	-o $sashimi_dir/rpl22_a_ko2_se \
	--group-info ../data/raw/sashimi_groups/rpl22_a_ko2.gf \
	--color $colors

rmats2sashimiplot \
	--b1 $zr751_gfp \
	--b2 $zr751_rpl22_1 \
	-t SE \
	-e ../data/raw/sashimi_events/rpl22_b_ko1_se.txt \
	--l1 control --l2 treatment \
	--exon_s $exon_scaledown \
	--intron_s $intron_scaledown \
	-o $sashimi_dir/rpl22_b_ko1_se \
	--group-info ../data/raw/sashimi_groups/rpl22_b_ko1.gf \
	--color $colors

rmats2sashimiplot \
	--b1 $zr751_gfp \
	--b2 $zr751_rpl22_2 \
	-t SE \
	-e ../data/raw/sashimi_events/rpl22_b_ko2_se.txt \
	--l1 control --l2 treatment \
	--exon_s $exon_scaledown \
	--intron_s $intron_scaledown \
	-o $sashimi_dir/rpl22_b_ko2_se \
	--group-info ../data/raw/sashimi_groups/rpl22_b_ko2.gf \
	--color $colors

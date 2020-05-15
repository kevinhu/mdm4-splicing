BIN_VERSION="0.10.0"

bcftools mpileup -Ob \
	--bam-list wgs_bams_1.txt \
	-f ../data/raw/reference/hg19.fa \
	-R MSI_exon_bounds.bed > ../data/intermediate/msi_exon_pileup_1.bcf.gz

bcftools mpileup -Ob \
	--bam-list wgs_bams_2.txt \
	-f ../data/raw/reference/hg19.fa \
	-R MSI_exon_bounds.bed > ../data/intermediate/msi_exon_pileup_2.bcf.gz

bcftools sort -Ob ../data/intermediate/msi_exon_pileup_1.bcf.gz -o ../data/intermediate/msi_exon_pileup_1_sorted.bcf.gz
bcftools sort -Ob ../data/intermediate/msi_exon_pileup_2.bcf.gz -o ../data/intermediate/msi_exon_pileup_2_sorted.bcf.gz

bcftools index ../data/intermediate/msi_exon_pileup_1_sorted.bcf.gz
bcftools index ../data/intermediate/msi_exon_pileup_2_sorted.bcf.gz

bcftools merge -Ob ../data/intermediate/msi_exon_pileup_1_sorted.bcf.gz ../data/intermediate/msi_exon_pileup_2_sorted.bcf.gz -o ../data/intermediate/msi_exon_pileup.bcf.gz
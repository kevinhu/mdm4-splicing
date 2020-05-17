# BIN_VERSION="0.10.0"

# bcftools mpileup -Ob \
# 	--bam-list wgs_bams_1.txt \
# 	-f ../data/raw/reference/hg19.fa \
# 	-R MSI_exon_bounds.bed > ../data/intermediate/msi_exon_pileup_1.bcf.gz

# bcftools mpileup -Ob \
# 	--bam-list wgs_bams_2.txt \
# 	-f ../data/raw/reference/hg19.fa \
# 	-R MSI_exon_bounds.bed > ../data/intermediate/msi_exon_pileup_2.bcf.gz

# bcftools sort -Ob ../data/intermediate/msi_exon_pileup_1.bcf.gz -o ../data/intermediate/msi_exon_pileup_1_sorted.bcf.gz
# bcftools sort -Ob ../data/intermediate/msi_exon_pileup_2.bcf.gz -o ../data/intermediate/msi_exon_pileup_2_sorted.bcf.gz

# bcftools index ../data/intermediate/msi_exon_pileup_1_sorted.bcf.gz
# bcftools index ../data/intermediate/msi_exon_pileup_2_sorted.bcf.gz

# bcftools merge -Ob ../data/intermediate/msi_exon_pileup_1_sorted.bcf.gz ../data/intermediate/msi_exon_pileup_2_sorted.bcf.gz -o ../data/intermediate/msi_exon_pileup.bcf.gz

# bcftools call -mv -Ob -o ../data/intermediate/msi_exon_calls.bcf ../data/intermediate/msi_exon_pileup.bcf.gz

# bcftools view -i '%QUAL>=20' ../data/intermediate/msi_exon_calls.bcf > ../data/intermediate/msi_exon_calls_filtered.vcf

# bcftools query ../data/intermediate/msi_exon_calls_filtered.vcf -f'[%CHROM\t%POS\t%SAMPLE\t%REF\t%ALT\t%GT\t\n]' -i'GT="alt"' > ../data/intermediate/msi_exon_calls_filtered.txt

WGS_IDS=wgs_ids.txt

printf "" > wgs_ids_samples.txt

for i in {1..329}
do
	sample_id=$(awk "NR==$i" $WGS_IDS)
	printf "$sample_id\t" >> wgs_ids_samples.txt
	samtools view -H "../data/raw/WGS_slices/${sample_id}.bam"  | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq >> wgs_ids_samples.txt
done
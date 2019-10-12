cd ../data

INDEX_PATH="raw/suppa_indices"
TPM_MATRIX="processed/tpm_matrix_suppa.txt"
IOI_OUTPUT_PATH="intermediate/suppa_transcript_psis"
IOE_OUTPUT_PATH="intermediate/suppa_exon_psis"
GTF_PATH="raw/kallisto_homo_sapiens/Homo_sapiens.GRCh38.96.gtf"

SPLICE_TYPES=("A3" "A5" "AF" "AL" "MX" "RI" "SE")

for s in ${SPLICE_TYPES[@]}; do

	suppa.py psiPerEvent --ioe-file "${INDEX_PATH}/suppa_index_ioe_${s}_strict.ioe" \
		--expression-file $TPM_MATRIX \
		-o "${IOE_OUTPUT_PATH}/${s}"

done

suppa.py psiPerIsoform -g $GTF_PATH \
	--expression-file $TPM_MATRIX \
	-o "intermediate/suppa_isoforms"

cd ../data

INDEX_PATH="raw/suppa_indices"
IOI_OUTPUT_PATH="intermediate/suppa_transcript_psis"
IOE_OUTPUT_PATH="intermediate/suppa_exon_psis"
TPM_PATH="processed/grouped_tpms"
PSI_PATH="intermediate/suppa_exon_psis_grouped"
DIFF_PATH="intermediate/suppa_diff_psis"

SPLICE_TYPES=("A3" "A5" "AF" "AL" "MX" "RI" "SE")

for s in ${SPLICE_TYPES[@]}; do

# 	suppa.py psiPerEvent --ioe-file "${INDEX_PATH}/suppa_index_ioe_${s}_strict.ioe" \
# 		--expression-file $TPM_MATRIX \
# 		-o "${IOE_OUTPUT_PATH}/${s}"

	echo "${INDEX_PATH}/suppa_index_ioe_${s}_strict.ioe"
	echo "${PSI_PATH}/${s}_RPL22_oe_c.psi"
	echo "${PSI_PATH}/${s}_RPL22_oe_t.psi"
	echo "${TPM_PATH}/RPL22_oe_c.txt"
	echo "${TPM_PATH}/RPL22_oe_t.txt"

	python3 /Users/khu/Desktop/SUPPA-master/suppa.py diffSplice --method empirical \
		--input "${INDEX_PATH}/suppa_index_ioe_${s}_strict.ioe" \
		--psi "${PSI_PATH}/${s}_RPL22_oe_c.psi" "${PSI_PATH}/${s}_RPL22_oe_t.psi" \
		--tpm "${TPM_PATH}/RPL22_oe_c.txt" "${TPM_PATH}/RPL22_oe_t.txt" \
		--alpha 0.05 \
		--lower-bound 0.01 \
		--tpm-threshold 0.1 \
		--nan-threshold 0 \
		--area 1000 \
		-o "${DIFF_PATH}/${s}_RPL22_oe"

done



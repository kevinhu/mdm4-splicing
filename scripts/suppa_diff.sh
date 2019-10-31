cd ../data

INDEX_PATH="raw/suppa_indices"
IOI_OUTPUT_PATH="intermediate/suppa_transcript_psis"
IOE_OUTPUT_PATH="intermediate/suppa_exon_psis"
TPM_PATH="processed/grouped_tpms"
PSI_PATH="intermediate/suppa_exon_psis_grouped"
ISOFORM_PSI_PATH="processed/grouped_transcript_psis"
DIFF_PATH="intermediate/suppa_diff_psis"
IOI_PATH="raw/suppa_indices/suppa_index_ioi.ioi"

SPLICE_TYPES=("A3" "A5" "AF" "AL" "MX" "RI" "SE")

for s in ${SPLICE_TYPES[@]}; do

	python3 ../SUPPA-master/suppa.py diffSplice --method empirical \
		--input "${INDEX_PATH}/suppa_index_ioe_${s}_strict.ioe" \
		--psi "${PSI_PATH}/${s}_RPL22_oe_c.psi" "${PSI_PATH}/${s}_RPL22_oe_t.psi" \
		--tpm "${TPM_PATH}/RPL22_oe_c.txt" "${TPM_PATH}/RPL22_oe_t.txt" \
		--lower-bound 0 \
		--tpm-threshold 1 \
		--tpm-threshold 0 \
		--area 1000 \
		--lower-bound 0.05 \
		-o "${DIFF_PATH}/${s}_RPL22_oe"

	# suppa.py diffSplice --method empirical \
	# 	--input "${INDEX_PATH}/suppa_index_ioe_${s}_strict.ioe" \
	# 	--psi "${PSI_PATH}/${s}_RPL22L1_oe_c.psi" "${PSI_PATH}/${s}_RPL22L1_oe_t.psi" \
	# 	--tpm "${TPM_PATH}/RPL22L1_oe_c.txt" "${TPM_PATH}/RPL22L1_oe_t.txt" \
	# 	--alpha 0.05 \
	# 	--lower-bound 0 \
	# 	--tpm-threshold 0 \
	# 	--nan-threshold 0 \
	# 	--area 1000 \
	# 	-o "${DIFF_PATH}/${s}_RPL22L1_oe"

	# suppa.py diffSplice --method empirical \
	# 	--input "${INDEX_PATH}/suppa_index_ioe_${s}_strict.ioe" \
	# 	--psi "${PSI_PATH}/${s}_shluc.psi" "${PSI_PATH}/${s}_sh704.psi" \
	# 	--tpm "${TPM_PATH}/shluc.txt" "${TPM_PATH}/sh704.txt" \
	# 	--alpha 0.05 \
	# 	--lower-bound 0 \
	# 	--tpm-threshold 0 \
	# 	--nan-threshold 0 \
	# 	--area 1000 \
	# 	-o "${DIFF_PATH}/${s}_sh704"

	# suppa.py diffSplice --method empirical \
	# 	--input "${INDEX_PATH}/suppa_index_ioe_${s}_strict.ioe" \
	# 	--psi "${PSI_PATH}/${s}_shluc.psi" "${PSI_PATH}/${s}_sh705.psi" \
	# 	--tpm "${TPM_PATH}/shluc.txt" "${TPM_PATH}/sh705.txt" \
	# 	--alpha 0.05 \
	# 	--lower-bound 0 \
	# 	--tpm-threshold 0 \
	# 	--nan-threshold 0 \
	# 	--area 1000 \
	# 	-o "${DIFF_PATH}/${s}_sh705"

done

# suppa.py diffSplice --method empirical \
# 	--input $IOI_PATH \
# 	--psi "${ISOFORM_PSI_PATH}/RPL22_oe_c.txt" "${ISOFORM_PSI_PATH}/RPL22_oe_t.txt" \
# 	--tpm "${TPM_PATH}/RPL22_oe_c.txt" "${TPM_PATH}/RPL22_oe_t.txt" \
# 	--alpha 0.05 \
# 	--lower-bound 0 \
# 	--tpm-threshold 0 \
# 	--nan-threshold 0 \
# 	--area 1000 \
# 	-o "${DIFF_PATH}/RPL22_oe"



cd ../data/raw/kallisto_homo_sapiens

suppa.py generateEvents -i Homo_sapiens.GRCh38.96.gtf \
	-o ../suppa_indices/suppa_index_ioe -f ioe \
	-e SE SS MX RI FL

suppa.py generateEvents -i Homo_sapiens.GRCh38.96.gtf \
	-o ../suppa_indices/suppa_index_ioi -f ioi
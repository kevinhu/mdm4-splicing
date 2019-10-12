FASTQ_PATH="../data/raw/zipped_fastqs"
OUTPUT_PATH="../data/intermediate/kallisto_quant"
INDEX="../data/raw/kallisto_homo_sapiens/transcriptome.idx"

# GFP overexpression in LNCaP
kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/LNCaP_GFP_1" \
	"${FASTQ_PATH}/LNCaP_GFP_1_R1.fq.gz" "${FASTQ_PATH}/LNCaP_GFP_1_R2.fq.gz"

kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/LNCaP_GFP_2" \
	"${FASTQ_PATH}/LNCaP_GFP_2_R1.fq.gz" "${FASTQ_PATH}/LNCaP_GFP_2_R2.fq.gz"

kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/LNCaP_GFP_3" \
	"${FASTQ_PATH}/LNCaP_GFP_3_R1.fq.gz" "${FASTQ_PATH}/LNCaP_GFP_3_R2.fq.gz"

# RPL22 overexpression in LNCaP
kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/LNCaP_RPL22_1" \
	"${FASTQ_PATH}/LNCaP_RPL22_1_R1.fq.gz" "${FASTQ_PATH}/LNCaP_RPL22_1_R2.fq.gz"

kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/LNCaP_RPL22_2" \
	"${FASTQ_PATH}/LNCaP_RPL22_2_R1.fq.gz" "${FASTQ_PATH}/LNCaP_RPL22_2_R2.fq.gz"

kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/LNCaP_RPL22_3" \
	"${FASTQ_PATH}/LNCaP_RPL22_3_R1.fq.gz" "${FASTQ_PATH}/LNCaP_RPL22_3_R2.fq.gz"

# GFP overexpression in CAL851
kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/CAL851_GFP_1" \
	"${FASTQ_PATH}/CAL851_GFP_1_R1.fq.gz" "${FASTQ_PATH}/CAL851_GFP_1_R2.fq.gz"

kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/CAL851_GFP_2" \
	"${FASTQ_PATH}/CAL851_GFP_2_R1.fq.gz" "${FASTQ_PATH}/CAL851_GFP_2_R2.fq.gz"

kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/CAL851_GFP_3" \
	"${FASTQ_PATH}/CAL851_GFP_3_R1.fq.gz" "${FASTQ_PATH}/CAL851_GFP_3_R2.fq.gz"

# RPL22L1 overexpression in CAL851
kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/CAL851_RPL22L1_1" \
	"${FASTQ_PATH}/CAL851_RPL22L1_1_R1.fq.gz" "${FASTQ_PATH}/CAL851_RPL22L1_1_R2.fq.gz"

kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/CAL851_RPL22L1_2" \
	"${FASTQ_PATH}/CAL851_RPL22L1_2_R1.fq.gz" "${FASTQ_PATH}/CAL851_RPL22L1_2_R2.fq.gz"

kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/CAL851_RPL22L1_3" \
	"${FASTQ_PATH}/CAL851_RPL22L1_3_R1.fq.gz" "${FASTQ_PATH}/CAL851_RPL22L1_3_R2.fq.gz"

# Luciferase knockdown in LNCaP
kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/LNCaP_shLuc_1" \
	"${FASTQ_PATH}/LNCaP_shLuc_1_R1.fq.gz" "${FASTQ_PATH}/LNCaP_shLuc_1_R2.fq.gz"

kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/LNCaP_shLuc_2" \
	"${FASTQ_PATH}/LNCaP_shLuc_2_R1.fq.gz" "${FASTQ_PATH}/LNCaP_shLuc_2_R2.fq.gz"

kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/LNCaP_shLuc_3" \
	"${FASTQ_PATH}/LNCaP_shLuc_3_R1.fq.gz" "${FASTQ_PATH}/LNCaP_shLuc_3_R2.fq.gz"

# sh704 RPL22L1 knockdown in LNCaP
kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/LNCaP_sh704_1" \
	"${FASTQ_PATH}/LNCaP_sh704_1_R1.fq.gz" "${FASTQ_PATH}/LNCaP_sh704_1_R2.fq.gz"

kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/LNCaP_sh704_2" \
	"${FASTQ_PATH}/LNCaP_sh704_2_R1.fq.gz" "${FASTQ_PATH}/LNCaP_sh704_2_R2.fq.gz"

kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/LNCaP_sh704_3" \
	"${FASTQ_PATH}/LNCaP_sh704_3_R1.fq.gz" "${FASTQ_PATH}/LNCaP_sh704_3_R2.fq.gz"

# sh705 RPL22L1 knockdown in LNCaP
kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/LNCaP_sh705_1" \
	"${FASTQ_PATH}/LNCaP_sh705_1_R1.fq.gz" "${FASTQ_PATH}/LNCaP_sh705_1_R2.fq.gz"

kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/LNCaP_sh705_2" \
	"${FASTQ_PATH}/LNCaP_sh705_2_R1.fq.gz" "${FASTQ_PATH}/LNCaP_sh705_2_R2.fq.gz"

kallisto quant -i $INDEX --bias -b 1000 -t 6 -o "${OUTPUT_PATH}/LNCaP_sh705_3" \
	"${FASTQ_PATH}/LNCaP_sh705_3_R1.fq.gz" "${FASTQ_PATH}/LNCaP_sh705_3_R2.fq.gz"
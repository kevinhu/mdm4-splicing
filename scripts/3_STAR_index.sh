STAR --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir STAR_index \
    --genomeFastaFiles \
    /broad/hptmp/khu/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa \
    --sjdbGTFfile  /broad/hptmp/khu/Homo_sapiens.GRCh37.75.gtf \
    --sjdbOverhang 100
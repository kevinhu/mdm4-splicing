## Accessing FASTQs

Zipped FASTQ files are all hosted on Pydio [here](https://distrib.dfci.harvard.edu/ws-mbcf/).

These are located under `/data/raw/zipped_fastqs/`.

To rename the files after downloading, run `rename.sh`.

## Transcript expression

### Kallisto

We used [kallisto](https://github.com/pachterlab/kallisto) for transcript expression quantification.

Hg37 indices for kallisto pseudoalignment were downloaded from [here](ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz). We use the **Ensembl transcriptomes v96** release for *Homo sapiens*.

To quantify transcript expression, run `kallisto_quant.sh`, which will output estimates to `/data/raw/kallisto_quant/`. 

### Sleuth

We used [Sleuth](https://github.com/pachterlab/sleuth) for differential expression analyses following estimation with kallisto. To run Sleuth, see `/notebooks/1_sleuth.R`.

## Splice quantification

### STAR alignment

We used the two-pass STAR alignment, which produces the BAM files for LeafCutter quantification.

First, we downloaded the hg19 FASTA reference from [here](https://console.cloud.google.com/storage/browser/_details/broad-references/hg19/v0/Homo_sapiens_assembly19.fasta), as well as the index from [here](https://console.cloud.google.com/storage/browser/_details/broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai)

Next, we build the STAR index as follows:

```bash
genomeDir=/broad/hptmp/khu/STAR_b37
fastaPath=/broad/hptmp/khu/broad_hg19/b37.fasta
gtfPath=/broad/hptmp/khu/broad_hg19/gencode.v19.genes.v7_model.patched_contigs.gtf
mkdir $genomeDir
STAR --runMode genomeGenerate \
	--genomeDir $genomeDir \
	--genomeFastaFiles $fastaPath \
	--sjdbGTFfile $gtfPath \
	--sjdbOverhang 75 \
	--runThreadN 4
```

We can now run the first-pass alignments:

```bash
fastqDir=/broad/hptmp/khu/zipped_fastqs
outputDir=/broad/hptmp/khu/STAR_pass1
genomeDir=/broad/hptmp/khu/STAR_b37

STAR --runMode alignReads \
	--genomeDir $genomeDir \
	--readFilesIn ${fastqDir}/${sample_id}_R1.fq.gz ${fastqDir}/${sample_id}_R2.fq.gz \
	--outFileNamePrefix ${outputDir}/${sample_id} \
	--outSAMtype BAM SortedByCoordinate \
	--runThreadN 4 \
```



### LeafCutter




## Experiments

To examine the role of RPL22 and RPL22L1 in splicing, we performed the following RNA-seq experiments, with triplicates for each setting:

1. Overexpression of RPL22 in LNCaP cells, which contain _RPL22_ with a p.K15fs frameshift mutation. A GFP overexpression construct was used as a control. This experiment is referenced using the `rpl22_oe` label.
2. Overexpression of RPL22L1 in CAL851 cells, which contain a wild-type RPL22 gene. A GFP overexpression construct was used as a control. This experiment is referenced using the `rpl22l1_oe` label.
3. RNA-interference (RNAi) knockdown of RPL22L1 in LNCaP cells with two different shRNAs (each with three replicates). A luciferase shRNA was used as a control. This experiment is referenced using the `rpl22l1_kd1` and `rpl22l1_kd2` labels.
4. CRISPR-Cas9 knockout of RPL22 in NCI-H2110 cells with two different sgRNAs (each with three replicates). A GFP sgRNA was used as a control. This experiment is referenced using the `rpl22_a_ko1` and `rpl22_a_ko2` labels.
5. CRISPR-Cas9 knockout of RPL22 in ZR75-1 cells with two different sgRNAs (each with three replicates, same as for NCIH-H2110). A GFP sgRNA was used as a control. This experiment is referenced using the `rpl22_b_ko1` and `rpl22_b_ko2` labels.

For programmatic processing of the experiments, we also summarized the setups in `/notebooks/experiments.json`.

## Accessing FASTQs

Zipped FASTQ files are all available on Google Drive at https://drive.google.com/drive/u/1/folders/1K4KKnPiGNBB-ZaLrUkoUue36axz1kmvs.

These should be downloaded to `/data/raw/zipped_fastqs/`. After downloading, verify that they have been transferred correctly using the MD5 sums provided in `/scripts/md5_sums.txt`. These 78 FASTQ files total about 100 GB.

## Transcript+gene expression

### Kallisto

We used [kallisto](https://github.com/pachterlab/kallisto) for transcript expression quantification.

GRCh37 cDNA sequences for constructing the Kallisto indices were downloaded from the Ensembl FTP archives at [ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz](ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz).

To construct the Kallisto index, place the cDNA FASTA file as-is into `/data/raw/kallisto_homo_sapiens`, and execute `1_kallisto_index.sh` from the `/scripts` directory.

Next, to quantify transcript expression, run `2_kallisto_quant.sh` from `/scripts/`, which will output estimates to `/data/raw/kallisto_quant/`. Kallisto is fast enough to process each paired-end run in about half an hour on most modern machines, so these can be run overnight on a laptop. We ran Kallisto with arguments `--bias -b 100 -t 6`.

### Sleuth

Kallisto outputs raw transcript quantifications and bootstrap samples for use with [Sleuth](https://github.com/pachterlab/sleuth) in differential expression analyses. To run Sleuth, see `/notebooks/p1_kallisto-sleuth.ipynb`. This will run differential expression analyses on each of the experiments, as well as annotate transcripts with matching genes and biotypes for downstream analyses.

## Splicing quantification

To compute splicing levels and perform differential splicing analyses, we aligned FASTQ files with [STAR](https://github.com/alexdobin/STAR), after which we used [rMATS](http://rnaseq-mats.sourceforge.net/) for splicing quantification.

### STAR alignment

First, we downloaded the Broad hg19 FASTA reference from [https://console.cloud.google.com/storage/browser/\_details/broad-references/hg19/v0/Homo_sapiens_assembly19.fasta](https://console.cloud.google.com/storage/browser/_details/broad-references/hg19/v0/Homo_sapiens_assembly19.fasta), and the associated index from [https://console.cloud.google.com/storage/browser/\_details/broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai](https://console.cloud.google.com/storage/browser/_details/broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai).

We also downloaded the Ensembl GTF file from [ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz](ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz)

STAR is a memory-heavy program, and building the index and alignments both require about 32-48 GB of RAM. An example script used for index construction on the Broad computing clusters is provided in `/scripts/3_STAR_index.sh`.

Once the index was built, we ran two-pass alignments on the Broad computing clusters using `4_STAR_2pass.sh`. Note the use of `/scripts/fastq_samples.txt`.

### rMATS

We used the rMATS Docker image described at http://rnaseq-mats.sourceforge.net/rmatsdockerbeta/. To run rMATS, change the `DATA_PATH` variable in `/scripts/5_run_rmats.sh` to the absolute path of the `/data` folder on your machine, and execute the script.

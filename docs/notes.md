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

### rMATS


---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.5
  kernelspec:
    display_name: R
    language: R
    name: ir
---

```R
# install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
BiocManager::install("devtools")
BiocManager::install("pachterlab/sleuth")

```

```R
# load packages
library(sleuth)
library("biomaRt")

# set working directory
setwd("../data/intermediate/kallisto_quant/")
```

```R
# define experiments
experiment_setup <- read.table(file.path("../experiment_setup.txt"), header = TRUE, stringsAsFactors=FALSE)
experiment_setup <- dplyr::mutate(experiment_setup, path = experiment_setup$sample)

# make experiment summaries
rpl22_oe_setup <- experiment_setup[c(1,2,3,4,5,6),]
rpl22l1_oe_setup <- experiment_setup[c(7,8,9,10,11,12),]
sh704_setup <- experiment_setup[c(13,14,15,16,17,18),]
sh705_setup <- experiment_setup[c(13,14,15,19,20,21),]
rpl22_a_ko1_setup <- experiment_setup[c(22,23,24,25,26,27),]
rpl22_a_ko2_setup <- experiment_setup[c(22,23,24,28,29,30),]
rpl22_b_ko1_setup <- experiment_setup[c(31,32,33,34,35,36),]
rpl22_b_ko2_setup <- experiment_setup[c(31,32,33,37,38,39),]
```

```R
# add Ensembl gene-to-transcript mappings
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         host="feb2014.archive.ensembl.org",
                         dataset = "hsapiens_gene_ensembl")

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", 
                                     "ensembl_gene_id",
                                     "ensembl_peptide_id",
                                     "hgnc_symbol",
                                     "entrezgene",
                                     "transcript_biotype"), 
                      mart = mart)

t2g <- dplyr::rename(t2g, 
                     target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, 
                     hgnc_gene = hgnc_symbol, 
                     entrez_gene = entrezgene)
t2g <- dplyr::select(t2g, c('target_id', 'ens_gene', 'hgnc_gene', 'entrez_gene', 'transcript_biotype'))

t2g$duplicate_transcript <- duplicated(t2g$target_id)

t2g <- t2g[which(!t2g$duplicate_transcript),]

write.table(t2g, file = "../sleuth_diff/ensembl_t2g.csv", sep=",",col.names=TRUE, row.names=TRUE)
```

```R
run_sleuth <- function(setup, condition_name, output_name){
  
  # make Sleuth object
  so <- sleuth_prep(setup, 
                    extra_bootstrap_summary = TRUE, 
                    read_bootstrap_tpm=TRUE,
                    target_mapping = t2g,
                    aggregation_column = 'entrez_gene'
  )
  
  # fit models
  so <- sleuth_fit(so, ~condition, 'full')
  so <- sleuth_fit(so, ~1, 'reduced')
  
  # run likelihood ratio tests
  so <- sleuth_wt(so, condition_name)
  
  # extract and save gene results
  table <- sleuth_results(so, condition_name, 'wt', show_all = TRUE)
  write.table(table, 
              file = paste("../sleuth_diff/",output_name,"_genes.csv",sep=""),
              sep=",",
              col.names=TRUE,
              row.names=TRUE
              )
  
  # extract and save transcript results
  table <- sleuth_results(so, condition_name, 'wt', show_all = TRUE, pval_aggregate=FALSE)
  write.table(table, 
              file = paste("../sleuth_diff/",output_name,"_transcripts.csv",sep=""),
              sep=",",
              col.names=TRUE,
              row.names=TRUE
              )
}
```

```R
# run sleuth
run_sleuth(rpl22_oe_setup,'conditionLNCaP_RPL22','rpl22_oe')
run_sleuth(rpl22l1_oe_setup,'conditionCAL851_RPL22L1','rpl22l1_oe')
run_sleuth(sh704_setup,'conditionLNCaP_shLuc','rpl22l1_kd1')
run_sleuth(sh705_setup,'conditionLNCaP_shLuc','rpl22l1_kd2')
run_sleuth(rpl22_a_ko1_setup,'conditionNCIH2110_RPL22-1A1','rpl22_a_ko1')
run_sleuth(rpl22_a_ko2_setup,'conditionNCIH2110_RPL22-4A1','rpl22_a_ko2')
run_sleuth(rpl22_b_ko1_setup,'conditionZR751_RPL22-1A1','rpl22_b_ko1')
run_sleuth(rpl22_b_ko2_setup,'conditionZR751_RPL22-4A1','rpl22_b_ko2')

```

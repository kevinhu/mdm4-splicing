# install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
BiocManager::install("devtools")
BiocManager::install("pachterlab/sleuth")

# load packages
library(sleuth)
library("biomaRt")

# set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../data/intermediate/kallisto_quant/")

# define experiments
experiment_setup <- read.table(file.path("../experiment_setup.txt"), header = TRUE, stringsAsFactors=FALSE)
experiment_setup <- dplyr::mutate(experiment_setup, path = experiment_setup$sample)

# make experiment summaries
rpl22_oe_setup <- experiment_setup[c(1,2,3,4,5,6),]
rpl22l1_oe_setup <- experiment_setup[c(7,8,9,10,11,12),]
sh704_setup <- experiment_setup[c(13,14,15,16,17,18),]
sh705_setup <- experiment_setup[c(13,14,15,19,20,21),]

# add Ensembl gene to transcript mappings
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         host="feb2014.archive.ensembl.org",
                         dataset = "hsapiens_gene_ensembl")

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", 
                                     "ensembl_gene_id",
                                     "ensembl_peptide_id",
                                     "hgnc_symbol"
), mart = mart)

t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = hgnc_symbol)
t2g <- dplyr::select(t2g, c('target_id', 'ens_gene', 'ext_gene'))

write.table(t2g, file = "../sleuth_diff/ensembl_t2g.csv", sep=",",col.names=TRUE, row.names=TRUE)

# import Sleuth objects
rpl22_oe_so <- sleuth_prep(rpl22_oe_setup, 
                           extra_bootstrap_summary = TRUE, 
                           read_bootstrap_tpm=TRUE,
                           target_mapping = t2g,
                           aggregation_column = 'ens_gene'
                           )
rpl22l1_oe_so <- sleuth_prep(rpl22l1_oe_setup, 
                             extra_bootstrap_summary = TRUE, 
                             read_bootstrap_tpm=TRUE,
                             target_mapping = t2g,
                             aggregation_column = 'ens_gene')
sh704_so <- sleuth_prep(sh704_setup, 
                        extra_bootstrap_summary = TRUE, 
                        read_bootstrap_tpm=TRUE,
                        target_mapping = t2g,
                        aggregation_column = 'ens_gene')
sh705_so <- sleuth_prep(sh705_setup, 
                        extra_bootstrap_summary = TRUE, 
                        read_bootstrap_tpm=TRUE,
                        target_mapping = t2g,
                        aggregation_column = 'ens_gene')

# likelihood ratio tests
rpl22_oe_so <- sleuth_fit(rpl22_oe_so, ~condition, 'full')
rpl22_oe_so <- sleuth_fit(rpl22_oe_so, ~1, 'reduced')
rpl22_oe_so <- sleuth_lrt(rpl22_oe_so, 'reduced', 'full')
rpl22_oe_table <- sleuth_results(rpl22_oe_so, 'reduced:full', 'lrt', show_all = TRUE)
write.table(rpl22_oe_table, file = "../kallisto_sleuth/rpl22_oe.csv", sep=",",col.names=TRUE, row.names=TRUE)

rpl22l1_oe_so <- sleuth_fit(rpl22l1_oe_so, ~condition, 'full')
rpl22l1_oe_so <- sleuth_fit(rpl22l1_oe_so, ~1, 'reduced')
rpl22l1_oe_so <- sleuth_lrt(rpl22l1_oe_so, 'reduced', 'full')
rpl22l1_oe_table <- sleuth_results(rpl22l1_oe_so, 'reduced:full', 'lrt', show_all = TRUE)
write.table(rpl22l1_oe_table, file = "../kallisto_sleuth/rpl22l1_oe.csv", sep=",",col.names=TRUE, row.names=TRUE)

sh704_so <- sleuth_fit(sh704_so, ~condition, 'full')
sh704_so <- sleuth_fit(sh704_so, ~1, 'reduced')
sh704_so <- sleuth_lrt(sh704_so, 'reduced', 'full')
sh704_table <- sleuth_results(sh704_so, 'reduced:full', 'lrt', show_all = TRUE)
write.table(sh704_table, file = "../kallisto_sleuth/sh704.csv", sep=",",col.names=TRUE, row.names=TRUE)

sh705_so <- sleuth_fit(sh705_so, ~condition, 'full')
sh705_so <- sleuth_fit(sh705_so, ~1, 'reduced')
sh705_so <- sleuth_lrt(sh705_so, 'reduced', 'full')
sh705_table <- sleuth_results(sh705_so, 'reduced:full', 'lrt', show_all = TRUE)
write.table(sh705_table, file = "../kallisto_sleuth/sh705.csv", sep=",",col.names=TRUE, row.names=TRUE)

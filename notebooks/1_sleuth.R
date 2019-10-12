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
RPL22_oe_setup <- experiment_setup[c(1,2,3,4,5,6),]
RPL22L1_oe_setup <- experiment_setup[c(7,8,9,10,11,12),]
sh704_setup <- experiment_setup[c(13,14,15,16,17,18),]
sh705_setup <- experiment_setup[c(13,14,15,19,20,21),]

# import Sleuth objects
RPL22_oe_so <- sleuth_prep(RPL22_oe_setup, extra_bootstrap_summary = TRUE, read_bootstrap_tpm=TRUE)
RPL22L1_oe_so <- sleuth_prep(RPL22L1_oe_setup, extra_bootstrap_summary = TRUE, read_bootstrap_tpm=TRUE)
sh704_so <- sleuth_prep(sh704_setup, extra_bootstrap_summary = TRUE, read_bootstrap_tpm=TRUE)
sh705_so <- sleuth_prep(sh705_setup, extra_bootstrap_summary = TRUE, read_bootstrap_tpm=TRUE)

# likelihood ratio tests
RPL22_oe_so <- sleuth_fit(RPL22_oe_so, ~condition, 'full')
RPL22_oe_so <- sleuth_fit(RPL22_oe_so, ~1, 'reduced')
RPL22_oe_so <- sleuth_lrt(RPL22_oe_so, 'reduced', 'full')
RPL22_oe_table <- sleuth_results(RPL22_oe_so, 'reduced:full', 'lrt', show_all = TRUE)
write.table(RPL22_oe_table, file = "../kallisto_sleuth/RPL22_oe.csv", sep=",",col.names=TRUE, row.names=TRUE)

RPL22L1_oe_so <- sleuth_fit(RPL22L1_oe_so, ~condition, 'full')
RPL22L1_oe_so <- sleuth_fit(RPL22L1_oe_so, ~1, 'reduced')
RPL22L1_oe_so <- sleuth_lrt(RPL22L1_oe_so, 'reduced', 'full')
RPL22L1_oe_table <- sleuth_results(RPL22L1_oe_so, 'reduced:full', 'lrt', show_all = TRUE)
write.table(RPL22L1_oe_table, file = "../kallisto_sleuth/RPL22L1_oe.csv", sep=",",col.names=TRUE, row.names=TRUE)

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

# add Ensembl gene to transcript mappings
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

write.table(t2g, file = "../kallisto_sleuth/ensembl_t2g.csv", sep=",",col.names=TRUE, row.names=TRUE)

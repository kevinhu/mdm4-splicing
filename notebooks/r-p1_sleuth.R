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
rpl22_1a1_setup <- experiment_setup[c(22,23,24,25,26,27),]
rpl22_4a1_setup <- experiment_setup[c(22,23,24,28,29,30),]

# add Ensembl gene to transcript mappings
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         host="feb2014.archive.ensembl.org",
                         dataset = "hsapiens_gene_ensembl")

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", 
                                     "ensembl_gene_id",
                                     "ensembl_peptide_id",
                                     "hgnc_symbol",
                                     "entrezgene"
), mart = mart)

t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, hgnc_gene = hgnc_symbol, entrez_gene = entrezgene)
t2g <- dplyr::select(t2g, c('target_id', 'ens_gene', 'hgnc_gene', 'entrez_gene'))

t2g$duplicate_transcript <- duplicated(t2g$target_id)

t2g <- t2g[which(!t2g$duplicate_transcript),]

write.table(t2g, file = "../sleuth_diff/ensembl_t2g.csv", sep=",",col.names=TRUE, row.names=TRUE)

# import Sleuth objects
rpl22_oe_so <- sleuth_prep(rpl22_oe_setup, 
                           extra_bootstrap_summary = TRUE, 
                           read_bootstrap_tpm=TRUE,
                           target_mapping = t2g,
                           aggregation_column = 'entrez_gene'
                           )

rpl22l1_oe_so <- sleuth_prep(rpl22l1_oe_setup, 
                             extra_bootstrap_summary = TRUE, 
                             read_bootstrap_tpm=TRUE,
                             target_mapping = t2g,
                             aggregation_column = 'entrez_gene')

sh704_so <- sleuth_prep(sh704_setup, 
                        extra_bootstrap_summary = TRUE, 
                        read_bootstrap_tpm=TRUE,
                        target_mapping = t2g,
                        aggregation_column = 'entrez_gene')

sh705_so <- sleuth_prep(sh705_setup, 
                        extra_bootstrap_summary = TRUE, 
                        read_bootstrap_tpm=TRUE,
                        target_mapping = t2g,
                        aggregation_column = 'entrez_gene')

rpl22_1a1_so <- sleuth_prep(rpl22_1a1_setup, 
                        extra_bootstrap_summary = TRUE, 
                        read_bootstrap_tpm=TRUE,
                        target_mapping = t2g,
                        aggregation_column = 'entrez_gene')

rpl22_4a1_so <- sleuth_prep(rpl22_4a1_setup, 
                            extra_bootstrap_summary = TRUE, 
                            read_bootstrap_tpm=TRUE,
                            target_mapping = t2g,
                            aggregation_column = 'entrez_gene')

# likelihood ratio tests
rpl22_oe_so <- sleuth_fit(rpl22_oe_so, ~condition, 'full')
rpl22_oe_so <- sleuth_fit(rpl22_oe_so, ~1, 'reduced')
rpl22_oe_so <- sleuth_wt(rpl22_oe_so, 'conditionLNCaP_RPL22')
rpl22_oe_table <- sleuth_results(rpl22_oe_so, 'conditionLNCaP_RPL22', 'wt', show_all = TRUE)
write.table(rpl22_oe_table, file = "../sleuth_diff/rpl22_oe_genes.csv", sep=",",col.names=TRUE, row.names=TRUE)
rpl22_oe_table <- sleuth_results(rpl22_oe_so, 'conditionLNCaP_RPL22', 'wt', show_all = TRUE, pval_aggregate=FALSE)
write.table(rpl22_oe_table, file = "../sleuth_diff/rpl22_oe_transcripts.csv", sep=",",col.names=TRUE, row.names=TRUE)

rpl22l1_oe_so <- sleuth_fit(rpl22l1_oe_so, ~condition, 'full')
rpl22l1_oe_so <- sleuth_fit(rpl22l1_oe_so, ~1, 'reduced')
rpl22l1_oe_so <- sleuth_wt(rpl22l1_oe_so, 'conditionCAL851_RPL22L1')
rpl22l1_oe_table <- sleuth_results(rpl22l1_oe_so, 'conditionCAL851_RPL22L1', 'wt', show_all = TRUE)
write.table(rpl22l1_oe_table, file = "../sleuth_diff/rpl22l1_oe_genes.csv", sep=",",col.names=TRUE, row.names=TRUE)
rpl22l1_oe_so <- sleuth_results(rpl22l1_oe_so, 'conditionCAL851_RPL22L1', 'wt', show_all = TRUE, pval_aggregate=FALSE)
write.table(rpl22l1_oe_so, file = "../sleuth_diff/rpl22l1_oe_transcripts.csv", sep=",",col.names=TRUE, row.names=TRUE)

sh704_so <- sleuth_fit(sh704_so, ~condition, 'full')
sh704_so <- sleuth_fit(sh704_so, ~1, 'reduced')
sh704_so <- sleuth_wt(sh704_so, 'conditionLNCaP_shLuc')
sh704_table <- sleuth_results(sh704_so, 'conditionLNCaP_shLuc', 'wt', show_all = TRUE)
write.table(sh704_table, file = "../sleuth_diff/rpl22l1_kd1_genes.csv", sep=",",col.names=TRUE, row.names=TRUE)
sh704_table <- sleuth_results(sh704_so, 'conditionLNCaP_shLuc', 'wt', show_all = TRUE, pval_aggregate=FALSE)
write.table(sh704_table, file = "../sleuth_diff/rpl22l1_kd1_transcripts.csv", sep=",",col.names=TRUE, row.names=TRUE)

sh705_so <- sleuth_fit(sh705_so, ~condition, 'full')
sh705_so <- sleuth_fit(sh705_so, ~1, 'reduced')
sh705_so <- sleuth_wt(sh705_so, 'conditionLNCaP_shLuc')
sh705_table <- sleuth_results(sh705_so, 'conditionLNCaP_shLuc', 'wt', show_all = TRUE)
write.table(sh705_table, file = "../sleuth_diff/rpl22l1_kd2_genes.csv", sep=",",col.names=TRUE, row.names=TRUE)
sh705_table <- sleuth_results(sh705_so, 'conditionLNCaP_shLuc', 'wt', show_all = TRUE, pval_aggregate=FALSE)
write.table(sh705_table, file = "../sleuth_diff/rpl22l1_kd2_transcripts.csv", sep=",",col.names=TRUE, row.names=TRUE)

rpl22_1a1_so <- sleuth_fit(rpl22_1a1_so, ~condition, 'full')
rpl22_1a1_so <- sleuth_fit(rpl22_1a1_so, ~1, 'reduced')
rpl22_1a1_so <- sleuth_wt(rpl22_1a1_so, 'conditionNCIH2110_RPL22-1A1')
rpl22_1a1_table <- sleuth_results(rpl22_1a1_so, 'conditionNCIH2110_RPL22-1A1', 'wt', show_all = TRUE)
write.table(rpl22_1a1_table, file = "../sleuth_diff/rpl22_ko1_genes.csv", sep=",",col.names=TRUE, row.names=TRUE)
rpl22_1a1_table <- sleuth_results(rpl22_1a1_so, 'conditionNCIH2110_RPL22-1A1', 'wt', show_all = TRUE, pval_aggregate=FALSE)
write.table(rpl22_1a1_table, file = "../sleuth_diff/rpl22_ko1_transcripts.csv", sep=",",col.names=TRUE, row.names=TRUE)

rpl22_4a1_so <- sleuth_fit(rpl22_4a1_so, ~condition, 'full')
rpl22_4a1_so <- sleuth_fit(rpl22_4a1_so, ~1, 'reduced')
rpl22_4a1_so <- sleuth_wt(rpl22_4a1_so, 'conditionNCIH2110_RPL22-4A1')
rpl22_4a1_table <- sleuth_results(rpl22_4a1_so, 'conditionNCIH2110_RPL22-4A1', 'wt', show_all = TRUE)
write.table(rpl22_4a1_table, file = "../sleuth_diff/rpl22_ko2_genes.csv", sep=",",col.names=TRUE, row.names=TRUE)
rpl22_4a1_table <- sleuth_results(rpl22_4a1_so, 'conditionNCIH2110_RPL22-4A1', 'wt', show_all = TRUE, pval_aggregate=FALSE)
write.table(rpl22_4a1_table, file = "../sleuth_diff/rpl22_ko2_transcripts.csv", sep=",",col.names=TRUE, row.names=TRUE)

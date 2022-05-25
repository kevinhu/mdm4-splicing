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

```R vscode={"languageId": "r"}
library(data.table)
library(fgsea)
library(ggplot2)

setwd("../data/processed/kallisto_sleuth_merge/")
```

```R vscode={"languageId": "r"}
pathways.hallmark <- gmtPathways("../../raw/msigdb/h.all.v7.0.entrez.gmt")
pathways.curated <- gmtPathways("../../raw/msigdb/c2.all.v7.0.entrez.gmt")
pathways.motif  <- gmtPathways("../../raw/msigdb/c3.all.v7.0.entrez.gmt")
pathways.go  <- gmtPathways("../../raw/msigdb/c5.all.v7.0.entrez.gmt")
pathways.oncogenic  <- gmtPathways("../../raw/msigdb/c6.all.v7.0.entrez.gmt")

# select pathways for analysis
all_pathways <- c(pathways.hallmark,pathways.go)
```

```R vscode={"languageId": "r"}
# load sleuth outputs
rpl22_oe <- read.table("rpl22_oe_genes.txt",sep="\t", header=TRUE)
rpl22l1_oe <- read.table("rpl22l1_oe_genes.txt",sep="\t", header=TRUE)
rpl22l1_kd1 <- read.table("rpl22l1_kd1_genes.txt",sep="\t", header=TRUE)
rpl22l1_kd2 <- read.table("rpl22l1_kd2_genes.txt",sep="\t", header=TRUE)
rpl22_a_ko1 <- read.table("rpl22_a_ko1_genes.txt",sep="\t", header=TRUE)
rpl22_a_ko2 <- read.table("rpl22_a_ko2_genes.txt",sep="\t", header=TRUE)
rpl22_b_ko1 <- read.table("rpl22_b_ko1_genes.txt",sep="\t", header=TRUE)
rpl22_b_ko2 <- read.table("rpl22_b_ko2_genes.txt",sep="\t", header=TRUE)

# drop Entrez duplicates, keep most significant
drop_entrez_duplicates <- function(sleuth_diff){
  sleuth_diff <- sleuth_diff[order(sleuth_diff$pval),]
  sleuth_diff <- sleuth_diff[!duplicated(sleuth_diff$target_id),]
  sleuth_diff <- sleuth_diff[order(-sleuth_diff$signed_pval),]
  
  return(sleuth_diff)
}

rpl22_oe <- drop_entrez_duplicates(rpl22_oe)
rpl22l1_oe <- drop_entrez_duplicates(rpl22l1_oe)
rpl22l1_kd1 <- drop_entrez_duplicates(rpl22l1_kd1)
rpl22l1_kd2 <- drop_entrez_duplicates(rpl22l1_kd2)
rpl22_a_ko1 <- drop_entrez_duplicates(rpl22_a_ko1)
rpl22_a_ko2 <- drop_entrez_duplicates(rpl22_a_ko2)
rpl22_b_ko1 <- drop_entrez_duplicates(rpl22_b_ko1)
rpl22_b_ko2 <- drop_entrez_duplicates(rpl22_b_ko2)
```

```R vscode={"languageId": "r"}
# put sets into a list for iteration
rank_sets = list(rpl22_oe,
                 rpl22l1_oe,
                 rpl22l1_kd1,
                 rpl22l1_kd2,
                 rpl22_a_ko1,
                 rpl22_a_ko2,
                 rpl22_b_ko1,
                 rpl22_b_ko2)

rank_set_names = list("rpl22_oe",
                      "rpl22l1_oe",
                      "rpl22l1_kd1",
                      "rpl22l1_kd2",
                      "rpl22_a_ko1",
                      "rpl22_a_ko2",
                      "rpl22_b_ko1",
                      "rpl22_b_ko2")
```

```R vscode={"languageId": "r"}
# helper function for running fgsea
run_fgsea <- function(rank_set, results_file){

  ranks <- setNames(rank_set$signed_pval, rank_set$target_id)
  fgseaRes <- fgseaMultilevel(pathways = all_pathways, 
                              stats = ranks,
                              minSize=5,
                              maxSize=500)
  
  fgseaRes$leadingEdge <- vapply(fgseaRes$leadingEdge, paste, collapse = ",", character(1L))
  
  write.table(fgseaRes, 
              file = paste("../../processed/fgsea_results/",results_file,sep=""), 
              sep="\t",
              col.names=TRUE, 
              row.names=TRUE)
}
```

```R vscode={"languageId": "r"}
# execute fgsea for each set
for(rank_set_idx in 1:length(rank_sets)){
    
  rank_set = rank_sets[[rank_set_idx]]
  rank_set_name = rank_set_names[[rank_set_idx]]
  
  run_fgsea(rank_set,paste(rank_set_name,".txt",sep=""))
    
}
```

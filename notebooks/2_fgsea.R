library(data.table)
library(fgsea)
library(ggplot2)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../data/processed/kallisto_sleuth_merge/")

pathways.hallmark <- gmtPathways("../../raw/msigdb/h.all.v7.0.entrez.gmt")
pathways.curated <- gmtPathways("../../raw/msigdb/c2.all.v7.0.entrez.gmt")
pathways.motif  <- gmtPathways("../../raw/msigdb/c3.all.v7.0.entrez.gmt")
pathways.go  <- gmtPathways("../../raw/msigdb/c5.all.v7.0.entrez.gmt")
pathways.oncogenic  <- gmtPathways("../../raw/msigdb/c6.all.v7.0.entrez.gmt")

all_pathways <- c(pathways.hallmark,pathways.curated,pathways.motif,pathways.go,pathways.oncogenic)

rpl22_oe <- read.table("rpl22_oe_genes.txt",sep="\t", header=TRUE)
rpl22l1_oe <- read.table("rpl22l1_oe_genes.txt",sep="\t", header=TRUE)
sh704 <- read.table("sh704_genes.txt",sep="\t", header=TRUE)
sh705 <- read.table("sh705_genes.txt",sep="\t", header=TRUE)

# drop Entrez duplicates, keep most significant
drop_entrez_duplicates <- function(sleuth_diff){
  sleuth_diff <- sleuth_diff[order(sleuth_diff$pval),]
  sleuth_diff <- sleuth_diff[!duplicated(sleuth_diff$target_id),]
  sleuth_diff <- sleuth_diff[order(-sleuth_diff$signed_pval),]
  
  return(sleuth_diff)
}

rpl22_oe <- drop_entrez_duplicates(rpl22_oe)
rpl22l1_oe <- drop_entrez_duplicates(rpl22l1_oe)
sh704 <- drop_entrez_duplicates(sh704)
sh705 <- drop_entrez_duplicates(sh705)

ranks <- setNames(sh705$signed_pval, sh705$target_id)

fgseaRes <- fgseaMultilevel(pathways = pathways.hallmark, 
                            stats = ranks,
                            minSize=15,
                            maxSize=1000)

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(all_pathways[topPathways], ranks, fgseaRes, 
              gseaParam = 0.5)

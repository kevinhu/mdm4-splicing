library(fgsea)
library(ggplot2)
library(tibble)
library(BiocParallel)
library(data.table)

setwd("~/Desktop/github/RPL22L1_analysis/data/gsea/")

pathways.hallmark <- gmtPathways("../msigdb/h.all.v6.2.symbols.gmt")
pathways.curated <- gmtPathways("../msigdb/c2.all.v6.2.symbols.gmt")
pathways.motif  <- gmtPathways("../msigdb/c3.all.v6.2.symbols.gmt")
pathways.go  <- gmtPathways("../msigdb/c5.all.v6.2.symbols.gmt")
pathways.oncogenic  <- gmtPathways("../msigdb/c6.all.v6.2.symbols.gmt")


all_pathways <- c(pathways.hallmark,pathways.curated,pathways.motif,pathways.go,pathways.oncogenic)
pathway_names <- list("hallmark","curated","motif","go","oncogenic")

rpl22_oe <- read.table("rpl22_oe_ranks.csv",sep=",", header=TRUE)
rpl22l1_oe <- read.table("rpl22l1_oe_ranks.csv",sep=",", header=TRUE)
sh704 <- read.table("sh704_ranks.csv",sep=",", header=TRUE)
sh705 <- read.table("sh705_ranks.csv",sep=",", header=TRUE)

all_ranks <- list(rpl22_oe,rpl22l1_oe,sh704,sh705)
ranked_names <- list("rpl22_oe","rpl22l1_oe","sh704","sh705")

for (rank_idx in c(1,2,3,4)){
  
  rank_idx = as.integer(rank_idx)
  
  pathways <- all_pathways[[pathway_idx]]
  
  ranks <- all_ranks[[rank_idx]]
  
  ranks <- setNames(ranks$b, ranks$hgnc_symbol)
  
  fgseaRes <- fgseaMultilevel(pathways = all_pathways, 
                              stats = ranks,
                              minSize=15,
                              maxSize=500)
  
  fwrite(fgseaRes, 
         file=paste("fgsea_results/",ranked_names[rank_idx],".tsv",sep=""),
         sep="\t", 
         sep2=c("", " ", ""))
}

ranks <- setNames(sh705$b, sh705$hgnc_symbol)
fgseaRes <- fgseaMultilevel(pathways = pathways.curated, 
                            stats = ranks,
                            minSize=15,
                            maxSize=500)


topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways.curated[topPathways], ranks, fgseaRes, 
              gseaParam = 0.5)

pdf("test.pdf", width=8, height=8)
plotEnrichment(all_pathways[["FISCHER_G2_M_CELL_CYCLE"]],
               ranks) + labs(title="G2-M checkpoint")
dev.off()




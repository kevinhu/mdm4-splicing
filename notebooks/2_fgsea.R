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

ranks <- setNames(rpl22l1_oe$signed_pval, rpl22l1_oe$target_id)

fgseaRes <- fgseaMultilevel(pathways = pathways.go, 
                            stats = ranks,
                            minSize=15,
                            maxSize=500)

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(all_pathways[topPathways], ranks, fgseaRes, 
              gseaParam = 0.5)

# plot fgsea results
plot_fgsea <- function(rank_set, gene_set, title, plot_file){
  ranks <- setNames(rank_set$signed_pval, rank_set$target_id)
  p <- plotEnrichment(all_pathways[[gene_set]],
                 ranks) + labs(title=title,
                               ticksSize=0.01
                               )
  pdf(paste("../../../plots/",plot_file,sep=""),width=3.2,height=2.5,paper='special')
  print(p)
  dev.off()
}

rank_sets = list(rpl22_oe,rpl22l1_oe,sh704,sh705)
rank_set_names = list("rpl22_oe","rpl22l1_oe","sh704","sh705")
gene_sets = list("HALLMARK_P53_PATHWAY",
              "HALLMARK_G2M_CHECKPOINT",
              "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
              "GO_MRNA_METABOLIC_PROCESS",
              "GO_RNA_SPLICING",
              "GO_EXTRACELLULAR_MATRIX",
              "GO_SPLICEOSOMAL_COMPLEX"
              )

for(rank_set_idx in list(1,2,3,4)){
  for(gene_set in gene_sets){
    
    rank_set = rank_sets[[rank_set_idx]]
    rank_set_name = rank_set_names[[rank_set_idx]]
    
    output_file = paste(rank_set_name,"_",gene_set,".pdf",sep="")
    
    plot_fgsea(rank_set,gene_set,gene_set,output_file)
    
  }
}
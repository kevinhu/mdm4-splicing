library(PathwaySplice)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

rpl22_oe <- read.table("../data/processed/rmats_merge/rpl22_oe.txt",sep="\t", header=TRUE, stringsAsFactors=FALSE)
rpl22l1_oe <- read.table("../data/processed/rmats_merge/rpl22l1_oe.txt",sep="\t", header=TRUE, stringsAsFactors=FALSE)
rpl22l1_kd1 <- read.table("../data/processed/rmats_merge/rpl22l1_kd1.txt",sep="\t", header=TRUE, stringsAsFactors=FALSE)
rpl22l1_kd2 <- read.table("../data/processed/rmats_merge/rpl22l1_kd2.txt",sep="\t", header=TRUE, stringsAsFactors=FALSE)
rpl22_a_ko1 <- read.table("../data/processed/rmats_merge/rpl22_a_ko1.txt",sep="\t", header=TRUE, stringsAsFactors=FALSE)
rpl22_a_ko2 <- read.table("../data/processed/rmats_merge/rpl22_a_ko2.txt",sep="\t", header=TRUE, stringsAsFactors=FALSE)
rpl22_b_ko1 <- read.table("../data/processed/rmats_merge/rpl22_b_ko1.txt",sep="\t", header=TRUE, stringsAsFactors=FALSE)
rpl22_b_ko2 <- read.table("../data/processed/rmats_merge/rpl22_b_ko2.txt",sep="\t", header=TRUE, stringsAsFactors=FALSE)

make_gene_table <- function(splicing_table){
  splicing_table$geneID = splicing_table$GeneID
  splicing_table$pvalue = splicing_table$PValue
  splicing_table$fdr = splicing_table$FDR
  
  return(makeGeneTable(splicing_table, stat = "fdr"))
}

rpl22_oe_reduced <- make_gene_table(rpl22_oe)
rpl22l1_oe_reduced <- make_gene_table(rpl22l1_oe)
rpl22l1_kd1_reduced <- make_gene_table(rpl22l1_kd1)
rpl22l1_kd2_reduced <- make_gene_table(rpl22l1_kd2)
rpl22_a_ko1_reduced <- make_gene_table(rpl22_a_ko1)
rpl22_a_ko2_reduced <- make_gene_table(rpl22_a_ko2)
rpl22_b_ko1_reduced <- make_gene_table(rpl22_b_ko1)
rpl22_b_ko2_reduced <- make_gene_table(rpl22_b_ko2)


dir.name <- system.file("extdata", package = "PathwaySplice")
hallmark.local.pathways <- file.path(dir.name, "h.all.v6.0.symbols.gmt.txt")
hlp <- gmtGene2Cat(hallmark.local.pathways, genomeID = "hg19")


run_splice <- function(gene_table){
  res <- runPathwaySplice(gene_table,
                          genome='hg19',
                          id='ensGene',
                          gene2cat = hlp, 
                          # test.cats=c("GO:CC", "GO:BP", "GO:MF"),
                          go.size.limit=c(10,500),
                          method='Wallenius',
                          use.genes.without.cat = TRUE)
  
  return(res)
}

rpl22_oe_res <- run_splice(rpl22_oe_reduced)
rpl22l1_oe_res <- run_splice(rpl22l1_oe_reduced)
rpl22l1_kd1_res <- run_splice(rpl22l1_kd1_reduced)
rpl22l1_kd2_res <- run_splice(rpl22l1_kd2_reduced)
rpl22_a_ko1_res <- run_splice(rpl22_a_ko1_reduced)
rpl22_a_ko2_res <- run_splice(rpl22_a_ko2_reduced)
rpl22_b_ko1_res <- run_splice(rpl22_b_ko1_reduced)
rpl22_b_ko2_res <- run_splice(rpl22_b_ko2_reduced)

enmap <- enrichmentMap(pathway.res = rpl22l1_oe_res,
                       n = 2,
                       output.file.dir = tempdir(),
                       similarity.threshold = 0.5, 
                       scaling.factor = 2)

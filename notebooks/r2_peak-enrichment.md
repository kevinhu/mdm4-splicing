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

See https://hbctraining.github.io/Intro-to-ChIPseq/lessons/12_functional_analysis.html

```R
BiocManager::install("ChIPseeker")
BiocManager::install("GO.db")
BiocManager::install("DO.db")
BiocManager::install("EnsDb.Hsapiens.v75")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
```

```R
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v75)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
```

```R
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakfile <- "../data/raw/eclip/RPL22-ZR751_clearCLIP.pool.tag.uniq.del.CIMS.fdr10.f10.bed"
peakAnno <- annotatePeak(peakfile, tssRegion=c(-1000, 1000), TxDb=txdb)
```

```R
pdf("../plots/peak_regions.pdf", width=6, height=3)
plotAnnoPie(peakAnno)
dev.off()
plotAnnoPie(peakAnno)
```

```R
peakAnnoDf <- data.frame(peakAnno@anno)
peakEntrez <- peakAnnoDf$geneId
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v75,
                                         keys = peakEntrez,
                                         columns = c("GENENAME","ENTREZID","GENEID"),
                                         keytype = "ENTREZID")
annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)

write.csv(peakAnnoDf,"../data/raw/peak_anno.csv", row.names = TRUE)
write.csv(annotations_edb,"../data/raw/peak_gene_ids.csv", row.names = TRUE)
```

```R
ego <- enrichGO(gene = peakEntrez, 
                    keyType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE)

pdf("../plots/go_clip_enrichment.pdf", width=6, height=4)
dotplot(ego, showCategory=10)+theme_light()
dev.off()
```

```R
ekegg <- enrichKEGG(gene = peakEntrez,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)

dotplot(ekegg)
```

```R
gene_matrix <- getTagMatrix(peak = peakAnno@anno, 
                           TxDb = txdb,
                           upstream = 1000,
                           downstream = 1000, 
                           by = "gene",
                           type = "end_site",
                           nbin=250
                           )

pdf("../plots/tts_enrichment.pdf", width=5, height=3) 
plotPeakProf(gene_matrix, xlab="Gene end (5'->3')", ylab = "Peak Frequency", conf = 0.95, resample = 1000) + theme_minimal()
dev.off()
```

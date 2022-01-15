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
library(biomaRt)
ens_human = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                    host="https://useast.ensembl.org",
                    path="/biomart/martservice",
                    dataset="hsapiens_gene_ensembl")

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

write.csv(listAttributes(mart=ens_human),"../data/raw/ensembl_attributes.csv")
```

```R
# all genes associated with "RNA splicing"

splicing_genes <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'go_id'),
                   filters = 'go', values = c('GO:0008380'), mart = ens_human)
write.csv(splicing_genes,"../data/raw/go_rna_splicing_genes.csv", row.names = TRUE)
```

```R
hgid <- getBM(attributes = "ensembl_gene_id",
              filters    = "with_hsapiens_paralog",
              values     = TRUE,
              mart       = ens_human)
```

```R
para <- getBM(attributes = c("ensembl_gene_id", 
                             "external_gene_name",
                             "entrezgene_id",
                             "hsapiens_paralog_ensembl_gene", 
                             "hsapiens_paralog_associated_gene_name"
                            ),
              filters    = "ensembl_gene_id",
              values     = hgid,
              mart       = ens_human)
```

```R
entrez <- getBM(attributes = c("ensembl_gene_id", 
                               "entrezgene_id",
                               "hgnc_id",
                               "hgnc_symbol"
                            ),
              filters    = "ensembl_gene_id",
              values     = hgid,
              mart       = ens_human)
```

```R
write.csv(para,"../data/raw/ensembl_paralogs.csv", row.names = TRUE)
write.csv(entrez,"../data/raw/ensembl_entrez.csv", row.names = TRUE)
```

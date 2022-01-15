---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.5
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

```python
import pandas as pd
```

```python
peak_anno = pd.read_csv("../data/raw/peak_anno.csv", index_col=0)
annotations_edb = pd.read_csv("../data/raw/peak_gene_ids.csv", index_col=0)

entrez_to_ensembl = dict(zip(annotations_edb["ENTREZID"], annotations_edb["GENEID"]))
entrez_to_name = dict(zip(annotations_edb["ENTREZID"], annotations_edb["GENENAME"]))

peak_anno["ensembl_gene"] = peak_anno["geneId"].apply(entrez_to_ensembl.get)
peak_anno["gene_name"] = peak_anno["geneId"].apply(entrez_to_name.get)

peak_ensembl_genes = set(peak_anno["ensembl_gene"])
```

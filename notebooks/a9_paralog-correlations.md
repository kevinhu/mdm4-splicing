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
import numpy
import many

import seaborn
import matplotlib.pyplot as plt
import cancer_data

from scipy.stats import pearsonr, spearmanr

from tqdm.notebook import tqdm_notebook
tqdm_notebook.pandas()
```

```python
ensembl_paralogs = pd.read_csv("../data/raw/ensembl_paralogs.csv", index_col=0)
```

```python
ccle_genex = cancer_data.load("ccle_gene_tpm")
ccle_exonusage = cancer_data.load("ccle_exonusage")
```

```python
genex_genes = pd.DataFrame(pd.Series(ccle_genex.columns, name="id"))
genex_genes["ensembl_id"] = (
    genex_genes["id"].apply(lambda x: x.split("_")[-1]).apply(lambda x: x.split(".")[0])
)

exonusage_genes = pd.DataFrame(pd.Series(ccle_exonusage.columns, name="id"))
exonusage_genes["ensembl_id"] = (
    exonusage_genes["id"].apply(lambda x: x.split("_")[-1]).apply(lambda x: x.split(".")[0])
)
```

```python
ensembl_to_genex = dict(zip(genex_genes["ensembl_id"],genex_genes["id"]))
ensembl_to_exonusage = exonusage_genes.groupby("ensembl_id")["id"].apply(list).to_dict()
```

```python
ensembl_paralogs["ccle_gene"] = ensembl_paralogs["ensembl_gene_id"].apply(
    ensembl_to_genex.get
)
ensembl_paralogs["ccle_paralog_exons"] = ensembl_paralogs[
    "hsapiens_paralog_ensembl_gene"
].apply(ensembl_to_exonusage.get)
ensembl_paralogs = ensembl_paralogs.dropna(subset=["ccle_gene", "ccle_paralog_exons"])
ensembl_paralogs = ensembl_paralogs.explode("ccle_paralog_exons")
ensembl_paralogs.reset_index(inplace=True)
```

```python
def compute_paralog_correlation(row):

    gene = row["ccle_gene"]
    exon = row["ccle_paralog_exons"]

    gene_genex = ccle_genex[gene].dropna()
    exon_exonusage = ccle_exonusage[exon].dropna()

    gene_genex, exon_exonusage = gene_genex.align(exon_exonusage, join="inner")

    intersection_size = gene_genex.size

    if intersection_size < 50:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)

    spearman_r, spearman_pval = spearmanr(gene_genex, exon_exonusage)
    pearson_r, pearson_pval = pearsonr(gene_genex, exon_exonusage)

    return (intersection_size, spearman_r, spearman_pval, pearson_r, pearson_pval)
```

```python
(
    ensembl_paralogs["num_samples"],
    ensembl_paralogs["spearman_r"],
    ensembl_paralogs["spearman_pval"],
    ensembl_paralogs["pearson_r"],
    ensembl_paralogs["pearson_pval"],
) = zip(*ensembl_paralogs.apply(compute_paralog_correlation, axis=1))
```

```python

```

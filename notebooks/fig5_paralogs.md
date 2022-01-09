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
import matplotlib.pyplot as plt
import numpy as np

import many
import cancer_data

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

from statsmodels.stats.multitest import multipletests
```

```python
depmap_damaging = cancer_data.load("depmap_damaging")
depmap_hotspot = cancer_data.load("depmap_hotspot")
ccle_exonusage = cancer_data.load("ccle_exonusage")
```

```python
ensembl_paralogs = pd.read_feather(
    "../data/intermediate/paralog_mutation_splicing_correlations.feather"
)

damaging_pvals = ensembl_paralogs["damaging_pval"].dropna()
hotspot_pvals = ensembl_paralogs["hotspot_pval"].dropna()

ensembl_paralogs["damaging_qval"] = pd.Series(
    multipletests(
        damaging_pvals,
        alpha=0.01,
        method="fdr_bh",
    )[1],
    index=damaging_pvals.index,
)

ensembl_paralogs["hotspot_qval"] = pd.Series(
    multipletests(
        hotspot_pvals,
        alpha=0.01,
        method="fdr_bh",
    )[1],
    index=hotspot_pvals.index,
)
```

```python
ensembl_paralogs["label"] = (
    ensembl_paralogs["external_gene_name"].fillna("")
    + " • "
    + ensembl_paralogs["hsapiens_paralog_associated_gene_name"].fillna("")
)
```

```python
ensembl_paralogs.columns
```

```python
ensembl_paralogs[ensembl_paralogs["damaging_qval"]<1e-10][["label","ccle_paralog_exons"]]
```

```python
sns.stripplot(
    depmap_damaging["KRT3"],
    ccle_exonusage["KRT76_3p_chr12_53165003_53164783_ENSG00000185069.2"],
)
```

```python
sns.stripplot(
    depmap_hotspot["GOLGA8R"],
    ccle_exonusage["GOLGA8S_5p_chr15_23609490_23609590_ENSG00000261739.2"],
)
```

```python
sns.boxplot(
    depmap_hotspot["GOLGA8R"],
    ccle_exonusage["GOLGA8S_3p_chr15_23609490_23609590_ENSG00000261739.2"],
    notch=True
)
```

```python
sns.boxplot(
    depmap_hotspot["RPL22"],
    ccle_exonusage["RPL22L1_5p_chr3_170585990_170585802_ENSG00000163584.13"],
    notch=True
)
```

```python
many.visuals.dense_plot(
    ensembl_paralogs["hotspot_biserial"],
    -np.log10(ensembl_paralogs["hotspot_qval"]),
    text_adjust=False,
    labels_mask=ensembl_paralogs["hotspot_qval"] < 1e-14,
    labels=ensembl_paralogs["label"],
    colormap=None,
)

plt.xlabel("Rank-biserial correlation")
plt.ylabel("-log10(Q value)")

plt.savefig(
    "../plots/hotspot_splicing_associations.pdf", dpi=512, bbox_inches="tight", transparent=True
)
```

```python
many.visuals.dense_plot(
    ensembl_paralogs["damaging_biserial"],
    -np.log10(ensembl_paralogs["damaging_qval"]),
    text_adjust=False,
    labels_mask=ensembl_paralogs["damaging_qval"] < 1e-15,
    labels=ensembl_paralogs["label"],
    colormap=None,
)

plt.xlabel("Rank-biserial correlation")
plt.ylabel("-log10(Q value)")

plt.savefig(
    "../plots/damaging_splicing_associations.pdf", dpi=512, bbox_inches="tight", transparent=True
)
```

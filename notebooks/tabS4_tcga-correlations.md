---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.2
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

```python
import numpy as np
import pandas as pd

from statsmodels.stats.multitest import multipletests

import gc

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

import cancer_data
import many

import config
config.config_visuals()

import stackprinter
stackprinter.set_excepthook(style='lightbg')
```

# Load data


## Merged TCGA info

```python
merged_tcga_info = pd.read_csv(
    "../data/supplementary/S2_merged-tcga-info.txt", sep="\t", index_col=0
)

tcga_msi = merged_tcga_info[merged_tcga_info["MSI"]==True]
```

## TCGA gene expression

```python
tcga_genex = cancer_data.load("tcga_normalized_gene_expression")

normal_genex = tcga_genex[tcga_genex.index.map(lambda x: x[-2:] == "11")]
tcga_genex = tcga_genex[tcga_genex.index.map(lambda x: x[-2:] != "11")]
```

## TCGA splicing

```python
def preprocess_splicing(df):
    
    df.index = df.index.map(lambda x: x[:15])
    df = df[~df.index.duplicated(keep="first")]
    
    # remove normals
    df = df[df.index.map(lambda x: x[-2:] != "11")]
    
    return df
```

```python
many.stats.mat_mwu_naive(
    tcga_se.iloc[:, :1000],
    merged_tcga_info[["RPL22_k15fs_mutation"]].dropna(),
    **corr_kwargs,
    pbar=True
)
```

```python
corr_kwargs = {"melt": True, "effect": "rank_biserial","pbar":True}

many.stats.mat_mwu_naive(
    tcga_a5ss,
    merged_tcga_info[["RPL22_k15fs_mutation"]],
    **corr_kwargs,
)
```

```python
tcga_splicing_sets = [
    "tcga_se",
    "tcga_a3ss",
    "tcga_a5ss",
    "tcga_ir",
    "tcga_mx_1",
    "tcga_mx_2",
    "tcga_mx_3",
    "tcga_mx_4",
    "tcga_mx_5",
]

mwu_kwargs = {"melt": True, "effect": "rank_biserial", "pbar": True}
corr_kwargs = {"melt": True, "method": "spearman", "pbar": True}


def splicing_vs_mutations(mut_series):

    merged = []

    for splicing_set in tcga_splicing_sets:

        print(f"Loading {splicing_set}")

        splicing = cancer_data.load(splicing_set)
        splicing = preprocess_splicing(splicing)

        correlations = many.stats.mat_mwu_naive(splicing, mut_series, **mwu_kwargs)
        merged.append(correlations)

        del splicing
        gc.collect()

    merged = pd.concat(merged)

    return merged


def recompute_qval(df):

    df["qval"] = multipletests(
        10 ** (-df["pval"]),
        alpha=0.01,
        method="fdr_bh",
    )[1]
    
    df["qval"] = -np.log10(df["qval"])

    return df


def corr_splicing(series):

    merged = []

    for splicing_set in tcga_splicing_sets:

        print(f"Loading {splicing_set}")

        splicing = cancer_data.load(splicing_set)
        splicing = preprocess_splicing(splicing)

        correlations = many.stats.mat_corr_naive(splicing, series, **corr_kwargs)
        merged.append(correlations)

        del splicing
        gc.collect()

    merged = pd.concat(merged)
    merged = merged[merged["n"] >= 100]
    merged = recompute_qval(merged)
    merged = merged.sort_values(by="qval", ascending=False)
    
    merged = merged.reset_index().rename(
        {
            "a_col": "first_splicing_event",
            "b_col": "second_splicing_event",
            "pval": "-log10(P value)",
            "qval": "-log10(Q value)",
        },
        axis=1
    )

    return merged
```

```python
tcga_mdm4_cosplicing = corr_splicing(merged_tcga_info["MDM4_exon_6_inclusion"])
tcga_rpl22l1_cosplicing = corr_splicing(merged_tcga_info["RPL22L1_exon_3A_inclusion"])
```

```python
tcga_mdm4_cosplicing.to_csv("../data/supplementary/S4-a_tcga-mdm4-cosplicing.txt", sep="\t")
tcga_rpl22l1_cosplicing.to_csv("../data/supplementary/S4-b_tcga-rpl22l1-cosplicing.txt", sep="\t")
```

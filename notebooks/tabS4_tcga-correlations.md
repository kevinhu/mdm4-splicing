---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.6.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

```python
import numpy as np
import pandas as pd

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

# subset
rpl22_cn_normal = np.abs(merged_tcga_info["RPL22_copy_number"])<0.25
rpl22_mut_wt = merged_tcga_info["RPL22_any_mutation"] == False

rpl22_wt_subset = merged_tcga_info[rpl22_cn_normal&rpl22_mut_wt]
```

## TCGA gene expression

```python
tcga_genex = cancer_data.load("tcga_normalized_gene_expression")

normal_genex = tcga_genex[tcga_genex.index.map(lambda x: x[-2:] == "11")]
tcga_genex = tcga_genex[tcga_genex.index.map(lambda x: x[-2:] != "11")]
```

## TCGA splicing

```python
tcga_se = cancer_data.load("tcga_se")
tcga_a3ss = cancer_data.load("tcga_a3ss")
tcga_a5ss = cancer_data.load("tcga_a5ss")
tcga_ir = cancer_data.load("tcga_ir")

def preprocess_splicing(df):
    
    df.index = df.index.map(lambda x: x[:15])
    df = df[~df.index.duplicated(keep="first")]
    
    # remove normals
    df = df[df.index.map(lambda x: x[-2:] != "11")]
    
    return df

tcga_se = preprocess_splicing(tcga_se)
tcga_a3ss = preprocess_splicing(tcga_a3ss)
tcga_a5ss = preprocess_splicing(tcga_a5ss)
tcga_ir = preprocess_splicing(tcga_ir)
```

# Compute correlations


## Splicing vs gene expression and splicing

```python
compare_sets = [tcga_genex, tcga_se, tcga_a3ss, tcga_a5ss, tcga_ir]

corr_kwargs = {"melt": True, "method": "spearman", "pbar": True}

(
    rpl22l1_3a_genex_overall_corrs,
    rpl22l1_3a_se_overall_corrs,
    rpl22l1_3a_a3ss_overall_corrs,
    rpl22l1_3a_a5ss_overall_corrs,
    rpl22l1_3a_ir_overall_corrs,
) = [
    many.stats.mat_corr_naive(
        merged_tcga_info["RPL22L1_exon_3A_inclusion"], df, **corr_kwargs
    )
    for df in compare_sets
]

(
    rpl22l1_3a_genex_wt_corrs,
    rpl22l1_3a_se_wt_corrs,
    rpl22l1_3a_a3ss_wt_corrs,
    rpl22l1_3a_a5ss_wt_corrs,
    rpl22l1_3a_ir_wt_corrs,
) = [
    many.stats.mat_corr_naive(
        rpl22_wt_subset["RPL22L1_exon_3A_inclusion"], df, **corr_kwargs
    )
    for df in compare_sets
]

(
    mdm4_6_genex_overall_corrs,
    mdm4_6_se_overall_corrs,
    mdm4_6_a3ss_overall_corrs,
    mdm4_6_a5ss_overall_corrs,
    mdm4_6_ir_overall_corrs,
) = [
    many.stats.mat_corr_naive(
        merged_tcga_info["MDM4_exon_6_inclusion"], df, **corr_kwargs
    )
    for df in compare_sets
]

(
    mdm4_6_genex_wt_corrs,
    mdm4_6_se_wt_corrs,
    mdm4_6_a3ss_wt_corrs,
    mdm4_6_a5ss_wt_corrs,
    mdm4_6_ir_wt_corrs,
) = [
    many.stats.mat_corr_naive(
        rpl22_wt_subset["MDM4_exon_6_inclusion"], df, **corr_kwargs
    )
    for df in compare_sets
]

(
    ubap2l_29_genex_overall_corrs,
    ubap2l_29_se_overall_corrs,
    ubap2l_29_a3ss_overall_corrs,
    ubap2l_29_a5ss_overall_corrs,
    ubap2l_29_ir_overall_corrs,
) = [
    many.stats.mat_corr_naive(
        merged_tcga_info["UBAP2L_exon_29_inclusion"], df, **corr_kwargs
    )
    for df in compare_sets
]

(
    ubap2l_29_genex_wt_corrs,
    ubap2l_29_se_wt_corrs,
    ubap2l_29_a3ss_wt_corrs,
    ubap2l_29_a5ss_wt_corrs,
    ubap2l_29_ir_wt_corrs,
) = [
    many.stats.mat_corr_naive(
        rpl22_wt_subset["UBAP2L_exon_29_inclusion"], df, **corr_kwargs
    )
    for df in compare_sets
]
```

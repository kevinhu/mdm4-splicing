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
import numpy as np
import pandas as pd

import cancer_data
import many

import config
config.config_visuals()

import stackprinter
stackprinter.set_excepthook(style='lightbg')
```

# Load data


## Merged CCLE info

```python
merged_ccle_info = pd.read_csv(
    "../data/supplementary/S1_merged-ccle-info.txt", sep="\t", index_col=0
)

merged_ccle_info["MDM4_6_mean"] = (
    merged_ccle_info["MDM4_3p_chr1_204506558_204506625_ENSG00000198625.8_exonusage"]
    + merged_ccle_info["MDM4_5p_chr1_204506558_204506625_ENSG00000198625.8_exonusage"]
)/2

# subset
rpl22_cn_normal = np.abs(1 - merged_ccle_info["RPL22_copynumber"]) < 0.125
rpl22_mut_wt = merged_ccle_info["RPL22_mutation_classification_collapsed"] == "WT"

rpl22_wt_subset = merged_ccle_info[rpl22_cn_normal & rpl22_mut_wt]
```

## CCLE gene expression and splicing

```python
ccle_genex = cancer_data.load("ccle_gene_tpm")
ccle_exonusage = cancer_data.load("ccle_exonusage")
```

# Compute correlations


## Splicing vs gene expression and splicing

```python
compare_sets = [ccle_genex, ccle_exonusage]

corr_kwargs = {"melt": True, "method": "spearman", "pbar": True}

(rpl22l1_3a_genex_overall_corrs, rpl22l1_3a_exonusage_overall_corrs,) = [
    many.stats.mat_corr_naive(
        merged_ccle_info[
            "RPL22L1_5p_chr3_170585990_170585802_ENSG00000163584.13_exonusage"
        ],
        df,
        **corr_kwargs
    )
    for df in compare_sets
]

(mdm4_6_genex_overall_corrs, mdm4_6_exonusage_overall_corrs,) = [
    many.stats.mat_corr_naive(merged_ccle_info["MDM4_6_mean"], df, **corr_kwargs)
    for df in compare_sets
]
```

```python
csv_kwargs = {"sep": "\t"}

outputs = [
    [rpl22l1_3a_genex_overall_corrs, "S3-a_rpl22l1-3a-genex-overall-corrs"],
    [rpl22l1_3a_exonusage_overall_corrs, "S3-b_rpl22l1-3a-exonusage-overall-corrs"],
]

for table, stem in outputs:

    table.to_csv(f"../data/supplementary/{stem}.txt", **csv_kwargs)
```

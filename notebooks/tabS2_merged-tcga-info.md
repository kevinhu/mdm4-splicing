---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.5.0
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

import sys
import os

sys.path.append(os.path.relpath("../../huygens"))
sys.path.append(os.path.relpath("../../galileo"))

import galileo as gal
import huygens as huy
```

```python
tcga_genex = pd.read_hdf("../../data/processed/TCGA/TCGA_genex_norm.h5",key="tcga_genex")

tcga_splicing = pd.read_hdf("../../data/processed/TCGA/merged.h5",key="tcga_splicing")
tcga_splicing.index = tcga_splicing.index.map(lambda x: x[:15])
tcga_splicing = tcga_splicing[~tcga_splicing.index.duplicated(keep="first")]

tcga_cn = pd.read_hdf("../../data/processed/TCGA/tcga_cn_whitelisted.hdf", key="tcga_cn")

tcga_mut_mat = pd.read_hdf("../../data/processed/TCGA/tcga_mut_mat.hdf", key="tcga_mut_mat")
tcga_msi = pd.read_hdf("../../data/processed/tcga/tcga_msi.h5",key="tcga_msi")
```

```python
rpl22_tcga = pd.read_csv("../data/raw/rpl22.tcga.data.csv")

rpl22_tcga = rpl22_tcga.dropna(subset=["sampleid"])
rpl22_tcga = rpl22_tcga.set_index("sampleid")
rpl22_tcga.index = rpl22_tcga.index.map(lambda x: x[:15])

rpl22_mut = rpl22_tcga["rpl22mut.mc3.k15"].dropna()

tcga_mut_mat["RPL22_chr1_6257785_6257785_T_-"] = rpl22_mut
```

```python
tcga_sample_info = pd.read_hdf("../../data/processed/TCGA/tcga_sample_info.hdf",key="tcga_sample_info")
```

# Aggregate attributes


## Tumor sample info

```python
select_sample_info = tcga_sample_info[["sample_type","_primary_disease","abbreviated_disease"]]

select_sample_info.columns = ["Sample_type","Primary_disease","Abbreviated_disease"]
```

## Mutations

```python
select_mutations = rpl22_tcga[["TP53mut","rpl22mut.mc3.k15"]]
select_mutations.columns = ["TP53_mutation_type","RPL22_k15fs_mutation"]
```

## MSI

```python
select_msi = tcga_msi[["MANTIS Score"]].copy()
select_msi["MSI"] = select_msi["MANTIS Score"] > 0.4

select_msi.columns = ["MANTIS_score","MSI"]
```

## Exonusage

```python
select_exons = [
    "MDM4_ENSG00000198625_ENSG00000198625.8_ES_1_204501318:204501374:204506557:204506625:204507336:204507436_204506557:204506625",
    "RPL22L1_ENSG00000163584_ENSG00000163584.13_A3_3_170586086:170586168:170585801:170585923:170585801:170585990_170585923:170585990",
    "UBAP2L_ENSG00000143569_ENSG00000143569.14_ES_1_154241382:154241430:154241837:154241888:154242675:154243040_154241837:154241888"
]

select_exonusage = tcga_splicing[select_exons]
select_exonusage.columns = ["MDM4_exon_6_inclusion","RPL22L1_exon_3A_inclusion","UBAP2L_exon_9_inclusion"]
```

## Gene expression

```python
select_genex_genes = [
    "MDM2_10743",
    "MDM4_10744",
    "RPL22_15208",
    "RPL22L1_15209"
]

select_genex = tcga_genex[select_genex_genes]
select_genex.columns = ["MDM2_mRNA","MDM4_mRNA","RPL22_mRNA","RPL22L1_mRNA"]
```

# Copy number

```python
select_copynumber_genes = [
    "MDM2_chr12_69201952_69244466",
    "MDM4_chr1_204485507_204527248",
    "RPL22_chr1_6241329_6260902",
    "RPL22L1_chr3_170582664_170588272"
]

select_copynumber = tcga_cn[select_copynumber_genes]
select_copynumber.columns = [
    "MDM2_copy_number",
    "MDM4_copy_number",
    "RPL22_copy_number",
    "RPL22L1_copy_number"
]
```

```python
merged_tcga_info = pd.concat([
    select_sample_info,
    select_mutations,
    select_msi,
    select_exonusage,
    select_genex,
    select_copynumber,
], join="outer", axis=1, sort=True)
```

```python
merged_tcga_info.to_csv("../data/supplementary/S2_merged-tcga-info.txt",sep="\t")
```

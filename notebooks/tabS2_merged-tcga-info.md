---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.9.1
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
```

```python
def process_tcga_splicing(tcga_splicing):
    
    # keep first four identifier fields
    tcga_splicing.index = tcga_splicing.index.map(lambda x: "-".join(x.split("-")[:4]))
    # remove last letter
    tcga_splicing.index = tcga_splicing.index.map(lambda x: x[:-1])
    tcga_splicing = tcga_splicing[tcga_splicing.index.map(lambda x: x[-2:] != "11")]

    tcga_splicing = tcga_splicing.loc[~tcga_splicing.index.duplicated(keep="first")]
    
    return tcga_splicing
```

```python
tcga_se = cancer_data.load("tcga_se")
tcga_a3ss = cancer_data.load("tcga_a3ss")

tcga_genex = cancer_data.load("tcga_normalized_gene_expression")
tcga_cn_continuous = cancer_data.load("tcga_cn_continuous")
tcga_cn_thresholded = cancer_data.load("tcga_cn_thresholded")
tcga_mutations = cancer_data.load("tcga_mutations")
tcga_msi = cancer_data.load("tcga_msi")

tcga_annotations = cancer_data.load("tcga_annotations")

tcga_se = process_tcga_splicing(tcga_se)
tcga_a3ss = process_tcga_splicing(tcga_a3ss)
```

```python
rpl22_tcga = pd.read_csv("../data/raw/rpl22.tcga.data.csv")

rpl22_tcga = rpl22_tcga.dropna(subset=["sampleid"])
rpl22_tcga = rpl22_tcga.set_index("sampleid")
rpl22_tcga.index = rpl22_tcga.index.map(lambda x: x[:15])

rpl22_mut = rpl22_tcga["rpl22mut.mc3.k15"].dropna()
```

```python
cbp_alterations = pd.read_csv("../data/external/cbioportal/tp53_rpl22.tsv",sep="\t")
cbp_alterations = cbp_alterations.set_index("Sample ID")
```

# Aggregate attributes


## Tumor sample info

```python
select_sample_info = tcga_annotations[
    ["sample_type", "_primary_disease", "abbreviated_disease"]
]

select_sample_info.columns = ["Sample_type", "Primary_disease", "Abbreviated_disease"]
```

## Mutations

```python
rpl22_tcga["TP53mut"] = rpl22_tcga["TP53mut"].fillna("WT")
select_mutations = rpl22_tcga[["TP53mut", "rpl22mut.mc3.k15", "rpl22mut.mc3.all"]]
select_mutations.columns = [
    "TP53_mutation_type",
    "RPL22_k15fs_mutation",
    "RPL22_any_mutation",
]
```

## cbioportal alterations

```python
select_cbp = cbp_alterations[
    ["RPL22: MUT", "RPL22: HOMDEL", "TP53: MUT", "TP53: HOMDEL"]
].replace({"not profiled":np.nan})
select_cbp.columns = [
    "RPL22_mutation_cbioportal",
    "RPL22_homdel_cbioportal",
    "TP53_mutation_cbioportal",
    "TP53_homdel_cbioportal",
]
```

## MSI

```python
select_msi = tcga_msi[["MANTIS Score"]].copy()
select_msi["MSI"] = select_msi["MANTIS Score"] > 0.4

select_msi.columns = ["MANTIS_score", "MSI"]
```

## Exonusage

```python
select_se = [
    "ENSG00000198625.8_ES_1_204501318:204501374:204506557:204506625:204507336:204507436_204506557:204506625",
    "ENSG00000143569.14_ES_1_154241382:154241430:154241837:154241888:154242675:154243040_154241837:154241888",
]

select_a3ss = [
    "ENSG00000163584.13_A3_3_170586086:170586168:170585801:170585923:170585801:170585990_170585923:170585990",
]

select_exonusage = pd.concat([tcga_se[select_se],tcga_a3ss[select_a3ss]],axis=1)
select_exonusage.columns = [
    "MDM4_exon_6_inclusion",
    "UBAP2L_exon_29_inclusion",
    "RPL22L1_exon_3A_inclusion",
]

select_exonusage["RPL22L1_exon_3A_inclusion"] = 1-select_exonusage["RPL22L1_exon_3A_inclusion"]
```

## Gene expression

```python
select_genex_genes = ["MDM2_10743", "MDM4_10744", "RPL22_15208", "RPL22L1_15209"]

select_genex = tcga_genex[select_genex_genes]
select_genex.columns = ["MDM2_mRNA", "MDM4_mRNA", "RPL22_mRNA", "RPL22L1_mRNA"]
```

# Copy number


## Continuous

```python
select_copynumber_genes = [
    "TP53",
    "MDM2",
    "MDM4",
    "RPL22",
    "RPL22L1",
]

select_copynumber = tcga_cn_continuous[select_copynumber_genes]
select_copynumber.columns = [
    "TP53_copy_number",
    "MDM2_copy_number",
    "MDM4_copy_number",
    "RPL22_copy_number",
    "RPL22L1_copy_number",
]
```

## Thresholded

```python
select_copynumber_thresholded_genes = ["TP53", "MDM2", "MDM4", "RPL22", "RPL22L1"]

select_copynumber_thresholded = tcga_cn_thresholded[select_copynumber_thresholded_genes]

select_copynumber_thresholded.columns = [
    "TP53_copy_number_thresholded",
    "MDM2_copy_number_thresholded",
    "MDM4_copy_number_thresholded",
    "RPL22_copy_number_thresholded",
    "RPL22L1_copy_number_thresholded",
]
```

# Merge

```python
merged_tcga_info = pd.concat(
    [
        select_sample_info,
        select_mutations,
        select_msi,
        select_exonusage,
        select_genex,
        select_copynumber,
        select_copynumber_thresholded,
        select_cbp,
    ],
    join="outer",
    axis=1,
    sort=True,
)
```

```python
merged_tcga_info.to_csv("../data/supplementary/S2_merged-tcga-info.txt", sep="\t")
```

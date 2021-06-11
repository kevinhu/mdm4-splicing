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

# Overview


# Aggregate attributes


## Cell line metadata

```python
# cell line characteristics
cell_line_info = cancer_data.load("depmap_annotations")
cell_line_info["COSMIC_ID"] = (
    cell_line_info["COSMICID"].astype(float).fillna(0).astype(int).astype(str).replace("0", "")
)

cell_line_info["CCLE_name"] = cell_line_info["CCLE_Name"]
cell_line_info["Sanger_ID"] = cell_line_info["Sanger_Model_ID"]
cell_line_info["disease"] = cell_line_info["primary_disease"]

select_cell_line_info = cell_line_info[
    ["CCLE_name", "COSMIC_ID", "Sanger_ID", "disease"]
]

select_cell_line_info.columns = [
    "CCLE_name",
    "COSMIC_ID",
    "Sanger_ID",
    "Primary_disease",
]
```

## Mutation calls

```python
mutation_calls = cancer_data.load("depmap_mutations")
mutation_calls["Change"] = "chr" + mutation_calls["Chromosome"].astype(str)
mutation_calls["Change"] = (
    mutation_calls["Change"] + "_" + mutation_calls["Start_position"].astype(str)
)
mutation_calls["Change"] = (
    mutation_calls["Change"] + "-" + mutation_calls["End_position"].astype(str)
)
mutation_calls["Change"] = (
    mutation_calls["Change"] + "_" + mutation_calls["Reference_Allele"].astype(str)
)
mutation_calls["Change"] = (
    mutation_calls["Change"] + ">" + mutation_calls["Tumor_Seq_Allele1"].astype(str)
)

# list of all cell lines with mutation data
all_mut_lines = sorted(list(set(mutation_calls["DepMap_ID"])))


def collapse_mutation(mutation_classes):
    if "damaging" in mutation_classes:
        return "damaging"
    if "other non-conserving" in mutation_classes:
        return "non-conserving"
    if "silent" in mutation_classes:
        return "silent"
    if "WT" in mutation_classes:
        return "wild-type"

    return np.nan


def get_mutations(gene):

    # fetch mutations for gene
    matched_muts = mutation_calls[mutation_calls["Hugo_Symbol"] == gene]

    # group by the cell line and convert to a list of mutations
    mutation_ids = matched_muts.groupby("DepMap_ID")["Change"].apply(list)
    mutation_ids = mutation_ids.apply(lambda x: ",".join([str(y) for y in x]))
    mutation_ids = mutation_ids.rename(gene + "_mutation")
    mutation_ids[mutation_ids.isna()] = ""

    # group by the type of mutation
    mutation_classes = matched_muts.groupby("DepMap_ID")["Variant_annotation"]

    # list the mutation classes
    mutation_classes = mutation_classes.apply(list)
    mutation_classes = mutation_classes.apply(lambda x: [str(y) for y in x])
    mutation_classes = mutation_classes.apply(lambda x: ",".join(x))
    mutation_classes = mutation_classes.rename(gene + "_mutation_classification")

    mutation_classes_collapsed = mutation_classes.apply(
        lambda x: collapse_mutation(x)
    ).rename(gene + "_mutation_classification_collapsed")

    # return mutation IDs and classes
    return pd.concat(
        [mutation_ids, mutation_classes, mutation_classes_collapsed], axis=1
    )
```

```python
mut_genes = ["TP53", "RPL22"]

select_muts = [get_mutations(gene) for gene in mut_genes]

merged_mutations = pd.concat(select_muts, axis=1, sort=False, join="outer")

merged_mutations, _ = merged_mutations.align(
    pd.Series(index=all_mut_lines), join="outer", axis=0
)

merged_mutations = merged_mutations.fillna("WT")
```

## MSI

```python
is_msi = cancer_data.load("ccle_msi")

select_msi = is_msi["MSI"]
```

## Splicing

```python
exonusage = cancer_data.load("ccle_exonusage")

select_exons = [
    "UBAP2L_5p_chr1_154242676_154243329_ENSG00000143569.14",
    "RPL22L1_5p_chr3_170585990_170585802_ENSG00000163584.13",
    "MDM4_3p_chr1_204506558_204506625_ENSG00000198625.8",
    "MDM4_5p_chr1_204506558_204506625_ENSG00000198625.8",
]

select_exonusage = exonusage[select_exons]
select_exonusage.columns = [x + "_exonusage" for x in select_exonusage.columns]

select_exonusage[
    "MDM4_mean_chr1_204506558_204506625_ENSG00000198625.8_exonusage"
] = exonusage[
    [
        "MDM4_3p_chr1_204506558_204506625_ENSG00000198625.8",
        "MDM4_5p_chr1_204506558_204506625_ENSG00000198625.8",
    ]
].mean(
    axis=1
)
```

## mRNA expression

```python
ccle_genex = cancer_data.load("ccle_gene_tpm")

select_genex_genes = [
    "MDM2_ENSG00000135679.17",
    "MDM4_ENSG00000198625.8",
    "RPL22_ENSG00000116251.5",
    "RPL22L1_ENSG00000163584.13",
]

select_genex = ccle_genex[select_genex_genes]

select_genex.columns = [x + "_mRNA" for x in select_genex.columns]
```

## Proteomics

```python
ms_prot = cancer_data.load("ccle_proteomics")
rppa = cancer_data.load("ccle_rppa")

select_ms_proteins = [
    "TP53_P04637",
    "MDM2_Q00987-11",
    "MDM4_O15151",
    "RPL22_P35268",
    "RPL22L1_Q6P5R6",
]

select_ms_prot = ms_prot[select_ms_proteins]
select_ms_prot.columns = [
    "TP53_MS_protein",
    "MDM2_MS_protein",
    "MDM4_MS_protein",
    "RPL22_MS_protein",
    "RPL22L1_MS_protein",
]

select_rppa_proteins = ["MDM2_MDM2_pS166", "MDM4_MDMX_MDM4(BetIHC-00108)_Caution"]

select_rppa = rppa[select_rppa_proteins]

select_rppa.columns = ["MDM2_RPPA_protein", "MDM4_RPPA_protein"]
```

## Copy number

```python
copynumber = cancer_data.load("depmap_copy_number")

select_copynumber_genes = ["TP53_7157", "RPL22_6146"]

select_copynumber = copynumber[select_copynumber_genes]
select_copynumber.columns = ["TP53_copynumber", "RPL22_copynumber"]
```

## Gene sensitivities

```python
avana = cancer_data.load("avana")
drive = cancer_data.load("drive")

select_avana_genes = [
    "TP53_7157",
    "MDM4_4194",
    "MDM2_4193",
    "RPL22_6146",
    "RPL22L1_200916",
]

select_drive_genes = [
    "TP53_7157",
    "MDM4_4194",
    "MDM2_4193",
    "RPL22_6146",
    "RPL22L1_200916",
]

select_avana = avana[select_avana_genes]
select_drive = drive[select_drive_genes]

select_avana.columns = [
    x.split("_")[0] + "_Avana_dependency" for x in select_avana.columns
]
select_drive.columns = [
    x.split("_")[0] + "_DRIVE_dependency" for x in select_drive.columns
]
```

## Drug sensitivities

```python
prism_primary_logfold = cancer_data.load("prism_primary_logfold")
prism_secondary_logfold = cancer_data.load("prism_secondary_logfold")

select_prism_primary = prism_primary_logfold[["nutlin-3_BRD-A12230535-001-06-7::2.5::HTS"]]

select_prism_primary.columns = ["nutlin-3_PRISM_primary_2.5"]
```

## Merge

```python
merged_ccle_info = pd.concat(
    [
        select_cell_line_info,
        merged_mutations,
        select_msi,
        select_exonusage,
        select_genex,
        select_ms_prot,
        select_rppa,
        select_copynumber,
        select_avana,
        select_drive,
        select_prism_primary
    ],
    join="outer",
    axis=1,
    sort=True,
)

merged_ccle_info = merged_ccle_info.dropna(how="all")
merged_ccle_info.index.name = "Achilles_ID"
```

```python
merged_ccle_info.to_csv("../data/supplementary/S1_merged-ccle-info.txt", sep="\t")
```

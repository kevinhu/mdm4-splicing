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
import numpy as np

import cancer_data

from scipy.stats import pearsonr, spearmanr, mannwhitneyu

from tqdm.notebook import tqdm_notebook
tqdm_notebook.pandas()
```

# Load data

```python
ensembl_paralogs = pd.read_csv("../data/raw/ensembl_paralogs.csv", index_col=0)
ensembl_entrez = pd.read_csv("../data/raw/ensembl_entrez.csv", index_col=0)
ensembl_entrez["entrezgene_id"] = ensembl_entrez["entrezgene_id"].astype("Int32")
```

```python
avana = cancer_data.load("avana")
depmap_damaging = cancer_data.load("depmap_damaging")
depmap_hotspot = cancer_data.load("depmap_hotspot")

ccle_genex = cancer_data.load("ccle_gene_tpm")
ccle_exonusage = cancer_data.load("ccle_exonusage")
depmap_cn = cancer_data.load("depmap_copy_number")
```

# Map IDs

```python
ensembl_to_entrez = dict(zip(ensembl_entrez["ensembl_gene_id"],ensembl_entrez["entrezgene_id"]))
ensembl_to_hgnc = dict(zip(ensembl_entrez["ensembl_gene_id"],ensembl_entrez["hgnc_symbol"]))
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

cn_genes = pd.DataFrame(pd.Series(depmap_cn.columns, name="id"))
cn_genes["entrez_id"] = (
    cn_genes["id"].apply(lambda x: x.split("_")[-1])
)

ensembl_to_genex = dict(zip(genex_genes["ensembl_id"],genex_genes["id"]))
ensembl_to_exonusage = exonusage_genes.groupby("ensembl_id")["id"].apply(list).to_dict()
entrez_to_cn = dict(zip(cn_genes["entrez_id"],cn_genes["id"]))
```

```python
ensembl_paralogs["entrez_id"] = ensembl_paralogs["ensembl_gene_id"].apply(
    ensembl_to_entrez.get
).astype(str)
ensembl_paralogs["paralog_entrez_id"] = ensembl_paralogs[
    "hsapiens_paralog_ensembl_gene"
].apply(ensembl_to_entrez.get)

ensembl_paralogs["hgnc_symbol"] = ensembl_paralogs["ensembl_gene_id"].apply(
    ensembl_to_hgnc.get
)
ensembl_paralogs["paralog_hgnc_symbol"] = ensembl_paralogs[
    "hsapiens_paralog_ensembl_gene"
].apply(ensembl_to_hgnc.get)
```

```python
ensembl_paralogs["ccle_gene"] = ensembl_paralogs["ensembl_gene_id"].apply(
    ensembl_to_genex.get
)
ensembl_paralogs["depmap_gene"] = ensembl_paralogs["entrez_id"].apply(
    entrez_to_cn.get
)
ensembl_paralogs["ccle_paralog_exons"] = ensembl_paralogs[
    "hsapiens_paralog_ensembl_gene"
].apply(ensembl_to_exonusage.get)
ensembl_paralogs = ensembl_paralogs.dropna(subset=["ccle_paralog_exons"])
ensembl_paralogs = ensembl_paralogs.explode("ccle_paralog_exons")
ensembl_paralogs.reset_index(inplace=True)
```

```python
def compute_paralog_cn_correlation(row):

    gene = row["depmap_gene"]
    exon = row["ccle_paralog_exons"]

    if gene is None:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)

    gene_cn = depmap_cn[gene].dropna()
    exon_exonusage = ccle_exonusage[exon].dropna()

    gene_cn, exon_exonusage = gene_cn.align(exon_exonusage, join="inner")

    intersection_size = gene_cn.size

    if intersection_size < 50:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)

    spearman_r, spearman_pval = spearmanr(gene_cn, exon_exonusage)
    pearson_r, pearson_pval = pearsonr(gene_cn, exon_exonusage)

    return (intersection_size, spearman_r, spearman_pval, pearson_r, pearson_pval)
```

```python
(
    ensembl_paralogs["cn_num_samples"],
    ensembl_paralogs["cn_spearman_r"],
    ensembl_paralogs["cn_spearman_pval"],
    ensembl_paralogs["cn_pearson_r"],
    ensembl_paralogs["cn_pearson_pval"],
) = zip(*ensembl_paralogs.progress_apply(compute_paralog_cn_correlation, axis=1))
```

```python
def compute_mwu(a, b):
    a, b = a.align(b, join="inner")

    pos = b[a == True]
    neg = b[a == False]

    if pos.size < 5 or neg.size < 5:
        return (np.nan, np.nan, np.nan, np.nan)

    U2, pval = mannwhitneyu(
        pos,
        neg,
        use_continuity=True,
        alternative="two-sided",
    )

    rank_biserial = 2 * U2 / (pos.size * neg.size) - 1

    return (pos.size, neg.size, rank_biserial, pval)


def compute_paralog_mut_correlation(row):

    gene = row["hgnc_symbol"]
    exon = row["ccle_paralog_exons"]

    if gene is None:
        return (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    exon_exonusage = ccle_exonusage[exon].dropna()

    if gene not in depmap_hotspot.columns:
        hotspot_pos_size = np.nan
        hotspot_neg_size = np.nan
        hotspot_biserial = np.nan
        hotspot_pval = np.nan
    else:
        gene_hotspot = depmap_hotspot[gene].dropna()
        (
            hotspot_pos_size,
            hotspot_neg_size,
            hotspot_biserial,
            hotspot_pval,
        ) = compute_mwu(gene_hotspot, exon_exonusage)

    if gene not in depmap_damaging.columns:
        damaging_pos_size = np.nan
        damaging_neg_size = np.nan
        damaging_biserial = np.nan
        damaging_pval = np.nan
    else:
        gene_damaging = depmap_damaging[gene].dropna()
        (
            damaging_pos_size,
            damaging_neg_size,
            damaging_biserial,
            damaging_pval,
        ) = compute_mwu(gene_damaging, exon_exonusage)

    return (
        hotspot_pos_size,
        hotspot_neg_size,
        hotspot_biserial,
        hotspot_pval,
        damaging_pos_size,
        damaging_neg_size,
        damaging_biserial,
        damaging_pval,
    )
```

```python
(
    ensembl_paralogs["hotspot_pos_size"],
    ensembl_paralogs["hotspot_neg_size"],
    ensembl_paralogs["hotspot_biserial"],
    ensembl_paralogs["hotspot_pval"],
    ensembl_paralogs["damaging_pos_size"],
    ensembl_paralogs["damaging_neg_size"],
    ensembl_paralogs["damaging_biserial"],
    ensembl_paralogs["damaging_pval"],
) = zip(*ensembl_paralogs.progress_apply(compute_paralog_mut_correlation, axis=1))
```

```python
ensembl_paralogs.to_feather("../data/intermediate/paralog_mutation_splicing_correlations.feather")
```

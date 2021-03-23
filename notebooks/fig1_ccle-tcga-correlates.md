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

from scipy.stats import kruskal
from scipy.stats import gaussian_kde

from statsmodels.stats.multitest import multipletests

import cancer_data 
import many

from adjustText import adjust_text

import config
config.config_visuals()
```

# Transcripts to genes

```python
t2g = pd.read_csv("../data/intermediate/sleuth_diff/ensembl_t2g.csv")
t2g["format_gene_id"] = t2g["hgnc_gene"].fillna("") + "_" + t2g["ens_gene"]

format_gene_map = dict(zip(t2g["ens_gene"], t2g["format_gene_id"]))
```

## RPL22 mutants

```python
rpl22_tcga = pd.read_csv("../data/raw/rpl22.tcga.data.csv")

rpl22_tcga = rpl22_tcga.dropna(subset=["sampleid"])
rpl22_tcga = rpl22_tcga.set_index("sampleid")
rpl22_tcga.index = rpl22_tcga.index.map(lambda x: x[:15])

rpl22_mut = rpl22_tcga["rpl22mut.mc3.k15"].dropna()
```

```python
rpl22_info = pd.concat(
    [rpl22_mut.rename("RPL22_k15"), tcga_cn_thresholded["RPL22"].rename("RPL22_cn")],
    axis=1,
    sort=True,
    join="inner",
)
```

# Boxplots

```python
def rpl22_status(row):
    if row["RPL22_k15"]:
        return "K15.fs"

    else:
        return "ΔCN=" + str(int(row["RPL22_cn"]))


rpl22_info["RPL22_status"] = rpl22_info.apply(rpl22_status, axis=1)

rpl22_order = ["K15.fs", "ΔCN=-2", "ΔCN=-1", "ΔCN=0", "ΔCN=1", "ΔCN=2"]

mut_hue = "#e23e57"
wt_hue = "#eaeaea"

rpl22_hues = dict(zip(rpl22_order, [mut_hue] + [wt_hue] * 5))
```

```python
def plot_rpl22(y, ylabel="y"):

    plt.figure(figsize=(4, 3))

    ax = plt.subplot(111)

    rpl22_info_join = rpl22_info.join(y.dropna(), how="inner")

    rpl22_counts = rpl22_info_join["RPL22_status"].value_counts()

    grouped_y = rpl22_info_join.groupby("RPL22_status")[y.name].apply(list)
    grouped_y = list(grouped_y)

    pval = kruskal(*grouped_y)[1]

    flierprops = dict(
        marker=".",
        markersize=8,
        markerfacecolor=(0, 0, 0, 0),
        markeredgecolor="black",
        markeredgewidth=0.5,
    )

    sns.boxplot(
        rpl22_info["RPL22_status"],
        y,
        notch=True,
        order=rpl22_order,
        bootstrap=1000,
        palette=rpl22_hues,
        flierprops=flierprops,
    )

    plt.xlabel("RPL22 status")
    plt.ylabel(ylabel)

    if pval == 0:
        plt.text(
            0.05, 1.025, "P < " + huy.as_si(10 ** (-320), 2), transform=ax.transAxes
        )

    else:
        plt.text(0.05, 1.025, "P = " + huy.as_si(pval, 2), transform=ax.transAxes)

    xticks = [x + "\n (" + str(int(rpl22_counts.loc[x])) + ")" for x in rpl22_order]

    ax.set_xticklabels(xticks)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["left"].set_position(("axes", -0.025))
    ax.spines["bottom"].set_position(("axes", -0.025))
```

```python
plot_rpl22(tcga_genex["MDM2_10743"], "MDM2 mRNA expression")
```

```python
plot_rpl22(tcga_genex["MDM4_10744"], "MDM4 mRNA expression")
```

```python
plot_rpl22(tcga_genex["RPL22L1_15209"], "RPL22L1 mRNA expression")
```

```python
plot_rpl22(
    1
    - tcga_splicing[
        "RPL22L1_ENSG00000163584_ENSG00000163584.13_A3_3_170586086:170586168:170585801:170585923:170585801:170585990_170585923:170585990"
    ],
    "RPL22L1 exon 3A inclusion",
)

plt.savefig(
    "../plots/rpl22l1-exon-3a_rpl22_tcga.pdf", bbox_inches="tight", transparent=True
)
```

```python
plot_rpl22(
    tcga_splicing[
        "MDM4_ENSG00000198625_ENSG00000198625.8_ES_1_204501318:204501374:204506557:204506625:204507336:204507436_204506557:204506625"
    ],
    "MDM4 exon 6 inclusion",
)

plt.savefig(
    "../plots/mdm4-exon-6_rpl22_tcga.pdf", bbox_inches="tight", transparent=True
)
```

```python
rpl22_mdm4 = pd.concat(
    [
        tcga_splicing[
            "RPL22L1_ENSG00000163584_ENSG00000163584.13_A3_3_170586086:170586168:170585801:170585923:170585801:170585990_170585923:170585990"
        ]
        .rename("RPL22L1_3A")
        .dropna(),
        tcga_cn_thresholded["RPL22"].dropna(),
    ],
    join="inner",
    axis=1,
)

rpl22_neutral = rpl22_mdm4[rpl22_mdm4["RPL22"] == 0]
```

```python
corrs = gal.mat_corrs_naive(rpl22_neutral["RPL22L1_3A"], tcga_splicing, pbar=True)
```

# Scatterplots

```python
def density_scatter(x, y):

    x = x.dropna()
    y = y.dropna()

    x, y = x.align(y, join="inner")

    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)

    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    plt.figure(figsize=(3, 3))
    ax = plt.subplot(111)

    ax.scatter(x, y, c=z, s=8, cmap="Blues", rasterized=True, vmin=min(z) - 1)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return ax
```

## RPL22L1 3A vs expression

```python
x = (
    1
    - tcga_splicing[
        "RPL22L1_ENSG00000163584_ENSG00000163584.13_A3_3_170586086:170586168:170585801:170585923:170585801:170585990_170585923:170585990"
    ]
)
y = tcga_genex["RPL22L1_15209"]

ax = density_scatter(x, y)

ax.set_xlabel("RPL22L1 exon 3A inclusion")
ax.set_ylabel("RPL22L1 mRNA expression")

plt.savefig(
    "../plots/rpl22l1_splicing_expression_density.pdf",
    dpi=2048,
    bbox_inches="tight",
    background="transparent",
)

plt.figure(figsize=(3, 3))

huy.regression(x, y)

plt.xlabel("RPL22L1 exon 3A inclusion")
plt.ylabel("RPL22L1 mRNA expression")

plt.savefig(
    "../plots/rpl22l1_splicing_expression_annotated.pdf",
    dpi=2048,
    bbox_inches="tight",
    background="transparent",
)
```

## RPL22L1 3A vs MDM4 6

```python
x = (
    1
    - tcga_splicing[
        "RPL22L1_ENSG00000163584_ENSG00000163584.13_A3_3_170586086:170586168:170585801:170585923:170585801:170585990_170585923:170585990"
    ]
)
y = tcga_splicing[
    "MDM4_ENSG00000198625_ENSG00000198625.8_ES_1_204501318:204501374:204506557:204506625:204507336:204507436_204506557:204506625"
]

ax = density_scatter(x, y)

ax.set_xlabel("RPL22L1 exon 3A inclusion")
ax.set_ylabel("MDM4 exon 6 inclusion")

plt.savefig(
    "../plots/rpl22l1_mdm4_co-splicing_density.pdf",
    dpi=2048,
    bbox_inches="tight",
    background="transparent",
)

plt.figure(figsize=(3, 3))

huy.regression(x, y)

plt.xlabel("RPL22L1 exon 3A inclusion")
plt.ylabel("MDM4 exon 6 inclusion")

plt.xlim(-0.05, 1.05)
plt.ylim(-0.05, 1.05)

plt.savefig(
    "../plots/rpl22l1_mdm4_co-splicing_annotated.pdf",
    dpi=2048,
    bbox_inches="tight",
    transparent=True,
)
```

# TCGA correlations

```python
MIN_SAMPLES = 100

mdm4_6_exonusage_overall_corrs = pd.read_csv(
    "../data/supplementary/S4-f_mdm4-6-exonusage-overall-corrs.txt",
    sep="\t",
    index_col=0,
)

mdm4_6_exonusage_overall_corrs["label"] = mdm4_6_exonusage_overall_corrs["b_col"].apply(
    lambda x: f"{format_gene_map[x.split('.')[0]].split('_')[0]}_{x.split('_')[-1]}"
)

mdm4_6_exonusage_overall_corrs = mdm4_6_exonusage_overall_corrs[
    mdm4_6_exonusage_overall_corrs["n"] >= MIN_SAMPLES
]
mdm4_6_exonusage_overall_corrs = mdm4_6_exonusage_overall_corrs[
    mdm4_6_exonusage_overall_corrs["spearman"] < 1
]
```

```python
plt.figure(figsize=(3, 4))

labels_mask = mdm4_6_exonusage_overall_corrs["qval"] > 70

labels_mask = labels_mask | (
    (mdm4_6_exonusage_overall_corrs["spearman"] < 0)
    & (mdm4_6_exonusage_overall_corrs["qval"] > 65)
)
labels_mask = labels_mask | (
    (mdm4_6_exonusage_overall_corrs["qval"] > 25)
    & (mdm4_6_exonusage_overall_corrs["spearman"] < -0.5)
)

many.visuals.dense_plot(
    mdm4_6_exonusage_overall_corrs["spearman"],
    mdm4_6_exonusage_overall_corrs["qval"],
    text_adjust=True,
    labels_mask=labels_mask,
    labels=mdm4_6_exonusage_overall_corrs["label"],
    colormap=None,
)
```

```python

```

```python
plt.scatter(mdm4_6_exonusage_overall_corrs["spearman"],mdm4_6_exonusage_overall_corrs["qval"])
```

```python
mdm4_6_exonusage_overall_corrs.
```

## Overall correlations

```python
def volcano(corrs):
    corrs = corrs.sort_values(by="pval")
    corrs["qval"] = multipletests(10 ** (-corrs["pval"]), alpha=0.01, method="fdr_bh")[
        1
    ]

    labels = pd.Series(corrs.index.map(lambda x: x.split("_")[0]), index=corrs.index)

    pos_ranks = corrs[corrs["corr"] > 0]["qval"].rank()
    neg_ranks = corrs[corrs["corr"] < 0]["qval"].rank()
    ranks = pd.concat([pos_ranks, neg_ranks])

    corrs["qval_rank"] = ranks

    pos_ranks = (-corrs[corrs["corr"] > 0]["corr"]).rank()
    neg_ranks = corrs[corrs["corr"] < 0]["corr"].rank()
    ranks = pd.concat([pos_ranks, neg_ranks])

    corrs["corr_rank"] = ranks

    ranks = pd.concat([pos_ranks, neg_ranks])

    huy.dense_plot(
        corrs["corr"],
        -np.log10(corrs["qval"]),
        labels_mask=(corrs["corr_rank"] < 6) | (corrs["qval_rank"] < 6),
        labels=labels,
        adjust=False,
        c="black",
    )

    plt.xlabel("Spearman correlation")
    plt.ylabel("-log10(q-value)")
```

```python
volcano(rpl22_cn_splicing)
plt.savefig(
    "../plots/RPL22_cn_vs_splicing.pdf", bbox_inches="tight", transparent=True, dpi=512
)
```

```python
volcano(mdm4_cosplicing)
plt.savefig(
    "../plots/MDM4_cosplicing.pdf", bbox_inches="tight", transparent=True, dpi=512
)
```

```python
volcano(rpl22l1_cosplicing)
plt.savefig(
    "../plots/RPL22L1_cosplicing.pdf", bbox_inches="tight", transparent=True, dpi=512
)
```

# RPL22 alterations by primary site

```python
rpl22_subtype_info = rpl22_info.join(
    tcga_sample_info["abbreviated_disease"], how="inner"
)

subtype_proportions = rpl22_subtype_info.groupby(["abbreviated_disease"])[
    "RPL22_status"
].value_counts()
subtype_proportions = subtype_proportions.unstack().fillna(0)

subtype_totals = subtype_proportions.sum(axis=1)

subtype_proportions = subtype_proportions.div(subtype_totals, axis=0) * 100
subtype_proportions["total"] = subtype_totals.astype(int)

altered_classes = [
    "ΔCN=-1",
    "ΔCN=-2",
    "K15.fs",
]

subtype_proportions["altered_proportion"] = subtype_proportions[altered_classes].sum(
    axis=1
)
subtype_proportions = subtype_proportions.sort_values(
    by="altered_proportion", ascending=False
)

subtype_proportions["display_disease"] = subtype_proportions.index
subtype_proportions["display_disease"] = (
    subtype_proportions["display_disease"]
    + " ("
    + subtype_proportions["total"].astype(str)
    + ")"
)

subtype_proportions = subtype_proportions[subtype_proportions["total"] >= 50]
```

```python
plt.figure(figsize=(7, 3))

ax = plt.subplot(111)

subtype_proportions.plot(
    x="display_disease",
    y=altered_classes,
    kind="bar",
    stacked=True,
    cmap=mpl.colors.ListedColormap(["#ebd5d5", "#fab57a", "#b83b5e"]),
    ax=ax,
    width=0.75,
)

plt.xlabel("Subtype disease")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_position(("axes", -0.025))

plt.xticks(rotation=45, ha="right")

plt.ylabel("% affected")

plt.savefig(
    "../plots/rpl22_subtype_distribution.pdf", bbox_inches="tight", transparent=True
)
```

# RPL22 boxplots

```python
rpl22_mdm4_merged = pd.concat(
    [
        rpl22_info,
        tcga_splicing[
            "MDM4_ENSG00000198625_ENSG00000198625.8_ES_1_204501318:204501374:204506557:204506625:204507336:204507436_204506557:204506625"
        ].rename("MDM4_exon_6"),
        tcga_sample_info,
    ],
    axis=1,
    sort=True,
)
```

```python
select_type = rpl22_mdm4_merged[rpl22_mdm4_merged["abbreviated_disease"] == "STAD"]
select_type = select_type.dropna(subset=["RPL22_k15", "MDM4_exon_6"], how="any")

huy.two_dists(select_type["MDM4_exon_6"], select_type["RPL22_k15"], summary_type="box")
```

# CCLE

```python
ccle_genex = pd.read_hdf(
    "../../data/processed/ccle/CCLE_RNAseq_rsem_genes_tpm_20180929.hdf",
    key="ccle_genex",
)
ccle_transcripts = pd.read_hdf(
    "../../data/processed/ccle/CCLE_RNAseq_rsem_transcripts_tpm_20180929.hdf",
    key="ccle_transcripts",
)
exonusage = pd.read_hdf(
    "../../data/processed/ccle/CCLE_RNAseq_ExonUsageRatio_20180929.hdf", key="exonusage"
)
ms_prot = pd.read_hdf("../../data/processed/ccle/ms_prot.h5", key="ms_prot")
rppa = pd.read_hdf("../../data/processed/ccle/CCLE_RPPA_20181003.hdf", key="rppa")

mdm4_exon_6 = exonusage[
    [
        "MDM4_3p_chr1_204506558_204506625_ENSG00000198625.8",
        "MDM4_5p_chr1_204506558_204506625_ENSG00000198625.8",
    ]
].mean(axis=1)

avana = pd.read_hdf("../../data/processed/depmap/avana.hdf", key="avana")
drive = pd.read_hdf(
    "../../data/processed/depmap/demeter2-drive_v12-gene-effect.hdf", key="drive"
)
```

```python
primary_logfold = pd.read_hdf(
    "../../data/processed/depmap/primary_logfold.h5", key="primary_logfold"
)
secondary_logfold = pd.read_hdf(
    "../../data/processed/depmap/secondary_logfold.h5", key="secondary_logfold"
)
```

```python
ubap2l_tx = ccle_transcripts[[x for x in ccle_transcripts.columns if "UBAP2L" in x]]
ubap2l_exons = exonusage[[x for x in exonusage.columns if "UBAP2L" in x]]
```

## RPPA correlations

```python
mdm4_rppa = pd.concat(
    [
        ccle_genex["MDM4_ENSG00000198625.8"].rename("MDM4_genex"),
        mdm4_exon_6.rename("MDM4_exon_6"),
        rppa["MDMX_MDM4(BetIHC-00108)_Caution"].rename("MDM4_protein"),
    ],
    axis=1,
    sort=True,
)
```

```python
fig, axes = plt.subplots(1, 3, figsize=(7, 7 / 3), sharey=True)

text_pos = (0.075, 0.925)

ax = axes[0]
huy.dense_regression(
    mdm4_rppa["MDM4_genex"], mdm4_rppa["MDM4_protein"], ax=ax, s=12, text_pos=text_pos
)
ax.set_ylabel("MDM4 protein (RPPA)")
ax.set_xlabel("MDM4 mRNA")

ax = axes[1]
huy.dense_regression(
    mdm4_rppa["MDM4_exon_6"], mdm4_rppa["MDM4_protein"], ax=ax, s=12, text_pos=text_pos
)
ax.set_xlabel("MDM4 exon 6 inclusion")

mult = np.log2(2 ** mdm4_rppa["MDM4_genex"] * mdm4_rppa["MDM4_exon_6"] + 1)

ax = axes[2]
huy.dense_regression(mult, mdm4_rppa["MDM4_protein"], ax=ax, s=12, text_pos=text_pos)

ax.set_xlabel("MDM4 mRNA × exon 6")

plt.savefig(
    "../plots/MDM4_RPPA_correlations.pdf",
    bbox_inches="tight",
    transparent=True,
    dpi=512,
)
```

## Mass-spec correlations

```python
mdm4_ms = pd.concat(
    [
        ccle_genex["MDM4_ENSG00000198625.8"].rename("MDM4_genex"),
        mdm4_exon_6.rename("MDM4_exon_6"),
        ms_prot["MDM4_HUMAN_O15151"].rename("MDM4_protein"),
    ],
    axis=1,
    sort=True,
)
```

```python
fig, axes = plt.subplots(1, 3, figsize=(7, 7 / 3), sharey=True)

text_pos = (0.075, 0.925)

ax = axes[0]
huy.regression(
    mdm4_ms["MDM4_genex"],
    mdm4_ms["MDM4_protein"],
    ax=ax,
    s=16,
    text_pos=text_pos,
    alpha=1,
    c="black",
)
ax.set_ylabel("MDM4 protein (MS)")
ax.set_xlabel("MDM4 mRNA")

ax = axes[1]
huy.regression(
    mdm4_ms["MDM4_exon_6"],
    mdm4_ms["MDM4_protein"],
    ax=ax,
    s=16,
    text_pos=text_pos,
    alpha=1,
    c="black",
)
ax.set_xlabel("MDM4 exon 6 inclusion")
ax.set_xlim(0, 1)
ax.set_xticks([0, 0.5, 1])

mult = np.log2(2 ** mdm4_ms["MDM4_genex"] * mdm4_ms["MDM4_exon_6"] + 1)

ax = axes[2]
huy.regression(
    mult, mdm4_ms["MDM4_protein"], ax=ax, s=16, text_pos=text_pos, alpha=1, c="black"
)

ax.set_xlabel("MDM4 mRNA × exon 6")

ax.set_ylim(ax.set_ylim()[0], ax.set_ylim()[1] * 1.75)

plt.savefig(
    "../plots/MDM4_MS_correlations.pdf", bbox_inches="tight", transparent=True, dpi=512
)
```

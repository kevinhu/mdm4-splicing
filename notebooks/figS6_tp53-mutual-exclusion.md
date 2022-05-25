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

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

import cancer_data
import many

import config

from scipy.stats import fisher_exact

config.config_visuals()
```

```python
def as_si(x: float, decimals: int) -> str:
    """
    Convert a number to scientific notation

    Parameters
    ----------
    x : float
        number to convert
    decimals: float
        number of decimal places

    Returns
    -------
    x_si : string
        x formatted in scientific notation
    """

    s = "{x:0.{ndp:d}e}".format(x=x, ndp=decimals)
    m, e = s.split("e")
    x_si = r"{m:s} × $10^{{{e:d}}}$".format(m=m, e=int(e))

    return x_si

def binary_contingency(a, b, ax=None, heatmap_kwargs={}):
    """
    Plot agreement between two binary variables, along with
    the odds ratio and Fisher's exact test p-value.

    Parameters
    ----------
    a : Boolean series
        boolean series of first variable
    y : Boolean series
        boolean series of second variable
    ax : MatPlotLib axis
        axis to plot in (will create new one if not provided)
    heatmap_kwargs : dictionary
        additional arguments to pass to seaborn.heatmap()

    Returns
    -------
    ax : MatPlotLib axis
        axis with plot data
    """

    # convert, align, and cast
    a, b = pd.Series(a), pd.Series(b)
    a, b = a.dropna(), b.dropna()
    a, b = a.astype(bool), b.astype(bool)

    a, b = a.align(b, join="inner")

    # store names before array conversion
    a_name = a.name
    b_name = b.name

    a, b = np.array(a), np.array(b)

    # compute contingency counts
    xx = np.sum(a & b)
    xy = np.sum(a & ~b)
    yx = np.sum(~a & b)
    yy = np.sum(~a & ~b)

    # create 2x2 contingency table
    contingency = pd.DataFrame(
        [[xx, xy], [yx, yy]],
        columns=["True", "False"],
        index=["True", "False"],
    )

    odds_ratio, p_val = fisher_exact([[xx, xy], [yx, yy]])

    print("Odds ratio:", odds_ratio)
    print("P-value:", p_val)

    # if axis not provided, make figure
    if ax is None:

        plt.figure(figsize=(4, 4))
        ax = plt.subplot(111)

    # plot the contingency table as heatmap
    g = sns.heatmap(
        contingency,
        fmt="d",
        annot=True,
        cbar=False,
        linewidths=2,
        ax=ax,
        **heatmap_kwargs,
    )

    plt.ylabel(a_name)
    plt.xlabel(b_name)

    g.xaxis.tick_top()
    g.xaxis.set_label_position("top")

    return ax, odds_ratio, p_val

```

```python
merged_ccle_info = pd.read_csv(
    "../data/supplementary/S1_merged-ccle-info.txt", sep="\t", index_col=0
)
merged_tcga_info = pd.read_csv(
    "../data/supplementary/S2_merged-tcga-info.txt", sep="\t", index_col=0
)
merged_tcga_info["TP53_mut"] = merged_tcga_info["TP53_mutation_cbioportal"].dropna() != "no alteration"
merged_tcga_info["TP53_del"] = merged_tcga_info["TP53_copy_number_thresholded"].dropna() <= -1
```

```python
binary_contingency(
        merged_tcga_info["TP53_del"],
        merged_tcga_info["RPL22_k15fs_mutation"].dropna(),
        heatmap_kwargs={"cmap":"Blues"},
    )
```

```python
def plot_subset_contingency(
    subset_type, first_col, second_col, first_col_labels, second_col_labels, title, ax
):

    subset = merged_tcga_info[merged_tcga_info["Sample_type"] == "Primary Tumor"].copy()
    subset = subset[subset["Primary_disease"] == subset_type]

    _, odds_ratio, pval = binary_contingency(
        subset[first_col].dropna(),
        subset[second_col].dropna(),
        ax=ax,
        heatmap_kwargs={"cmap": "Blues"},
    )

    ax.set_title(title)

    ax.set_xticklabels(second_col_labels)
    ax.set_yticklabels(first_col_labels)

    ax.set_xlabel(f"Odds={odds_ratio:.2f}, P={as_si(pval,1)}")
    ax.set_ylabel(None)

    ax.tick_params(
        axis="x",
        which="both",
        bottom=False,
        top=False,
    )

    ax.tick_params(
        axis="y",
        which="both",
        left=False,
        right=False,
    )


fig, axes = plt.subplots(2, 3, figsize=(9, 6))

plot_subset_contingency(
    "colon adenocarcinoma",
    "TP53_mut",
    "RPL22_k15fs_mutation",
    ["TP53 mut", "TP53 WT"],
    ["RPL22 k15fs", "RPL22 WT"],
    "COAD",
    axes[0][0],
)
plot_subset_contingency(
    "stomach adenocarcinoma",
    "TP53_mut",
    "RPL22_k15fs_mutation",
    ["TP53 mut", "TP53 WT"],
    ["RPL22 k15fs", "RPL22 WT"],
    "STAD",
    axes[0][1],
)
plot_subset_contingency(
    "uterine corpus endometrioid carcinoma",
    "TP53_mut",
    "RPL22_k15fs_mutation",
    ["TP53 mut", "TP53 WT"],
    ["RPL22 k15fs", "RPL22 WT"],
    "UCEC",
    axes[0][2],
)

plot_subset_contingency(
    "colon adenocarcinoma",
    "TP53_del",
    "RPL22_k15fs_mutation",
    ["TP53 del", "TP53 WT"],
    ["RPL22 k15fs", "RPL22 WT"],
    "COAD",
    axes[1][0],
)
plot_subset_contingency(
    "stomach adenocarcinoma",
    "TP53_del",
    "RPL22_k15fs_mutation",
    ["TP53 del", "TP53 WT"],
    ["RPL22 k15fs", "RPL22 WT"],
    "STAD",
    axes[1][1],
)
plot_subset_contingency(
    "uterine corpus endometrioid carcinoma",
    "TP53_del",
    "RPL22_k15fs_mutation",
    ["TP53 del", "TP53 WT"],
    ["RPL22 k15fs", "RPL22 WT"],
    "UCEC",
    axes[1][2],
)

plt.subplots_adjust(wspace=0.6, hspace=0.6)

plt.savefig("../plots/fig-s5_TP53-RPL22k15-MSI_contingency.pdf", bbox_inches="tight")
```

```python
def plot_subset(subset_type=None, ax=None, heatmap_kwargs={}):

    if ax is None:

        plt.figure(figsize=(12, 2))
        ax = plt.subplot(111)

    if subset_type is None:

        subset = merged_tcga_info[
            merged_tcga_info["Sample_type"] == "Primary Tumor"
        ].copy()

    else:

        subset = merged_tcga_info[
            (merged_tcga_info["Abbreviated_disease"] == subset_type)
            & (merged_tcga_info["Sample_type"] == "Primary Tumor")
        ].copy()

    # define TP53 altered
    subset["TP53_mut"] = subset["TP53_mutation_cbioportal"].dropna() != "no alteration"
    subset["TP53_del"] = subset["TP53_copy_number_thresholded"].dropna() <= -1
#     subset["RPL22_del"] = subset["RPL22_copy_number_thresholded"].dropna() <= -1
    
    # drop missing annotations
    subset = subset[["TP53_del", "TP53_mut", "RPL22_k15fs_mutation", "MSI"]].dropna()

    # sort before plotting
    subset = subset.sort_values(
        ["TP53_mut", "TP53_del", "RPL22_k15fs_mutation", "MSI"], ascending=False
    )

    subset["RPL22_k15fs_mutation"] *= 2
    subset["MSI"] *= 3

    sns.heatmap(
        subset[["TP53_mut", "TP53_del", "RPL22_k15fs_mutation", "MSI"]].astype(int).T,
        cbar=False,
        xticklabels=False,
        cmap=sns.color_palette(["whitesmoke", "black", "red", "red"]),
        ax=ax,
        **heatmap_kwargs,
    )

    # set tick labels
    ax.set_yticklabels(["TP53 mutation", "TP53 deletion", "RPL22 k15fs", "MSI"])
    ax.tick_params(axis="both", which="both", length=0)

    if subset_type:

        ax.set_title(f"{subset_type} (n={len(subset)})")


fig, axes = plt.subplots(3, 1, figsize=(12, 6))

heatmap_kwargs = {"linewidth": 0.25}

plot_subset("COAD", ax=axes[0], heatmap_kwargs=heatmap_kwargs)
plot_subset("STAD", ax=axes[1], heatmap_kwargs=heatmap_kwargs)
plot_subset(
    "UCEC", ax=axes[2], heatmap_kwargs=heatmap_kwargs
)

plt.subplots_adjust(hspace=0.5)

plt.savefig("../plots/fig-s5_TP53-RPL22k15.pdf", bbox_inches="tight")
```

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
import ujson
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib_venn import venn2

import config
config.config_visuals()
```

```python
rpl22_oe_rmats = pd.read_hdf("../data/processed/rmats_merge/rpl22_oe.h5", key="rmats")
rpl22l1_oe_rmats = pd.read_hdf(
    "../data/processed/rmats_merge/rpl22l1_oe.h5", key="rmats"
)
rpl22l1_kd1_rmats = pd.read_hdf(
    "../data/processed/rmats_merge/rpl22l1_kd1.h5", key="rmats"
)
rpl22l1_kd2_rmats = pd.read_hdf(
    "../data/processed/rmats_merge/rpl22l1_kd2.h5", key="rmats"
)
rpl22_a_ko1_rmats = pd.read_hdf(
    "../data/processed/rmats_merge/rpl22_a_ko1.h5", key="rmats"
)
rpl22_a_ko2_rmats = pd.read_hdf(
    "../data/processed/rmats_merge/rpl22_a_ko2.h5", key="rmats"
)
rpl22_b_ko1_rmats = pd.read_hdf(
    "../data/processed/rmats_merge/rpl22_b_ko1.h5", key="rmats"
)
rpl22_b_ko2_rmats = pd.read_hdf(
    "../data/processed/rmats_merge/rpl22_b_ko2.h5", key="rmats"
)

rpl22_oe_rmats = rpl22_oe_rmats.rename({"PValue": "pval", "FDR": "qval"}, axis=1)
rpl22l1_oe_rmats = rpl22l1_oe_rmats.rename({"PValue": "pval", "FDR": "qval"}, axis=1)
rpl22l1_kd1_rmats = rpl22l1_kd1_rmats.rename({"PValue": "pval", "FDR": "qval"}, axis=1)
rpl22l1_kd2_rmats = rpl22l1_kd2_rmats.rename({"PValue": "pval", "FDR": "qval"}, axis=1)
rpl22_a_ko1_rmats = rpl22_a_ko1_rmats.rename({"PValue": "pval", "FDR": "qval"}, axis=1)
rpl22_a_ko2_rmats = rpl22_a_ko2_rmats.rename({"PValue": "pval", "FDR": "qval"}, axis=1)
rpl22_b_ko1_rmats = rpl22_b_ko1_rmats.rename({"PValue": "pval", "FDR": "qval"}, axis=1)
rpl22_b_ko2_rmats = rpl22_b_ko2_rmats.rename({"PValue": "pval", "FDR": "qval"}, axis=1)
```

```python
from functools import reduce
def concat_cols(df, cols, delim):
    cols_str = [df[x].astype(str) for x in cols]

    return reduce(lambda a, b: a + delim + b, cols_str)
def load_se(experiment):
    se = pd.read_csv("../data/raw/rmats_output/"+ experiment +
                     "/SE.MATS.JC.txt", sep="\t", index_col=0)
    
    se["exon_length"] = se["exonEnd"]-se["exonStart_0base"]
    
    se["exon_id"] = concat_cols(
        se, ['chr','exonStart_0base', 'exonEnd',
       'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'], "_")
    
    se = se.set_index("exon_id")
    se = se.sort_values(by="FDR")
    
    return se

rpl22_oe_se = load_se("rpl22_oe")
rpl22l1_oe_se = load_se("rpl22l1_oe")
rpl22l1_kd1_se = load_se("rpl22l1_kd1")
rpl22l1_kd2_se = load_se("rpl22l1_kd2")
rpl22_a_ko1_se = load_se("rpl22_a_ko1")
rpl22_a_ko2_se = load_se("rpl22_a_ko2")
rpl22_b_ko1_se = load_se("rpl22_b_ko1")
rpl22_b_ko2_se = load_se("rpl22_b_ko2")


```

```python
diff_exons = rpl22l1_oe_se["FDR"] < 0.01
# diff_exons = (rpl22l1_oe_se["FDR"] < 0.01) & ((rpl22_a_ko2_se["FDR"] < 0.01))

padding = 100
# max_len = 10000

exon_set = rpl22l1_oe_se[
    ["chr", "exonStart_0base", "exonEnd", "strand", "exon_length", "geneSymbol"]
].copy(deep=True)
# exon_set = exon_set[exon_set["exon_length"] <= max_len]

exon_set["chr"] = exon_set["chr"].apply(lambda x: x[3:])
exon_set["start"] = exon_set["exonStart_0base"] - padding
exon_set["end"] = exon_set["exonEnd"] + padding
exon_set["id"] = range(len(exon_set))
exon_set["id"] = exon_set["geneSymbol"] + "_" + exon_set["id"].astype(str)

pos_exons = exon_set[diff_exons]
neg_exons = exon_set[~diff_exons]

print(len(pos_exons),len(neg_exons))

pos_exons[["chr", "start", "end", "id", "exon_length", "strand"]].to_csv(
    "../data/intermediate/diff_exons_pos.bed", sep="\t", header=False, index=False
)
neg_exons[["chr", "start", "end", "id", "exon_length", "strand"]].to_csv(
    "../data/intermediate/diff_exons_neg.bed", sep="\t", header=False, index=False
)
```

```python
!bedtools getfasta -s -name \
    -fi ../data/raw/reference/hg19.fa -bed \
    ../data/intermediate/diff_exons_pos.bed \
    > ../data/intermediate/diff_exons_pos.fasta

!bedtools getfasta -s -name \
    -fi ../data/raw/reference/hg19.fa -bed \
    ../data/intermediate/diff_exons_neg.bed \
    > ../data/intermediate/diff_exons_neg.fasta
```

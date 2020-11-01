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

# Load splicing results

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

# Extract exon intervals

```python
def output_diff_exons(exon_set, filter_set, output_name, start_col, end_col, padding):

    exon_set = exon_set.copy(deep=True)

    exon_set["chr"] = exon_set["chr"].apply(lambda x: x[3:])
    exon_set["start"] = exon_set[start_col] - padding
    exon_set["end"] = exon_set[end_col] + padding
    exon_set["id"] = range(len(exon_set))
    exon_set["id"] = exon_set["geneSymbol"] + "_" + exon_set["id"].astype(str)

    pos_exons = exon_set[filter_set]
    neg_exons = exon_set[~filter_set]

    print(f"Positive exons: {len(pos_exons)}")
    print(f"Negative exons: {len(neg_exons)}")

    pos_exons[["chr", "start", "end", "id", "exon_length", "strand"]].to_csv(
        f"../data/intermediate/extracted_sequences/{output_name}_pos.bed",
        sep="\t",
        header=False,
        index=False,
    )
    neg_exons[["chr", "start", "end", "id", "exon_length", "strand"]].to_csv(
        f"../data/intermediate/extracted_sequences/{output_name}_neg.bed",
        sep="\t",
        header=False,
        index=False,
    )


output_diff_exons(
    rpl22l1_oe_se,
    rpl22l1_oe_se["FDR"] < 0.01,
    "RPL22L1_OE_SE",
    "exonStart_0base",
    "exonEnd",
    padding=100,
)
output_diff_exons(
    rpl22l1_kd1_se,
    (rpl22l1_kd1_se["FDR"] < 0.01) & (rpl22l1_kd2_se["FDR"] < 0.01),
    "RPL22L1_KD_SE",
    "exonStart_0base",
    "exonEnd",
    padding=100,
)
output_diff_exons(
    rpl22_a_ko1_se,
    (rpl22_a_ko1_se["FDR"] < 0.01) & (rpl22_a_ko2_se["FDR"] < 0.01),
    "RPL22_A_KO_SE",
    "exonStart_0base",
    "exonEnd",
    padding=100,
)
output_diff_exons(
    rpl22_b_ko1_se,
    (rpl22_b_ko1_se["FDR"] < 0.01) & (rpl22_b_ko2_se["FDR"] < 0.01),
    "RPL22_B_KO_SE",
    "exonStart_0base",
    "exonEnd",
    padding=100,
)
```

# Extract output sequences

```python
output_names = ["RPL22L1_OE_SE","RPL22L1_KD_SE","RPL22_A_KO_SE","RPL22_B_KO_SE"]

for name in output_names:
    !bedtools getfasta -s -name \
    -fi ../data/raw/reference/hg19.fa -bed \
    ../data/intermediate/extracted_sequences/{name}_pos.bed \
    > ../data/intermediate/extracted_sequences/{name}_pos.fasta

    !bedtools getfasta -s -name \
        -fi ../data/raw/reference/hg19.fa -bed \
        ../data/intermediate/extracted_sequences/{name}_neg.bed \
        > ../data/intermediate/extracted_sequences/{name}_neg.fasta
```

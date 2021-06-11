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
import pandas as pd

import cancer_data
import many
```

```python
ccle_exonusage = cancer_data.load("ccle_exonusage")
ccle_proteomics = cancer_data.load("ccle_proteomics")

merged_ccle_info = pd.read_csv("../data/supplementary/S1_merged-ccle-info.txt",sep="\t",index_col=0)
```

# Compute differences

```python
msi_prot_diffs = many.stats.mat_mwu_naive(
    ccle_proteomics, merged_ccle_info["MSI"], effect="rank_biserial", melt=True, pbar=True
)
msi_prot_diffs = msi_prot_diffs.reset_index()
msi_prot_diffs = msi_prot_diffs.rename({"a_col":"protein"},axis=1).drop("b_col",axis=1)
msi_prot_diffs.to_csv("../data/supplementary/S5-a__msi-prot-diffs.txt", sep="\t")

msi_exon_diffs = many.stats.mat_mwu_naive(
    ccle_exonusage, merged_ccle_info["MSI"], effect="rank_biserial", melt=True, pbar=True
)
msi_exon_diffs = msi_exon_diffs.reset_index()
msi_exon_diffs = msi_exon_diffs.rename({"a_col":"exon"},axis=1).drop("b_col",axis=1)
msi_exon_diffs.to_csv("../data/intermediate/S5-b_msi-exon-diffs.txt", sep="\t")
```

# CCLE BAMs

```python
ccle_samples = pd.read_csv("../data/raw/fullccle_samples.csv")
wgs_samples = ccle_samples[ccle_samples["datatype"] == "wgs"]
wgs_samples = wgs_samples.sort_values(by=["arxspan_id", "version", "size"])
wgs_samples = wgs_samples.drop_duplicates(subset=["arxspan_id"], keep="first")
```

```python
msi_exons = msi_exon_diffs.copy()[msi_exon_diffs["qval"] >= 4]

msi_exons["gene"] = msi_exons["exon"].map(lambda x: x.split("_")[-1])
msi_exons["exon_chrom"] = msi_exons["exon"].map(lambda x: x.split("_")[-4][3:])
msi_exons["exon_start"] = msi_exons["exon"].map(lambda x: x.split("_")[-3])
msi_exons["exon_end"] = msi_exons["exon"].map(lambda x: x.split("_")[-2])

msi_exons["exon"] = (
    msi_exons["exon_chrom"]
    + "_"
    + msi_exons["exon_start"]
    + "_"
    + msi_exons["exon_end"]
)

msi_exons["exon_start"] = msi_exons["exon_start"].astype(int)
msi_exons["exon_end"] = msi_exons["exon_end"].astype(int)

msi_exons = msi_exons.drop_duplicates(subset=["exon"])
```

```python
def get_exon_bounds(row, padding=1000):

    if row["exon_start"] <= row["exon_end"]:
        row["bound_start"] = row["exon_start"] - padding
        row["bound_end"] = row["exon_end"] + padding

    elif row["exon_start"] >= row["exon_end"]:
        row["bound_start"] = row["exon_end"] - padding
        row["bound_end"] = row["exon_start"] + padding

    return row


msi_exons = msi_exons.apply(get_exon_bounds, axis=1)
```

```python
msi_exons[["exon_chrom", "bound_start", "bound_end"]].to_csv(
    "../scripts/MSI_exon_bounds.bed", sep="\t", header=False, index=False
)
msi_exons[["exon_chrom", "bound_start", "bound_end"]].to_csv(
    "../data/raw/MSI_exon_bounds.bed", sep="\t", header=False, index=False
)
```

```python
with open("../scripts/wgs_paths.txt", "w") as f:
    for bam_path in list(wgs_samples["internal_bam_filepath"]):
        f.write(bam_path + "\n")

with open("../scripts/wgs_ids.txt", "w") as f:
    for bam_path in list(wgs_samples["arxspan_id"]):
        f.write(bam_path + "\n")
```

```python
with open("../scripts/7_fetch-msi-slices.sh", "w") as f:
    for bam_path, ach_id in zip(
        list(wgs_samples["internal_bam_filepath"]), list(wgs_samples["arxspan_id"])
    ):

        f.write(
            "GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token) samtools view -b -h -M -L ../data/raw/MSI_exon_bounds.bed {} > ../data/raw/WGS_slices/{}.bam\n".format(
                bam_path, ach_id
            )
        )
```

# Read mutations

```python
msi_exon_calls = pd.read_csv(
    "../data/intermediate/msi_exon_calls/msi_exon_calls_filtered.txt",
    sep="\t",
    names=["chrom", "pos", "sample", "ref", "alt", "genotype"],
)
```

```python
msi_exon_calls["ref_len"] = msi_exon_calls[""]
```

```python
msi_exon_calls["truncated"] = msi_exon_calls[msi_exon_calls["alt"]]
```

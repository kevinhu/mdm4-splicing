---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.2
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

```python
import pandas as pd
import numpy as np
import json
```

# Make TPM matrix


## Experiments list

```python
experiments = pd.read_csv("../data/intermediate/experiment_setup.txt", sep="\t")

experiments["kallisto_path"] = (
    "../data/intermediate/kallisto_quant/" + experiments["sample"] + "/abundance.tsv"
)


def load_tpms(exp_path, exp_name):
    abundances = pd.read_csv(exp_path, sep="\t", index_col=0)

    tpms = abundances["tpm"]
    tpms = tpms.rename(exp_name)
    tpms = tpms.astype(np.float64)

    return tpms
```

## Merge transcript TPMs

```python
exps = zip(experiments["kallisto_path"], experiments["sample"])

# read in TPMs for each experiment
transcript_tpms = [load_tpms(x[0], x[1]) for x in exps]
transcript_tpms = pd.concat(transcript_tpms, axis=1)

transcript_tpms.to_csv("../data/processed/transcript_tpms.txt", sep="\t")
```

## Merge gene TPMs

```python
t2g = pd.read_csv("../data/intermediate/sleuth_diff/ensembl_t2g.csv")

_, align_ensembl_genes = transcript_tpms.align(
    t2g.set_index("target_id")["ens_gene"], axis=0, join="inner"
)

# group transcripts by gene and sum
gene_tpms = transcript_tpms.groupby(align_ensembl_genes).sum()

gene_tpms.to_csv("../data/processed/gene_tpms.txt", sep="\t")
```

# Sleuth outputs

```python
with open("experiments.json", "r") as f:
    exp = json.load(f)

    experiments = exp["experiments"]
    experiment_ids = exp["experiment_ids"]
```

## Helper functions

```python
def compute_medians(sleuth_diff, experiment):

    # compute the medians in control/treatment groups
    sleuth_diff["control_median"] = sleuth_diff[experiments[experiment][0]].median(
        axis=1
    )
    sleuth_diff["treatment_median"] = sleuth_diff[experiments[experiment][1]].median(
        axis=1
    )
    sleuth_diff["median_foldchange"] = (
        sleuth_diff["treatment_median"] / sleuth_diff["control_median"]
    )


def signed_p_rank(sleuth_diff):

    # compute -log10 P-values
    sleuth_diff["-log_pval"] = -np.log10(sleuth_diff["pval"])

    # replace ultra-low P-values
    sleuth_diff["-log_pval"] = sleuth_diff["-log_pval"].replace(np.inf, 320)

    sleuth_diff["treatment_increase"] = -1 + 2 * (sleuth_diff["median_foldchange"] > 1)

    sleuth_diff["signed_pval"] = (
        sleuth_diff["-log_pval"] * sleuth_diff["treatment_increase"]
    )

    sleuth_diff = sleuth_diff.sort_values(by="signed_pval")
```

## Process transcripts

```python
def process_sleuth_transcripts(experiment):
    sleuth_diff = pd.read_csv(
        "../data/intermediate/sleuth_diff/" + experiment + "_transcripts.csv",
        index_col=1,
    )

    controls = experiments[experiment][0]
    treatments = experiments[experiment][1]

    sleuth_tpms = transcript_tpms.loc[sleuth_diff.index, controls + treatments]

    sleuth_diff = pd.concat([sleuth_diff, sleuth_tpms], axis=1)

    # compute medians and P-values
    compute_medians(sleuth_diff, experiment)
    signed_p_rank(sleuth_diff)

    # remove missing values
    sleuth_diff = sleuth_diff.dropna(
        subset=["pval", "median_foldchange", "entrez_gene"], how="any"
    )

    # format Entrez ID as string
    sleuth_diff["entrez_gene"] = sleuth_diff["entrez_gene"].astype(int).astype(str)

    sleuth_diff.to_csv(
        "../data/processed/kallisto_sleuth_merge/" + experiment + "_transcripts.txt",
        sep="\t",
    )
    sleuth_diff.to_hdf(
        "../data/processed/kallisto_sleuth_merge/" + experiment + "_transcripts.h5",
        key="sleuth_diff",
        mode="w",
    )
```

## Process genes

```python
def process_sleuth_genes(experiment):
    sleuth_diff = pd.read_csv(
        "../data/intermediate/sleuth_diff/" + experiment + "_genes.csv", index_col=1
    )

    controls = experiments[experiment][0]
    treatments = experiments[experiment][1]

    sleuth_tpms = gene_tpms.loc[sleuth_diff.index, controls + treatments]

    sleuth_diff = pd.concat([sleuth_diff, sleuth_tpms], axis=1)

    # compute medians and P-values
    compute_medians(sleuth_diff, experiment)
    signed_p_rank(sleuth_diff)

    # remove missing values
    sleuth_diff = sleuth_diff.dropna(
        subset=["pval", "median_foldchange", "target_id"], how="any"
    )

    # drop biotype column
    sleuth_diff = sleuth_diff.drop(["transcript_biotype"], axis=1)
    sleuth_diff = sleuth_diff[~sleuth_diff.index.duplicated(keep="first")]

    # format Entrez ID as string
    sleuth_diff["target_id"] = sleuth_diff["target_id"].astype(int).astype(str)

    sleuth_diff.to_csv(
        "../data/processed/kallisto_sleuth_merge/" + experiment + "_genes.txt", sep="\t"
    )
    sleuth_diff.to_hdf(
        "../data/processed/kallisto_sleuth_merge/" + experiment + "_genes.h5",
        key="sleuth_diff",
        mode="w",
    )
```

## Apply over experiments

```python
for exp_id in experiment_ids:

    print(exp_id)

    process_sleuth_genes(exp_id)
    process_sleuth_transcripts(exp_id)
```

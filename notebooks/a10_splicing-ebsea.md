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
from functools import reduce

import ujson
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib_venn import venn2

import config
config.config_visuals()

FDR_CUTOFF = 0.01
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
```

```python
def get_coords(rmats_set, name):

    filtered_set = rmats_set.copy(deep=True)

    filtered_set = filtered_set[filtered_set["FDR"] < FDR_CUTOFF]

    filtered_set["chrom"] = filtered_set.index.map(lambda x: x.split("_")[-7])

    filtered_set["exon_start"] = filtered_set.index.map(lambda x: x.split("_")[-6])
    filtered_set["exon_end"] = filtered_set.index.map(lambda x: x.split("_")[-5])
    filtered_set["dpsi"] = (
        filtered_set["treatment_median"] - filtered_set["control_median"]
    )

    filtered_set[["chrom", "exon_start", "exon_end", "dpsi"]].to_csv(
        f"../data/processed/rmats_dease_hg19/{name}.bed",
        index=False,
        sep="\t",
        header=False,
    )


get_coords(rpl22_oe_rmats, "rpl22_oe")
```

```python

```

```python
rpl22_oe_rmats
```

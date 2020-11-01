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
import pandas as pd
```

```python
merged_ccle_info = pd.read_csv("../data/supplementary/S1_merged-ccle-info.txt", sep="\t", index_col=0)
```

```python
merged_ccle_info
```

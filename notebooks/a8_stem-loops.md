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
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
```

```python
peak_anno = pd.read_csv("../data/intermediate/RPL22-ZR751_clearCLIP.pool.tag.uniq.del.CIMS.fdr10.f10.bed.anno.csv",
                        index_col=0)
```

```python
def get_distance_to_tss(row):
    if row["geneStrand"] == 1:
        
        start_diff = row["start"] - row["geneStart"]
        end_diff = row["end"] - row["geneStart"]
    else:
        start_diff = row["geneEnd"] - row["start"]
        end_diff = row["geneEnd"]-row["end"]
        
    if abs(start_diff) >= abs(end_diff):
        return end_diff
    else:
        return start_diff
    
def get_distance_to_tts(row):
    if row["geneStrand"] == 1:
        
        start_diff = row["start"] - row["geneEnd"]
        end_diff = row["end"] - row["geneEnd"]
    else:
        start_diff = row["geneStart"] - row["start"]
        end_diff = row["geneStart"]-row["end"]
        
    if abs(start_diff) >= abs(end_diff):
        return end_diff
    else:
        return start_diff
    
peak_anno["distance_to_tss"] = peak_anno.apply(get_distance_to_tss, axis=1)
peak_anno["distance_to_tts"] = peak_anno.apply(get_distance_to_tts, axis=1)
```

```python
sns.distplot(peak_anno["distanceToTSS"],bins=np.linspace(-1000,1000,100))
plt.xlim(-1000,1000)
```

```python
peak_anno[(peak_anno["distance_to_tts"]<0)&(peak_anno["distance_to_tts"]>-250)]
```

```python
sns.distplot(peak_anno["distance_to_tts"],bins=np.linspace(-1000,1000,100))
plt.xlim(-1000,1000)
```

```python

```

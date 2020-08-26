# MDM4-splicing

Code, figures, and prior literature for *RPL22 mutations drive splicing-mediated activation of MDM4*.

## Overview

This repository contains the scripts for generating all the figures and supplementary data in the manuscript. The raw data used here is quite large and is therefore not included, but may be found at [insert link here].

## Requirements

The `cancer_data` and `many` repositories contain supplementary methods for external datasets as well as several statistical pipelines and plots used here. These may be found at their relevant GitHub repositories as well as on PyPI. Additional packages such as `pandas`, `numpy`, `matplotlib`, and `seaborn` are also required

## Structure

The analysis scripts are contained within `/scripts` and `/notebooks`. To reproduce the processed data and plots of this paper, use the following steps:

1. Make sure raw files have been downloaded from [insert link here]
2. Run the scripts, which process the RNAseq data.
3. Run the notebooks for generating the supplementary tables (in order). These are identified as having a name that begins with `tabS[n]` .
4. Run the notebooks for generating the main figures. These are identified as having a name that begins with `fig[n]`.
5. Run the notebooks for generating the supplementary figures. These are identified as having a name that begins with `figS[n]`.

The figures are generated as subplots and output to `/plots`. These subplots were then manually assembled and edited to produce the final figures found in `/figures`.


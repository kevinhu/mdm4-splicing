# MDM4-splicing

Code, figures, and prior literature for [title].

## Overview

This repository contains the scripts for generating all the figures and supplementary data in the manuscript. The raw data used here is quite and is therefore not included, but may be found [insert link here].

## Requirements

The `galileo` and `huygens` repositories contain supplementary methods for several statistical pipelines and plots used here. These can be found at the relevant GitHub repositories on my profile. Additional packages such as `pandas`, `numpy`, `matplotlib`, and `seaborn` are also referenced.

## Structure

The analysis scripts are all contained within `/notebooks`. To reproduce the processed data and plots of this paper, use the following steps:

1. Make sure raw files have been downloaded from [insert link here]
2. Run the notebooks for generating the supplementary tables (in order). These are identified as having a name that begins with `tabS[n]` .
3. Run the notebooks for generating the main figures. These are identified as having a name that begins with `fig[n]`.
4. Run the notebooks for generating the supplementary figures. These are identified as having a name that begins with `figS[n]`.

The figures are generated as subplots and output to `/plots`. These subplots were then manually assembled and edited to produce the final figures found in `/figures`.


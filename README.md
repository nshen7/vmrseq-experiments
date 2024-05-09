# Scripts for experiments in vmrseq paper

This repository contains the code to reproduce the analyses and figures in "vmrseq: Probabilistic Modeling of Single-cell Methylation Heterogeneity". These scripts rely on the package vmrseq from [github](https://github.com/nshen7/vmrseq) (v0.99.0).

The directory contains three folders:

1. **code** - scripts for preprocessing steps, the analyses with `vmrseq` and other tools used in the comparisons.
2. **plots** - original plots generated from the scripts in **code** folder.
3. **manuscript_related** - publication-level figures and tables in the manuscript along with the scripts generated them.

To reproduce all results from the paper, run the code in the scripts in **code/raw_data_process** first to dowload and pre-process datasets used in the analyses, then scripts in **code/sim_studies** and **code/case_studies** can be run independently. Finally, run the scripts in **manuscript_related** to reproduce the results published in paper.

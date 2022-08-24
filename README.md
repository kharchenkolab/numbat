# Numbat

<!-- badges: start -->

[![<kharchenkolab>](https://circleci.com/gh/kharchenkolab/numbat.svg?style=svg)](https://app.circleci.com/pipelines/github/kharchenkolab/numbat)
  
<!-- badges: end -->

<img src="logo.png" align="right" width="150">

Numbat is a haplotype-aware CNV caller from single-cell transcriptomics data. It integrates signals from gene expression, allelic ratio, and population-derived haplotype information to accurately infer allele-specific CNVs in single cells and reconstruct their lineage relationship. 

Numbat can be used to:
 1. Detect allele-specific copy number variations from scRNA-seq 
 2. Differentiate tumor versus normal cells in the tumor microenvironment 
 3. Infer the clonal architecture and evolutionary history of profiled tumors. 

![image](https://user-images.githubusercontent.com/13375875/153020818-2e782689-09db-427f-ad98-2c175021a936.png)

Numbat does not require paired DNA or genotype data and operates solely on the donor scRNA-data data (for example, 10x Cell Ranger output). For details of the method, please checkout our preprint:

[Teng Gao, Ruslan Soldatov, Hirak Sarkar, et al. Haplotype-enhanced inference of somatic copy number profiles from single-cell transcriptomes. bioRxiv 2022.](https://www.biorxiv.org/content/10.1101/2022.02.07.479314v1)

# User Guide
For a complete guide, please see [Numbat User Guide](https://kharchenkolab.github.io/numbat/).

# Questions?
We appreciate your feedback! Please raise a github issue or [email](mailto:tgaoteng@gmail.com) us.

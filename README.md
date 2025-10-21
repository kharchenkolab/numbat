# Numbat

<!-- badges: start -->
[![<kharchenkolab>](https://circleci.com/gh/kharchenkolab/numbat.svg?style=svg)](https://app.circleci.com/pipelines/github/kharchenkolab/numbat)
[![CRAN status](https://www.r-pkg.org/badges/version/numbat)](https://cran.r-project.org/package=numbat)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/numbat)](https://cran.r-project.org/package=numbat)
<!-- badges: end -->

<img src="logo.png" align="right" width="150">

Numbat is a haplotype-aware CNV caller from single-cell and spatial transcriptomics data. It integrates signals from gene expression, allelic ratio, and population-derived haplotype information to accurately infer allele-specific CNVs in single cells and reconstruct their lineage relationship. 

Numbat can be used to:
 1. Detect allele-specific copy number variations from scRNA-seq and spatial transcriptomics
 2. Differentiate tumor versus normal cells in the tumor microenvironment 
 3. Infer the clonal architecture and evolutionary history of profiled tumors. 

![image](https://user-images.githubusercontent.com/13375875/153020818-2e782689-09db-427f-ad98-2c175021a936.png)

Numbat does not require paired DNA or genotype data and operates solely on the donor scRNA-seq data (for example, 10x Cell Ranger output). For details of the method, please checkout our paper:

> [Teng Gao, Ruslan Soldatov, Hirak Sarkar, Adam Kurkiewicz, Evan Biederstedt, Po-Ru Loh, Peter Kharchenko. Haplotype-aware analysis of somatic copy number variations from single-cell transcriptomes. Nature Biotechnology (2022).](https://www.nature.com/articles/s41587-022-01468-y)

## Numbat-multiome
Numbat was later extended to multi-modality (single-cell RNA and ATAC) data. Check out the [vignette](https://kharchenkolab.github.io/numbat/articles/numbat-multiome.html) and paper below:

> [Ruitong Li, Jean-Baptiste Alberge, Tina Keshavarzian, Junko Tsuji, Johan Gustafsson, Mahshid Rahmat, Elizabeth D Lightbody, Stephanie L Deng, Santiago Riviero, Mendy Miller, F Naz Cemre Kalayci, Adrian Wiestner, Clare Sun, Mathieu Lupien, Irene Ghobrial, Erin Parry, Teng Gao, Gad Getz. Numbat-multiome: inferring copy number variations by combining RNA and chromatin accessibility information from single-cell data. Briefings in Bioinformatics (2025).](https://academic.oup.com/bib/article/26/5/bbaf516/8290422)


# User Guide
For a complete guide, please see [Numbat User Guide](https://kharchenkolab.github.io/numbat/).

# Questions?
We appreciate your feedback! Please raise a github [issue](https://github.com/kharchenkolab/numbat/issues) for bugs, questions and new feature requests. For bug reports, please attach full log, error message, input parameters, and ideally a reproducible example (if possible).


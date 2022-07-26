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

# Docker

We've provided Dockerfiles for both base R and RStudio. 

## Ready-to-run Docker image

To start a docker container for numbat, first install docker on your platform and then start the numbat for RStudio  with the following command in the shell:

```
docker run -p 8080:8080 -e PASSWORD=pass pkharchenkolab/numbat-rstudio:latest
```

The first time you run this command, it will download several large images so make sure that you have fast internet access setup. You can then point your browser to http://localhost:8080/ to get an Rstudio environment with `numbat` installed (please log in using credentials username=`rstudio`, password=`pass`). Explore the [docker --mount option](https://docs.docker.com/storage/volumes/) to allow access of the docker image to your local files.

The base R image can be run as follows:

```
docker run -it pkharchenkolab/numbat-rbase:latest /bin/bash
```

**Note:** If you already downloaded the docker image and want to update it, please pull the latest image with: 
```
docker pull pkharchenkolab/numbat-rbase:latest
```
or 
```
docker pull pkharchenkolab/numbat-rstudio:latest
```

The docker images can be found on Dockerhub for [R Studio here](https://hub.docker.com/r/pkharchenkolab/numbat-rstudio) and [base R here](https://hub.docker.com/r/pkharchenkolab/numbat-rbase).


## Building Docker image from the Dockerfile

The Dockerfiles are located in the subfolder `/docker`

For R studio, the command to build the Docker image from scratch is:

```
docker build -f Dockerfile.rstudio -t numbat-rstudio .
```

For base R, the analogous command is:


```
docker build -f Dockerfile.rbase -t numbat-rbase .
```


# Questions?
We appreciate your feedback! Please raise a github issue or [email](mailto:tgaoteng@gmail.com) us.

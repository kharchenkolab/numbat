# Preprocess allele data

Preprocess allele data

## Usage

``` r
preprocess_allele(sample, vcf_pu, vcf_phased, AD, DP, barcodes, gtf, gmap)
```

## Arguments

- sample:

  character Sample label

- vcf_pu:

  dataframe Pileup VCF from cell-snp-lite

- vcf_phased:

  dataframe Phased VCF from eagle2

- AD:

  dgTMatrix Alt allele depth matrix from pileup

- DP:

  dgTMatrix Total allele depth matrix from pileup

- barcodes:

  vector List of barcodes from pileup

- gtf:

  dataframe Transcript GTF

- gmap:

  dataframe Genetic map

## Value

dataframe Allele counts by cell

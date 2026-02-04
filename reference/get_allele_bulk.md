# Aggregate into pseudobulk alelle profile

Aggregate into pseudobulk alelle profile

## Usage

``` r
get_allele_bulk(df_allele, nu = 1, min_depth = 0)
```

## Arguments

- df_allele:

  dataframe Single-cell allele counts

- nu:

  numeric Phase switch rate

- min_depth:

  integer Minimum coverage to filter SNPs

## Value

dataframe Pseudobulk allele profile

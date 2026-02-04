# Aggregate single-cell data into combined bulk expression and allele profile

Aggregate single-cell data into combined bulk expression and allele
profile

## Usage

``` r
get_bulk(
  count_mat,
  lambdas_ref,
  df_allele,
  gtf,
  subset = NULL,
  min_depth = 0,
  nu = 1,
  segs_loh = NULL,
  verbose = TRUE
)
```

## Arguments

- count_mat:

  dgCMatrix Gene expression counts

- lambdas_ref:

  matrix Reference expression profiles

- df_allele:

  dataframe Single-cell allele counts

- gtf:

  dataframe Transcript gtf

- subset:

  vector Subset of cells to aggregate

- min_depth:

  integer Minimum coverage to filter SNPs

- nu:

  numeric Phase switch rate

- segs_loh:

  dataframe Segments with clonal LOH to be excluded

- verbose:

  logical Verbosity

## Value

dataframe Pseudobulk gene expression and allele profile

## Examples

``` r
bulk_example = get_bulk(
    count_mat = count_mat_example,
    lambdas_ref = ref_hca,
    df_allele = df_allele_example,
    gtf = gtf_hg38)
```

# Make a group of pseudobulks

Make a group of pseudobulks

## Usage

``` r
make_group_bulks(
  groups,
  count_mat,
  df_allele,
  lambdas_ref,
  gtf,
  min_depth = 0,
  nu = 1,
  segs_loh = NULL,
  ncores = NULL
)
```

## Arguments

- groups:

  list Contains fields named "sample", "cells", "size", "members"

- count_mat:

  dgCMatrix Gene counts

- df_allele:

  dataframe Alelle counts

- lambdas_ref:

  matrix Reference expression profiles

- gtf:

  dataframe Transcript GTF

- min_depth:

  integer Minimum allele depth to include

- segs_loh:

  dataframe Segments with clonal LOH to be excluded

- ncores:

  integer Number of cores

## Value

dataframe Pseudobulk profiles

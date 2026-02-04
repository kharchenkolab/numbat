# Run smoothed expression-based hclust

Run smoothed expression-based hclust

## Usage

``` r
exp_hclust(
  count_mat,
  lambdas_ref,
  gtf,
  sc_refs = NULL,
  window = 101,
  ncores = 1,
  verbose = TRUE
)
```

## Arguments

- count_mat:

  dgCMatrix Gene counts

- lambdas_ref:

  matrix Reference expression profiles

- gtf:

  dataframe Transcript GTF

- sc_refs:

  named list Reference choices for single cells

- window:

  integer Sliding window size

- ncores:

  integer Number of cores

- verbose:

  logical Verbosity

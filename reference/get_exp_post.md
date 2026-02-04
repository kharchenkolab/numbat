# compute single-cell expression posteriors

compute single-cell expression posteriors

## Usage

``` r
get_exp_post(
  segs_consensus,
  count_mat,
  gtf,
  lambdas_ref,
  sc_refs = NULL,
  diploid_chroms = NULL,
  use_loh = NULL,
  segs_loh = NULL,
  ncores = 30,
  verbose = TRUE,
  debug = FALSE
)
```

## Arguments

- segs_consensus:

  dataframe Consensus segments

- count_mat:

  dgCMatrix gene expression count matrix

- gtf:

  dataframe transcript gtf

- lambdas_ref:

  matrix Reference expression profiles

## Value

dataframe Expression posteriors

# get the single cell expression likelihoods

get the single cell expression likelihoods

## Usage

``` r
get_exp_likelihoods(
  exp_counts,
  diploid_chroms = NULL,
  use_loh = FALSE,
  depth_obs = NULL,
  mu = NULL,
  sigma = NULL
)
```

## Arguments

- exp_counts:

  dataframe Single-cell expression counts (CHROM, seg, cnv_state, gene,
  Y_obs, lambda_ref)

- diploid_chroms:

  character vector Known diploid chromosomes

- use_loh:

  logical Whether to include CNLOH regions in baseline

## Value

dataframe Single-cell CNV likelihood scores

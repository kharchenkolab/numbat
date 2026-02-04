# retest consensus segments on pseudobulks

retest consensus segments on pseudobulks

## Usage

``` r
retest_bulks(
  bulks,
  segs_consensus = NULL,
  t = 1e-05,
  min_genes = 10,
  gamma = 20,
  nu = 1,
  use_loh = FALSE,
  diploid_chroms = NULL,
  ncores = 1,
  exclude_neu = TRUE,
  min_LLR = 5
)
```

## Arguments

- bulks:

  dataframe Pseudobulk profiles

- segs_consensus:

  dataframe Consensus segments

- use_loh:

  logical Whether to use loh in the baseline

- diploid_chroms:

  vector User-provided diploid chromosomes

## Value

dataframe Retested pseudobulks

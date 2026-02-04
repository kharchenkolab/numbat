# Run multiple HMMs

Run multiple HMMs

## Usage

``` r
run_group_hmms(
  bulks,
  t = 1e-04,
  gamma = 20,
  alpha = 1e-04,
  min_genes = 10,
  nu = 1,
  common_diploid = TRUE,
  diploid_chroms = NULL,
  allele_only = FALSE,
  retest = TRUE,
  run_hmm = TRUE,
  exclude_neu = TRUE,
  ncores = 1,
  verbose = FALSE,
  debug = FALSE
)
```

## Arguments

- bulks:

  dataframe Pseudobulk profiles

- t:

  numeric Transition probability

- gamma:

  numeric Dispersion parameter for the Beta-Binomial allele model

- alpha:

  numeric P value cut-off to determine segment clusters in find_diploid

- common_diploid:

  logical Whether to find common diploid regions between pseudobulks

- diploid_chroms:

  character vector Known diploid chromosomes to use as baseline

- allele_only:

  logical Whether only use allele data to run HMM

- retest:

  logcial Whether to retest CNVs

- run_hmm:

  logical Whether to run HMM segments or just retest

- ncores:

  integer Number of cores

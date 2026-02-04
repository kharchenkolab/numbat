# retest CNVs in a pseudobulk

retest CNVs in a pseudobulk

## Usage

``` r
retest_cnv(
  bulk,
  theta_min = 0.08,
  logphi_min = 0.25,
  gamma = 20,
  allele_only = FALSE,
  exclude_neu = TRUE
)
```

## Arguments

- bulk:

  pesudobulk dataframe

- gamma:

  numeric Dispersion parameter for the Beta-Binomial allele model

- allele_only:

  whether to retest only using allele data

## Value

a dataframe of segments with CNV posterior information

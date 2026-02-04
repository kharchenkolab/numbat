# Get an allele HMM

Get an allele HMM

## Usage

``` r
get_allele_hmm(pAD, DP, p_s, theta, gamma = 20)
```

## Arguments

- pAD:

  integer vector Paternal allele counts

- DP:

  integer vector Total alelle counts

- p_s:

  numeric vector Phase switch probabilities

- theta:

  numeric Haplotype imbalance

- gamma:

  numeric Overdispersion in the allele-specific expression

## Value

HMM object

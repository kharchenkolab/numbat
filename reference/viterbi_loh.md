# Viterbi for clonal LOH detection

Viterbi for clonal LOH detection

## Usage

``` r
viterbi_loh(hmm, ...)
```

## Arguments

- hmm:

  HMM object; expect variables x (SNP count), snp_sig (snp rate standard
  deviation), pm (snp density for ref and loh states), pn (gene
  lengths), d (total expression depth), y (expression count),
  lambda_star (reference expression rate), mu (global expression mean),
  sig (global expression standard deviation), Pi (transition prob
  matrix), delta (prior for each state), phi (expression fold change for
  each state)

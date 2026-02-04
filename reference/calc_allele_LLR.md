# Calculate LLR for an allele HMM

Calculate LLR for an allele HMM

## Usage

``` r
calc_allele_LLR(pAD, DP, p_s, theta_mle, theta_0 = 0, gamma = 20)
```

## Arguments

- pAD:

  numeric vector Phased allele depth

- DP:

  numeric vector Total allele depth

- p_s:

  numeric vector Phase switch probabilities

- theta_mle:

  numeric MLE of imbalance level theta (alternative hypothesis)

- theta_0:

  numeric Imbalance level in the null hypothesis

- gamma:

  numeric Dispersion parameter for the Beta-Binomial allele model

## Value

numeric Log-likelihood ratio

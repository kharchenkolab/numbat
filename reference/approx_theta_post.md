# Laplace approximation of the posterior of allelic imbalance theta

Laplace approximation of the posterior of allelic imbalance theta

## Usage

``` r
approx_theta_post(
  pAD,
  DP,
  p_s,
  lower = 0.001,
  upper = 0.499,
  start = 0.25,
  gamma = 20
)
```

## Arguments

- pAD:

  numeric vector Variant allele depth

- DP:

  numeric vector Total allele depth

- p_s:

  numeric vector Variant allele frequency

- lower:

  numeric Lower bound of theta

- upper:

  numeric Upper bound of theta

- start:

  numeric Starting value of theta

- gamma:

  numeric Gamma parameter of the beta-binomial distribution

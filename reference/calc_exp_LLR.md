# Calculate LLR for an expression HMM

Calculate LLR for an expression HMM

## Usage

``` r
calc_exp_LLR(
  Y_obs,
  lambda_ref,
  d,
  phi_mle,
  mu = NULL,
  sig = NULL,
  alpha = NULL,
  beta = NULL
)
```

## Arguments

- Y_obs:

  numeric vector Gene expression counts

- lambda_ref:

  numeric vector Reference expression levels

- d:

  numeric vector Total library size

- phi_mle:

  numeric MLE of expression fold change phi (alternative hypothesis)

- mu:

  numeric Mean parameter for the PLN expression model

- sig:

  numeric Dispersion parameter for the PLN expression model

- alpha:

  numeric Hyperparameter for the gamma poisson model (not used)

- beta:

  numeric Hyperparameter for the gamma poisson model (not used)

## Value

numeric Log-likelihood ratio

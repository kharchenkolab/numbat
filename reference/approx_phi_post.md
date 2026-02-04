# Laplace approximation of the posterior of expression fold change phi

Laplace approximation of the posterior of expression fold change phi

## Usage

``` r
approx_phi_post(
  Y_obs,
  lambda_ref,
  d,
  alpha = NULL,
  beta = NULL,
  mu = NULL,
  sig = NULL,
  lower = 0.2,
  upper = 10,
  start = 1
)
```

## Arguments

- Y_obs:

  numeric vector Gene expression counts

- lambda_ref:

  numeric vector Reference expression levels

- d:

  numeric Total library size

- alpha:

  numeric Shape parameter of the gamma distribution

- beta:

  numeric Rate parameter of the gamma distribution

- mu:

  numeric Mean of the normal distribution

- sig:

  numeric Standard deviation of the normal distribution

- lower:

  numeric Lower bound of phi

- upper:

  numeric Upper bound of phi

- start:

  numeric Starting value of phi

## Value

numeric MLE of phi and its standard deviation

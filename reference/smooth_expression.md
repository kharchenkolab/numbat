# filtering, normalization and capping

filtering, normalization and capping

## Usage

``` r
smooth_expression(count_mat, lambdas_ref, gtf, window = 101, verbose = FALSE)
```

## Arguments

- count_mat:

  dgCMatrix Gene expression counts

- lambdas_ref:

  matrix Reference expression profiles

- gtf:

  dataframe Transcript gtf

## Value

dataframe Log(x+1) transformed normalized expression values for single
cells

# filter for mutually expressed genes

filter for mutually expressed genes

## Usage

``` r
filter_genes(count_mat, lambdas_ref, gtf, verbose = FALSE)
```

## Arguments

- count_mat:

  dgCMatrix Gene expression counts

- lambdas_ref:

  named numeric vector A reference expression profile

- gtf:

  dataframe Transcript gtf

## Value

vector Genes that are kept after filtering

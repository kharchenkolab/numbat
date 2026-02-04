# choose beest reference for each cell based on correlation

choose beest reference for each cell based on correlation

## Usage

``` r
choose_ref_cor(count_mat, lambdas_ref, gtf)
```

## Arguments

- count_mat:

  dgCMatrix Gene expression counts

- lambdas_ref:

  matrix Reference expression profiles

- gtf:

  dataframe Transcript gtf

## Value

named vector Best references for each cell

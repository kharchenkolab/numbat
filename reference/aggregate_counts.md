# Utility function to make reference gene expression profiles

Utility function to make reference gene expression profiles

## Usage

``` r
aggregate_counts(count_mat, annot, normalized = TRUE, verbose = TRUE)
```

## Arguments

- count_mat:

  matrix/dgCMatrix Gene expression counts

- annot:

  dataframe Cell annotation with columns "cell" and "group"

- normalized:

  logical Whether to return normalized expression values

- verbose:

  logical Verbosity

## Value

matrix Reference gene expression levels

## Examples

``` r
ref_custom = aggregate_counts(count_mat_ref, annot_ref, verbose = FALSE)
```

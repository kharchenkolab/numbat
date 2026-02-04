# Plot a group of pseudobulk HMM profiles

Plot a group of pseudobulk HMM profiles

## Usage

``` r
plot_bulks(bulks, ..., ncol = 1, title = TRUE, title_size = 8)
```

## Arguments

- bulks:

  dataframe Pseudobulk profiles annotated with "sample" column

- ...:

  additional parameters passed to plot_psbulk()

- ncol:

  integer Number of columns

- title:

  logical Whether to add titles to individual plots

- title_size:

  numeric Size of titles

## Value

a ggplot object

## Examples

``` r
p = plot_bulks(bulk_example)
#> Warning: Arguments in `...` must be used.
#> ✖ Problematic argument:
#> • na.rm = TRUE
#> ℹ Did you misspell an argument name?
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the numbat package.
#>   Please report the issue to the authors.
```

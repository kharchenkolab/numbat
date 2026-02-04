# Plot single-cell smoothed expression magnitude heatmap

Plot single-cell smoothed expression magnitude heatmap

## Usage

``` r
plot_sc_tree(
  gtree,
  label_mut = TRUE,
  label_size = 3,
  dot_size = 2,
  branch_width = 0.5,
  tip = TRUE,
  tip_length = 0.5,
  pal_clone = NULL
)
```

## Arguments

- gtree:

  tbl_graph The single-cell phylogeny

- label_mut:

  logical Whether to label mutations

- label_size:

  numeric Size of mutation labels

- dot_size:

  numeric Size of mutation nodes

- branch_width:

  numeric Width of branches in tree

- tip:

  logical Whether to plot tip point

- tip_length:

  numeric Length of the tips

- pal_clone:

  named vector Clone colors

## Value

ggplot A single-cell phylogeny with mutation history labeled

## Examples

``` r
p = plot_sc_tree(phylogeny_example)
#> Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the numbat package.
#>   Please report the issue to the authors.
```

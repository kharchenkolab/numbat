# Plot single-cell smoothed expression magnitude heatmap

Plot single-cell smoothed expression magnitude heatmap

## Usage

``` r
plot_exp_roll(
  gexp_roll_wide,
  hc,
  k,
  gtf,
  lim = 0.8,
  n_sample = 300,
  reverse = TRUE,
  plot_tree = TRUE
)
```

## Arguments

- gexp_roll_wide:

  matrix Cell x gene smoothed expression magnitudes

- hc:

  hclust Hierarchical clustring result

- k:

  integer Number of clusters

- gtf:

  dataframe Transcript GTF

- lim:

  numeric Limit for expression magnitudes

- n_sample:

  integer Number of cells to subsample

- reverse:

  logical Whether to reverse the cell order

- plot_tree:

  logical Whether to plot the dendrogram

## Value

ggplot A single-cell heatmap of window-smoothed expression CNV signals

## Examples

``` r
p = plot_exp_roll(gexp_roll_example, gtf = gtf_hg38, hc = hc_example, k = 3)
```

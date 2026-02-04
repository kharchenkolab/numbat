# Plot CNV heatmap

Plot CNV heatmap

## Usage

``` r
cnv_heatmap(
  segs,
  var = "group",
  label_group = TRUE,
  legend = TRUE,
  exclude_gap = TRUE,
  genome = "hg38"
)
```

## Arguments

- segs:

  dataframe Segments to plot. Need columns "seg_start", "seg_end",
  "cnv_state"

- var:

  character Column to facet by

- label_group:

  logical Label the groups

- legend:

  logical Display the legend

- exclude_gap:

  logical Whether to mark gap regions

- genome:

  character Genome build, either 'hg38' or 'hg19'

## Value

ggplot Heatmap of CNVs along the genome

## Examples

``` r
p = cnv_heatmap(segs_example)
#> Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the numbat package.
#>   Please report the issue to the authors.
```

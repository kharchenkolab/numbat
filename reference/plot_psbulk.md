# Plot a pseudobulk HMM profile

Plot a pseudobulk HMM profile

## Usage

``` r
plot_psbulk(
  bulk,
  use_pos = TRUE,
  allele_only = FALSE,
  min_LLR = 5,
  min_depth = 8,
  exp_limit = 2,
  phi_mle = TRUE,
  theta_roll = FALSE,
  dot_size = 0.8,
  dot_alpha = 0.5,
  legend = TRUE,
  exclude_gap = TRUE,
  genome = "hg38",
  text_size = 10,
  raster = FALSE
)
```

## Arguments

- bulk:

  dataframe Pseudobulk profile

- use_pos:

  logical Use marker position instead of index as x coordinate

- allele_only:

  logical Only plot alleles

- min_LLR:

  numeric LLR threshold for event filtering

- min_depth:

  numeric Minimum coverage depth for a SNP to be plotted

- exp_limit:

  numeric Expression logFC axis limit

- phi_mle:

  logical Whether to plot estimates of segmental expression fold change

- theta_roll:

  logical Whether to plot rolling estimates of allele imbalance

- dot_size:

  numeric Size of marker dots

- dot_alpha:

  numeric Transparency of the marker dots

- legend:

  logical Whether to show legend

- exclude_gap:

  logical Whether to mark gap regions and centromeres

- genome:

  character Genome build, either 'hg38' or 'hg19'

- text_size:

  numeric Size of text in the plot

- raster:

  logical Whether to raster images

## Value

ggplot Plot of pseudobulk HMM profile

## Examples

``` r
p = plot_psbulk(bulk_example)
#> Warning: Arguments in `...` must be used.
#> ✖ Problematic argument:
#> • na.rm = TRUE
#> ℹ Did you misspell an argument name?
```

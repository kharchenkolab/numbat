# Find the common diploid region in a group of pseudobulks

Find the common diploid region in a group of pseudobulks

## Usage

``` r
find_common_diploid(
  bulks,
  grouping = "clique",
  gamma = 20,
  theta_min = 0.08,
  t = 1e-05,
  fc_min = 2^0.25,
  alpha = 1e-04,
  min_genes = 10,
  ncores = 1,
  debug = FALSE,
  verbose = TRUE
)
```

## Arguments

- bulks:

  dataframe Pseudobulk profiles (differentiated by "sample" column)

- grouping:

  logical Whether to use cliques or components in the graph to find
  dipoid cluster

- gamma:

  numeric Dispersion parameter for the Beta-Binomial allele model

- theta_min:

  numeric Minimum imbalance threshold

- t:

  numeric Transition probability

- fc_min:

  numeric Minimum fold change to call quadruploid cluster

- alpha:

  numeric FDR cut-off for q values to determine edges

- ncores:

  integer Number of cores to use

## Value

list Ploidy information

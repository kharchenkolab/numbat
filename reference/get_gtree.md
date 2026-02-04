# Get a tidygraph tree with simplified mutational history.

Specify either max_cost or n_cut. max_cost works similarly as h and
n_cut works similarly as k in stats::cutree. The top-level normal
diploid clone is always included.

## Usage

``` r
get_gtree(tree, P, n_cut = 0, max_cost = 0)
```

## Arguments

- tree:

  phylo Single-cell phylogenetic tree

- P:

  matrix Genotype probability matrix

- n_cut:

  integer Number of cuts on the phylogeny to define subclones

- max_cost:

  numeric Likelihood threshold to collapse internal branches

## Value

tbl_graph Phylogeny annotated with branch lengths and mutation events

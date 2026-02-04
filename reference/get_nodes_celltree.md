# Get the internal nodes of a dendrogram and the leafs in each subtree

Get the internal nodes of a dendrogram and the leafs in each subtree

## Usage

``` r
get_nodes_celltree(hc, clusters)
```

## Arguments

- hc:

  hclust Clustering results

- clusters:

  named vector Cutree output specifying the terminal clusters

## Value

list Interal node subtrees with leaf memberships

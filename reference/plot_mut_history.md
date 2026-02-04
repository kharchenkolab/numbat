# Plot mutational history

Plot mutational history

## Usage

``` r
plot_mut_history(
  G,
  clone_post = NULL,
  edge_label_size = 4,
  node_label_size = 6,
  node_size = 10,
  arrow_size = 2,
  show_clone_size = TRUE,
  show_distance = TRUE,
  legend = TRUE,
  edge_label = TRUE,
  node_label = TRUE,
  horizontal = TRUE,
  pal = NULL
)
```

## Arguments

- G:

  igraph Mutation history graph

- clone_post:

  dataframe Clone assignment posteriors

- edge_label_size:

  numeric Size of edge label

- node_label_size:

  numeric Size of node label

- node_size:

  numeric Size of nodes

- arrow_size:

  numeric Size of arrows

- show_clone_size:

  logical Whether to show clone size

- show_distance:

  logical Whether to show evolutionary distance between clones

- legend:

  logical Whether to show legend

- edge_label:

  logical Whether to label edges

- node_label:

  logical Whether to label nodes

- horizontal:

  logical Whether to use horizontal layout

- pal:

  named vector Node colors

## Value

ggplot object

## Examples

``` r
p = plot_mut_history(mut_graph_example)
```

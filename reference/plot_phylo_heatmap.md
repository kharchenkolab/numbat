# Plot single-cell CNV calls along with the clonal phylogeny

Plot single-cell CNV calls along with the clonal phylogeny

## Usage

``` r
plot_phylo_heatmap(
  gtree,
  joint_post,
  segs_consensus,
  clone_post = NULL,
  p_min = 0.9,
  annot = NULL,
  pal_annot = NULL,
  annot_title = "Annotation",
  annot_scale = NULL,
  clone_dict = NULL,
  clone_bar = TRUE,
  clone_stack = TRUE,
  pal_clone = NULL,
  clone_title = "Genotype",
  clone_legend = TRUE,
  line_width = 0.1,
  tree_height = 1,
  branch_width = 0.2,
  tip_length = 0.2,
  branch_length = TRUE,
  annot_bar_width = 0.25,
  clone_bar_width = 0.25,
  bar_label_size = 7,
  tvn_line = TRUE,
  clone_line = FALSE,
  exclude_gap = FALSE,
  root_edge = TRUE,
  raster = FALSE,
  show_phylo = TRUE
)
```

## Arguments

- gtree:

  tbl_graph The single-cell phylogeny

- joint_post:

  dataframe Joint single cell CNV posteriors

- segs_consensus:

  datatframe Consensus segment dataframe

- clone_post:

  dataframe Clone assignment posteriors

- p_min:

  numeric Probability threshold to display CNV calls

- annot:

  dataframe Cell annotations, dataframe with 'cell' and additional
  annotation columns

- pal_annot:

  named vector Colors for cell annotations

- annot_title:

  character Legend title for the annotation bar

- annot_scale:

  ggplot scale Color scale for the annotation bar

- clone_dict:

  named vector Clone annotations, mapping from cell name to clones

- clone_bar:

  logical Whether to display clone bar plot

- clone_stack:

  character Whether to plot clone assignment probabilities as stacked
  bar

- pal_clone:

  named vector Clone colors

- clone_title:

  character Legend title for the clone bar

- clone_legend:

  logical Whether to display the clone legend

- line_width:

  numeric Line width for CNV heatmap

- tree_height:

  numeric Relative height of the phylogeny plot

- branch_width:

  numeric Line width in the phylogeny

- tip_length:

  numeric Length of tips in the phylogeny

- branch_length:

  logical Whether to use branch length in the phylogeny

- annot_bar_width:

  numeric Width of annotation bar

- clone_bar_width:

  numeric Width of clone genotype bar

- bar_label_size:

  numeric Size of sidebar text labels

- tvn_line:

  logical Whether to draw line separating tumor and normal cells

- clone_line:

  logical Whether to display borders for clones in the heatmap

- exclude_gap:

  logical Whether to mark gap regions

- root_edge:

  logical Whether to plot root edge

- raster:

  logical Whether to raster images

- show_phylo:

  logical Whether to display phylogeny on y axis

## Value

ggplot panel

## Examples

``` r
p = plot_phylo_heatmap(
    gtree = phylogeny_example,
    joint_post = joint_post_example,
    segs_consensus = segs_example)
```

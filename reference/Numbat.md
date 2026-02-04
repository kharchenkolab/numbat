# Numbat R6 class

Used to allow users to plot results

## Value

a new 'Numbat' object

## Public fields

- `label`:

  character Sample name

- `gtf`:

  dataframe Transcript annotation

- `joint_post`:

  dataframe Joint posterior

- `exp_post`:

  dataframe Expression posterior

- `allele_post`:

  dataframe Allele posetrior

- `bulk_subtrees`:

  dataframe Bulk profiles of lineage subtrees

- `bulk_clones`:

  dataframe Bulk profiles of clones

- `segs_consensus`:

  dataframe Consensus segments

- `tree_post`:

  list Tree posterior

- `mut_graph`:

  igraph Mutation history graph

- `gtree`:

  tbl_graph Single-cell phylogeny

- `clone_post`:

  dataframe Clone posteriors

- `gexp_roll_wide`:

  matrix Smoothed expression of single cells

- `P`:

  matrix Genotype probability matrix

- `treeML`:

  matrix Maximum likelihood tree as phylo object

- `hc`:

  hclust Initial hierarchical clustering

## Methods

### Public methods

- [`Numbat$new()`](#method-Numbat-new)

- [`Numbat$plot_phylo_heatmap()`](#method-Numbat-plot_phylo_heatmap)

- [`Numbat$plot_exp_roll()`](#method-Numbat-plot_exp_roll)

- [`Numbat$plot_mut_history()`](#method-Numbat-plot_mut_history)

- [`Numbat$plot_sc_tree()`](#method-Numbat-plot_sc_tree)

- [`Numbat$plot_consensus()`](#method-Numbat-plot_consensus)

- [`Numbat$plot_clone_profile()`](#method-Numbat-plot_clone_profile)

- [`Numbat$cutree()`](#method-Numbat-cutree)

- [`Numbat$clone()`](#method-Numbat-clone)

------------------------------------------------------------------------

### Method `new()`

initialize Numbat class

#### Usage

    Numbat$new(out_dir, i = 2, gtf = gtf_hg38, verbose = TRUE)

#### Arguments

- `out_dir`:

  character string Output directory

- `i`:

  integer Get results from which iteration (default=2)

- `gtf`:

  dataframe Transcript gtf (default=gtf_hg38)

- `verbose`:

  logical Whether to output verbose results (default=TRUE)

#### Returns

a new 'Numbat' object

------------------------------------------------------------------------

### Method [`plot_phylo_heatmap()`](https://kharchenkolab.github.io/numbat/reference/plot_phylo_heatmap.md)

Plot the single-cell CNV calls in a heatmap and the corresponding
phylogeny

#### Usage

    Numbat$plot_phylo_heatmap(...)

#### Arguments

- `...`:

  additional parameters passed to plot_phylo_heatmap()

------------------------------------------------------------------------

### Method [`plot_exp_roll()`](https://kharchenkolab.github.io/numbat/reference/plot_exp_roll.md)

Plot window-smoothed expression profiles

#### Usage

    Numbat$plot_exp_roll(k = 3, n_sample = 300, ...)

#### Arguments

- `k`:

  integer Number of clusters

- `n_sample`:

  integer Number of cells to subsample

- `...`:

  additional parameters passed to plot_exp_roll()

------------------------------------------------------------------------

### Method [`plot_mut_history()`](https://kharchenkolab.github.io/numbat/reference/plot_mut_history.md)

Plot the mutation history of the tumor

#### Usage

    Numbat$plot_mut_history(...)

#### Arguments

- `...`:

  additional parameters passed to plot_mut_history()

------------------------------------------------------------------------

### Method [`plot_sc_tree()`](https://kharchenkolab.github.io/numbat/reference/plot_sc_tree.md)

Plot the single cell phylogeny

#### Usage

    Numbat$plot_sc_tree(...)

#### Arguments

- `...`:

  additional parameters passed to plot_sc_tree()

------------------------------------------------------------------------

### Method [`plot_consensus()`](https://kharchenkolab.github.io/numbat/reference/plot_consensus.md)

Plot consensus segments

#### Usage

    Numbat$plot_consensus(...)

#### Arguments

- `...`:

  additional parameters passed to plot_sc_tree()

------------------------------------------------------------------------

### Method `plot_clone_profile()`

Plot clone cnv profiles

#### Usage

    Numbat$plot_clone_profile(...)

#### Arguments

- `...`:

  additional parameters passed to plot_clone_profile()

------------------------------------------------------------------------

### Method [`cutree()`](https://rdrr.io/r/stats/cutree.html)

Re-define subclones on the phylogeny.

#### Usage

    Numbat$cutree(max_cost = 0, n_cut = 0)

#### Arguments

- `max_cost`:

  numeric Likelihood threshold to collapse internal branches

- `n_cut`:

  integer Number of cuts on the phylogeny to define subclones

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Numbat$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

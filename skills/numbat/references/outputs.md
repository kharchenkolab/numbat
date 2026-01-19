# Outputs and Interpretation

## Key output files
- `bulk_subtrees_{i}.tsv.gz`: subtree pseudobulk profiles
- `segs_consensus_{i}.tsv.gz`: consensus CNV segments
- `bulk_clones_{i}.tsv.gz`, `bulk_clones_final.tsv.gz`: clone-level pseudobulk profiles
- `exp_post_{i}.tsv`: expression-based posteriors per segment per cell
- `allele_post_{i}.tsv`: allele-based posteriors per segment per cell
- `joint_post_{i}.tsv`: joint posteriors per segment per cell
- `clone_post_{i}.tsv`: clone assignment and tumor-vs-normal posteriors
- `panel_{i}.png`: integrated phylogeny + CNV landscape
- `exp_roll_clust.png`: expression-based clustering view
- `tree_final_{i}.rds`: final lineage tree
- `log.txt`: run log

## Numbat object helpers
- Instantiate: `nb <- Numbat$new(out_dir)`
- Useful plots: `plot_phylo_heatmap()`, `plot_consensus()`, `plot_bulks()`, `plot_sc_tree()`, `plot_mut_history()`
- Key tables: `nb$joint_post`, `nb$clone_post`, `nb$bulk_clones`, `nb$segs_consensus`
- Refine clones: `nb$cutree(n_cut=...)` or `nb$cutree(max_cost=...)`

## Interpretation notes
- `clone_post$p_cnv` is the tumor-vs-normal probability.
- `joint_post$p_cnv` is event-level posterior per segment.
- Use `min_LLR` and `max_entropy` thresholds when deciding which events to emphasize.

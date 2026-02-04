# Package index

## Preparing data

- [`aggregate_counts()`](https://kharchenkolab.github.io/numbat/reference/aggregate_counts.md)
  : Utility function to make reference gene expression profiles

## CNV detection

- [`get_bulk()`](https://kharchenkolab.github.io/numbat/reference/get_bulk.md)
  : Aggregate single-cell data into combined bulk expression and allele
  profile
- [`analyze_bulk()`](https://kharchenkolab.github.io/numbat/reference/analyze_bulk.md)
  : Call CNVs in a pseudobulk profile using the Numbat joint HMM
- [`detect_clonal_loh()`](https://kharchenkolab.github.io/numbat/reference/detect_clonal_loh.md)
  : Call clonal LOH using SNP density. Rcommended for cell lines or
  tumor samples with no normal cells.

## Single-cell clonal decomposition

- [`run_numbat()`](https://kharchenkolab.github.io/numbat/reference/run_numbat.md)
  : Run workflow to decompose tumor subclones

## Visualizing results

- [`Numbat`](https://kharchenkolab.github.io/numbat/reference/Numbat.md)
  : Numbat R6 class
- [`get_gtree()`](https://kharchenkolab.github.io/numbat/reference/get_gtree.md)
  : Get a tidygraph tree with simplified mutational history.
- [`plot_bulks()`](https://kharchenkolab.github.io/numbat/reference/plot_bulks.md)
  : Plot a group of pseudobulk HMM profiles
- [`plot_consensus()`](https://kharchenkolab.github.io/numbat/reference/plot_consensus.md)
  : Plot consensus CNVs
- [`plot_exp_roll()`](https://kharchenkolab.github.io/numbat/reference/plot_exp_roll.md)
  : Plot single-cell smoothed expression magnitude heatmap
- [`plot_mut_history()`](https://kharchenkolab.github.io/numbat/reference/plot_mut_history.md)
  : Plot mutational history
- [`plot_phylo_heatmap()`](https://kharchenkolab.github.io/numbat/reference/plot_phylo_heatmap.md)
  : Plot single-cell CNV calls along with the clonal phylogeny
- [`plot_psbulk()`](https://kharchenkolab.github.io/numbat/reference/plot_psbulk.md)
  : Plot a pseudobulk HMM profile
- [`plot_sc_tree()`](https://kharchenkolab.github.io/numbat/reference/plot_sc_tree.md)
  : Plot single-cell smoothed expression magnitude heatmap
- [`cnv_heatmap()`](https://kharchenkolab.github.io/numbat/reference/cnv_heatmap.md)
  : Plot CNV heatmap

## Example data

- [`df_allele_example`](https://kharchenkolab.github.io/numbat/reference/df_allele_example.md)
  : example allele count dataframe
- [`count_mat_example`](https://kharchenkolab.github.io/numbat/reference/count_mat_example.md)
  : example gene expression count matrix
- [`ref_hca`](https://kharchenkolab.github.io/numbat/reference/ref_hca.md)
  : reference expression magnitudes from HCA
- [`ref_hca_counts`](https://kharchenkolab.github.io/numbat/reference/ref_hca_counts.md)
  : reference expression counts from HCA
- [`bulk_example`](https://kharchenkolab.github.io/numbat/reference/bulk_example.md)
  : example pseudobulk dataframe
- [`annot_ref`](https://kharchenkolab.github.io/numbat/reference/annot_ref.md)
  : example reference cell annotation
- [`count_mat_ref`](https://kharchenkolab.github.io/numbat/reference/count_mat_ref.md)
  : example reference count matrix
- [`gexp_roll_example`](https://kharchenkolab.github.io/numbat/reference/gexp_roll_example.md)
  : example smoothed gene expression dataframe
- [`hc_example`](https://kharchenkolab.github.io/numbat/reference/hc_example.md)
  : example hclust tree
- [`joint_post_example`](https://kharchenkolab.github.io/numbat/reference/joint_post_example.md)
  : example joint single-cell cnv posterior dataframe
- [`mut_graph_example`](https://kharchenkolab.github.io/numbat/reference/mut_graph_example.md)
  : example mutation graph
- [`phylogeny_example`](https://kharchenkolab.github.io/numbat/reference/phylogeny_example.md)
  : example single-cell phylogeny
- [`pre_likelihood_hmm`](https://kharchenkolab.github.io/numbat/reference/pre_likelihood_hmm.md)
  : HMM object for unit tests
- [`segs_example`](https://kharchenkolab.github.io/numbat/reference/segs_example.md)
  : example CNV segments dataframe

## Genome annotations

- [`acen_hg19`](https://kharchenkolab.github.io/numbat/reference/acen_hg19.md)
  : centromere regions (hg19)
- [`acen_hg38`](https://kharchenkolab.github.io/numbat/reference/acen_hg38.md)
  : centromere regions (hg38)
- [`chrom_sizes_hg19`](https://kharchenkolab.github.io/numbat/reference/chrom_sizes_hg19.md)
  : chromosome sizes (hg19)
- [`chrom_sizes_hg38`](https://kharchenkolab.github.io/numbat/reference/chrom_sizes_hg38.md)
  : chromosome sizes (hg38)
- [`gaps_hg19`](https://kharchenkolab.github.io/numbat/reference/gaps_hg19.md)
  : genome gap regions (hg19)
- [`gaps_hg38`](https://kharchenkolab.github.io/numbat/reference/gaps_hg38.md)
  : genome gap regions (hg38)
- [`gtf_hg19`](https://kharchenkolab.github.io/numbat/reference/gtf_hg19.md)
  : gene model (hg19)
- [`gtf_hg38`](https://kharchenkolab.github.io/numbat/reference/gtf_hg38.md)
  : gene model (hg38)
- [`gtf_mm10`](https://kharchenkolab.github.io/numbat/reference/gtf_mm10.md)
  : gene model (mm10)
- [`vcf_meta`](https://kharchenkolab.github.io/numbat/reference/vcf_meta.md)
  : example VCF header

## Misc

- [`annotate_genes()`](https://kharchenkolab.github.io/numbat/reference/annotate_genes.md)
  : Annotate genes on allele dataframe
- [`upgma()`](https://kharchenkolab.github.io/numbat/reference/upgma.md)
  : UPGMA and WPGMA clustering

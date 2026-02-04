# Call CNVs in a pseudobulk profile using the Numbat joint HMM

Call CNVs in a pseudobulk profile using the Numbat joint HMM

## Usage

``` r
analyze_bulk(
  bulk,
  t = 1e-05,
  gamma = 20,
  theta_min = 0.08,
  logphi_min = 0.25,
  nu = 1,
  min_genes = 10,
  exp_only = FALSE,
  allele_only = FALSE,
  bal_cnv = TRUE,
  retest = TRUE,
  find_diploid = TRUE,
  diploid_chroms = NULL,
  classify_allele = FALSE,
  run_hmm = TRUE,
  prior = NULL,
  exclude_neu = TRUE,
  phasing = TRUE,
  verbose = TRUE
)
```

## Arguments

- bulk:

  dataframe Pesudobulk profile

- t:

  numeric Transition probability

- gamma:

  numeric Dispersion parameter for the Beta-Binomial allele model

- theta_min:

  numeric Minimum imbalance threshold

- logphi_min:

  numeric Minimum log expression deviation threshold

- nu:

  numeric Phase switch rate

- min_genes:

  integer Minimum number of genes to call an event

- exp_only:

  logical Whether to run expression-only HMM

- allele_only:

  logical Whether to run allele-only HMM

- bal_cnv:

  logical Whether to call balanced amplifications/deletions

- retest:

  logical Whether to retest CNVs after Viterbi decoding

- find_diploid:

  logical Whether to run diploid region identification routine

- diploid_chroms:

  character vector User-given chromosomes that are known to be in
  diploid state

- classify_allele:

  logical Whether to only classify allele (internal use only)

- run_hmm:

  logical Whether to run HMM (internal use only)

- prior:

  numeric vector Prior probabilities of states (internal use only)

- exclude_neu:

  logical Whether to exclude neutral segments from retesting (internal
  use only)

- phasing:

  logical Whether to use phasing information (internal use only)

- verbose:

  logical Verbosity

## Value

a pseudobulk profile dataframe with called CNV information

## Examples

``` r
bulk_analyzed = analyze_bulk(bulk_example, t = 1e-5, find_diploid = FALSE, retest = FALSE)
```

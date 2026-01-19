---
name: numbat
description: Guide to run Numbat analysis for allele-specific CNV calling from single-cell or spatial transcriptomics data. Use when users need to prepare inputs, run pileup and phasing, execute run_numbat, tune key parameters, or interpret Numbat outputs in this repo.
---

# Numbat Analysis Skill

## Scope
- Use this skill to help users run Numbat on scRNA-seq, spatial, or multiome data using this repo.
- Must-use sources: docs site in `docs/`, GitHub issues (triage patterns), codebase (scripts and functions), and the primary paper.
- Prefer repo docs: `README.md`, `vignettes/numbat.Rmd`, `vignettes/numbat-multiome.Rmd`, `vignettes/spatial-rna.Rmd`, `vignettes/mouse.Rmd`, `vignettes/results.Rmd`.
- Track sources consulted in `skills/numbat/references/sources.md`.

## Quick Start
1. Confirm prerequisites and references.
2. Generate allele counts and phased SNPs with `inst/bin/pileup_and_phase.R`.
3. Provide gene x cell count matrix and expression reference.
4. Run `run_numbat` and review outputs.

## Prerequisites
- External tools in `PATH`: `cellsnp-lite`, `eagle2`, `samtools`.
- SNP VCF and 1000G phasing panel; use the URLs in `vignettes/numbat.Rmd` if needed.
- For container runs: `docker run -v /work:/mnt/mydata -it pkharchenkolab/numbat-rbase:latest /bin/bash`.

## Authoritative Sources
- Docs site (repo): use `docs/index.html` and `docs/articles/*.html` for canonical guidance.
- GitHub issues: use for common errors, edge cases, and user-facing troubleshooting patterns.
- Codebase: prefer `inst/bin/*.R` for CLI interfaces and `R/` for defaults and core logic.
- Paper: use to justify methodological assumptions and explain outputs at a high level. Full text in `skills/numbat/references/paper-full-article.html` and `skills/numbat/references/paper-full-clean.txt`, summary in `skills/numbat/references/paper-summary.md`.

## Prepare Inputs
1. Allele data: run `inst/bin/pileup_and_phase.R` (details in `skills/numbat/references/cli.md`).
   - Output: `{sample}_allele_counts.tsv.gz` in `outdir` plus `phasing/` and `pileup/`.
2. Expression matrix: gene x cell integer UMI counts from upstream pipeline (e.g. CellRanger).
3. Expression reference:
   - Use `ref_hca` for quick start, or compute with `aggregate_counts(count_mat, cell_annot)`.
   - The reference regresses out cell-type expression and reduces noise in CNV detection.

## Run Numbat (Core)
Use `run_numbat` with the count matrix, expression reference, and allele dataframe.

Example:
```r
library(numbat)

out <- run_numbat(
  count_mat,
  ref_hca,
  df_allele,
  genome = "hg38",
  t = 1e-5,
  ncores = 4,
  plot = TRUE,
  out_dir = "./numbat_out"
)
```

## Function Reference (Key Signatures)
- `run_numbat(count_mat, lambdas_ref, df_allele, genome = "hg38", out_dir = ..., max_iter = 2, t = 1e-5, gamma = 20, min_LLR = 5, max_entropy = 0.5, init_k = 3, min_cells = 50, tau = 0.3, nu = 1, n_cut = 0, plot_min_depth = 8, multi_allelic = TRUE, ...)`
- `aggregate_counts(count_mat, annot, normalized = TRUE)` where `annot` has `cell` and `group`.
- `detect_clonal_loh(bulk, t = 1e-5, snp_rate_loh = 5, min_depth = 0)` for clonal LOH in high-purity samples.

## Key Parameters to Tune
- `t`: HMM transition probability. Lower for complex CNV landscapes, higher to reduce false positives.
- `gamma`: allele overdispersion. Default 20 for 10x; use smaller (e.g. 5) for noisier protocols.
- `min_cells`: minimum cells per pseudobulk HMM; increase if per-cell coverage is sparse.
- `multi_allelic`: enable multiallelic CNV calling.
- `min_LLR`, `max_entropy`: CNV filtering before phylogeny.
- `n_cut`, `max_cost`, `tau`: control clone number and phylogeny simplification.
- `init_k`, `max_iter`, `check_convergence`: iterative optimization settings.
- `call_clonal_loh`: identify clonal LOH in high-purity samples.
- `segs_consensus_fix`: fix CNV boundaries and states using external CNV calls.
  
## Phasing Options
- Default phasing uses 1000G panel; larger panels can improve power (gnomAD HGDP + 1KG, TOPMed).
- If you have DNA-derived phased VCFs, use those to improve SNP density and phasing accuracy.

## Common Variants
- Multiome: see `skills/numbat/references/cli.md` and `vignettes/numbat-multiome.Rmd`.
- Spatial RNA: see `skills/numbat/references/faq.md` and `vignettes/spatial-rna.Rmd`.
- Mouse: see `skills/numbat/references/faq.md` and `vignettes/mouse.Rmd`.

## Outputs and Interpretation
Review plots/tables in `out_dir` and use `skills/numbat/references/outputs.md` for file/column descriptions and interpretation patterns. For detailed walkthroughs, use `vignettes/results.Rmd` and `vignettes/descriptions.Rmd`.

## Troubleshooting via Issues
- Map the user's error to similar GitHub issues before suggesting changes.
- Always ask for: full error, parameters used, genome build, and minimal repro if possible.
Track recurring failures in `skills/numbat/references/issues.md` as the skill grows.

## Troubleshooting Guidance
Keep the core questions in mind and use `skills/numbat/references/faq.md` for detailed patterns and fixes.

## Common Misconceptions

- **`p_cnv` is for visualization, not calling**: The `p_cnv` threshold controls which CNV events are displayed in plots (e.g., `plot_phylo_heatmap` uses `p_cnv > 0.9` by default). It does not affect the CNV calling itself - use `min_LLR` and `max_entropy` for that.
- **Clone 1 is the diploid baseline**: Clone 1 representing "normal" cells is expected behavior. Numbat uses diploid cells as the CNV calling baseline. Tumor cells without detectable CNVs may also appear in Clone 1.
- **High MSE requires fixing the reference, not thresholds**: When MSE is high, the core problem is reference-sample mismatch. Adjusting downstream thresholds won't fix this - instead fix the reference (cell type matching, use `aggregate_counts()`, exclude tumor cells from reference).

## Multi-Sample Analysis Rules

- **Different individuals**: Must always run separately (different germline genotypes). Compare results post-hoc using `joint_post` outputs.
- **Same individual**: Run `pileup_and_phase.R` jointly with comma-delimited BAM list, concatenate `allele_df` and `count_mat`, then run single `run_numbat()` for integrated phylogeny.

## Making CNV Calling More Stringent

To reduce false positive CNV calls:
- Increase `min_LLR` (default 5) - stricter pseudobulk CNV filtering
- Decrease `max_entropy` (default 0.5) - stricter single-cell CNV filtering

To increase sensitivity for subclonal CNVs:
- Increase `init_k` - higher clustering resolution
- Increase `max_entropy` - more permissive filtering

## Useful Post-Processing Features

- **Adjust clone number**: Use `nb$cutree(n_cut = N)` to redefine clones after running (N cuts → N+1 clones)
- **Visualize sample metadata**: `plot_phylo_heatmap()` supports annotation bars via the `annot` parameter for overlaying conditions, timepoints, etc. See [function reference](https://kharchenkolab.github.io/numbat/reference/plot_phylo_heatmap.html)

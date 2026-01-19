# Numbat GitHub Issues: Quick Reference

Condensed maintainer guidance from closed issues. Link format: `#NNN` → github.com/kharchenkolab/numbat/issues/NNN

---

## Platform & Data Compatibility

### Unsupported Data Types
- **Probe-based/FFPE chemistry** (10x Flex, Visium FFPE): Not supported - no SNP info captured. Use InferCNV/CopyKAT instead. (#48, #95, #113, #163, #167)
- **Inbred mouse strains** (pure C57BL/6J): Lacks heterozygous SNPs required for allele analysis. (#96, #122, #198)
- **Expression-only analysis**: Full workflow doesn't support it; only pseudobulk HMM via internal hahmmr function. (#187)

### Supported Data Types
- **ATAC-seq**: Same workflow as before. (#244)
- **Mouse hybrid crosses**: Works if parental haplotypes known. See [mouse tutorial](https://kharchenkolab.github.io/numbat/articles/mouse.html). (#15, #174)
- **Smart-seq**: Use `--smartseq` mode with bam file list. (#22, #60, #61)
- **Visium fresh frozen**: Works like scRNA-seq. (#113)
- **Bulk RNA-seq**: Use `get_bulk` → `analyze_bulk`, not `run_numbat`. (#53)
- **Other species with phased genome**: See mouse tutorial; may need to hack code for custom GTF. (#174)

---

## Input Data Requirements

### Count Matrix
- Use **raw counts**, not normalized. (#205, #206)
- Include only cells from ONE individual per run. (#111)
- Filter out reference cells before running. (#31)

### Reference Expression
- Must be **gene × cell-type** matrix from `numbat::aggregate_counts`. (#97, #98)
- Input should be raw counts + cell type annotation. (#29)
- Don't include tumor cells in diploid reference. (#62)
- External ref (e.g., `ref_hca`) works if cell types match. (#31)
- Match intron settings between reference and target. (#155)

### Allele Data
- Include only **heterozygous SNPs**. (#122)
- Check column data types match example (`df_allele_ATC2`), especially `DP` column. (#63)
- Ensure chr prefix matches between BAM and SNP file. (#167)

---

## Multi-Sample Analysis

### Same Individual
- Run `pileup_and_phase.R` jointly with comma-delimited BAM list. (#32, #129, #130, #176)
- Concatenate allele_df and count_mat, then feed to single `run_numbat`. (#130)
- No barcode renaming needed; no normalization needed. (#176, #206)

### Different Individuals
- Always run separately (different genotypes). (#111)
- Compare results post-hoc using `joint_post` CNV calls. (#111)

---

## Common Errors & Solutions

### "No CNV remains after filtering"
| Filter | Solution |
|--------|----------|
| LLR in pseudobulks | Check pseudobulk profiles - may have no CNVs; try `init_k` higher for subclonal detection. (#91, #199) |
| Entropy in single cells | Increase `max_entropy` (0.5→0.8 or higher). (#74, #141) |
| General | Sample may genuinely lack CNVs; balanced BAF = no CNV. (#134) |

### "non-finite value supplied by optim" / curly brace error
- Usually: reference cells included in analysis. Exclude them. (#110, #234)

### Memory Issues
- `mclapply` duplicates workspace: reduce `ncores`. (#70)
- Use `ncores_nni` separately for phylogeny step. (#70)
- v0.1.3+ uses RcppParallel for phylogeny (10-20x speedup, constant memory). (#34)

### "C stack usage too close to limit"
- `ape::ladderize` bug when only 1 event passes filter. Raise `max_entropy`. (#30, #79)

### Missing Files / Empty Output
- Output folder created at start; check incrementally. (#55)
- Check gene name overlap between ref and count_mat. (#73)
- Check for 0-coverage cells. (#80)

---

## Output Interpretation

### CNV States
- `cnv_state_post` is final assignment (use this). (#124)
- Rough integer mapping: del=(1,0), neu=(1,1), amp=(2,1), bamp=(2,2). (#100)
- Compare to WGS using `cnv_state` column (states, not absolute CN). (#194)

### Clone Naming
- `clone_post` letters (1a, 1b): chr + segment order (a=first, b=second). (#235)
- Segment names can be any unique identifier. (#152)

### Probabilities
- `p_cnv > 0.5-1` reasonable for single-cell CNV comparison. (#94)
- `p_*` columns in `joint_post` = single-cell posteriors. (#146)
- `prior_*` = consensus from pseudobulk HMM. (#146)
- `LLR` for pseudobulk; `p_cnv` for single cells. (#94, #194)

### Plots
- Top track (expression DE) is noisy; try y-axis scale -5 to 5. (#215)
- Circle vs square: low vs high allelic imbalance. (#215)
- logFC ≈ logR; pHF ≈ phased BAF. (#69)

---

## Parameter Tuning

### When to Adjust
| Symptom | Parameter | Action |
|---------|-----------|--------|
| Noisy expression | Custom reference | Use matched cell-type reference |
| Focal CNVs filtered | `max_entropy` | Increase (0.8+) |
| Subclonal CNVs missed | `init_k` | Increase clustering resolution |
| High ploidy / WGD | `diploid_chroms` | Supply quadruploid regions as baseline (#144) |

### Using External CNV Calls
- v1.3.0+: Use `segs_consensus_fix` to supply pre-defined segments. (#45, #102, #108)
- Include diploid regions, not just CNV regions. (#108)

---

## Post-Processing

### Adjust Clone Number
- Use `nb$cutree(n_cut = N)` for N+1 clones. (#103, #129)

### Custom Iteration
- Set `i` parameter in `Numbat$new` for specific iteration results. (#117)

### Recreate Clone Pseudobulks
```r
nb$cutree(n_cut = 1)
clones = nb$clone_post %>% split(.$clone_opt) %>%
  purrr::map(function(c){list(...)})
# See #220 for full code
```

### Custom Cell Order in Heatmap
- Provide custom tree via `gtree` argument (tbl_graph). (#192)

---

## Data Downloads

### 1000G Reference Panel
- URL: `http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip` (#227)
- Alternative: Copy from Docker image `/data/1000G_hg38`. (#101, #109)

### Eagle Genetic Map
- Copy from Docker or run inside Docker. (#109)

---

## Key Documentation Links
- Output descriptions: https://kharchenkolab.github.io/numbat/articles/descriptions.html
- Refine subclones: https://kharchenkolab.github.io/numbat/articles/results.html#refine-subclones-on-the-phylogeny
- Using existing CNV calls: https://kharchenkolab.github.io/numbat/articles/numbat.html#using-existing-cnv-calls
- Pre-phased SNP profiles: https://kharchenkolab.github.io/numbat/articles/numbat.html#using-prephased-snp-profiles

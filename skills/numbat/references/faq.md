# FAQ and Troubleshooting Patterns

## Expression reference and noise
Q: What is the expression reference and why does it matter?
A: It regresses out cell-type expression differences unrelated to CNV. A closely matched normal reference reduces logFC noise and improves CNV detection.

## Expected CNV not detected
Potential causes and checks:
1) Not called in pseudobulk HMM
   - Check `bulk_subtrees_*.png` and `exp_roll_clust.png`.
   - Try a better reference or adjust `t`, `gamma`, or `init_k`.
2) Called in pseudobulk but filtered
   - Lower `min_LLR` to keep weak events.
3) Present in consensus but missing from phylogeny
   - Increase `max_entropy`; check `joint_post_*` for `avg_entropy`.

## Spatial transcriptomics guidance
- Spots may be mixed; interpret `clone_post$p_cnv` as a probability.
- Use Visium/normal tissue references when possible.
- Increase `max_entropy` (e.g. 0.8) for sparse allele coverage.

## Mouse data caveats
- Best for F1 hybrids with known parental haplotypes.
- Use `genome = "mm10"` and `nu = 0` when phasing is perfect.
- If admixture is not 50/50, filter to het SNPs and keep default `nu`.

## Common preflight checks
- Genome build consistency across VCF, GTF, and `--gmap`.
- Verify barcodes and sample names align with BAMs.
- Ensure `paneldir` has `chr{1..22}.genotypes.bcf`.
- Confirm `cellsnp-lite`, `eagle2`, and `samtools` are in `PATH`.

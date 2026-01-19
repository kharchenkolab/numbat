# CLI and Input Preparation

## pileup_and_phase.R
Purpose: run SNP pileup with cellsnp-lite, phase with Eagle2, and generate per-sample allele count tables.

Required args:
- `--bams`, `--snpvcf`, `--outdir`, `--paneldir`, `--gmap` always required.
- 10x mode: also `--samples` and `--barcodes`.

Modes:
- 10x (default): `--bams`, `--barcodes`, `--UMItag` (Auto/UB), `--cellTAG` (CB).
- Mixed 10x: multiple `--UMItag` values per sample (comma-delimited).
- SMART-seq: `--smartseq` with `--bams` and `--barcodes` as text files.
- Bulk: `--bulk` with `--bams` and `--samples` only (no barcodes).

UMI tags:
- `Auto` or `UB` for 10x scRNA-seq
- `None` for scATAC or bulk
- `XM` for Slide-seq

Outputs:
- `{outdir}/{sample}_allele_counts.tsv.gz`
- `{outdir}/pileup/`, `{outdir}/phasing/`
- Logs and scripts: `run_pileup.sh`, `pileup.log`, `run_phasing.sh`, `phasing.log`

Notes:
- Genome is inferred from `--gmap` path (hg19 vs hg38); phasing runs Eagle2 for chr1-22.
- If input is multiplexed, run one individual at a time with matching barcodes.

## Multiome prep scripts
- `get_gene_binned_intersections.R`: gene-to-bin mapping; supports `hg38`, `hg19`, `mm10` or a custom GTF.
- `get_binned_rna.R`: bins RNA counts (Seurat, matrix, or 10x h5).
- `get_binned_atac.R`: bins ATAC fragments into `GRanges` bins.
- `input_prep.R`: `binCnt()` to combine RNA+ATAC bins and `agg_refs()` to aggregate references.
- `run_numbat_multiome.R`: runs combined-bin, RNA-bin, or ATAC-bin analysis.

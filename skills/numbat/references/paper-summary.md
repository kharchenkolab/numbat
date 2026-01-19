# Paper Summary (Gao et al., Nat Biotechnol 2022)

Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC10289836/

## Short summary
Numbat performs haplotype-aware, allele-specific CNV inference from single-cell transcriptomes by integrating expression shifts, allelic imbalance, and population-based phasing. It detects CNVs and classifies tumor vs normal cells without paired DNA, and reconstructs tumor subclonal structure and phylogeny from scRNA-seq data. Key methodological assumptions include reliable SNP pileups, adequate expression references for regressing cell-type effects, and phased genotypes to link allele counts across loci.

## When to consult the full text
- Method details (model assumptions, likelihoods, priors)
- Benchmarking data and sensitivity analyses
- Limitations and failure modes

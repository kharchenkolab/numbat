# Genotyping main function

Genotyping main function

## Usage

``` r
genotype(label, samples, vcfs, outdir, het_only = FALSE, chr_prefix = TRUE)
```

## Arguments

- label:

  character Individual/sample label

- samples:

  vector Sample names

- vcfs:

  list of vcfR VCFs from cellsnp-lite pileup

- outdir:

  character Output directory

- het_only:

  logical Whether to only use heterozygous SNPs

- chr_prefix:

  logical Whether to add chr prefix

## Value

integer Status code

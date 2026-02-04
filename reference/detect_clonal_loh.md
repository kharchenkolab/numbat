# Call clonal LOH using SNP density. Rcommended for cell lines or tumor samples with no normal cells.

Call clonal LOH using SNP density. Rcommended for cell lines or tumor
samples with no normal cells.

## Usage

``` r
detect_clonal_loh(bulk, t = 1e-05, snp_rate_loh = 5, min_depth = 0)
```

## Arguments

- bulk:

  dataframe Pseudobulk profile

- t:

  numeric Transition probability

- snp_rate_loh:

  numeric The assumed SNP density in clonal LOH regions

- min_depth:

  integer Minimum coverage to filter SNPs

## Value

dataframe LOH segments

## Examples

``` r
segs_loh = detect_clonal_loh(bulk_example)
```

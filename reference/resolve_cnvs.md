# Get unique CNVs from set of segments

Get unique CNVs from set of segments

## Usage

``` r
resolve_cnvs(segs_all, min_overlap = 0.5, debug = FALSE)
```

## Arguments

- segs_all:

  dataframe CNV segments from multiple samples

- min_overlap:

  numeric scalar Minimum overlap fraction to determine count two events
  as as overlapping

## Value

dataframe Consensus CNV segments

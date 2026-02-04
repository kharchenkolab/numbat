# Extract consensus CNV segments

Extract consensus CNV segments

## Usage

``` r
get_segs_consensus(bulks, min_LLR = 5, min_overlap = 0.45, retest = TRUE)
```

## Arguments

- bulks:

  dataframe Pseudobulks

- min_LLR:

  numeric LLR threshold to filter CNVs

- min_overlap:

  numeric Minimum overlap fraction to determine count two events as as
  overlapping

## Value

dataframe Consensus segments

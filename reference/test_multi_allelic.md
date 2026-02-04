# test for multi-allelic CNVs

test for multi-allelic CNVs

## Usage

``` r
test_multi_allelic(bulks, segs_consensus, min_LLR = 5, p_min = 0.999)
```

## Arguments

- bulks:

  dataframe Pseudobulk profiles

- segs_consensus:

  dataframe Consensus segments

- min_LLR:

  numeric CNV LLR threshold to filter events

- p_min:

  numeric Probability threshold to call multi-allelic events

## Value

dataframe Consensus segments annotated with multi-allelic events

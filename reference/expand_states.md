# expand multi-allelic CNVs into separate entries in the single-cell posterior dataframe

expand multi-allelic CNVs into separate entries in the single-cell
posterior dataframe

## Usage

``` r
expand_states(sc_post, segs_consensus)
```

## Arguments

- sc_post:

  dataframe Single-cell posteriors

- segs_consensus:

  dataframe Consensus segments

## Value

dataframe Single-cell posteriors with multi-allelic CNVs split into
different entries

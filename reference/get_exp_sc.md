# get the single cell expression dataframe

get the single cell expression dataframe

## Usage

``` r
get_exp_sc(segs_consensus, count_mat, gtf, segs_loh = NULL)
```

## Arguments

- segs_consensus:

  dataframe Consensus segments

- count_mat:

  dgCMatrix gene expression count matrix

- gtf:

  dataframe Transcript gtf

## Value

dataframe single cell expression counts annotated with segments

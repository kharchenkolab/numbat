# Check and format the GTF input

This function validates the input GTF object, ensuring it is either a
properly structured \`GRanges\` or \`data.frame\` that matches the
expected Numbat GTF format. If a \`GRanges\` is provided, it will be
processed into a standard GTF-like \`data.frame\` with gene coordinates
and annotations.

## Usage

``` r
check_gtf_input(gtf_obj, reference_gtf_cols = colnames(gtf_hg38))
```

## Arguments

- gtf_obj:

  GRanges or data.frame; input GTF object

- reference_gtf_cols:

  character vector; expected column names (default from
  numbat::gtf_hg38)

## Value

data.frame Formatted GTF object with required structure

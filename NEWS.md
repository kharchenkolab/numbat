# numbat 1.0.5 - 11/27/2022

* Fixing bugs #65, #66, #67

* Retire dependency on `reshape2`

# numbat 1.0.4 - 11/20/2022

* Improving error handling and removing python dependency (`argparse`) in `pileup_and_phase.R`

* Allow plotting of mutliple annotations in `plot_phylo_heatmap` (thanks to @whtns)

* Adding diagnostic messages

# numbat 1.0.3 - 10/09/2022

* Fail gracefully when no CNV remains after `retest_bulks`

* Passing `gamma` parameter to `retest_bulks`

# numbat 1.0.2 - 09/07/2022

* Conform to CRAN guidelines

* Removed ATC2 examples from package data - users can download from lab server link instead

* New option to specify genome version (`genome = 'hg38' or 'hg19'`). Support plotting of centromeres and gap regions for hg19.

* Removed genetic maps from package data and they are no longer provided as input to `run_numbat`. Annotation of genetic distance is performed in `pileup_and_phase.R` script instead, using the genetic map included in Eagle2.

# numbat 1.0.0 - 08/12/2022

* Archival version for the paper

# numbat 0.1.3 - 07/02/2022

* Speed up of NNI using RcppParallel (#34). 10x faster and much more memory efficient (memory requirement is constant with respect to the number of threads).

* Speed up of expression single-cell testing using roptim. Approximately 2x speedup.

* New LLR metric for CNV filtering that is not inflated (default: 5).

* Only keep heterozygous SNPs in alelle dataframe to reduce memory usage

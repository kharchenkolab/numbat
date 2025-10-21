# numbat 1.5.1 - 10/03/2025

* Allow users to input their own GTF (#174)

* Allow CNV probabilities in `segs_consensus_fix` parameter (#207)

* Add multiome analysis helper functions and vignette

* Improve input validation and error handling

# numbat 1.4.2 - 09/18/2023

* Fix pseudobulk plotting legend (#182)

* Requirement for dplyr and tidyr versions (#189, #190)

* Fix Numbat$plot_exp_roll (#169)

* Fix CNV states reporting when `segs_loh` is provided (#183)

* Fix `n_states` reporting (#178)

* Improve error handling in `pileup_and_phase` (#179)

# numbat 1.4.0 - 02/23/2023

* Integration with hahmmr 

* Better input checking for pileup_and_phase

* Fix compatibility with igraph v2.0+ and tidygraph v1.3+ (#150)

* Fix multiallelic CNV state probability reporting (#146)

# numbat 1.3.3 - 08/15/2023

* Fix plotting issue #135

* Fix CRAN check compilation issues

# numbat 1.3.2 - 06/05/2023

* Adding better checks for input files

* Improve error handling (#122, #127)

# numbat 1.3.1 - 04/14/2023

* Fixing bug #68 - this may cause slight changes in the results for runs with `segs_loh`/`call_segs_loh` enabled.

# numbat 1.3.0 - 03/31/2023

* Allows users to supply existing CNV profiles (e.g. from bulk WGS/WES analysis) via `segs_consensus_fix` parameter

* Adding `call_clonal_loh` option to call clonal LOH events within `run_numbat`

* Fixing bug #81

* Fixing oversegmentation issue in `find_common_diploid` caused by `annot_segs`

# numbat 1.2.2 - 02/13/2023

* Introduce `n_cut` parameter to specify the number of clones to define from the phylogeny 

* Allows users to redefine subclones from the phylogeny via `nb$cutree`

# numbat 1.2.1 - 01/11/2023

* Fixing bugs #30, #79, #89

# numbat 1.2.0 - 12/26/2022

* Numbat now works for F1 hybrid mice! Check out the new tutorial under `Articles`.

* Fixing bugs #80, #82

* Offers stacked clone bars in `plot_phylo_heatmap`

# numbat 1.1.0 - 11/28/2022

* Externalize phylogeny module as separate package (`scistreer`)

* Prepare for new CRAN version

* Better CNV state legends for `plot_bulks`

# numbat 1.0.5 - 11/27/2022

* Fixing bugs #65, #66, #67

* Retire dependency on `reshape2`

# numbat 1.0.4 - 11/20/2022

* Improving error handling and removing python dependency (`argparse`) in `pileup_and_phase.R`

* Allows plotting of mutliple annotations in `plot_phylo_heatmap` (thanks to @whtns)

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

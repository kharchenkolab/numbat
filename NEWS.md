# numbat 0.1.3 - 07/02/2022

* Speed up of NNI using RcppParallel (#34). 10x faster and much more memory efficient (memory requirement is constant with respect to the number of threads).

* Speed up of expression single-cell testing using roptim. Approximately 2x speedup.

* New LLR metric for CNV filtering that is not inflated (default: 5).

* Only keep heterozygous SNPs in alelle dataframe to reduce memory usage

# Numbat

<!-- badges: start -->

[![<kharchenkolab>](https://circleci.com/gh/kharchenkolab/Numbat.svg?style=svg)](https://app.circleci.com/pipelines/github/kharchenkolab/Numbat)
  
<!-- badges: end -->

<img src="logo.png" align="right" width="150">

Numbat is a haplotype-enhanced CNV caller from single-cell transcriptomics data. It integrates signals from gene expression, allelic ratio, and population-derived haplotype information to accurately infer allele-specific CNVs in single cells and reconstruct their lineage relationship. 

Numbat can be used to 1. detect allele-specific copy number variations from scRNA-seq 2. differentiate tumor versus normal cells in the tumor microenvironment 3. infer the clonal architecture and evolutionary history of profiled tumors. 

  <img width="867" alt="Screen Shot 2022-02-03 at 12 05 47 PM" src="https://user-images.githubusercontent.com/13375875/152392495-024f89c6-0010-4e53-8375-9eb8a64dbe73.png">

Numbat does not require paired DNA or genotype data and operates solely on the donor scRNA-data data (for example, 10x Cell Ranger output).

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Preparing data](#preparing-data)
- [Running Numbat](#running-numbat)
- [Understanding results](#understanding-results)

# Prerequisites
Numbat uses cellsnp-lite for generating SNP pileup data and eagle2 for phasing. Please follow their installation instructions and make sure their binary executables can be found in your $PATH.

1. [cellsnp-lite](https://github.com/single-cell-genetics/cellsnp-lite)
2. [eagle2](https://alkesgroup.broadinstitute.org/Eagle/)

Additionally, Numbat needs a common SNP VCF and phasing reference panel. You can use the 1000 Genome reference below:

3. 1000G SNP reference file 
```
# hg38
wget https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
# hg19
wget https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz
```
4. 1000G Reference Panel
```
# hg38
wget http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip
# hg19
wget http://pklab.med.harvard.edu/teng/data/1000G_hg19.zip
```

# Installation
Please first install the below dependencies via `BiocManager`:
```
BiocManager::install("GenomicRanges")
BiocManager::install("Rsamtools")
BiocManager::install("ggtree")
```
and make sure that `samtools` is already installed.
Then install the Numbat R package via:
```
devtools::install_github("https://github.com/kharchenkolab/Numbat")
```

# Preparing data
1. Prepare the allele data. Run the preprocessing script (`pileup_and_phase.R`) to count alleles and phase SNPs
```
usage: pileup_and_phase.R [-h] --label LABEL --samples SAMPLES --bams BAMS
                          --barcodes BARCODES --gmap GMAP --snpvcf SNPVCF
                          --paneldir PANELDIR --outdir OUTDIR --ncores NCORES
                          [--UMItag UMITAG] [--cellTAG CELLTAG]

Run SNP pileup and phasing with 1000G

optional arguments:
  -h, --help           show this help message and exit
  --label LABEL        Individual label
  --samples SAMPLES    Sample names, comma delimited
  --bams BAMS          BAM files, one per sample, comma delimited
  --barcodes BARCODES  Cell barcodes, one per sample, comma delimited
  --gmap GMAP          Path to genetic map provided by Eagle2
  --snpvcf SNPVCF      SNP VCF for pileup
  --paneldir PANELDIR  Directory to phasing reference panel (BCF files)
  --outdir OUTDIR      Output directory
  --ncores NCORES      Number of cores
  --UMItag UMITAG      UMI tag in bam. Should be Auto for 10x and XM for
                       Slide-seq
  --cellTAG CELLTAG    Cell tag in bam. Should be CB for 10x and XC for Slide-
                       seq
```
This will produce a file named `{sample}_allele_counts.tsv.gz` under the specified output directory, which contains cell-level phased allele counts. If multiple samples from the same individual was provided, there will be one allele count file for each sample. Other outputs include phased vcfs under `phasing/` folder and raw pileup counts under `pileup/`.

2. Prepare the expression data. Numbat takes a gene by cell integer UMI count matrix as input. You can directly use results from upstream transcriptome quantification pipelines such as 10x CellRanger.
  
3. Prepare the expression reference, which is a gene by cell type matrix of normalized expression values (not raw counts!). For a quick start, you may use a our HCA collection (`ref_hca`) that ships with the package. If you have matched normal cells (ideally, of various cell type) from the same patient or dataset and would like to make your own references, you may use this utility function:
```
# count_mat is a gene x cell raw count matrices
# cell_annot is a dataframe with columns "gene" and "cell_type"
ref_internal = aggregate_counts(count_mat, cell_annot)$exp_mat
```
  
# Running Numbat
  
In this example (ATC2 from [Gao et al](https://www.nature.com/articles/s41587-020-00795-2)), the gene expression count matrix and allele dataframe are already prepared for you.
```
library(numbat)

# run
out = run_numbat(
    count_mat_ATC2, # gene x cell integer UMI count matrix 
    ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
    df_allele_ATC2, # allele dataframe generated by pileup_and_phase script
    gtf_hg38, # provided upon loading the package
    genetic_map_hg38, # provided upon loading the package
    min_cells = 20,
    t = 1e-3,
    max_iter = 2,
    min_LLR = 50,
    init_k = 3,
    ncores = 20,
    plot = TRUE,
    out_dir = './test'
)
```
## Run parameters
There are a few parameters you can consider tuning to a specific dataset. 
- `t`: the transition probability used in the HMM. A lower `t` is more appropriate for tumors with more complex copy number landscapes (from which you can expect more breakpoints) and is sometimes better for detecting subclonal events. A higher `t` is more effective for controlling false-positive rates of CNV calls.
- `max_iter`: maximum number of iterations. Numbat iteratively optimizes the phylogeny and copy number estimations. In practice, we find that results after 2 iterations are usually stable.  
- `min_LLR`: minimum log-likelihood ratio threshold to filter CNVs by. To ensure quality of phylogeny inference, we only use confident CNVs to reconstruct the phylogeny. By default, this threshold is 50.
- `max_entropy`: another criteria that we use to filter CNVs before phylogeny construction. The entropy of the binary posterior quantifies the uncertainty of an event across single cells. The value can be from 0 to 1 where 1 is the least stringent.
- `min_cells`: minimum number of cells for which an pseudobulk HMM will be run. If the allele coverage per cell is very sparse for a dataset, then I would consider setting this threshold to be higher.
- `init_k`: initial number of subclusters to use for the `hclust` initialization. Numbat by default uses hierarchical clustering (`hclust`) of smoothed expression values to approximate an initial phylogeny. This will cut the initial tree into k clusters. More clusters means more resolution at the initial stage for subclonal CNV detection. By default, we set init_k to be 3.
- `skip_nj`: if the number of cells is extremely large (>100k), the initial NJ tree construction may take a long time. You can set skip_nj to be TRUE to only use the faster UPGMA as seed for maximum likelihood tree search.
  
# Understanding results
A detailed vignette on how to interpret and visualize Numbat results is available:  
- [Interpreting Numbat results](https://kharchenkolab.github.io/Numbat)
  
Numbat generates a number of files in the output folder. The file names are post-fixed with the `i`th iteration of phylogeny optimization. Here is a detailed list:
- `gexp_roll_wide.tsv.gz`: window-smoothed normalized expression profiles of single cells
- `hc.rds`: hierarchical clustering result based on smoothed expression
- `bulk_subtrees_{i}.tsv.gz`: pseudobulk HMM profiles based on subtrees defined by current cell lineage tree
- `segs_consensus_{i}.tsv.gz`: consensus segments from subtree pseudobulk HMMs
- `bulk_clones_{i}.tsv.gz`: pseudobulk HMM profiles based on clones defined by current cell lineage tree
- `bulk_clones_{i}.png`: visualization of clone pseudobulk HMM profiles
- `exp_sc_{i}.tsv.gz`: single-cell expression profiles used for single-cell CNV testing
- `exp_post_{i}.tsv`: single-cell expression posteriors 
- `allele_post_{i}.tsv`: single-cell allele posteriors 
- `joint_post_{i}.tsv`: single-cell joint posteriors 
- `treeUPGMA_{i}.rds`: UPGMA tree
- `treeNJ_{i}.rds`: NJ tree
- `tree_list_{i}.rds`: list of candidate phylogeneies in the maximum likelihood tree search
- `tree_final_{i}.rds`: final tree after simplification
- `mut_graph_{i}.rds`: final mutation history
- `clone_post_{i}.rds`: clone assignment and tumor versus normal classification posteriors
- `bulk_subtrees_{i}.png`: visualization of subtree pseudobulk HMM profiles 
- `bulk_clones_{i}.png`: visualization of clone pseudobulk HMM profiles 
- `panel_{i}.png`: visualization of combined phylogeny and CNV heatmap

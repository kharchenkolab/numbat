# Numbat
Haplotype-aware detection of copy number variations in scRNA-Seq

![image](https://user-images.githubusercontent.com/13375875/136429050-609ee367-8d5d-4a63-8fa8-a87171aff01c.png)

# Installation
Install the Numbat R package via:
```
devtools::install_github("kharchenkolab/Numbat")
```
## Other prerequisites
1. [cellsnp-lite](https://github.com/single-cell-genetics/cellsnp-lite)
2. [ScisTree](https://github.com/kharchenkolab/ScisTree)
```
git clone https://github.com/kharchenkolab/ScisTree
cd ScisTree/ScisTree-ver1.2.0.6-src
make
./scistree
```
Please make sure this binary executable can be found in $PATH.

3. 1000 Genome SNP reference file 
```
wget https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
```
If you are doing phasing locally instead of on the cloud [server](https://imputation.biodatacatalyst.nhlbi.nih.gov), please make sure you also have:

4. [eagle2](https://alkesgroup.broadinstitute.org/Eagle/)
```
wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz
tar -xvf Eagle_v2.4.1.tar.gz
./Eagle_v2.4.1/eagle
```
5. 1000 Genome Reference Panel
```
mkdir ref
wget -r --no-parent -A "*vcf*" http://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes
mv hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/* ./ref
```

# Usage
1. Run the preprocessing script: collect allele data and phase SNPs
```
Rscript pileup_and_phase.r \
      --label T200 \
      --sample T200 \
      --bams /home/tenggao/external/NB_Dong/SRR12148217/SRR12148217_transcriptome/outs/possorted_genome_bam.bam \
      --barcodes /home/tenggao/external/NB_Dong/SRR12148217/SRR12148217_transcriptome/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
      --gmap /home/tenggao/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
      --snpvcf /home/tenggao/ref/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
      --paneldir /home/tenggao/ref/1000G \
      --outdir /home/tenggao/external/NB_Dong/processed/T200 \
      --ncores 25
```

## Running Numbat
```
library(numbat)

# Generate reference expression profile(s) using a known reference dataset
lambdas_ref = make_psbulk(count_mat_ref, cell_annot)

# run
out = numbat_subclone(
      count_mat_obs,
      lambdas_ref,
      df,
      gtf_transcript,
      t = 1e-5,
      out_dir = glue('~/results/{sample}')
)
```

# Visualizations of results
1. Collect results
```res = fetch_results(out_dir, i = 2)```

2. Evolutionary history and cell-cnv heatmap
```
plot_clone_panel(res[[sample]], ratio = 1)
```
![image](https://user-images.githubusercontent.com/13375875/136427928-ed7f67ed-4bd1-4f24-9b9e-f381b5920f54.png)

3. Aggregated clone CNV profile
```
bulk_clones %>% plot_bulks(ncol = 1)
```
![image](https://user-images.githubusercontent.com/13375875/136428374-06100e23-1527-4e35-b945-a1528dae93b3.png)

4. Original single-cell phylogeny from ScisTree
```
tree_heatmap2(
    res$joint_post %>% filter(seg %in% rownames(res$geno)),
    res$gtree
)
```
![image](https://user-images.githubusercontent.com/13375875/136428423-9f92b303-5577-482d-8214-f4bbe2115b50.png)


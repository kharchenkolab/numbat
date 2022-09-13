library(dplyr)
library(data.table)
library(stringr)
library(glue)
library(pagoda2)
library(conos)

## gtf_hg19.rda ##
# https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz
gtf = fread('~/ref/hg38.refGene.gtf')
cols = c('CHROM', 'source', 'biotype', 'start', 'end', 'a', 'strand', 'b', 'info')
colnames(gtf) = cols

gtf = gtf %>% filter(CHROM %in% paste0('chr', 1:22))

gtf_transcript = gtf %>% filter(biotype == 'transcript') %>% distinct(start, end, `.keep_all` = TRUE) %>%
    mutate(gene = str_extract(info, '(?<=gene_name \\").*(?=\\";)')) %>%
    mutate(CHROM = as.integer(str_remove(CHROM, 'chr'))) %>%
    select(gene, start, end, CHROM) %>%
    mutate(region = paste0('chr', CHROM, ':', start, '-', end)) %>%
    mutate(gene_length = end-start) %>%
    arrange(CHROM, gene, -gene_length) %>%
    distinct(gene, .keep_all = TRUE) %>%
    rename(gene_start = start, gene_end = end) %>%
    arrange(CHROM, gene_start)

## gtf_hg38.rda ##
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
gtf = fread('~/ref/hg19.refGene.gtf')
cols = c('CHROM', 'source', 'biotype', 'start', 'end', 'a', 'strand', 'b', 'info')
colnames(gtf) = cols

gtf = gtf %>% filter(CHROM %in% paste0('chr', 1:22))

gtf_transcript = gtf %>% filter(biotype == 'transcript') %>% distinct(start, end, `.keep_all` = TRUE) %>%
    mutate(gene = str_extract(info, '(?<=gene_name \\").*(?=\\";)')) %>%
    mutate(CHROM = as.integer(str_remove(CHROM, 'chr'))) %>%
    select(gene, start, end, CHROM) %>%
    mutate(region = paste0('chr', CHROM, ':', start, '-', end)) %>%
    mutate(gene_length = end-start) %>%
    arrange(CHROM, gene, -gene_length) %>%
    distinct(gene, .keep_all = TRUE) %>%
    rename(gene_start = start, gene_end = end) %>%
    arrange(CHROM, gene_start)

## chrom_sizes_hg38.rda ##
chrom_sizes_hg38 = fread('~/ref/hg38.chrom.sizes.txt') %>% 
    setNames(c('CHROM', 'size')) %>%
    mutate(CHROM = str_remove(CHROM, 'chr')) %>%
    filter(CHROM %in% 1:22) %>%
    mutate(CHROM = factor(as.integer(CHROM)))

## gaps_hg38.rda ##
# https://github.com/hartwigmedical/hmftools/blob/master/purple/src/main/resources/circos/gaps.38.txt
gaps_hg38 = fread('~/ref/gaps.38.txt') %>%
        setNames(c('CHROM', 'start', 'end')) %>%
        mutate(CHROM = str_remove(CHROM, 'hs')) %>%
        filter(CHROM %in% 1:22) %>%
        mutate(CHROM = as.factor(as.integer(CHROM))) %>%
        arrange(CHROM)

## acen_hg38.rda ##
acen_hg38 = fread('~/ref/chromosome.band.hg38.txt') %>%
    rename(CHROM = `#chrom`) %>%
    filter(gieStain == 'acen') %>%
    mutate(start = chromStart, end = chromEnd) %>%
    mutate(CHROM = str_remove(CHROM, 'chr')) %>%
    filter(CHROM %in% 1:22) %>%
    mutate(CHROM = factor(CHROM, 1:22)) %>%
    group_by(CHROM) %>%
    summarise(start = min(start), end = max(end))

## chrom_sizes_hg19.rda ##
chrom_sizes_hg19 = fread('~/ref/hg19.chrom.sizes.txt') %>% 
    setNames(c('CHROM', 'size')) %>%
    mutate(CHROM = str_remove(CHROM, 'chr')) %>%
    filter(CHROM %in% 1:22) %>%
    mutate(CHROM = factor(as.integer(CHROM)))

## acen_hg19.rda ##
acen_hg19 = fread('~/ref/chromosome.band.hg19.txt') %>%
    setNames(c('CHROM', 'start', 'end', 'band', 'gieStain')) %>%
    filter(gieStain == 'acen') %>%
    mutate(CHROM = str_remove(CHROM, 'chr')) %>%
    filter(CHROM %in% 1:22) %>%
    mutate(CHROM = factor(CHROM, 1:22)) %>%
    group_by(CHROM) %>%
    summarise(start = min(start), end = max(end))

## gaps_hg19.rda ##
# https://github.com/hartwigmedical/hmftools/blob/master/purple/src/main/resources/circos/gaps.19.txt
gaps_hg19 = fread('~/ref/gaps.19.txt') %>%
        setNames(c('CHROM', 'start', 'end')) %>%
        mutate(CHROM = str_remove(CHROM, 'hs')) %>%
        filter(CHROM %in% 1:22) %>%
        mutate(CHROM = as.factor(as.integer(CHROM))) %>%
        arrange(CHROM)

## vcf_meta.rda ##
# wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz
# https://github.com/picrin/pysccnv/blob/master/Eagle2_benchmark.ipynb
vcf_template = read.vcfR('~/bench/reverse_jewish_son_chr1.vcf.gz', verbose = F)
vcf_meta = vcf_template@meta

## ATC2 subsampled ##
con = readRDS(glue('~/paper_data/conos_objects/conos_ATC.rds'))

cell_annot = fread(glue('~/paper_data/cell_annotations/cell_annot_MDA.tsv')) %>% 
    mutate(annot = copykat.pred) %>% 
    split(.$sample)

set.seed(0)
tumor_cells = cell_annot[['ATC2']] %>% filter(annot == 'T') %>% pull(cell)
normal_cells = cell_annot[['ATC2']] %>% filter(annot == 'N') %>% pull(cell) %>% sample(50)
cells = c(tumor_cells, normal_cells)

count_mat = as.matrix(t(con$samples[['ATC2']]$misc$rawCounts))
df = fread(glue('~/paper_data/processed/ATC2_allele_counts.tsv.gz'), sep = '\t')

count_mat_ATC2 = count_mat[,cells]
df_allele_ATC2 = df %>% filter(cell %in% cells)

## ref_hca.rda ##
library(stringr)
library(dplyr)
library(magrittr)
library(Seurat)
library(biomaRt)

# Seurat object from Synapse under the ID syn21041850
seu = readRDS('~/external/HCA/krasnow.rds')
count_mat = GetAssayData(object = seu, slot = "counts")

# get gene names
gids = count_mat %>% rownames
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
res = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'),
      filters = 'ensembl_gene_id',
      values = gids, 
      mart = ensembl)

# create the pseudobulks
gene_dict = res %>% {setNames(.$external_gene_name, .$ensembl_gene_id)}
genes = sapply(rownames(count_mat), function(x){unname(gene_dict[x])})
rownames(count_mat) = genes

cell_annot = data.frame(
    cell = names(seu$cell_type),
    cell_type = unname(seu$cell_type)
)

ref = numbat::aggregate_counts(as.matrix(count_mat), cell_annot)

# rename cell types
rename_dict = c('natural killer cell' = 'NK', 
  'fibroblast' = 'Fibroblast', 
  'macrophage' = 'Macrophage', 
  'effector memory CD4-positive, alpha-beta T cell' = 'CD4+T', 
  'effector memory CD8-positive, alpha-beta T cell' = 'CD8+T',
  'endothelial cell' = 'Endothelial',
  'myeloid leukocyte' = 'Myeloid',
  'classical monocyte' = 'Monocyte',
  'dendritic cell' = 'Dendritic',
  'plasma cell' = 'Plasma',
  'B cell' = 'B',
  'epithelial cell' = 'Epithelial')

ref_hca = ref$exp_mat[,names(rename_dict)]
colnames(ref_hca) = unname(sapply(colnames(ref_hca), function(x){rename_dict[x]}))

ref_hca_counts = ref$count_mat[,names(rename_dict)]
colnames(ref_hca_counts) = unname(sapply(colnames(ref_hca_counts), function(x){rename_dict[x]}))


## count_mat_example.rda ##
genes = gtf_hg38 %>% filter(CHROM == 5) %>% pull(gene)
genes = intersect(genes, rownames(count_mat_ATC2))
count_mat_example = count_mat_ATC2[genes,]

## df_allele_example.rda ##
df_allele_example = df_allele_ATC2 %>% filter(CHROM == 5)

## count_mat_ref.rda ##
genes_common = intersect(gtf_hg38$gene, rownames(count_mat_ATC2))
count_mat_ref = count_mat_ATC2[genes_common[1:1000], normal_cells]

## annot_ref.rda ##
annot_ref = fread('~/paper_data/cell_annotations/cell_annot_MDA.tsv') %>% 
    filter(cell %in% colnames(count_mat_ATC2)) %>%
    filter(copykat.pred == 'N') %>%
    mutate(group = c(rep('Immune', 25), rep('Endothelial', 20), rep('Fibroblast', 5))) %>%
    select(cell, group) %>%
    as.data.frame()

## mut_graph_example.rda ##
mut_graph_example = readRDS('~/ATC2_test/mut_graph_1.rds')

## segs_example.rda ##
segs_example = fread('~/ATC2_test/segs_consensus_1.tsv') %>% mutate(group = 'ATC2') %>% relevel_chrom()

## phylogeny_example.rda ##
phylogeny_example = readRDS('~/ATC2_test/tree_final_1.rds')

## gexp_roll_example.rda ##
gexp_roll_example = nb$gexp_roll_wide[1:10,sample(1:ncol(nb$gexp_roll_wide), 2000)]

## hc_example.rda ##
hc_example = hclust(dist(gexp_roll_example))

## joint_post_example.rda ##
joint_post_example = fread('~/ATC2_test/joint_post_1.tsv')

## Data used for unit tests

### pre_likelihood_hmm
### Input to likelihood_allele() for unit tests

### Created by running the example with the following line changed
### in the function calc_allele_lik():
###
### calc_allele_lik = function (pAD, DP, p_s, theta, gamma = 20) {
###     hmm = get_allele_hmm(pAD, DP, p_s, theta, gamma)
###     saveRDS(hmm, "pre_likelihood_hmm.rds")
###     LL = likelihood_allele(hmm)
###     return(LL)
### }

### The example run to create this RDS file was the following:
### """
### library(numbat)
### 
### bulk = get_bulk(
###      count_mat = count_mat_ATC2,
###      df_allele = df_allele_ATC2,
###      lambdas_ref = ref_hca,
###      gtf_transcript = gtf_hg38,
###      genetic_map = genetic_map_hg38)
### 
### bulk = analyze_bulk(bulk, t = 1e-5)
### """

library(dplyr)
library(data.table)
library(stringr)
library(glue)
library(pagoda2)
library(conos)

## genetic_map_hg19.rda ##
genetic_map = fread('~/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz') %>% 
    setNames(c('CHROM', 'POS', 'rate', 'cM')) %>%
    group_by(CHROM) %>%
    mutate(
        start = POS,
        end = c(POS[2:length(POS)], POS[length(POS)])
    ) %>%
    ungroup()

## genetic_map_hg38.rda ##
genetic_map = fread('~/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz') %>% 
    setNames(c('CHROM', 'POS', 'rate', 'cM')) %>%
    group_by(CHROM) %>%
    mutate(
        start = POS,
        end = c(POS[2:length(POS)], POS[length(POS)])
    ) %>%
    ungroup()

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
#!/usr/bin/env Rscript

library(argparse, quietly = T)
library(ggplot2, quietly = T)
library(dplyr, quietly = T)
library(Matrix, quietly = T)
library(data.table, quietly = T)
library(stringr, quietly = T)
library(patchwork, quietly = T)
library(glue, quietly = T)
library(ggforce, quietly = T)
library(vcfR, quietly = T)
library(ggpubr, quietly = T)
library(extraDistr, quietly = T)

parser <- ArgumentParser(description='Process some integers')
parser$add_argument('--sample', type = "character", help = "sample name")
parser$add_argument('--vcf', type = "character", help = "pileup VCF")
parser$add_argument('--AD', type = "character", help = "allele count matrix")
parser$add_argument('--DP', type = "character", help = "total depth matrix")
parser$add_argument('--barcodes', type = "character", help = "total depth matrix")
parser$add_argument('--pvcf', type = "character", help = "phased VCFs directory")
parser$add_argument('--out', type = "character", help = "output directory")

# parser$add_argument('gtf', type = "character")
args <- parser$parse_args()
print(args)

dir.create(args$out)

cat('Analyzing', args$sample, '\n')

sample = args$sample


######## Preprocessing #########
chr_lengths = fread('~/ref/hg38.chrom.sizes.txt') %>% setNames(c('CHROM', 'CHROM_SIZE')) %>%
    mutate(CHROM = str_remove(CHROM, 'chr'))

# read in GTF
gtf = fread('~/ref/hg38.refGene.gtf')
cols = c('CHROM', 'source', 'biotype', 'start', 'end', 'a', 'strand', 'b', 'info')
colnames(gtf) = cols

gtf = gtf %>% filter(CHROM %in% paste0('chr', 1:22))

transcript_regions = gtf %>% filter(biotype == 'transcript') %>%
    mutate(region = paste0(CHROM, ':', start, '-', end)) %>%
    pull(region) %>%
    unique %>% bedr::bedr.sort.region(verbose = F)

cat('Reading in snp count matrix..', '\n')

vcf = fread(args$vcf) %>% rename(CHROM = `#CHROM`)

vcf = vcf %>%
    mutate(INFO = str_remove_all(INFO, '[:alpha:]|=')) %>%
    tidyr::separate(col = 'INFO', into = c('AD', 'DP', 'OTH'), sep = ';') %>%
    mutate_at(c('AD', 'DP', 'OTH'), as.integer)

vcf = vcf %>% mutate(snp_id = paste(CHROM, POS, REF, ALT, sep = '_'))

AD = readMM(args$AD)
DP = readMM(args$DP)
cell_barcodes = fread(args$barcodes, header = F) %>% pull(V1)

DP = as.data.frame(summary(DP)) %>%
    mutate(
        cell = cell_barcodes[j],
        snp_id = vcf$snp_id[i]
    ) %>%
    select(-i, -j) %>%
    rename(DP = x) %>%
    select(cell, snp_id, DP)

AD = as.data.frame(summary(AD)) %>%
    mutate(
        cell = cell_barcodes[j],
        snp_id = vcf$snp_id[i]
    ) %>%
    select(-i, -j) %>%
    rename(AD = x) %>%
    select(cell, snp_id, AD)

df = DP %>% left_join(AD, by = c("cell", "snp_id")) %>%
    mutate(AD = ifelse(is.na(AD), 0, AD))

df = df %>% left_join(
    vcf %>% rename(AD_all = AD, DP_all = DP, OTH_all = OTH),
    by = 'snp_id')

df = df %>% mutate(
        AR = AD/DP,
        AR_all = AD_all/DP_all
    )

df = df %>% filter(DP_all > 1 & OTH_all == 0)

df = df %>% mutate(
    snp_index = as.integer(factor(snp_id, unique(snp_id))),
    cell_index = as.integer(factor(cell, sample(unique(cell))))
)

# read. inphased VCF
cat('Reading in phased VCF..', '\n')

vcf_phased = lapply(1:22, function(chr) {
    fread(glue('{args$pvcf}/{sample}_chr{chr}_phased.vcf.gz')) %>%
    rename(CHROM = `#CHROM`) %>%
    mutate(CHROM = str_remove(CHROM, 'chr'))   
}) %>% Reduce(rbind, .) %>%
mutate(CHROM = factor(CHROM, unique(CHROM)))

vcf_phased = vcf_phased %>% mutate(INFO = str_remove_all(INFO, '[:alpha:]|=')) %>%
    tidyr::separate(col = 'INFO', into = c('AD', 'DP', 'OTH'), sep = ';') %>%
    mutate_at(c('AD', 'DP', 'OTH'), as.integer)

vcf_phased = vcf_phased %>% mutate(snp_id = paste(CHROM, POS, REF, ALT, sep = '_')) %>%
    mutate(GT = get(sample))

vcf_phased = vcf_phased %>% mutate(region = paste0('chr', CHROM, ':', POS, '-', format(POS+1, scientific = F, trim = T)))

# intersect with gene model
vcf_phased_regions = vcf_phased %>%
    pull(region) %>%
    bedr::bedr.sort.region(verbose = F)

overlap_transcript = bedr::bedr(
    input = list(a = vcf_phased_regions, b = transcript_regions), 
    method = "intersect", 
    params = "-loj -sorted",
    verbose = F
) %>% filter(V4 != '.')

# annotate SNP by gene
cat("Annotating SNPs by gene..", '\n')

vcf_phased = vcf_phased %>% left_join(
    overlap_transcript %>%
        rename(start = V5, end = V6) %>%
        mutate(start = as.integer(start), end = as.integer(end)) %>%
        left_join(
            gtf %>% filter(biotype == 'transcript') %>% distinct(start, end, `.keep_all` = TRUE) %>%
            mutate(gene = str_extract(info, '(?<=gene_name \\").*(?=\\";)')) %>%
            select(start, end, gene),
            by = c('start', 'end')
        ) %>%
    arrange(index, gene) %>%
    distinct(index, `.keep_all` = T),
    by = c('region' = 'index')
) %>% rename(gene_start = start, gene_end = end)


# add annotation to cell counts and pseudobulk
cat('Aggregating into pseudobulk..', '\n')

df = df %>% left_join(vcf_phased %>% select(snp_id, gene, gene_start, gene_end, GT), by = 'snp_id') 
df = df %>% mutate(CHROM = factor(CHROM, unique(CHROM)))

pseudobulk = df %>% distinct(snp_index, .keep_all = TRUE) %>% 
    mutate(CHROM = factor(CHROM, unique(CHROM))) %>%
    select(-cell, -AD, -DP, -AR) %>% 
    rename(AD = AD_all, DP = DP_all, OTH = OTH_all, AR = AR_all) %>%
    arrange(CHROM, POS) %>%
    mutate(snp_index = 1:n())

# plot SNP densities
D = vcf_phased %>% 
    mutate(genotype = ifelse(GT == '1|1', '1/1', '0/1')) %>%
    count(CHROM, genotype) %>%
    dcast(CHROM ~ genotype) %>%
    mutate(total = `0/1` + `1/1`) %>%
    left_join(chr_lengths) %>%
    mutate(
        het = `0/1`/CHROM_SIZE * 1e6,
        hom = `1/1`/CHROM_SIZE * 1e6
    ) %>%
    melt(
        measure.vars = c('het', 'hom'),
        variable.name = 'genotype',
        value.name = 'snp_density'
    ) %>%
    mutate(CHROM = factor(CHROM, unique(CHROM)))

p = ggplot(
        D,
        aes(x = CHROM, y = snp_density, fill = genotype)
    ) +
    geom_col() +
    theme_bw() +
    ylab('SNPs/Mb')

ggsave(glue('{args$out}/snp_density_{sample}.png'), p, width = 7, height = 3, dpi = 200)

# plot phased BAF
D = pseudobulk %>%
    filter(DP >= 20) %>%
    filter(AR > 0.1 & AR < 0.9) %>%
    filter(!is.na(GT)) %>%
    mutate(CHROM = factor(CHROM, unique(CHROM))) %>%
    arrange(CHROM, POS) %>%
    mutate(snp_index = 1:n()) 

p = ggplot(
    D,
    aes(x = snp_index, y = AR, color = GT)
) +
geom_point(size = 1, alpha = 0.5) +
theme_classic() +
theme(
    legend.position = 'none'
) +
facet_wrap(~CHROM, nrow = 3, scale = 'free_x') 

ggsave(glue('{args$out}/BAF_{sample}.png'), p, width = 12, height = 5, dpi = 300)

######## All cell HMM #########
cat('Running HMM on all cells', '\n')
Obs = pseudobulk %>%
    filter(DP >= 10) %>%
    filter(GT %in% c('1|0', '0|1')) %>%
    mutate(pBAF = ifelse(GT == '1|0', AR, 1-AR)) %>%
    mutate(pAD = ifelse(GT == '1|0', AD, DP - AD)) %>%
    mutate(CHROM = factor(CHROM, unique(CHROM))) %>%
    arrange(CHROM, POS) %>%
    group_by(CHROM) %>%
    mutate(snp_index = 1:n()) %>%
    ungroup

# states
N = c("theta_up", "neu", "theta_down")

# probabilty to stay in the same state
p_0 = 1-1e-5
# probability of phase switch
p_s = 0.1

# transition matrix
A <- matrix(
    c(p_0 * (1-p_s), 1 - p_0, p_0 * p_s, 
     (1-p_0)/2, p_0, (1-p_0)/2,
     p_0 * p_s, 1-p_0, p_0 * (1 - p_s)),
    ncol = length(N),
    byrow = TRUE
)

# intitial probabilities
prior = rep(1/length(N), length(N))

Obs = Obs %>% 
    group_by(CHROM) %>%
    mutate(
        state = HiddenMarkov::Viterbi(HiddenMarkov::dthmm(
            x = pAD, 
            Pi = A, 
            delta = prior, 
            distn = "bbinom",
            pm = list(alpha=c(10,10,6), beta=c(6,10,10)),
            pn = list(size = DP),
            discrete=TRUE))
    ) %>%
    mutate(state = c('1' = 'theta_up', '2' = 'neutral', '3' = 'theta_down')[state])

# plot result
p1 = ggplot(
        Obs %>% filter(DP >= 10) %>% mutate(cut_2 = 'All'),
        aes(x = snp_index, y = AR, color = GT)
    ) +
    geom_point(
        size = 0.5, alpha = 0.5
    ) +
    theme_classic() +
    theme(
        legend.position = 'none',
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
        strip.text.x = element_text(size = 6)
    ) +
    facet_grid(cut_2~CHROM, scale = 'free_x', space = 'free_x')

p2 = ggplot(
        Obs %>% filter(DP >= 10) %>% mutate(cut_2 = 'All'),
        aes(x = snp_index, y = pBAF, color = state)
    ) +
    geom_point(
        size = 0.5, alpha = 0.5
    ) +
    theme_classic() +
    theme(
        legend.position = 'none',
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
        strip.text.x = element_text(size = 6)
    ) +
    facet_grid(cut_2~CHROM, scale = 'free_x', space = 'free_x') +
    scale_color_manual(values = c('gray', 'red', 'darkred'))

# ggsave(glue('{args$out}/HMM_all_{sample}.png'), p, width = 12, height = 3.5, dpi = 300)

######## Clustering #########
df = df %>% 
    select(-any_of('GT_MAF')) %>%
    left_join(
    pseudobulk %>%
    filter(GT %in% c('1|0', '0|1')) %>%
    filter(DP >= 10) %>%
    mutate(GT_MAF = ifelse(AR >= 0.5, '1|0', '0|1')) %>%
    select(snp_id, GT_MAF),
    by = "snp_id"
) %>%
mutate(MAF = ifelse(GT_MAF == '1|0', AR, 1-AR))

D = df %>%
    filter(!is.na(MAF)) %>%
    filter(DP_all >= 30) %>%
    group_by(CHROM) %>%
    mutate(
        snp_index = as.integer(factor(snp_id, unique(snp_id)))
    ) %>%
    ungroup() %>%
    mutate(cell = factor(cell, unique(cell)))

mat = df %>% filter(!is.na(MAF)) %>%
    filter(DP_all >= 30) %>%
    reshape2::dcast(snp_id ~ cell, value.var = 'MAF') %>%
    tibble::column_to_rownames('snp_id')

mat.smooth = apply(mat, 2, caTools::runmean, k = 100)

row.names(mat.smooth) = row.names(mat)

D_smooth = mat.smooth %>% reshape2::melt() %>%
    setNames(c('snp_id', 'cell', 'MAF')) %>%
    mutate(snp_id = as.character(snp_id), cell = as.character(cell)) %>%
    filter(!is.na(MAF)) %>% 
    left_join(
        D %>% distinct(snp_id, CHROM, snp_index),
        by = "snp_id"
    )

dist_mat_file = glue('{args$out}/dist_mat_{sample}.tsv')

if (file.exists(dist_mat_file)) {
    cat('Using existing distance matrix..', '\n')
    d = fread(dist_mat_file)
    d = as.dist(as.matrix(d))
} else {
    cat('Computing distance matrix..', '\n')
    d = amap::Dist(t(mat.smooth), nbproc = 20)
    fwrite(as.matrix(d), dist_mat_file, sep = '\t')
}

d[is.na(d)] <- 0
d[is.nan(d)] <- 0
d[is.infinite(d)] <- 0

hc <- hclust(d, method="ward.D2")

cluster_cuts = data.frame(
    sapply(1:5, function(h){cutree(hc, k = h)})
) %>%
setNames(str_replace(names(.), 'X', 'cut_')) %>% 
tibble::rownames_to_column('cell')

df = df %>% 
    left_join(
        cluster_cuts,
        by = "cell"
    )

png(file = glue("{args$out}/clusters_{sample}.png"), width = 800, height = 500)
plot(hc, hang = -1, cex = 0.1)
dev.off()

######## cluster-specific HMM ##########
pseudobulk_cut2 = df %>%
    filter(GT %in% c('1|0', '0|1')) %>%
    group_by(snp_id, cut_2) %>%
    summarise(
        AD = sum(AD),
        DP = sum(DP),
        AR = AD/DP,
        .groups = 'drop'
    ) %>%
    left_join(
        vcf_phased %>% select(CHROM, POS, REF, ALT, snp_id, GT, gene, gene_start, gene_end),
        by = "snp_id"
    ) %>%
    arrange(CHROM, POS) %>%
    mutate(snp_index = as.integer(factor(snp_id, unique(snp_id)))) %>%
    ungroup() %>%
    mutate(pBAF = ifelse(GT == '1|0', AR, 1-AR)) %>%
    mutate(pAD = ifelse(GT == '1|0', AD, DP - AD))

p3 = ggplot(
        pseudobulk_cut2 %>% filter(DP >= 10),
        aes(x = snp_index, y = AR, color = GT)
    ) +
    geom_point(size = 1, alpha = 0.5) +
    theme_classic() +
    theme(
        legend.position = 'none',
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
        strip.text.x = element_text(size = 6)
    ) +
    facet_grid(cut_2~CHROM, scale = 'free_x', space = 'free_x') 

# run HMM
Obs = pseudobulk_cut2 %>%
    filter(DP >= 10)

Obs = Obs %>% 
    group_by(CHROM, cut_2) %>%
    mutate(
        state = HiddenMarkov::Viterbi(HiddenMarkov::dthmm(
            x = pAD, 
            Pi = A, 
            delta = prior, 
            distn = "bbinom",
            pm = list(alpha=c(10,10,6), beta=c(6,10,10)),
            pn = list(size = DP),
            discrete=TRUE))
    ) %>%
    mutate(state = c('1' = 'theta_up', '2' = 'neutral', '3' = 'theta_down')[state])

p4 = ggplot(
        Obs %>% filter(DP >= 10),
        aes(x = snp_index, y = pBAF, color = state)
    ) +
    geom_point(
        size = 0.5, alpha = 0.5
    ) +
    theme_classic() +
    theme(
        legend.position = 'none',
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
        strip.text.x = element_text(size = 6)
    ) +
    facet_grid(cut_2~CHROM, scale = 'free_x', space = 'free_x') +
    scale_color_manual(values = c('gray', 'red', 'darkred'))

## level 3
pseudobulk_cut3 = df %>%
    filter(GT %in% c('1|0', '0|1')) %>%
    rename(cut = cut_3) %>%
    group_by(snp_id, cut) %>%
    summarise(
        AD = sum(AD),
        DP = sum(DP),
        AR = AD/DP,
        .groups = 'drop'
    ) %>%
    left_join(
        vcf_phased %>% select(CHROM, POS, REF, ALT, snp_id, GT, gene, gene_start, gene_end),
        by = "snp_id"
    ) %>%
    arrange(CHROM, POS) %>%
    mutate(snp_index = as.integer(factor(snp_id, unique(snp_id)))) %>%
    ungroup() %>%
    mutate(pBAF = ifelse(GT == '1|0', AR, 1-AR)) %>%
    mutate(pAD = ifelse(GT == '1|0', AD, DP - AD))

p5 = ggplot(
        pseudobulk_cut3 %>% filter(DP >= 10),
        aes(x = snp_index, y = AR, color = GT)
    ) +
    geom_point(size = 1, alpha = 0.5) +
    theme_classic() +
    theme(
        legend.position = 'none',
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
        strip.text.x = element_text(size = 6)
    ) +
    facet_grid(cut~CHROM, scale = 'free_x', space = 'free_x') 

# run HMM
Obs = pseudobulk_cut3 %>%
    filter(DP >= 10) %>%
    group_by(CHROM, cut) %>%
    mutate(
        state = HiddenMarkov::Viterbi(HiddenMarkov::dthmm(
            x = pAD, 
            Pi = A, 
            delta = prior, 
            distn = "bbinom",
            pm = list(alpha=c(10,10,6), beta=c(6,10,10)),
            pn = list(size = DP),
            discrete=TRUE))
    ) %>%
    mutate(state = c('1' = 'theta_up', '2' = 'neutral', '3' = 'theta_down')[state])

p6 = ggplot(
        Obs %>% filter(DP >= 10),
        aes(x = snp_index, y = pBAF, color = state)
    ) +
    geom_point(
        size = 0.5, alpha = 0.5
    ) +
    theme_classic() +
    theme(
        legend.position = 'none',
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
        strip.text.x = element_text(size = 6)
    ) +
    facet_grid(cut~CHROM, scale = 'free_x', space = 'free_x') +
    scale_color_manual(values = c('gray', 'red', 'darkred'))

panel = p1/p2/p3/p4/p5/p6 + plot_layout(heights = c(1,1,2,2,3,3))

ggsave(glue('{args$out}/BAF_cut2_{sample}.png'), panel, width = 18, height = 15, dpi = 300)

## cell annotations
cell_annot = readRDS('/home/meisl/neuroblastoma/conos/cell.annotation.rds')

cell_annot = data.frame(
        cell = names(cell_annot),
        cell_type = unname(cell_annot)
    ) %>%
    filter(str_detect(cell, sample)) %>%
    mutate(cell = str_remove(cell, paste0(sample, '_')))

df = df %>% left_join(
    cell_annot,
    by = "cell"
) 

p = ggplot(
    df %>% distinct(cell, `.keep_all` = TRUE) %>%
        count(cut_2, cell_type) %>%
        group_by(cut_2) %>%
        mutate(prop = n/sum(n)) %>%
        ungroup,
    aes(x = cell_type, y = factor(cut_2, c(2, 1)), fill = prop, label = n)
) +
geom_tile() +
geom_text(size = 2) +
theme_classic() +
scale_fill_gradient(low = 'white', high = 'red') +
theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = 'top'
) +
ylab('Cluster')

ggsave(glue('{args$out}/celltypes_cut2_{sample}.png'), p, width = 6, height = 2.5, dpi = 300)

## save cell-level dataframe
fwrite(df, glue('{args$out}/df_{sample}.tsv'), sep = '\t')


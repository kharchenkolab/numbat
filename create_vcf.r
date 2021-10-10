#!/usr/bin/env Rscript
library(dplyr)
library(Matrix)
library(data.table)
library(stringr)
library(glue)
library(vcfR)
library(argparse)
source('genotyping.r')

parser <- ArgumentParser(description='Create VCF from scRNA dataset for phasing')
parser$add_argument('--samples', type = "character", help = "sample names")
parser$add_argument('--label', type = "character", help = "output sample label")

args <- parser$parse_args()

samples = str_split(args$samples, ',')[[1]]
label = args$label

###### Main ######

vcf_template = read.vcfR('/home/tenggao/bench/reverse_jewish_son_chr1.vcf.gz', verbose = F)

if (length(samples) == 1) {
    
    sample = samples

    vcf_original = read.vcfR(glue('/home/tenggao/pileup/{sample}/cellSNP.base.vcf'), verbose = F)

    snps = get_snps(sample)

    for (chr in 1:22) {
        create_vcf(chr, snps, vcf_original, vcf_template, sample)
    }

} else {
    
    vcf_original = read.vcfR(glue('/home/tenggao/pileup/{samples[1]}/cellSNP.base.vcf'), verbose = F)

    snps = lapply(samples, function(sample){get_snps(sample)}) %>%
        bind_rows() %>%
        group_by(CHROM, POS, REF, ALT, snp_id) %>% 
        summarise(
            AD = sum(AD),
            DP = sum(DP),
            OTH = sum(OTH),
            .groups = 'drop'
        ) %>%
        mutate(AR = AD/DP) %>%
        arrange(CHROM, POS)

    vcf_original@fix = snps %>% 
        arrange(CHROM, POS) %>%
        mutate(ID = NA, QUAL = NA, FILTER = 'PASS', INFO = paste0('AD=', AD, ';', 'DP=', DP, ';', 'OTH=', OTH)) %>%
        select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) %>%
        # to prevent as.matrix from adding leading whitespace
        mutate(POS = as.character(POS)) %>%
        as.matrix

    for (chr in 1:22) {
        create_vcf(chr, snps, vcf_original, vcf_template, label, het_only = FALSE)
    }
}


  



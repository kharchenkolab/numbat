#!/usr/bin/env Rscript
library(dplyr)
library(Matrix)
library(data.table)
library(stringr)
library(glue)
library(vcfR)
library(argparse)

parser <- ArgumentParser(description='Create VCF from scRNA dataset for phasing')
parser$add_argument('--samples', type = "character", help = "sample names")
parser$add_argument('--label', type = "character", help = "output sample label")

args <- parser$parse_args()

samples = str_split(args$samples, ',')[[1]]
label = args$label

## per chrom function
create_vcf = function(chr, snps, vcf_original, vcf_template, sample, chr_prefix = TRUE) {
    
    chr_snps = snps %>%
        filter(CHROM == chr) %>%
        mutate(
            het = 0.1 <= AR & AR <= 0.9,
            hom_alt = AR == 1 & DP >= 10,
            hom_ref = AR == 0 & DP >= 10
        ) %>% 
        filter(het | hom_alt) %>%
        mutate(
            GT = case_when(
                het ~ '0/1',
                hom_alt ~ '1/1',
                hom_ref ~ '0/0'
            )
        ) %>%
    distinct(snp_id, `.keep_all` = T)
    
    vcf_chr = vcf_original[vcf_original@fix[,1] == chr, ]

    vcf_chr_fix = data.frame(vcf_chr@fix) %>%
        mutate(POS = as.integer(POS)) %>%
        mutate(snp_id = paste(CHROM, POS, REF, ALT, sep = '_'))
    
    # sometimes the vcf from cell-snp-lite has duplicated records
    vcf_chr = vcf_chr[vcf_chr_fix$snp_id %in% chr_snps$snp_id & (!duplicated(vcf_chr_fix$snp_id)), ]

    vcf_chr@meta = vcf_template@meta

    # make sure they are sorted the same way
    vcf_chr = vcf_chr[order(as.integer(vcf_chr@fix[,2])),]
    chr_snps = chr_snps %>% arrange(POS)

    if (any(chr_snps$POS != as.integer(vcf_chr@fix[,2]))) {
        stop('Not sorted')
    }

    vcf_chr@gt = data.frame(
        FORMAT = rep('GT', nrow(vcf_chr)),
        GT = chr_snps$GT
    ) %>% 
    setNames(c('FORMAT', sample)) %>% 
    as.matrix
    
    if (chr_prefix) {
        vcf_chr@fix[,1] = paste0('chr', vcf_chr@fix[,1])
    }

    file_name = glue('~/phasing/{sample}_chr{chr}.vcf.gz')
    
    write.vcf(vcf_chr, file_name)

    vcfR_file_fix <- function(file) {

      out <- R.utils::gunzip(file)
      out2 <- Rsamtools::bgzip(out, dest=sprintf("%s.gz", sub("\\.gz$", "", out)),
                     overwrite = TRUE)
      file.remove(out)

    }

    # compress using bgzip
    vcfR_file_fix(file_name)

    # index the vcf
    cmd = paste0('tabix ', file_name)

    system(cmd)
    
    return(vcf_chr)
    
}

# read in VCFs

get_snps = function(sample) {
    
    snps = fread(glue('/home/tenggao/pileup/{sample}/cellSNP.base.vcf')) %>% rename(CHROM = `#CHROM`)

    snps = snps %>%
        mutate(INFO = str_remove_all(INFO, '[:alpha:]|=')) %>%
        tidyr::separate(col = 'INFO', into = c('AD', 'DP', 'OTH'), sep = ';') %>%
        mutate_at(c('AD', 'DP', 'OTH'), as.integer)

    snps = snps %>% mutate(snp_id = paste(CHROM, POS, REF, ALT, sep = '_'))

    snps = snps %>% arrange(CHROM, POS) %>%
        mutate(CHROM = factor(CHROM, unique(CHROM))) %>%
        mutate(snp_index = 1:n()) %>%
        mutate(AR = AD/DP) %>%
        filter(!is.na(AR))
    
    return(snps)
}

merge_vcfs = function(samples) {

    

    return(list('vcf' = vcf, 'snps' = snps))
}

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
        as.matrix

    for (chr in 1:22) {
        create_vcf(chr, snps, vcf_original, vcf_template, label)
    }
}


  



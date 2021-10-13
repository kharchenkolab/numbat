#' Main function
#' @param label character string 
#' @param samples list of strings
#' @param vcfs vcfR objects
#' @return a status code
#' @export
genotype = function(label, samples, vcfs, outdir) {

    snps = lapply(
            vcfs, function(vcf){get_snps(vcf)}
        ) %>%
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

    vcf_original = vcfs[[1]]

    vcf_original@fix = snps %>% 
        arrange(CHROM, POS) %>%
        mutate(ID = NA, QUAL = NA, FILTER = 'PASS', INFO = paste0('AD=', AD, ';', 'DP=', DP, ';', 'OTH=', OTH)) %>%
        select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) %>%
        # to prevent as.matrix from adding leading whitespace
        mutate(POS = as.character(POS)) %>%
        as.matrix

    for (chr in 1:22) {
        make_vcf_chr(chr, snps, vcf_original, label, outdir, het_only = FALSE)
    }

    return(0)
}

#' process VCFs into SNP dataframe
#' @param vcf vcfR object
get_snps = function(vcf) {

    snps = as.data.frame(vcf@fix) %>%
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

## per chrom function
make_vcf_chr = function(chr, snps, vcf_original, vcf_template, label, outdir = '~/phasing', het_only = TRUE, chr_prefix = TRUE) {
    
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

    if (het_only) {
        chr_snps = chr_snps %>% filter(het)
    }
    
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
    setNames(c('FORMAT', label)) %>% 
    as.matrix
    
    if (chr_prefix) {
        vcf_chr@fix[,1] = paste0('chr', vcf_chr@fix[,1])
    }

    file_name = glue('{outdir}/{label}_chr{chr}.vcf.gz')
    
    write.vcf(vcf_chr, file_name)

    # compress using bgzip
    vcfR_file_fix(file_name)

    # index the vcf
    cmd = paste0('tabix ', file_name)

    system(cmd)
    
    return(vcf_chr)
    
}

vcfR_file_fix <- function(file) {

    out <- R.utils::gunzip(file)
    out2 <- Rsamtools::bgzip(out, dest=sprintf("%s.gz", sub("\\.gz$", "", out)),
                    overwrite = TRUE)
    file.remove(out)

}

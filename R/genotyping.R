#' Genotyping main function
#'
#' @param label character Individual/sample label 
#' @param samples vector Sample names
#' @param vcfs list of vcfR VCFs from cellsnp-lite pileup
#' @param outdir character Output directory
#' @param het_only logical Whether to only use heterozygous SNPs
#' @param chr_prefix logical Whether to add chr prefix
#' @return integer Status code
#' @keywords internal
genotype = function(label, samples, vcfs, outdir, het_only = FALSE, chr_prefix = TRUE) {

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
        make_vcf_chr(chr, snps, vcf_original, label, outdir, het_only = het_only, chr_prefix = chr_prefix)
    }

    return(0)
}

#' process VCFs into SNP dataframe
#' @param vcf vcfR object
#' @return dataframe SNP information
#' @keywords internal 
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
#' @keywords internal 
make_vcf_chr = function(chr, snps, vcf_original, label, outdir, het_only = FALSE, chr_prefix = TRUE) {
    
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
        distinct(snp_id, `.keep_all` = TRUE)

    if (nrow(chr_snps) == 0) {
        message(glue('No SNPs left for chr{chr}!'))
        return(NULL)
    }

    if (het_only) {
        chr_snps = chr_snps %>% filter(het)
    }
    
    vcf_chr = vcf_original[vcf_original@fix[,1] == chr, ]

    vcf_chr_fix = data.frame(vcf_chr@fix) %>%
        mutate(POS = as.integer(POS)) %>%
        mutate(snp_id = paste(CHROM, POS, REF, ALT, sep = '_'))
    
    # sometimes the vcf from cell-snp-lite has duplicated records
    vcf_chr = vcf_chr[vcf_chr_fix$snp_id %in% chr_snps$snp_id & (!duplicated(vcf_chr_fix$snp_id)), ]

    vcf_chr@meta = vcf_meta

    # make sure they are sorted the same way
    vcf_chr = vcf_chr[order(as.integer(vcf_chr@fix[,2])),]
    chr_snps = chr_snps %>% arrange(as.integer(POS))

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
    system(paste0('gunzip ', file_name))
    system(paste0('bgzip ', str_remove(file_name, '.gz')))

    # index the vcf
    system(paste0('tabix ', file_name))
    
    return(vcf_chr)
    
}

#' Preprocess allele data
#' 
#' @param sample character Sample label
#' @param vcf_pu dataframe Pileup VCF from cell-snp-lite
#' @param vcf_phased dataframe Phased VCF from eagle2
#' @param AD dgTMatrix Alt allele depth matrix from pileup
#' @param DP dgTMatrix Total allele depth matrix from pileup
#' @param barcodes vector List of barcodes from pileup
#' @param gtf dataframe Transcript GTF
#' @param gmap dataframe Genetic map
#' @return dataframe Allele counts by cell
#' @keywords internal
preprocess_allele = function(
    sample,
    vcf_pu,
    vcf_phased,
    AD,
    DP,
    barcodes,
    gtf,
    gmap
) {

    # pileup VCF
    vcf_pu = vcf_pu %>%
        mutate(INFO = str_remove_all(INFO, '[:alpha:]|=')) %>%
        tidyr::separate(col = 'INFO', into = c('AD', 'DP', 'OTH'), sep = ';') %>%
        mutate_at(c('AD', 'DP', 'OTH'), as.integer) %>%
        mutate(snp_id = paste(CHROM, POS, REF, ALT, sep = '_'))

    # pileup count matrices
    DP = as.data.frame(Matrix::summary(DP)) %>%
        mutate(
            cell = barcodes[j],
            snp_id = vcf_pu$snp_id[i]
        ) %>%
        select(-i, -j) %>%
        rename(DP = x) %>%
        select(cell, snp_id, DP)

    AD = as.data.frame(Matrix::summary(AD)) %>%
        mutate(
            cell = barcodes[j],
            snp_id = vcf_pu$snp_id[i]
        ) %>%
        select(-i, -j) %>%
        rename(AD = x) %>%
        select(cell, snp_id, AD)

    df = DP %>% left_join(AD, by = c("cell", "snp_id")) %>%
        mutate(AD = ifelse(is.na(AD), 0, AD))

    df = df %>% left_join(
        vcf_pu %>% rename(AD_all = AD, DP_all = DP, OTH_all = OTH),
        by = 'snp_id')

    df = df %>% mutate(
            AR = AD/DP,
            AR_all = AD_all/DP_all
        )

    df = df %>% dplyr::filter(DP_all > 1 & OTH_all == 0)

    # vcf has duplicated records sometimes
    df = df %>% distinct() 

    df = df %>% mutate(
        snp_index = as.integer(factor(snp_id, unique(snp_id))),
        cell_index = as.integer(factor(cell, sample(unique(cell))))
    )
    
    # phased VCF
    vcf_phased = vcf_phased %>% 
        mutate(snp_id = paste(CHROM, POS, REF, ALT, sep = '_')) %>%
        mutate(GT = get(sample))

    # annotate SNP by gene
    overlap_transcript = GenomicRanges::findOverlaps(
            vcf_phased %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$POS,
                    end = .$POS)
            )},
            gtf %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$gene_start,
                    end = .$gene_end)
            )}
        ) %>%
        as.data.frame() %>%
        setNames(c('snp_index', 'gene_index')) %>%
        left_join(
            vcf_phased %>% mutate(snp_index = 1:n()) %>%
                select(snp_index, snp_id),
            by = c('snp_index')
        ) %>%
        left_join(
            gtf %>% mutate(gene_index = 1:n()),
            by = c('gene_index')
        ) %>%
        arrange(snp_index, gene) %>%
        distinct(snp_index, `.keep_all` = TRUE)

    vcf_phased = vcf_phased %>%
        left_join(
            overlap_transcript %>% select(snp_id, gene, gene_start, gene_end),
            by = c('snp_id')
        )

    # annotate SNPs by genetic map
    marker_map = GenomicRanges::findOverlaps(
            vcf_phased %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$POS,
                    end = .$POS)
            )},
            gmap %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$start,
                    end = .$end)
            )}
        ) %>%
        as.data.frame() %>%
        setNames(c('marker_index', 'map_index')) %>%
        left_join(
            vcf_phased %>% mutate(marker_index = 1:n()) %>%
                select(marker_index, snp_id),
            by = c('marker_index')
        ) %>%
        left_join(
            gmap %>% mutate(map_index = 1:n()),
            by = c('map_index')
        ) %>%
        arrange(marker_index, -start) %>%
        distinct(marker_index, `.keep_all` = TRUE) %>%
        select(snp_id, cM)

    vcf_phased = vcf_phased %>% 
        left_join(marker_map, by = 'snp_id')

    # add annotation to cell counts
    df = df %>% left_join(vcf_phased %>% select(snp_id, gene, GT, cM), by = 'snp_id')
    df = df %>% mutate(CHROM = factor(CHROM, unique(CHROM)))
    df = df %>% select(cell, snp_id, CHROM, POS, cM, REF, ALT, AD, DP, GT, gene)

    # only keep hets
    df = df %>% filter(GT %in% c('1|0', '0|1'))
    
    return(df)
}


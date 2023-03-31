log_mem = function() {
    m = pryr::mem_used()
    msg = paste0('Mem used: ', signif(m/1e9, 3), 'Gb')
    log_message(msg)
}

log_message = function(msg, verbose = TRUE) {
    log_info(msg)
    if (verbose) {
        message(msg)
    }
}

#' check the format of a count matrix
#' @keywords internal
check_matrix = function(count_mat) {

    if ('dgCMatrix' %in% class(count_mat)) {
        return(count_mat)
    } else if ('matrix' %in% class(count_mat)) {
        count_mat <- as(Matrix(count_mat, sparse=TRUE), "dgCMatrix")
        return(count_mat)
    } else {

        if (!('dgCMatrix' %in% class(count_mat))) {
            msg = "count_mat is not of class dgCMatrix or matrix"
        } else if (!is.numeric(count_mat@x)) {
            msg = "The parameter 'count_mat' must be of type 'integer'. Please fix."
        } else if (all(count_mat@x != as.integer(count_mat@x))) {
            msg = "The parameter 'count_mat' must be of type 'integer'. Please fix."
        } else if (any(duplicated(rownames(count_mat)))) {
            msg = "Please remove duplicated genes in count matrix"
        }

        log_error(msg)
        stop(msg)
    }
}

#' check the format of a allele dataframe
#' @keywords internal
check_allele_df = function(df) {

    expected_colnames = c('cell', 'snp_id', 'CHROM', 'POS', 'cM', 'REF', 'ALT', 'AD', 'DP', 'GT', 'gene')

    potential_missing_columns = return_missing_columns(df, expected_colnames)
    
    if (!is.null(potential_missing_columns)) {
        stop(paste0("The allele count dataframe appears to be malformed; expected column names: ", potential_missing_columns, ". Please fix."))
    }

    snps = df %>% 
        filter(GT != '') %>% 
        group_by(snp_id) %>%
        summarise(
            n = length(unique(GT))
        )
    
    if (any(snps$n > 1)) {
        msg = 'Inconsistent SNP genotypes; Are cells from two different individuals mixed together?'
        log_error(msg)
        stop(msg)
    }

    df = df %>% mutate(CHROM = factor(CHROM, 1:22))

    return(df)

}

#' Annotate genes on allele dataframe
#' @param df dataframe Allele count dataframe 
#' @param gtf dataframe Gene gtf
#' @return dataframe Allele dataframe with gene column 
#' @export
annotate_genes = function(df, gtf) {

    snps = df %>% distinct(snp_id, CHROM, POS)

    hits = GenomicRanges::findOverlaps(
            snps %>% {GenomicRanges::GRanges(
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
            snps %>% mutate(snp_index = 1:n()) %>%
                select(snp_index, snp_id),
            by = c('snp_index')
        ) %>%
        left_join(
            gtf %>% mutate(gene_index = 1:n()),
            by = c('gene_index')
        ) %>%
        arrange(snp_index, gene) %>%
        distinct(snp_index, `.keep_all` = TRUE)

    snps = snps %>%
        left_join(
            hits %>% select(snp_id, gene),
            by = c('snp_id')
        )
    
    df = df %>% select(-any_of(c('gene', 'gene_start', 'gene_end'))) %>%
        left_join(snps, by = c('snp_id', 'CHROM', 'POS'))

    return(df)
}

#' check the format of lambdas_ref
#' @keywords internal
check_exp_ref = function(lambdas_ref) {

    if (!is.matrix(lambdas_ref)) {
        lambdas_ref = as.matrix(lambdas_ref) %>% magrittr::set_colnames('ref')
    }

    return(lambdas_ref)

}

#' check inter-individual contamination
#' @param bulk dataframe Pseudobulk profile
#' @keywords internal
check_contam = function(bulk) {

    hom_rate = bulk %>% filter(DP >= 8) %>%
        {mean(na.omit(.$AR == 0 | .$AR == 1))}

    if (hom_rate > 0.4) {
        msg = paste0(
            'High SNP contamination detected ',
            '(', round(hom_rate*100, 1), '%)',
            '. Please make sure that cells from only one individual are included in genotyping step.')
        message(msg)
        log_warn(msg)
    }

}

#' check noise level
#' @param bulk dataframe Pseudobulk profile
#' @keywords internal
check_exp_noise = function(bulk) {

    sig = unique(na.omit(bulk$sig))

    if (sig > 1) {
        noise_level  = 'high'
        noise_msg = 'Consider using a custom expression reference profile.'
    } else if (sig > 0.5) {
        noise_level = 'medium'
        noise_msg = ''
    } else {
        noise_level = 'low'
        noise_msg = ''
    }

    msg = paste0(
        'Expression noise level: ',
        noise_level,
        ' (', signif(sig, 2), '). ',
        noise_msg)

    # message(msg)
    log_message(msg)
}

# check the format of a given consensus segment dataframe
check_segs_fix = function(segs_consensus_fix) {

    if (is.null(segs_consensus_fix)) {
        return(NULL)
    }

    if (!all(c('cnv_state', 'seg', 'seg_start', 'seg_end') %in% colnames(segs_consensus_fix))) {
        stop('The consensus segment dataframe appears to be malformed. Please fix.')
    }

    segs_consensus_fix = segs_consensus_fix %>% 
        mutate(
            cnv_state_post = cnv_state,
            seg_cons = seg,
            p_amp = ifelse(cnv_state == 'amp', 1, 0),
            p_del = ifelse(cnv_state == 'del', 1, 0),
            p_loh = ifelse(cnv_state == 'loh', 1, 0),
            p_bamp = ifelse(cnv_state == 'bamp', 1, 0),
            p_bdel = ifelse(cnv_state == 'bdel', 1, 0),
            seg_length = seg_end - seg_start,
            LLR = Inf
        ) %>%
        relevel_chrom() %>%
        as.data.frame()

    return(segs_consensus_fix)
    
}

#' @keywords internal
return_missing_columns = function(file, expected_colnames = NULL) {
    ## if user sets expected_colnames = NULL, return NULL
    if (is.null(expected_colnames)) {
        return(NULL)
    }
    if (!is.vector(expected_colnames) || !is.character(expected_colnames)) {
        stop("The parameter 'expected_colnames' needs to be a character vector")
    }
    '%ni%' <- Negate('%in%')
    if (any(expected_colnames %ni% colnames(file))) {
        missing_columns = expected_colnames[!(expected_colnames %in% colnames(file))]
        if (length(missing_columns) == 0) {
            stop("Some mismatch exists between the expected columns and the columns in the file. This error shouldn't happen. Check and fix.")
        }
        return(missing_columns)
    } else {
        return(NULL)
    }
}


relevel_chrom = function(df) {
    if (!is.null(df)) {
        df = df %>% mutate(CHROM = factor(CHROM, 1:22))
    }
    return(df)
}

#' @keywords internal
check_fread_works = function(input) {
    tryCatch({
        return(data.table::fread(input))
    },
    error = function(e){
        message(paste0("Could not read the input file ", input, " with data.table::fread(). Please check that the file is valid."))
        return(NULL)
    })
}

#' @keywords internal
check_rds_works = function(input) {
    tryCatch({
        return(readRDS(input))
    },
    error = function(e){
        message(paste0("Could not read the input file ", input, " with readRDS(). Please check that the file is valid."))
        return(NULL)
    })
}


#' @keywords internal
read_file = function(inputfile, expected_colnames = NULL, filetype="tsv") {
    if (filetype == "tsv") {
        file = check_fread_works(inputfile)
    } else if (filetype == "rds") {
        file = check_rds_works(inputfile)       
    } else {
        stop("The parameter 'filetype' must be either 'tsv' or 'rds'. Please fix.")
    }
    potential_missing_columns = return_missing_columns(file, expected_colnames)
    if (!is.null(potential_missing_columns)) {
        stop(paste0("The file ", inputfile, " appears to be malformed; expected column names: ", potential_missing_columns, ". Please fix."))
    } else {
        return(file)
    }
}

#' @keywords internal
read_hc_rds = function(inputfile) {
    file = check_rds_works(inputfile)
    if (!is.list(file)) {
        stop(paste0("The file: ", inputfile, " is malformed; should be a list. Please fix."))
    }        
    hc_colnames = c("merge", "height", "order", "labels", "method", "call", "dist.method")
    '%ni%' <- Negate('%in%')
    if (any(hc_colnames %ni% names(file))) {
        missing_columns = expected_colnames[!(expected_colnames %in% names(file))]
        stop(paste0("The file ", inputfile, " appears to be malformed; expected column names: ", potential_missing_columns, ". Please fix."))
    } else {
        return(file)
    }
}


            


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

    snps = df %>% 
        filter(GT != '') %>% 
        group_by(snp_id) %>%
        summarise(
            n = length(unique(GT))
        )
    
    if (any(snps$n > 1)) {
        log_err('Inconsistent SNP genotypes; Are cells from two different individuals mixed together?')
    }

    df = df %>% mutate(CHROM = factor(CHROM, 1:22))

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

    if (hom_rate > 0.15) {
        msg = paste0(
            'High SNP contamination detected ',
            '(', round(hom_rate*100, 1), '%)',
            '. Are cells from only one individual included in genotyping step?')
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
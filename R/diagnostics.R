
#' Log memory usage
#' @return NULL
#' @keywords internal
log_mem = function() {
    m = pryr::mem_used()
    msg = paste0('Mem used: ', signif(m/1e9, 3), 'Gb')
    log_message(msg)
}

#' Log a message
#' @param msg string Message to log
#' @param verbose boolean Whether to print message to console
#' @return NULL
#' @keywords internal
log_message = function(msg, verbose = TRUE) {
    log_info(msg)
    if (verbose) {
        message(msg)
    }
}

#' Check the format of a count matrix
#' @param count_mat matrix Count matrix
#' @return matrix Count matrix
#' @keywords internal
check_matrix = function(count_mat) {

    # Make sure that the count matrix is of class dgCMatrix
    if ('matrix' %in% class(count_mat)) {
        count_mat <- as(Matrix(count_mat, sparse=TRUE), "dgCMatrix")
    } else if (!('dgCMatrix' %in% class(count_mat))) {
        msg = "count_mat should be of class dgCMatrix or matrix"
        log_error(msg)
        stop(msg)
    }

    # Make sure that the count matrix is of type integer
    if (!is.numeric(count_mat@x)) {
        msg = "The parameter 'count_mat' should be of type 'integer'. Please fix."
        log_error(msg)
        stop(msg)
    } else if (all(count_mat@x != as.integer(count_mat@x))) {
        msg = "The parameter 'count_mat' should be of type 'integer'. Please fix."
        log_error(msg)
        stop(msg)
    } else if (any(duplicated(rownames(count_mat)))) {
        msg = "Please remove duplicated genes in count matrix"
        log_error(msg)
        stop(msg)
    }

    return(count_mat)
}

#' Check the format of a allele dataframe
#' @param df dataframe Allele dataframe
#' @return dataframe Allele dataframe
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
    
    # check chrom prefix 
    if (any(str_detect(df$CHROM[1], '^chr'))) {
        df = df %>% mutate(CHROM = str_remove(CHROM, 'chr'))
    } 

    df = df %>% 
        filter(CHROM %in% 1:22) %>%
        mutate(CHROM = factor(CHROM, 1:22))

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
#' @param lambdas_ref matrix Expression reference profile
#' @return matrix Expression reference profile
#' @keywords internal
check_exp_ref = function(lambdas_ref) {

    if (!is.matrix(lambdas_ref)) {
        lambdas_ref = as.matrix(lambdas_ref) %>% magrittr::set_colnames('ref')
    }

    if (any(is.na(lambdas_ref))) {
        msg = "The reference expression matrix 'lambdas_ref' should not contain any NA values."
        log_error(msg)
        stop(msg)
    }
    
    # check if all entries in the reference profile are integers
    if (all(lambdas_ref == as.integer(lambdas_ref))) {
        msg = "The reference expression matrix 'lambdas_ref' should be normalized gene expression magnitudes. Please use aggregate_counts() function to prepare the reference profile from raw counts."
        log_error(msg)
        stop(msg)
    } else if (any(duplicated(rownames(lambdas_ref)))) {
        msg = "Please remove duplicated genes in reference profile"
        log_error(msg)
        stop(msg)
    }


    return(lambdas_ref)

}



#' check inter-individual contamination
#' @param bulk dataframe Pseudobulk profile
#' @return NULL
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
#' @return NULL
#' @keywords internal
check_exp_noise = function(bulk) {

    mse = unique(na.omit(bulk$mse))

    if (mse > 1.5) {
        noise_level  = 'high'
        noise_msg = 'Consider using a custom expression reference profile.'
    } else if (mse > 0.5) {
        noise_level = 'medium'
        noise_msg = ''
    } else {
        noise_level = 'low'
        noise_msg = ''
    }

    msg = paste0(
        'Expression noise level (MSE): ',
        noise_level,
        ' (', signif(mse, 2), '). ',
        noise_msg)

    # message(msg)
    log_message(msg)
}

#' Check the format of a given clonal LOH segment dataframe
#' @param segs_loh dataframe Clonal LOH segment dataframe
#' @return dataframe Clonal LOH segment dataframe
#' @keywords internal
check_segs_loh = function(segs_loh) {
    
        if (is.null(segs_loh)) {
            return(NULL)
        }
    
        if (!all(c('CHROM', 'seg', 'seg_start', 'seg_end') %in% colnames(segs_loh))) {
            stop('The clonal LOH segment dataframe appears to be malformed. Please fix.')
        }

        if (is.integer(segs_loh$seg)) {
            segs_loh = segs_loh %>% mutate(seg = paste0(CHROM, '_', seg))
        }

        segs_loh = segs_loh %>% 
            mutate(loh = TRUE) %>%
            relevel_chrom() %>%
            arrange(CHROM, seg_start)
    
        return(segs_loh)

}

#' check the format of a given consensus segment dataframe
#' @param segs_consensus_fix dataframe Consensus segment dataframe
#' @return dataframe Consensus segment dataframe
#' @keywords internal
check_segs_fix = function(segs_consensus_fix) {

    if (is.null(segs_consensus_fix)) {
        return(NULL)
    }

    if (!all(c('CHROM', 'seg', 'seg_start', 'seg_end', 'cnv_state') %in% colnames(segs_consensus_fix))) {
        stop('The consensus segment dataframe appears to be malformed. Please fix.')
    }

    segs_consensus_fix = segs_consensus_fix %>% 
        relevel_chrom() %>%
        arrange(CHROM, seg_start)

    if (is.integer(segs_consensus_fix$seg)) {
        segs_consensus_fix = segs_consensus_fix %>% mutate(seg = paste0(CHROM, '_', seg))
    }

    segs_consensus_fix = segs_consensus_fix %>%
        arrange(CHROM) %>%
        mutate(
            cnv_state_post = cnv_state,
            seg_cons = seg,
            p_amp = ifelse(cnv_state == 'amp', 1, 0),
            p_del = ifelse(cnv_state == 'del', 1, 0),
            p_loh = ifelse(cnv_state == 'loh', 1, 0),
            p_bamp = ifelse(cnv_state == 'bamp', 1, 0),
            p_bdel = ifelse(cnv_state == 'bdel', 1, 0),
            seg_length = seg_end - seg_start,
            LLR = ifelse(cnv_state == 'neu', NA, Inf)
        ) %>%
        as.data.frame()

    return(segs_consensus_fix)
    
}

#' Check the format of a given file
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

#' Relevel chromosome column
#' @param df dataframe Dataframe with chromosome column
#' @keywords internal 
relevel_chrom = function(df) {
    if (!is.null(df)) {
        df = df %>% mutate(CHROM = factor(CHROM, 1:22))
    }
    return(df)
}

#' @keywords internal
check_fread_works = function(input, ...) {
    tryCatch({
        return(data.table::fread(input, ...))
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
read_file = function(inputfile, expected_colnames = NULL, filetype="tsv", ...) {
    if (filetype == "tsv") {
        file = check_fread_works(inputfile, ...)
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


            


########################### Data processing ############################

#' choose beest reference for each cell based on correlation
#' @param count_mat dgCMatrix Gene expression counts
#' @param lambdas_ref matrix Reference expression profiles
#' @param gtf dataframe Transcript gtf
#' @return named vector Best references for each cell
#' @keywords internal
choose_ref_cor = function(count_mat, lambdas_ref, gtf) {

    if (any(duplicated(rownames(lambdas_ref)))) {
        log_err('Duplicated genes in lambdas_ref')
    }

    if (ncol(lambdas_ref) == 1) {
        ref_name = colnames(lambdas_ref)
        cells = colnames(count_mat)
        best_refs = setNames(rep(ref_name, length(cells)), cells)
        return(best_refs)
    }
    
    genes_annotated = gtf %>% 
        pull(gene) %>% 
        intersect(rownames(count_mat)) %>%
        intersect(rownames(lambdas_ref))

    count_mat = count_mat[genes_annotated,,drop=FALSE]
    lambdas_ref = lambdas_ref[genes_annotated,,drop=FALSE]
    
    exp_mat = scale_counts(count_mat)
    # keep highly expressed genes in at least one of the references
    exp_mat = exp_mat[rowSums(lambdas_ref * 1e6 > 2) > 0,,drop=FALSE]
    
    cors = cor(as.matrix(log(exp_mat * 1e6 + 1)), log(lambdas_ref * 1e6 + 1)[rownames(exp_mat),])
    best_refs = apply(cors, 1, function(x) {colnames(cors)[which.max(x)]})
    
    return(best_refs)
}

scale_counts = function(count_mat) {
    count_mat@x <- count_mat@x/rep.int(colSums(count_mat), diff(count_mat@p))
    return(count_mat)
}

#' Utility function to make reference gene expression profiles
#' @param count_mat matrix/dgCMatrix Gene expression counts
#' @param annot dataframe Cell annotation with columns "cell" and "group"
#' @param normalized logical Whether to return normalized expression values
#' @param verbose logical Verbosity
#' @return matrix Reference gene expression levels
#' @export
aggregate_counts = function(count_mat, annot, normalized = TRUE, verbose = TRUE) {

    annot = annot %>% mutate(group = factor(group))

    cells = intersect(colnames(count_mat), annot$cell)
    
    count_mat = count_mat[,cells]

    annot = annot %>% filter(cell %in% cells)
    
    cell_dict = annot %>% {setNames(.$group, .$cell)} %>% droplevels

    cell_dict = cell_dict[cells]

    if (verbose) {
        print(table(cell_dict))
    }
    
    if (length(levels(cell_dict)) == 1) {
        count_mat_clust = count_mat %>% rowSums() %>% as.matrix %>% magrittr::set_colnames(levels(cell_dict))
        exp_mat_clust = count_mat_clust/sum(count_mat_clust)
    } else {
        M = model.matrix(~ 0 + cell_dict) %>% magrittr::set_colnames(levels(cell_dict))
        count_mat_clust = count_mat %*% M
        exp_mat_clust = count_mat_clust %*% diag(1/colSums(count_mat_clust)) %>% magrittr::set_colnames(colnames(count_mat_clust))
    }

    if (normalized) {
        return(exp_mat_clust)
    } else {
        return(count_mat_clust)
    }

}

#' filtering, normalization and capping
#' @param count_mat dgCMatrix Gene expression counts
#' @param lambdas_ref matrix Reference expression profiles
#' @param gtf dataframe Transcript gtf
#' @return dataframe Log(x+1) transformed normalized expression values for single cells
#' @keywords internal
smooth_expression = function(count_mat, lambdas_ref, gtf, window = 101, verbose = FALSE) {

    mut_expressed = filter_genes(count_mat, lambdas_ref, gtf, verbose = verbose)
    count_mat = count_mat[mut_expressed,,drop=FALSE]
    lambdas_ref = lambdas_ref[mut_expressed]

    exp_mat = scale_counts(count_mat)

    exp_mat_norm = log2(exp_mat * 1e6 + 1) - log2(lambdas_ref * 1e6 + 1)

    # capping 
    cap = 3
    exp_mat_norm[exp_mat_norm > cap] = cap
    exp_mat_norm[exp_mat_norm < -cap] = -cap

    # centering by cell
    row_mean = apply(exp_mat_norm, 2, function(x) { mean(x, na.rm=TRUE) })
    exp_mat_norm = t(apply(exp_mat_norm, 1, "-", row_mean))

    # smoothing
    exp_mat_smooth = caTools::runmean(exp_mat_norm, k = window, align = "center")
    rownames(exp_mat_smooth) = rownames(exp_mat_norm)
    colnames(exp_mat_smooth) = colnames(exp_mat_norm)
    
    return(exp_mat_smooth)
}


#' filter for mutually expressed genes
#' @param count_mat dgCMatrix Gene expression counts
#' @param lambdas_ref named numeric vector A reference expression profile
#' @param gtf dataframe Transcript gtf
#' @return vector Genes that are kept after filtering
#' @keywords internal
filter_genes = function(count_mat, lambdas_ref, gtf, verbose = FALSE) {

    genes_keep = gtf$gene %>% 
        intersect(rownames(count_mat)) %>%
        intersect(names(lambdas_ref))

    genes_exclude = gtf %>%
        filter(CHROM == 6 & gene_start < 33480577 & gene_end > 28510120) %>%
        pull(gene)

    genes_keep = genes_keep[!genes_keep %in% genes_exclude]

    count_mat = count_mat[genes_keep,,drop=FALSE]
    lambdas_ref = lambdas_ref[genes_keep]
    lambdas_obs = rowSums(count_mat)/sum(count_mat)

    min_both = 2

    mut_expressed = ((lambdas_ref * 1e6 > min_both & lambdas_obs * 1e6 > min_both) |
        (lambdas_ref > mean(lambdas_ref[lambdas_ref != 0])) |
        (lambdas_obs > mean(lambdas_obs[lambdas_obs != 0]))) &
        (lambdas_ref > 0)

    retained = names(mut_expressed)[mut_expressed]

    if (verbose) {
        message(glue('number of genes left: {length(retained)}'))
    }

    return(retained)
}

#' Aggregate into bulk expression profile
#' @param count_mat dgCMatrix Gene expression counts
#' @param lambdas_ref matrix Reference expression profiles
#' @param gtf dataframe Transcript gtf
#' @return dataframe Pseudobulk gene expression profile
#' @keywords internal
get_exp_bulk = function(count_mat, lambdas_ref, gtf, verbose = FALSE) {

    depth_obs = sum(count_mat)
    
    mut_expressed = filter_genes(count_mat, lambdas_ref, gtf)
    count_mat = count_mat[mut_expressed,,drop=FALSE]
    lambdas_ref = lambdas_ref[mut_expressed]

    bulk_obs = count_mat %>%
        rowSums() %>%
        data.frame() %>%
        setNames('Y_obs') %>%
        tibble::rownames_to_column('gene') %>%
        mutate(lambda_obs = (Y_obs/depth_obs)) %>%
        mutate(lambda_ref = lambdas_ref[gene]) %>%
        mutate(d_obs = depth_obs) %>%
        left_join(gtf, by = "gene") 
    
    # annotate using GTF
    bulk_obs = bulk_obs %>%
        mutate(gene = droplevels(factor(gene, gtf$gene))) %>%
        mutate(gene_index = as.integer(gene)) %>%
        arrange(gene) %>%
        mutate(CHROM = factor(CHROM)) %>%
        mutate(
            logFC = log2(lambda_obs) - log2(lambda_ref),
            lnFC = log(lambda_obs) - log(lambda_ref),
            logFC = ifelse(is.infinite(logFC), NA, logFC),
            lnFC = ifelse(is.infinite(lnFC), NA, lnFC)
        )
    
    return(bulk_obs)
}


#' Aggregate into pseudobulk alelle profile
#' @param df_allele dataframe Single-cell allele counts
#' @param genetic_map dataframe Genetic map
#' @param lambda numeric Phase switch rate
#' @param min_depth integer Minimum coverage to filter SNPs
#' @return dataframe Pseudobulk allele profile
#' @keywords internal
get_allele_bulk = function(df_allele, genetic_map, lambda = 1, min_depth = 0) {
    df_allele %>%
        filter(GT %in% c('1|0', '0|1')) %>%
        group_by(snp_id, CHROM, POS, REF, ALT, GT, gene) %>%
        summarise(
            AD = sum(AD),
            DP = sum(DP),
            AR = AD/DP,
            .groups = 'drop'
        ) %>%
        arrange(CHROM, POS) %>%
        group_by(CHROM) %>%
        mutate(snp_index = as.integer(factor(snp_id, unique(snp_id)))) %>%
        ungroup() %>%
        filter(DP >= min_depth) %>%
        mutate(pBAF = ifelse(GT == '1|0', AR, 1-AR)) %>%
        mutate(pAD = ifelse(GT == '1|0', AD, DP - AD)) %>%
        mutate(CHROM = factor(CHROM, unique(CHROM))) %>%
        arrange(CHROM, POS) %>%
        annot_cm(genetic_map) %>%
        group_by(CHROM) %>%
        filter(n() > 1) %>%
        mutate(
            inter_snp_cm = c(NA, cM[2:length(cM)] - cM[1:(length(cM)-1)]),
            p_s = switch_prob_cm(inter_snp_cm, lambda = lambda)
        ) %>%
        ungroup()
}


#' Combine allele and expression pseudobulks
#' @param allele_bulk dataframe Bulk allele profile
#' @param genetic_map dataframe Genetic map
#' @return dataframe Pseudobulk allele and expression profile
#' @keywords internal
combine_bulk = function(allele_bulk, exp_bulk) {
    
    bulk = allele_bulk %>% 
        full_join(
            exp_bulk,
            by = c("CHROM", "gene")
        ) %>%
        mutate(
            snp_id = ifelse(is.na(snp_id), gene, snp_id),
            # levels will be missing if not expressed
            gene = factor(gene, levels(exp_bulk$gene)),
            POS = ifelse(is.na(POS), gene_start, POS),
            # phase switch is forbidden if not heteroSNP
            p_s = ifelse(is.na(p_s), 0, p_s)
        ) %>%
        arrange(CHROM, POS) %>%
        filter(!(CHROM == 6 & POS < 33480577 & POS > 28510120)) %>%
        group_by(CHROM) %>%
        mutate(snp_index = as.integer(factor(snp_id, unique(snp_id)))) %>%
        ungroup()

    # get rid of duplicate gene expression values
    bulk = bulk %>% 
        group_by(CHROM, gene) %>%
        mutate(
            Y_obs = ifelse(
                !is.na(gene) & n() > 1,
                c(unique(Y_obs), rep(NA, n()-1)), Y_obs
            )
        ) %>%
        ungroup() %>% 
        mutate(
            lambda_obs = Y_obs/d_obs,
            logFC = log2(lambda_obs/lambda_ref),
            lnFC = log(lambda_obs/lambda_ref)
        ) %>%
        mutate_at(
            c('logFC', 'lnFC'),
            function(x) ifelse(is.infinite(x), NA, x)
        )
    
    return(bulk)
    
}

#' Get average reference expressio profile based on single-cell ref choices
#' @keywords internal
get_lambdas_bar = function(lambdas_ref, sc_refs, verbose = TRUE) {

    w = sapply(colnames(lambdas_ref), function(ref){sum(sc_refs == ref)})
    w = w/sum(w)
    lambdas_bar = lambdas_ref %*% w %>% {setNames(as.vector(.), rownames(.))}
    
    if (verbose) {
        log_message('Fitted reference proportions:')
        log_message(paste0(paste0(names(w), ':', signif(w, 2)), collapse = ','))
    }

    return(lambdas_bar)
}

#' Aggregate single-cell data into combined bulk expression and allele profile
#'
#' @param count_mat dgCMatrix Gene expression counts
#' @param lambdas_ref matrix Reference expression profiles
#' @param df_allele dataframe Single-cell allele counts
#' @param gtf dataframe Transcript gtf
#' @param genetic_map dataframe Genetic map
#' @return dataframe Pseudobulk gene expression and allele profile
#' @export
get_bulk = function(count_mat, lambdas_ref, df_allele, gtf, genetic_map, min_depth = 0, lambda = 1, verbose = TRUE) {

    count_mat = check_matrix(count_mat)

    fit = fit_ref_sse(rowSums(count_mat), lambdas_ref, gtf)

    exp_bulk = get_exp_bulk(
            count_mat,
            fit$lambdas_bar,
            gtf,
            verbose = verbose
        ) %>%
        filter((logFC < 5 & logFC > -5) | Y_obs == 0) %>%
        mutate(sse = fit$sse)

    allele_bulk = get_allele_bulk(
        df_allele,
        genetic_map,
        lambda = lambda,
        min_depth = min_depth)
            
    bulk = combine_bulk(
        allele_bulk = allele_bulk,
        exp_bulk = exp_bulk)

    if (any(duplicated(bulk$snp_id))) {
        stop('duplicated SNPs, check genotypes')
    }

    # doesn't work with 0s in the ref
    bulk = bulk %>% filter(lambda_ref != 0 | is.na(gene))

    bulk = bulk %>%
        mutate(CHROM = as.character(CHROM)) %>%
        mutate(CHROM = ifelse(CHROM == 'X', 23, CHROM)) %>%
        mutate(CHROM = factor(as.integer(CHROM)))
}

#' Fit a reference profile from multiple references using constrained least square
#' @param Y_obs vector 
#' @param lambdas_ref named vector 
#' @param gtf dataframe 
#' @return fitted expression profile
#' @keywords internal
fit_ref_sse = function(Y_obs, lambdas_ref, gtf, min_lambda = 2e-6, verbose = FALSE) {

    if (any(duplicated(rownames(lambdas_ref)))) {
        log_err('Duplicated genes in lambdas_ref')
    }

    if (length(dim(lambdas_ref)) == 1 | is.null(dim(lambdas_ref))) {
        return(list('w' = 1, 'lambdas_bar' = lambdas_ref))
    }
    
    Y_obs = Y_obs[Y_obs > 0]

    # take the union of expressed genes across cell type
    genes_common = gtf$gene %>% 
        intersect(names(Y_obs)) %>%
        intersect(rownames(lambdas_ref)[rowMeans(lambdas_ref) > min_lambda])

    if (verbose) {
        log_info(glue('{length(genes_common)} genes common in reference and observation'))
    }

    Y_obs = Y_obs[genes_common]
    lambdas_obs = Y_obs/sum(Y_obs)
    lambdas_ref = lambdas_ref[genes_common,,drop=FALSE]

    n_ref = ncol(lambdas_ref)
    
    fit = optim(
        fn = function(w) {
            w = w/sum(w)
            sum(log(lambdas_obs/as.vector(lambdas_ref %*% w))^2)
        },
        method = 'L-BFGS-B',
        par = rep(1/n_ref, n_ref),
        lower = rep(1e-6, n_ref)
    )

    w = fit$par
    w = w/sum(w)
    w = setNames(w, colnames(lambdas_ref))

    lambdas_bar = lambdas_ref %*% w %>% {setNames(as.vector(.), rownames(.))}

    return(list('w' = w, 'lambdas_bar' = lambdas_bar, 'sse' = fit$value/length(Y_obs)))
}

#' Annotate genetic distance between markers
#' @param bulk dataframe Pseudobulk profile
#' @param genetic_map dataframe Genetic map
#' @return dataframe Annotated pseudobulk profile
#' @keywords internal
annot_cm = function(bulk, genetic_map) {

    bulk = bulk %>% ungroup()
    
    marker_map = GenomicRanges::findOverlaps(
            bulk %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$POS,
                    end = .$POS)
            )},
            genetic_map %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$start,
                    end = .$end)
            )}
        ) %>%
        as.data.frame() %>%
        setNames(c('marker_index', 'map_index')) %>%
        left_join(
            bulk %>% mutate(marker_index = 1:n()) %>%
                select(marker_index, snp_id),
            by = c('marker_index')
        ) %>%
        left_join(
            genetic_map %>% mutate(map_index = 1:n()),
            by = c('map_index')
        ) %>%
        arrange(marker_index, -start) %>%
        distinct(marker_index, `.keep_all` = T) %>%
        select(snp_id, cM)

    bulk = bulk %>% left_join(marker_map, by = 'snp_id') %>%
        filter(!is.na(cM))

    return(bulk)
    
}

#' predict phase switch probablity as a function of genetic distance
#' @param d numeric vector Genetic distance in cM
#' @param lambda numeric Phase switch rate
#' @param min_p numeric Minimum phase switch probability 
#' @return numeric vector Phase switch probability
#' @keywords internal
switch_prob_cm = function(d, lambda = 1, min_p = 1e-10) {
    p = (1-exp(-2*lambda*d))/2
    p = pmax(p, min_p)
    p = ifelse(is.na(d), 0, p)
    return(p)
}

########################### Analysis ############################

#' Call CNVs in a pseudobulk profile using the Numbat joint HMM
#' @param bulk dataframe Pesudobulk profile
#' @param t numeric Transition probability
#' @param gamma numeric Dispersion parameter for the Beta-Binomial allele model
#' @param theta_min numeric Minimum imbalance threshold
#' @param logphi_min numeric Minimum log expression deviation threshold
#' @param lambda numeric Phase switch rate
#' @param min_genes integer Minimum number of genes to call an event
#' @param exp_only logical Whether to run expression-only HMM
#' @param allele_only logical Whether to run allele-only HMM
#' @param bal_cnv logical Whether to call balanced amplifications/deletions
#' @param retest whether to retest CNVs after Viterbi decoding
#' @param find_diploid whether to run diploid region identification routine
#' @param diploid_chroms user-given chromsomes that are known to be in diploid state
#' @return a dataframe of segments with CNV posterior information
#' @export
analyze_bulk = function(
    bulk, t = 1e-5, gamma = 20, theta_min = 0.08, logphi_min = 0.25,
    lambda = 1, min_genes = 10,
    exp_only = FALSE, allele_only = FALSE, bal_cnv = TRUE, retest = TRUE, 
    find_diploid = TRUE, diploid_chroms = NULL, segs_loh = NULL,
    classify_allele = FALSE, run_hmm = TRUE, prior = NULL, exclude_neu = TRUE,
    phasing = TRUE, verbose = TRUE
) {
    
    if (!is.numeric(t)) {
        stop('transition probability is not numeric')
    }

    if ('gamma' %in% colnames(bulk)) {
        bulk = bulk %>% select(-gamma)
    }

    # update transition probablity
    bulk = bulk %>% mutate(p_s = switch_prob_cm(inter_snp_cm, lambda = UQ(lambda)))

    if (exp_only | allele_only) {
        bulk$diploid = TRUE
    } else if (!is.null(diploid_chroms)) {
        log_info(glue('Using diploid chromosomes given: {paste0(diploid_chroms, collapse = ",")}'))
        bulk = bulk %>% mutate(diploid = CHROM %in% diploid_chroms)
    } else if (find_diploid) {
        bulk = find_common_diploid(
            bulk, gamma = gamma, t = t, theta_min = theta_min, 
            min_genes = min_genes, fc_min = 2^logphi_min, segs_loh = segs_loh)
    } else if (!'diploid' %in% colnames(bulk)) {
        stop('Must define diploid region if not given')
    }

    if (is.null(segs_loh)) {
        bulk = bulk %>% mutate(loh = FALSE)
    }

    if (!allele_only) {
        # fit expression baseline
        fit = bulk %>%
            filter(!is.na(Y_obs)) %>%
            filter(logFC < 8 & logFC > -8) %>%
            filter(diploid) %>%
            {fit_lnpois(.$Y_obs, .$lambda_ref, unique(.$d_obs))}
            
        bulk = bulk %>% mutate(mu = fit@coef[1], sig = fit@coef[2])
    } else {
        bulk = bulk %>% mutate(mu = NA, sig = NA)
    }

    if (run_hmm) {
        bulk = bulk %>% 
            group_by(CHROM) %>%
            mutate(state = 
                run_joint_hmm(
                    pAD = pAD,
                    DP = DP, 
                    p_s = p_s,
                    Y_obs = Y_obs, 
                    lambda_ref = lambda_ref, 
                    d_total = na.omit(unique(d_obs)),
                    phi_amp = 2^(logphi_min),
                    phi_del = 2^(-logphi_min),
                    mu = mu,
                    sig = sig,
                    t = t,
                    gamma = unique(gamma),
                    theta_min = theta_min,
                    prior = prior,
                    bal_cnv = bal_cnv,
                    exp_only = exp_only,
                    allele_only = allele_only,
                    classify_allele = classify_allele,
                    phasing = phasing
                )
            ) %>% 
            mutate(state = ifelse(loh, 'del_2_up', state)) %>%
            mutate(cnv_state = str_remove(state, '_down|_up')) %>%
            annot_segs(var = 'cnv_state') %>%
            smooth_segs(min_genes = min_genes) %>%
            annot_segs(var = 'cnv_state')
    }
    
    # retest CNVs
    if (retest & (!exp_only)) {
        
        if (verbose) {
            message('Retesting CNVs..')
        }

        segs_post = retest_cnv(
            bulk, gamma = gamma, theta_min = theta_min,
            logphi_min = logphi_min, exclude_neu = exclude_neu,
            allele_only = allele_only)
        
        bulk = bulk %>% 
            select(-any_of(colnames(segs_post)[!colnames(segs_post) %in% c('seg', 'CHROM', 'seg_start', 'seg_end')])) %>%
            left_join(
                segs_post, by = c('seg', 'CHROM', 'seg_start', 'seg_end')
            ) %>%
            mutate(
                cnv_state_post = ifelse(is.na(cnv_state_post), 'neu', cnv_state_post),
                cnv_state = ifelse(is.na(cnv_state), 'neu', cnv_state)
            ) %>%
            mutate(state_post = ifelse(
                cnv_state_post %in% c('amp', 'del', 'loh') & (!cnv_state %in% c('bamp', 'bdel')),
                paste0(cnv_state_post, '_', str_extract(state, 'up_1|down_1|up_2|down_2|up|down|1_up|2_up|1_down|2_down')),
                cnv_state_post
            )) %>%
            mutate(state_post = str_remove(state_post, '_NA'))

        # note that this uses retest cnv states
        bulk = bulk %>% classify_alleles()

        bulk = annot_theta_roll(bulk)
        
    } else {
        bulk = bulk %>% mutate(state_post = state, cnv_state_post = cnv_state)
    }

    if (!allele_only) {

        # annotate phi MLE for all segments
        bulk = bulk %>%
            group_by(seg) %>%
            mutate(
                approx_phi_post(
                    Y_obs[!is.na(Y_obs)], lambda_ref[!is.na(Y_obs)], unique(na.omit(d_obs)),
                    mu = mu[!is.na(Y_obs)],
                    sig = sig[!is.na(Y_obs)]
                )
            ) %>%
            ungroup()

        # rolling phi estimates
        bulk = bulk %>% 
            select(-any_of('phi_mle_roll')) %>%
            left_join(
                bulk %>% 
                    group_by(CHROM) %>%
                    filter(!is.na(Y_obs)) %>%
                    filter(Y_obs > 0) %>%
                    mutate(
                        phi_mle_roll = phi_hat_roll(Y_obs, lambda_ref, d_obs, mu, sig, h = 50)
                    ) %>%
                    select(phi_mle_roll, CHROM, gene),
            by = c('CHROM', 'gene')
        ) %>%
        group_by(CHROM) %>%
        mutate(phi_mle_roll = zoo::na.locf(phi_mle_roll, na.rm=FALSE)) %>%
        ungroup()

    }

    # store these info here
    bulk$lambda = lambda 
    bulk$gamma = gamma

    return(bulk)
}


#' retest CNVs in a pseudobulk
#' @param bulk pesudobulk dataframe
#' @param gamma numeric Dispersion parameter for the Beta-Binomial allele model
#' @param allele_only whether to retest only using allele data
#' @return a dataframe of segments with CNV posterior information
#' @keywords internal
retest_cnv = function(bulk, theta_min = 0.08, logphi_min = 0.25, gamma = 20, allele_only = FALSE, exclude_neu = TRUE) {

    # rolling theta estimates
    bulk = annot_theta_roll(bulk)

    if (exclude_neu) {
        bulk = bulk %>% filter(cnv_state != 'neu')
    }
    
    G = c('20' = 1/5, '10' = 1/5, '21' = 1/10, '31' = 1/10, '22' = 1/5, '00' = 1/5)
    
    if (allele_only) {
        segs_post = bulk %>% 
            group_by(CHROM, seg, cnv_state) %>%
            summarise(
                n_genes = length(na.omit(unique(gene))),
                n_snps = sum(!is.na(pAD)),
                seg_start = min(POS),
                seg_end = max(POS),
                theta_hat = theta_hat_seg(major_count[!is.na(major_count)], minor_count[!is.na(minor_count)]),
                approx_theta_post(pAD[!is.na(pAD)], DP[!is.na(pAD)], p_s[!is.na(pAD)], gamma = unique(gamma), start = theta_hat),
                p_loh = 1, p_amp = 0, p_del = 0, p_bamp = 0, p_bdel = 0,
                LLR_x = 0,
                LLR_y = calc_allele_LLR(pAD[!is.na(pAD)], DP[!is.na(pAD)], p_s[!is.na(pAD)], theta_mle, gamma = unique(gamma)),
                LLR = LLR_x + LLR_y,
                phi_mle = 1,
                phi_sigma = 0,
                .groups = 'drop'
            ) %>%
            rowwise() %>%
            mutate(cnv_state_post = 'loh') %>%
            ungroup()
    } else {
        segs_post = bulk %>% 
            group_by(CHROM, seg, seg_start, seg_end, cnv_state) %>%
            summarise(
                n_genes = length(na.omit(unique(gene))),
                n_snps = sum(!is.na(pAD)),
                theta_hat = theta_hat_seg(major_count[!is.na(major_count)], minor_count[!is.na(minor_count)]),
                approx_theta_post(pAD[!is.na(pAD)], DP[!is.na(pAD)], p_s[!is.na(pAD)], gamma = unique(gamma), start = theta_hat),
                P_y_n = pnorm.range(0, theta_min, theta_mle, theta_sigma),
                P_y_d = pnorm.range(theta_min, 0.499, theta_mle, theta_sigma),
                P_y_a = pnorm.range(theta_min, 0.375, theta_mle, theta_sigma),
                approx_phi_post(
                    Y_obs[!is.na(Y_obs)], lambda_ref[!is.na(Y_obs)], unique(na.omit(d_obs)),
                    alpha = alpha[!is.na(Y_obs)],
                    beta = beta[!is.na(Y_obs)],
                    mu = mu[!is.na(Y_obs)],
                    sig = sig[!is.na(Y_obs)]
                ),
                P_x_n = pnorm.range(2^(-logphi_min), 2^logphi_min, phi_mle, phi_sigma),
                P_x_d = pnorm.range(0.1, 2^(-logphi_min), phi_mle, phi_sigma),
                P_x_a = pnorm.range(2^logphi_min, 3, phi_mle, phi_sigma),
                Z_cnv = sum(G['20'] * P_x_n * P_y_d,
                        G['10'] * P_x_d * P_y_d,
                        G['21'] * P_x_a * P_y_a,
                        G['31'] * P_x_a * P_y_a,
                        G['22'] * P_x_a * P_y_n, 
                        G['00'] * P_x_d * P_y_n),
                Z = Z_cnv/2 + 0.5 * P_x_n * P_y_n,
                p_neu = 0.5 * P_x_n * P_y_n/Z,
                p_loh = (G['20'] * P_x_n * P_y_d)/Z_cnv,
                p_amp = ((G['31'] + G['21']) * P_x_a * P_y_a)/Z_cnv,
                p_del = (G['10'] * P_x_d * P_y_d)/Z_cnv,
                p_bamp = (G['22'] * P_x_a * P_y_n)/Z_cnv,
                p_bdel = (G['00'] * P_x_d * P_y_n)/Z_cnv,
                LLR_x = calc_exp_LLR(
                    Y_obs[!is.na(Y_obs)],
                    lambda_ref[!is.na(Y_obs)], 
                    unique(na.omit(d_obs)),
                    phi_mle,
                    mu = mu[!is.na(Y_obs)],
                    sig = sig[!is.na(Y_obs)]
                ),
                LLR_y = calc_allele_LLR(pAD[!is.na(pAD)], DP[!is.na(pAD)], p_s[!is.na(pAD)], theta_mle, gamma = unique(gamma)),
                LLR = log(1-p_neu) - log(p_neu),
                .groups = 'drop'
            ) %>%
            mutate_at(
                c('p_loh', 'p_amp', 'p_del', 'p_bamp', 'p_bdel'),
                function(x) {ifelse(is.na(x), 0, x)}
            ) %>%
            rowwise() %>%
            mutate(cnv_state_post = c('loh', 'amp', 'del', 'bamp', 'bdel')[
                which.max(c(p_loh, p_amp, p_del, p_bamp, p_bdel))
            ]) %>%
            mutate(
                cnv_state_post = ifelse(p_neu >= 0.5, 'neu', cnv_state_post)
            )
    }

    return(segs_post)
}

#' classify alleles using viterbi and forward-backward
#' @param bulk dataframe Pesudobulk profile
#' @return dataframe Pesudobulk profile
#' @keywords internal
classify_alleles = function(bulk) {

    if (all(bulk$cnv_state_post %in% c('neu'))) {
        return(bulk)
    }
    
    allele_post = bulk %>%
        filter(!cnv_state_post %in% c('neu')) %>%
        filter(!is.na(AD)) %>%
        group_by(CHROM, seg) %>%
        filter(n() > 1) %>%
        mutate(
            p_up = forward_back_allele(get_allele_hmm(pAD, DP, p_s, theta = unique(theta_mle), gamma = 20))[,1],
            haplo_post = case_when(
                p_up >= 0.5 & GT == '1|0' ~ 'major',
                p_up >= 0.5 & GT == '0|1' ~ 'minor',
                p_up < 0.5 & GT == '1|0' ~ 'minor',
                p_up < 0.5 & GT == '0|1' ~ 'major'
            ),
            haplo_naive = ifelse(AR < 0.5, 'minor', 'major')
        ) %>%
        ungroup() %>%
        select(snp_id, p_up, haplo_post, haplo_naive)

    bulk = bulk %>% 
            select(-any_of(colnames(allele_post)[!colnames(allele_post) %in% c('snp_id')])) %>%
            left_join(
                allele_post,
            by = c('snp_id')
        )

    return(bulk)
}

#' Annotate rolling estimate of imbalance level theta
#' @param bulk a pseudobulk dataframe
#' @return a pseudobulk dataframe
#' @keywords internal
annot_theta_roll = function(bulk) {

    bulk = bulk %>% 
        mutate(haplo_theta_min = case_when(
            str_detect(state, 'up') ~ 'major',
            str_detect(state, 'down') ~ 'minor',
            TRUE ~ ifelse(pBAF > 0.5, 'major', 'minor')
        )) %>% 
        mutate(
            major_count = ifelse(haplo_theta_min == 'major', pAD, DP - pAD),
            minor_count = DP - major_count
        )

    bulk = bulk %>%
        select(-any_of('theta_hat_roll')) %>%
        left_join(
            bulk %>%
                group_by(CHROM) %>%
                filter(!is.na(pAD)) %>%
                mutate(theta_hat_roll = theta_hat_roll(major_count, minor_count, h = 100)) %>%
                select(theta_hat_roll, CHROM, snp_id),
            by = c('CHROM', 'snp_id')
        ) %>%
        group_by(CHROM) %>%
        mutate(theta_hat_roll = zoo::na.locf(theta_hat_roll, na.rm=FALSE)) %>%
        ungroup()

    return(bulk)
}

#' Estimate of imbalance level theta in a segment
#' @param major_count vector of major allele count
#' @param minor_count vector of minor allele count
#' @return estimate of theta
#' @keywords internal
theta_hat_seg = function(major_count, minor_count) {
    major_total = sum(major_count)
    minor_total = sum(minor_count)
    MAF = major_total/(major_total + minor_total)
    return(MAF - 0.5)
}

#' Rolling estimate of imbalance level theta
#' @param major_count vector of major allele count
#' @param minor_count vector of minor allele count
#' @param h window size
#' @return rolling estimate of theta
#' @keywords internal
theta_hat_roll = function(major_count, minor_count, h) {
    n = length(major_count)
    sapply(
        1:n,
        function(c) {
            slice = max(c - h - 1, 1):min(c + h, n)
            theta_hat_seg(major_count[slice], minor_count[slice])
        }
    )
}

#' Estimate of expression fold change phi in a segment
#' @keywords internal
phi_hat_seg = function(Y_obs, lambda_ref, d, mu, sig) {
    logFC = log((Y_obs/d)/lambda_ref)
    exp(mean(logFC - mu))
}

#' Rolling estimate of expression fold change phi
#' @keywords internal
phi_hat_roll = function(Y_obs, lambda_ref, d_obs, mu, sig, h) {
    n = length(Y_obs)
    
    if (length(mu) == 1 & length(sig) == 1) {
        mu = rep(mu, n)
        sig = rep(sig, n)
    }
    
    sapply(
        1:n,
        function(c) {
            slice = max(c - h - 1, 1):min(c + h, n)
            phi_hat_seg(Y_obs[slice], lambda_ref[slice], unique(d_obs), mu[slice], sig[slice])
        }
    )
}

#' @export
letters_all = c(letters, paste0(letters, letters), paste0(letters, letters, letters))

#' Annotate copy number segments after HMM decoding 
#' @param bulk dataframe Pseudobulk profile
#' @return a pseudobulk dataframe
#' @keywords internal
annot_segs = function(bulk, var = 'cnv_state') {

    bulk = bulk %>% 
            group_by(CHROM) %>%
            arrange(CHROM, snp_index) %>%
            mutate(boundary = c(0, get(var)[2:length(get(var))] != get(var)[1:(length(get(var))-1)])) %>%
            group_by(CHROM) %>%
            mutate(seg = paste0(CHROM, letters_all[cumsum(boundary)+1])) %>%
            arrange(CHROM) %>%
            mutate(seg = factor(seg, unique(seg))) %>%
            ungroup() %>%
            group_by(seg) %>%
            mutate(
                seg_start = min(POS),
                seg_end = max(POS),
                seg_start_index = min(snp_index),
                seg_end_index = max(snp_index),
                n_genes = length(na.omit(unique(gene))),
                n_snps = sum(!is.na(pAD))
            ) %>%
            ungroup()

    return(bulk)
}

#' Annotate haplotype segments after HMM decoding 
#' @keywords internal
annot_haplo_segs = function(bulk) {

    bulk = bulk %>% 
            group_by(CHROM) %>%
            arrange(CHROM, snp_index) %>%
            mutate(boundary = c(0, state[2:length(state)] != state[1:(length(state)-1)])) %>%
            group_by(CHROM) %>%
            mutate(seg = paste0(CHROM, '_', cumsum(boundary))) %>%
            arrange(CHROM) %>%
            mutate(seg = factor(seg, unique(seg))) %>%
            ungroup() %>%
            group_by(seg) %>%
            mutate(
                seg_start_index = min(snp_index),
                seg_end_index = max(snp_index),
                n_genes = length(na.omit(unique(gene))),
                n_snps = sum(!is.na(pAD))
            ) %>%
            ungroup()

    return(bulk)
}

#' Smooth the segments after HMM decoding 
#' @param bulk dataframe Pseudobulk profile
#' @param min_genes integer Minimum number of genes to call a segment
#' @return dataframe Pseudobulk profile
#' @keywords internal
smooth_segs = function(bulk, min_genes = 10) {
    bulk %>% group_by(seg) %>%
        mutate(
            cnv_state = ifelse(n_genes <= min_genes, NA, cnv_state)
        ) %>%
        ungroup() %>%
        group_by(CHROM) %>%
        mutate(cnv_state = zoo::na.locf(cnv_state, fromLast = FALSE, na.rm=FALSE)) %>%
        mutate(cnv_state = zoo::na.locf(cnv_state, fromLast = TRUE, na.rm=FALSE)) %>%
        ungroup()
}

#' Annotate a consensus segments on a pseudobulk dataframe
#' @param bulk dataframe Pseudobulk profile
#' @param segs_consensus datatframe Consensus segment dataframe
#' @return dataframe Pseudobulk profile
#' @keywords internal
annot_consensus = function(bulk, segs_consensus, join_mode = 'inner') {

    if (join_mode == 'inner') {
        join = inner_join
    } else {
        join = left_join
    }

    if (!'seg_cons' %in% colnames(segs_consensus)) {
        segs_consensus = segs_consensus %>% mutate(seg_cons = seg)
    }

    bulk = bulk %>% ungroup()
    
    marker_seg = GenomicRanges::findOverlaps(
            bulk %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$POS,
                    end = .$POS)
            )}, 
            segs_consensus %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$seg_start,
                    end = .$seg_end)
            )}
        ) %>%
        as.data.frame() %>%
        setNames(c('marker_index', 'seg_index')) %>%
        left_join(
            bulk %>% mutate(marker_index = 1:n()) %>%
                select(marker_index, snp_id),
            by = c('marker_index')
        ) %>%
        left_join(
            segs_consensus %>% mutate(seg_index = 1:n()),
            by = c('seg_index')
        ) %>%
        distinct(snp_id, `.keep_all` = TRUE) %>%
        select(-any_of(c('sample', 'marker_index', 'seg_index')))
    
    bulk = bulk %>% 
        select(-any_of(colnames(marker_seg)[!colnames(marker_seg) %in% c('snp_id', 'CHROM')])) %>% 
        join(marker_seg, by = c('snp_id', 'CHROM'))  %>%
        mutate(seg = seg_cons) %>%
        mutate(seg = factor(seg, gtools::mixedsort(unique(seg)))) 

    return(bulk)
}

#' Find the common diploid region in a group of pseudobulks
#' @param bulks dataframe Pseudobulk profiles (differentiated by "sample" column)
#' @param grouping logical Whether to use cliques or components in the graph to find dipoid cluster
#' @param gamma numeric Dispersion parameter for the Beta-Binomial allele model
#' @param theta_min numeric Minimum imbalance threshold
#' @param t numeric Transition probability
#' @param fc_min numeric Minimum fold change to call quadruploid cluster
#' @param alpha numeric FDR cut-off for q values to determine edges
#' @param ncores integer Number of cores to use  
#' @return list Ploidy information
#' @keywords internal
find_common_diploid = function(
    bulks, grouping = 'clique', gamma = 20, theta_min = 0.08, t = 1e-5, fc_min = 2^0.25, alpha = 1e-4, min_genes = 10, 
    ncores = 1, segs_loh = NULL, debug = FALSE, verbose = TRUE) {

    if (!'sample' %in% colnames(bulks)) {
        bulks$sample = '1'
    }

    # define balanced regions in each sample
    bulks = mclapply(
        bulks %>% split(.$sample),
        mc.cores = ncores,
        function(bulk) {
            
            bulk %>% 
                group_by(CHROM) %>%
                mutate(state = 
                    run_allele_hmm(
                        pAD = pAD,
                        DP = DP, 
                        p_s = p_s,
                        t = t,
                        theta_min = theta_min,
                        gamma = gamma
                    )
                ) %>% ungroup() %>%
                mutate(cnv_state = str_remove(state, '_down|_up')) %>%
                annot_segs() %>%
                smooth_segs(min_genes = min_genes) %>%
                annot_segs()

        }) %>%
        bind_rows()

    # annotate clonal LOH regions
    if (!is.null(segs_loh)) {
        bulks = bulks %>% 
            annot_consensus(segs_loh, join_mode = 'left') %>%
            mutate(loh = ifelse(is.na(loh), FALSE, TRUE)) %>%
            mutate(cnv_state = ifelse(loh, 'loh', cnv_state)) %>%
            annot_segs(var = 'cnv_state')
    } else {
        bulks = bulks %>% mutate(loh = FALSE)
    }

    # unionize imbalanced segs
    segs_imbal = bulks %>% 
        filter(cnv_state != 'neu') %>%
        distinct(sample, CHROM, seg, seg_start, seg_end) %>%
        {GenomicRanges::GRanges(
            seqnames = .$CHROM,
            IRanges::IRanges(start = .$seg_start,
                end = .$seg_end)
        )} %>%
        GenomicRanges::reduce() %>%
        as.data.frame() %>%
        select(CHROM = seqnames, seg_start = start, seg_end = end) %>%
        mutate(seg_length = seg_end - seg_start) %>%
        group_by(CHROM) %>%
        mutate(
            seg = paste0(CHROM, '_', 1:max(n(),1)),
            cnv_state = 'theta_1'
        ) %>%
        ungroup()

    segs_consensus = fill_neu_segs(segs_imbal, get_segs_neu(bulks)) %>% mutate(seg = seg_cons)
    bulks = bulks %>% annot_consensus(segs_consensus)

    # only keep segments with sufficient SNPs and genes
    segs_bal = bulks %>% 
        filter(cnv_state == 'neu') %>%
        group_by(seg, sample) %>%
        summarise(
            n_snps = sum(DP >= 5, na.rm = TRUE),
            n_genes = sum(!is.na(gene)),
            .groups = 'drop'
        ) %>%
        group_by(seg) %>%
        filter(any(n_genes > 50 & n_snps > 50)) %>%
        pull(seg) %>% 
        as.character %>% 
        unique()

    bulks_bal = bulks %>% filter(seg %in% segs_bal) %>% filter(!is.na(lnFC))

    if (length(segs_bal) == 0) {
        msg = 'No balanced segments, using all segments as baseline'
        log_warn(msg)
        diploid_segs = bulks %>% pull(seg) %>% unique
        bamp = TRUE
    } else if (length(segs_bal) == 1) {
        diploid_segs = segs_bal
        bamp = FALSE
    } else {
        test_dat = bulks_bal %>%
            select(gene, seg, lnFC, sample) %>%
            reshape2::dcast(seg+gene ~ sample, value.var = 'lnFC') %>%
            na.omit() %>%
            mutate(seg = as.character(seg)) %>%
            select(-gene)

        # sime's p
        samples = unique(bulks_bal$sample)

        tests = lapply(
                samples,
                function(sample) {
                    t(combn(segs_bal, 2)) %>% 
                        as.data.frame() %>%
                        setNames(c('i', 'j')) %>%
                        mutate(s = sample)
                }
            ) %>% 
            bind_rows() %>%
            rowwise() %>%
            mutate(
                p = t.test(
                    x = bulks_bal$lnFC[bulks_bal$seg == i & bulks_bal$sample == s],
                    y = bulks_bal$lnFC[bulks_bal$seg == j & bulks_bal$sample == s]
                )$p.value,
                lnFC_i = mean(bulks_bal$lnFC[bulks_bal$seg == i & bulks_bal$sample == s]),
                lnFC_j = mean(bulks_bal$lnFC[bulks_bal$seg == j & bulks_bal$sample == s])
            ) %>%
            group_by(i,j) %>%
            summarise(
                p = simes_p(p, length(samples)),
                lnFC_max_i = max(lnFC_i),
                lnFC_max_j = max(lnFC_j),
                delta_max = abs(lnFC_max_i - lnFC_max_j),
                .groups = 'drop'
            ) %>%
            ungroup() %>%
            mutate(q = p.adjust(p))

        # build graph
        V = segs_bal

        E = tests %>% filter(q > alpha)

        G = igraph::graph_from_data_frame(d=E, vertices=V, directed=FALSE)

        if (grouping == 'component') {

            components = igraph::components(G)

            FC = bulks_bal %>%
                mutate(component = components$membership[as.character(seg)]) %>%
                filter(!is.na(component)) %>%
                group_by(component, sample) %>%
                summarise(
                    lnFC = mean(lnFC, na.rm = T),
                    .groups = 'drop'
                ) %>%
                reshape2::dcast(sample ~ component, value.var = 'lnFC') %>%
                tibble::column_to_rownames('sample')

        } else {

            cliques = igraph::maximal.cliques(G)

            FC = lapply(
                1:length(cliques),
                function(i) {
                    bulks_bal %>% 
                    filter(seg %in% names(cliques[[i]])) %>%
                    group_by(sample) %>%
                    summarise(
                        lnFC = mean(lnFC, na.rm = T),
                        .groups = 'drop'
                    ) %>%
                    setNames(c('sample', i))
            }) %>%
            Reduce(x = ., f = function(x,y){full_join(x, y, by = 'sample')}) %>%
            tibble::remove_rownames() %>%
            tibble::column_to_rownames('sample')

        }

        # choose diploid clique based on min FC in case there's a tie
        diploid_cluster = as.integer(names(sort(apply(FC[,Modes(apply(FC, 1, which.min)),drop=FALSE], 2, min))))[1]

        fc_max = sapply(colnames(FC), function(c) {FC[,c] - FC[,diploid_cluster]}) %>% max %>% exp
        
        # seperation needs to be clear (2 vs 4 copy)
        if (fc_max > fc_min) {
            log_info('quadruploid state enabled')

            if (grouping == 'component') {
                diploid_segs = components$membership %>% keep(function(x){x == diploid_cluster}) %>% names
            } else {
                diploid_segs = names(cliques[[diploid_cluster]])
            }
            
            bamp = TRUE
        } else {
            diploid_segs = segs_bal
            bamp = FALSE
        }
    }

    log_info(glue('diploid regions: {paste0(gtools::mixedsort(diploid_segs), collapse = ",")}'))
    
    bulks = bulks %>% mutate(diploid = seg %in% diploid_segs)
    
    if (debug) {
        res = list(
            'bamp' = bamp, 'bulks' = bulks, 'diploid_segs' = diploid_segs,
            'segs_consensus' = segs_consensus, 'G' = G, 'tests' = tests,
            'test_dat' = test_dat, 'FC' = FC, 'cliques' = cliques, 'segs_bal' = segs_bal
        )
        return(res)
    }

    return(bulks)
}

#' get neutral segments from multiple pseudobulks
#' @keywords internal
get_segs_neu = function(bulks) {
    segs_neu = bulks %>% filter(cnv_state == "neu") %>%
        group_by(sample, seg, CHROM) %>% 
        mutate(seg_start = min(POS), seg_end = max(POS)) %>% 
        distinct(sample, CHROM, seg, seg_start, seg_end)

    segs_neu = segs_neu %>% arrange(CHROM) %>% {
            GenomicRanges::GRanges(seqnames = .$CHROM, IRanges::IRanges(start = .$seg_start, end = .$seg_end))
        } %>% 
        GenomicRanges::reduce() %>% 
        as.data.frame() %>% 
        select(CHROM = seqnames, seg_start = start, seg_end = end) %>% 
        mutate(seg_length = seg_end - seg_start)
    return(segs_neu)
}

#' calculate joint likelihood of allele data
#' @param AD numeric vector Variant allele depth
#' @param DP numeric vector Total allele depth
#' @param alpha numeric Alpha parameter of Beta-Binomial distribution
#' @param beta numeric Beta parameter of Beta-Binomial distribution
#' @return numeric Joint log likelihood
#' @keywords internal
l_bbinom = function(AD, DP, alpha, beta) {
    sum(dbbinom(AD, DP, alpha, beta, log = TRUE))
}

#' fit a Beta-Binomial model by maximum likelihood
#' @param AD numeric vector Variant allele depth
#' @param DP numeric vector Total allele depth
#' @return MLE of alpha and beta 
#' @keywords internal
fit_bbinom = function(AD, DP) {

    fit = stats4::mle(
        minuslogl = function(alpha, beta) {
            -l_bbinom(AD, DP, alpha, beta)
        },
        start = c(5, 5),
        lower = c(0.0001, 0.0001)
    )

    alpha = fit@coef[1]
    beta = fit@coef[2]
    
    return(fit)
}

#' fit gamma maximum likelihood
#' @param AD numeric vector Variant allele depth
#' @param DP numeric vector Total allele depth
#' @return a fit
#' @keywords internal
fit_gamma = function(AD, DP, start = 20) {

    fit = stats4::mle(
        minuslogl = function(gamma) {
            -l_bbinom(AD, DP, gamma/2, gamma/2)
        },
        start = start,
        lower = 0.0001
    )

    gamma = fit@coef[1]

    return(gamma)
}

#' calculate joint likelihood of a gamma-poisson model
#' @param Y_obs numeric vector Gene expression counts
#' @param lambda_ref numeric vector Reference expression levels
#' @param d numeric Total library size
#' @param alpha numeric Alpha parameter of Gamma-Poisson distribution
#' @param beta numeric Beta parameter of Gamma-Poisson distribution
#' @return numeric Joint log likelihood
#' @keywords internal
l_gpois = function(Y_obs, lambda_ref, d, alpha, beta, phi = 1) {
    sum(dgpois(Y_obs, shape = alpha, rate = beta/(phi * d * lambda_ref), log = TRUE))
}

#' fit a Gamma-Poisson model by maximum likelihood
#' @param Y_obs numeric vector Gene expression counts
#' @param lambda_ref numeric vector Reference expression levels
#' @param d numeric Total library size
#' @return numeric MLE of alpha and beta
#' @keywords internal
fit_gpois = function(Y_obs, lambda_ref, d) {

    Y_obs = Y_obs[lambda_ref > 0]
    lambda_ref = lambda_ref[lambda_ref > 0]
    
    fit = stats4::mle(
        minuslogl = function(alpha, beta) {
            -l_gpois(Y_obs, lambda_ref, d, alpha, beta)
        },
        start = c(1, 1),
        lower = c(0.01, 0.01),
        control = list('trace' = FALSE)
    )

    alpha = fit@coef[1]
    beta = fit@coef[2]
    
    return(fit)
}

#' ccalculate joint likelihood of a PLN model
#' @param Y_obs numeric vector Gene expression counts
#' @param lambda_ref numeric vector Reference expression levels
#' @param d numeric Total library size
#' @return numeric Joint log likelihood
#' @keywords internal
l_lnpois = function(Y_obs, lambda_ref, d, mu, sig, phi = 1) {
    if (any(sig <= 0)) {stop(glue('negative sigma. value: {sig}'))}
    if (length(sig) == 1) {sig = rep(sig, length(Y_obs))}
    sum(log(dpoilog(Y_obs, mu + log(phi * d * lambda_ref), sig)))
}

#' fit a PLN model by maximum likelihood
#' @param Y_obs numeric vector Gene expression counts
#' @param lambda_ref numeric vector Reference expression levels
#' @param d numeric Total library size
#' @return numeric MLE of mu and sig
#' @keywords internal
fit_lnpois = function(Y_obs, lambda_ref, d) {

    Y_obs = Y_obs[lambda_ref > 0]
    lambda_ref = lambda_ref[lambda_ref > 0]
    
    fit = stats4::mle(
        minuslogl = function(mu, sig) {
            if (sig < 0) {stop('optim is trying negative sigma')}
            -l_lnpois(Y_obs, lambda_ref, d, mu, sig)
        },
        start = c(0, 1),
        lower = c(-Inf, 0.01),
        control = list('trace' = FALSE)
    )

    return(fit)
}

#' Calculate the MLE of expression fold change phi
#' @keywords internal
calc_phi_mle_lnpois = function (Y_obs, lambda_ref, d, mu, sig, lower = 0.1, upper = 10) {
    
    if (length(Y_obs) == 0) {
        return(1)
    }
    
    start = max(min(1, upper), lower)
    
    res = stats4::mle(minuslogl = function(phi) {
        -l_lnpois(Y_obs, lambda_ref, d, mu, sig, phi)
    }, start = start, lower = lower, upper = upper)
    
    return(res@coef)
    
}

#' Laplace approximation of the posterior of allelic imbalance theta
#' @keywords internal
approx_theta_post = function(pAD, DP, p_s, lower = 0.001, upper = 0.499, start = 0.25, gamma = 20) {

    gamma = unique(gamma)

    if (length(gamma) > 1) {
        stop('gamma has to be a single value')
    }
    
    if (length(pAD) <= 10) {
        return(tibble('theta_mle' = 0, 'theta_sigma' = 0))
    }

    fit = optim(
        start, 
        function(theta) {-calc_allele_lik(pAD, DP, p_s, theta, gamma)},
        method = 'L-BFGS-B',
        lower = lower,
        upper = upper,
        hessian = TRUE
    )
    
    mu = fit$par
    sigma = sqrt(as.numeric(1/(fit$hessian)))

    if (is.na(sigma)) {
        sigma = 0
    }
    
    return(tibble('theta_mle' = mu, 'theta_sigma' = sigma))
}

#' Laplace approximation of the posterior of expression fold change phi
#' @keywords internal
approx_phi_post = function(Y_obs, lambda_ref, d, alpha = NULL, beta = NULL, mu = NULL, sig = NULL, lower = 0.2, upper = 10, start = 1) {
    
    if (length(Y_obs) == 0) {
        return(tibble('phi_mle' = 1, 'phi_sigma' = 0))
    }
    
    start = max(min(1, upper), lower)

    l = function(phi) {l_lnpois(Y_obs, lambda_ref, d, mu, sig, phi = phi)}

    fit = optim(
        start,
        function(phi) {
            -l(phi)
        },
        method = 'L-BFGS-B',
        lower = lower,
        upper = upper,
        hessian = TRUE
    )

    mean = fit$par
    sd = sqrt(as.numeric(1/(fit$hessian)))

    if (is.na(sd)) {
        mean = 1
        sd = 0
    }

    return(tibble('phi_mle' = mean, 'phi_sigma' = sd))
}

#' Helper function to get the internal nodes of a dendrogram and the leafs in each subtree 
#' @keywords internal
get_internal_nodes = function(den, node, labels) {

    membership = data.frame(
        cell = dendextend::get_leaves_attr(den, attribute = 'label'),
        node = node
    )
    
    is_leaf = length(unique(labels[dendextend::get_leaves_attr(den, attribute = 'label')])) == 1
    
    if (is_leaf) {
        return(data.frame())
    }
                
    sub_dens = dendextend::get_subdendrograms(den, k = 2)
    
    membership_l = get_internal_nodes(
        sub_dens[[1]],
        paste0(node, '.', 1),
        labels
    )
    
    membership_r = get_internal_nodes(
        sub_dens[[2]],
        paste0(node, '.', 2),
        labels
    )
        
    return(rbind(membership, membership_l, membership_r))
}

#' Get the internal nodes of a dendrogram and the leafs in each subtree 
#' @param hc hclust Clustering results
#' @param clusters named vector Cutree output specifying the terminal clusters
#' @return list Interal node subtrees with leaf memberships
#' @keywords internal
get_nodes_celltree = function(hc, clusters) {
        
    # internal nodes
    nodes = get_internal_nodes(as.dendrogram(hc), '0', clusters)
    
    # add terminal nodes
    nodes = nodes %>% mutate(cluster = clusters[cell]) %>%
        rbind(
            data.frame(
                cell = names(clusters),
                node = unname(clusters),
                cluster = unname(clusters)
            )
        )
    
    # convert to list
    nodes = nodes %>%
        split(.$node) %>%
        purrr::map(function(node){
            list(
                sample = unique(node$node),
                members = unique(node$cluster),
                cells = node$cell,
                size = length(node$cell)
            )
        })
    
    return(nodes)
    
}

#' Calculate expression distance matrix between cell populatoins
#' @param count_mat dgCMatrix Gene expression counts
#' @param cell_annot dataframe specifying the cell ID and cluster memberships
#' @return a distance matrix
#' @keywords internal
calc_cluster_dist = function(count_mat, cell_annot) {

    if (!is(count_mat, 'matrix')) {
        count_mat = as.matrix(count_mat)
    }

    cells = intersect(colnames(count_mat), cell_annot$cell)
    
    count_mat = count_mat %>% extract(,cells)
    cell_annot = cell_annot %>% filter(cell %in% cells)
    
    cell_dict = cell_annot %>% {setNames(.$cluster, .$cell)} %>% droplevels()
    cell_dict = cell_dict[cells]
    
    M = model.matrix(~ 0 + cell_dict) %>% set_colnames(levels(cell_dict))

    count_mat_clust = count_mat %*% M
    count_mat_clust = count_mat_clust[rowSums(count_mat_clust) > 0,]
    exp_mat_clust = count_mat_clust %*% diag(1/colSums(count_mat_clust)) %>% set_colnames(colnames(count_mat_clust))
    exp_mat_clust = log10(exp_mat_clust * 1e6 + 1)

    dist_mat = 1-cor(exp_mat_clust)
    
    return(dist_mat)
    
}

#' Calculate LLR for an allele HMM
#' @param pAD numeric vector Phased allele depth
#' @param DP numeric vector Total allele depth
#' @param p_s numeric vector Phase switch probabilities
#' @param theta_mle numeric MLE of imbalance level theta (alternative hypothesis)
#' @param theta_0 numeric Imbalance level in the null hypothesis 
#' @param gamma numeric Dispersion parameter for the Beta-Binomial allele model
#' @return numeric Log-likelihood ratio
#' @keywords internal
calc_allele_LLR = function(pAD, DP, p_s, theta_mle, theta_0 = 0, gamma = 20) {
    if (length(pAD) <= 1) {
        return(0)
    }
    l_1 = calc_allele_lik(pAD, DP, p_s, theta = theta_mle, gamma = gamma) 
    l_0 = l_bbinom(pAD, DP, gamma*0.5, gamma*0.5)
    return(l_1 - l_0)
}

#' Calculate LLR for an expression HMM
#' @param Y_obs numeric vector Gene expression counts
#' @param lambda_ref numeric vector Reference expression levels 
#' @param d numeric vector Total library size
#' @param phi_mle numeric MLE of expression fold change phi (alternative hypothesis)
#' @param mu numeric Mean parameter for the PLN expression model
#' @param sig numeric Dispersion parameter for the PLN expression model
#' @param alpha numeric Hyperparameter for the gamma poisson model (not used)
#' @param beta numeric Hyperparameter for the gamma poisson model (not used)
#' @return numeric Log-likelihood ratio
#' @keywords internal
calc_exp_LLR = function(Y_obs, lambda_ref, d, phi_mle, mu = NULL, sig = NULL, alpha = NULL, beta = NULL) {
    if (length(Y_obs) == 0) {
        return(0)
    }

    l_1 = l_lnpois(Y_obs, lambda_ref, d, mu, sig, phi = phi_mle)
    l_0 = l_lnpois(Y_obs, lambda_ref, d, mu, sig, phi = 1)
    
    return(l_1 - l_0)
}

########################### Misc ############################

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

log_err = function(msg) {
    log_error(msg)
    stop(msg)
}

#' check the format of a count matrix
#' @keywords internal
check_matrix = function(count_mat) {
    if ('matrix' %in% class(count_mat)) {
        count_mat <- as(Matrix(count_mat, sparse=TRUE), "dgCMatrix")
    }
    if (!('dgCMatrix' %in% class(count_mat))) {
        log_err("count_mat is not of class dgCMatrix or matrix")
    }
    if (any(duplicated(rownames(count_mat)))) {
        log_err("Please remove duplicated genes in count matrix")
    }
    return(count_mat)
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

#' Calculate simes' p
#' @keywords internal
simes_p = function(p.vals, n_dim) {
    n_dim * min(sort(p.vals)/seq_along(p.vals))
}

#' Get the total probability from a region of a normal pdf
#' @keywords internal
pnorm.range = function(lower, upper, mu, sd) {
    if (sd == 0) {
        return(1)
    }
    pnorm(upper, mu, sd) - pnorm(lower, mu, sd)
}

#' Get the modes of a vector
#' @keywords internal
Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

#' calculate entropy for a binary variable 
#' @keywords internal
binary_entropy = function(p) {
    H = -p*log2(p)-(1-p)*log2(1-p)
    H[is.na(H)] = 0
    return(H)
}


fit_switch_prob = function(y, d) {
    
    eta = function(d, lambda, min_p = 1e-10) {
        switch_prob_cm(d, lambda)
    }

    l_lambda = function(y, d, lambda) {
        sum(log(eta(d[y == 1], lambda))) + sum(log(1-eta(d[y == 0], lambda)))
    }

    fit = stats4::mle(
        minuslogl = function(lambda) {
            -l_lambda(y, d, lambda)
        },
        start = c(5),
        lower = c(1e-7)
    )
    
    return(fit@coef)
}

switch_prob_mle = function(pAD, DP, d, theta, gamma = 20) {

    fit = optim(
            par = 1, 
            function(lambda) {
                p_s = switch_prob_cm(d, lambda)
                -calc_allele_lik(pAD, DP, p_s, theta, gamma)
            },
            method = 'L-BFGS-B',
            lower = 0,
            hessian = TRUE
        )
    
    return(fit)
    
}

switch_prob_mle2 = function(pAD, DP, d, gamma = 20) {

    fit = optim(
            par = c(1, 0.4), 
            function(params) {
                lambda = params[1]
                theta = params[2]
                p_s = switch_prob_cm(d, lambda)
                -calc_allele_lik(pAD, DP, p_s, theta, gamma)
            },
            method = 'L-BFGS-B',
            lower = c(0,0),
            upper = c(NA,0.499),
            hessian = TRUE
        )
    
    return(fit)
    
}

read_if_exist = function(f) {

    if (file.exists(f)) {
        if (str_detect(f, 'tsv|txt|csv')) {
            return(fread(f))
        } else if (str_detect(f, 'rds')) {
            return(readRDS(f))
        } else {
            return(NULL)
        }
    }
}


l_geno_2d = function(g, theta_mle, phi_mle, theta_sigma, phi_sigma) {

    if (g[[1]] == 2 & g[[2]] == 2) {
        integrand = function(f) {
            pnorm.range(0, 0.04, theta_mle, theta_sigma) * 
            dnorm(2^logphi_f(f, g[[1]], g[[2]]), phi_mle, phi_sigma) 
        }
    } else {
        integrand = function(f) {
            dnorm(theta_f(f, g[[1]], g[[2]])-0.5, theta_mle, theta_sigma) * 
            dnorm(2^logphi_f(f, g[[1]], g[[2]]), phi_mle, phi_sigma) 
        }
    }

    integrate(
        integrand,
        lower = 0, 
        upper = 1
    )$value

}

p_geno_2d = function(theta_mle, phi_mle, theta_sigma, phi_sigma) {
    
    ps = sapply(
        list(c(3,1), c(2,1), c(2,2), c(1,0), c(2,0)),
        function(g) {
            l_geno_2d(g, theta_mle, phi_mle, theta_sigma, phi_sigma)
        }
    )
    
    ps = ps/sum(ps)

    tibble('p_amp' = ps[1] + ps[2], 'p_bamp' = ps[3], 'p_del' = ps[4], 'p_loh' = ps[5])
    
}

logphi_f = function(f, m = 0, p = 1) {
  log2((2 * (1 - f) + (m + p) * f) / 2)
}

theta_f = function(f, m = 0, p = 1) {
    baf = (m * f + 1 - f)/((m * f + 1 - f) + (p * f + 1 - f))
    baf
}

snp_rate_roll = function(gene_snps, gene_length, h) {
    n = length(gene_snps)
    sapply(
        1:n,
        function(c) {
            slice = max(c - h - 1, 1):min(c + h, n)
            fit_snp_rate(gene_snps[slice], gene_length[slice])[1]
        }
    )
}

fit_snp_rate_pois = function(gene_snps, gene_length) {
    
    fit = optim(
        par = 1,
        fn = function(lambda) {
            -sum(dpois(x = gene_snps, lambda = lambda * gene_length/1e6, log = TRUE))
        },
        method = 'L-BFGS-B',
        lower = 1e-10
    )
    
    return(fit$par)
}

fit_snp_rate_poilog = function(gene_snps, gene_length) {
    
    n = length(gene_snps)

    fit = optim(
        par = c(30, 1),
        fn = function(params) {
            lambda = params[1]
            sig = params[2]
            -sum(dpoilog(x = gene_snps, mu = log(lambda * gene_length/1e6), sig = rep(sig, n), log = TRUE))
        },
        method = 'L-BFGS-B',
        lower = c(1e-10, 1e-10)
    )
    
    return(fit$par)
}

# negative binomial model
fit_snp_rate = function(gene_snps, gene_length) {
    
    n = length(gene_snps)

    fit = optim(
        par = c(10, 1),
        fn = function(params) {
            lambda = params[1]
            sig = params[2]
            
            -sum(dnbinom(x = gene_snps, mu = lambda * gene_length/1e6, size = sig, log = TRUE))
        },
        method = 'L-BFGS-B',
        lower = c(1e-10, 1e-10)
    )
    
    return(fit$par)
}

# detect clonal LOH
detect_loh = function(bulk, min_depth = 0, t = 1e-5, mu = NULL, sig = NULL) {

    # bulk = find_common_diploid(bulk)

    # fit = bulk %>%
    #     filter(!is.na(Y_obs)) %>%
    #     filter(logFC < 8 & logFC > -8) %>%
    #     filter(diploid) %>%
    #     {fit_lnpois(.$Y_obs, .$lambda_ref, unique(.$d_obs))}

    bulk_snps = bulk %>% 
        filter(!is.na(gene)) %>%
        group_by(CHROM, gene, gene_start, gene_end) %>%
        summarise(
            gene_snps = sum(!is.na(AD[DP > min_depth]), na.rm = TRUE),
            Y_obs = unique(na.omit(Y_obs)),
            lambda_ref = unique(na.omit(lambda_ref)),
            logFC = unique(na.omit(logFC)),
            d_obs = unique(na.omit(d_obs)),
            .groups = 'drop'
        ) %>%
        mutate(gene_length = gene_end - gene_start)

    if (is.null(mu) | is.null(sig)) {

        fit = bulk_snps %>%
            filter(logFC < 8 & logFC > -8) %>%
            {fit_lnpois(.$Y_obs, .$lambda_ref, unique(.$d_obs))}

        mu = fit@coef[1]
        sig = fit@coef[2]

    }

    snp_fit = fit_snp_rate(bulk_snps$gene_snps, bulk_snps$gene_length)
    snp_ref = snp_fit[1]
    snp_sig = snp_fit[2]

    n = nrow(bulk_snps)
    A = matrix(c(1-t, t, t, 1-t), nrow = 2)
    As = array(rep(A,n), dim = c(2,2,n))

    hmm = list(
        x = bulk_snps$gene_snps, 
        Pi = As, 
        delta = c(1-t,t), 
        pm = c(snp_ref, 5),
        pn = bulk_snps$gene_length/1e6,
        snp_sig = snp_sig,
        y = bulk_snps$Y_obs,
        phi = c(1, 0.5),
        lambda_star = bulk_snps$lambda_ref,
        d = unique(bulk_snps$d_obs),
        mu = mu,
        sig = sig,
        states = c('neu', 'loh')
    )

    bulk_snps = bulk_snps %>% mutate(cnv_state = viterbi_loh(hmm))

    segs_loh = bulk_snps %>%
        mutate(snp_index = 1:n(), POS = gene_start, pAD = 1) %>%
        annot_segs() %>%
        group_by(CHROM, seg, seg_start, seg_end, cnv_state) %>%
        summarise(
            snp_rate = fit_snp_rate(gene_snps, gene_length)[1],
            .groups = 'drop'
        ) %>%
        ungroup() %>%
        filter(cnv_state == 'loh') %>%
        mutate(loh = TRUE) %>% 
        select(CHROM, seg, seg_start, seg_end, snp_rate, loh)

    if (nrow(segs_loh) == 0) {
        segs_loh = NULL
    }
    
    return(segs_loh)
}

########################### Experimental ############################

test_branch = function(count_mat, df_allele, gtree, segs_consensus, from_node, to_node, ncores, ...) {
    
    segs_consensus = segs_consensus %>% mutate(CHROM = factor(CHROM))
    
    gtree = gtree %>% activate(nodes)
    
    # from root
    if (is.na(from_node)) {
        from_node = gtree %>% filter(root) %>% pull(name)
    }
    
    from_node = gtree %>% filter(name == from_node) %>% pull(id)
    to_node = gtree %>% filter(name == to_node) %>% pull(id)
    
    cells_1 = get_desc(gtree, from_node)
    cells_2 = get_desc(gtree, to_node)

    cells_1 = cells_1[!cells_1 %in% cells_2]
    
    groups = list(
        list(sample = '1', cells = cells_1, size = length(cells_1)), 
        list(sample = '2', cells = cells_2, size = length(cells_2))
    )
    
    bulks = make_group_bulks(
        groups = groups,
        count_mat = count_mat,
        df_allele = df_allele, 
        lambdas_ref = lambdas_ref,
        gtf = gtf,
        genetic_map = genetic_map,
        min_depth = min_depth,
        ncores = ncores
    ) %>%
    annot_consensus(segs_consensus)
        
    bulks = lapply(
        bulks %>% split(.$sample),
        function(bulk) {

            fit = bulk %>%
                filter(cnv_state == 'neu') %>%
                filter(!is.na(Y_obs)) %>%
                filter(logFC < 8 & logFC > -8) %>%
                {fit_lnpois(.$Y_obs, .$lambda_ref, unique(.$d_obs))}

            bulk %>% mutate(mu = fit@coef[1], sig = fit@coef[2])

        }
    ) %>% bind_rows()
    
    return(bulks)
}

compare_bulks = function(bulks, segs_consensus, t = 1e-5, gamma = 20, ncores = 1) {
    
    segs_consensus = segs_consensus %>% mutate(CHROM = factor(CHROM))
    
    bulks = bulks %>%
        arrange(sample) %>%
        mutate(sample = as.integer(factor(sample, unique(sample)))) %>%
        select(-phi_mle) %>%
        arrange(CHROM, POS) %>%
        mutate(snp_index = as.integer(factor(snp_id, unique(snp_id)))) %>%
        group_by(CHROM) %>%
        ungroup()

    bulks_unified = bulks %>% select(sample, CHROM, POS, snp_id, snp_index, cM, pAD, DP, Y_obs, p_s, mu, sig, d_obs, lambda_ref) %>%
        tidyr::pivot_wider(names_from = "sample", "values_from" = c("pAD", "DP", "Y_obs", "mu", "sig", "d_obs", "lambda_ref")) %>%
        mutate(
            inter_snp_cm = c(NA, cM[2:length(cM)] - cM[1:(length(cM)-1)]),
            p_s = switch_prob_cm(inter_snp_cm)
        )
    
    bulks_unified = bulks_unified %>% annot_consensus(segs_consensus)
    
    LLRs = mclapply(
            mc.cores = ncores,
            bulks_unified %>%
                filter(cnv_state != 'neu') %>%
                split(droplevels(.$seg)),
            function(bulks) {

                logphi_min = log(unique(bulks$phi_mle))
                theta_min = unique(bulks$theta_mle)

                hmm1 = bulks %>%
                    {get_joint_hmm(
                        pAD = .$pAD_1,
                        DP = .$DP_1, 
                        Y_obs = .$Y_obs_1, 
                        d_total = na.omit(unique(.$d_obs_1)),
                        p_s = .$p_s,
                        lambda_ref = .$lambda_ref_1, 
                        phi_amp = 2^(logphi_min),
                        phi_del = 2^(-logphi_min),
                        mu = .$mu_1,
                        sig = .$sig_1,
                        t = t,
                        gamma = gamma,
                        theta_min = theta_min
                    )}

                hmm2 = bulks %>%
                    {get_joint_hmm(
                        pAD = .$pAD_2,
                        DP = .$DP_2, 
                        Y_obs = .$Y_obs_2, 
                        d_total = na.omit(unique(.$d_obs_2)),
                        p_s = .$p_s,
                        lambda_ref = .$lambda_ref_2, 
                        phi_amp = 2^(logphi_min),
                        phi_del = 2^(-logphi_min),
                        mu = .$mu_2,
                        sig = .$sig_2,
                        t = t,
                        gamma = gamma,
                        theta_min = theta_min
                    )}

                LL = compare_hmm(hmm1, hmm2)

                return(LL)

            }
        ) %>% unlist()
    
    return(LLRs)
    
}

# evaluate pseduobulks against one set of segmentations
evaluate_segs = function(bulks, segs, t = 1e-5, ncores = 1) {

    M = length(unique(bulks$sample))

    bulks %>% 
    annot_consensus(segs) %>%
    split(droplevels(.$CHROM)) %>%
    mclapply(
        mc.cores = ncores,
        function(bulks){
            bulks %>%
            group_by(CHROM, sample, seg) %>%
            mutate(
                approx_theta_post(
                    pAD[!is.na(pAD)], DP[!is.na(pAD)], p_s[!is.na(pAD)], start = 0.25
                ),
                approx_phi_post(
                    Y_obs[!is.na(Y_obs)], 
                    lambda_ref[!is.na(Y_obs)],
                    unique(na.omit(d_obs)),
                    mu = mu[!is.na(Y_obs)],
                    sig = sig[!is.na(Y_obs)]
                )
            ) %>%
            ungroup() %>%
            group_by(CHROM, sample) %>%
            summarise(
                l_x = l_lnpois(
                    Y_obs[!is.na(Y_obs)], lambda_ref[!is.na(Y_obs)], unique(na.omit(d_obs)),
                    mu[!is.na(Y_obs)], sig[!is.na(Y_obs)], phi = phi_mle[!is.na(Y_obs)]
                ),
                l_y = calc_allele_lik(pAD[!is.na(pAD)], DP[!is.na(pAD)], p_s[!is.na(pAD)], theta = theta_mle[!is.na(pAD)]),
                n_x = sum(!is.na(Y_obs)),
                n_y = sum(!is.na(pAD)),
                bkp = length(unique(seg)) - 1,
                l = l_x + l_y,
                .groups = 'drop'
            ) %>%
            group_by(CHROM) %>%
            summarise(
                l = sum(l),
                l_x = sum(l_x),
                l_y = sum(l_y),
                n_x = sum(n_x),
                n_y = sum(n_y),
                bkp = min(bkp),
                l_adj = l + bkp * log(t),
                .groups = 'drop'
            )
    }) %>%
    bind_rows()
}

# get optimal segmentation from subtree and clone bulk results
get_segs_optimal = function(bulk_subtrees, bulk_clones, t = 1e-5, min_LLR = 10, ncores = 1) {
    
    samples = unique(bulk_subtrees$sample)

    segs_all = lapply(
            samples,
            function(sample) {
                bulk_subtrees %>% 
                filter(sample == UQ(sample)) %>%
                get_segs_consensus(min_LLR = min_LLR) %>%
                mutate(sample = UQ(sample))
            }
        )

    scores = lapply(
            segs_all,
            function(segs) {
                bulk_clones %>% 
                    evaluate_segs(segs, t = t, ncores = ncores) %>%
                    mutate(sample = unique(segs$sample))
            }
        ) %>% bind_rows()

    segs_optimal = inner_join(
            bind_rows(segs_all),
            scores %>%
                arrange(CHROM, -l_adj) %>%
                distinct(CHROM, .keep_all = TRUE) %>%
                select(CHROM, sample),
            by = c('sample', 'CHROM')
        ) %>%
        arrange(CHROM) %>%
        group_by(CHROM) %>%
        mutate(seg_cons = paste0(CHROM, letters_all[1:n()])) %>%
        ungroup() %>%
        mutate(seg = seg_cons)
    
    return(list(segs_optimal = segs_optimal, scores = scores))
    
}
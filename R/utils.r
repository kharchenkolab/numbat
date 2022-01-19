########################### Data processing ############################

# take in a reference matrix, decide which one is best 
#' @export
choose_ref = function(count_mat_obs, lambda_mat_ref, gtf_transcript, debug = F) {

    d = sum(count_mat_obs)

    Y_obs = rowSums(count_mat_obs)

    genes_common = gtf_transcript$gene %>% 
        intersect(rownames(count_mat_obs)) %>%
        intersect(rownames(lambda_mat_ref)[rowSums(lambda_mat_ref == 0) == 0])
    
    Y_obs = Y_obs[genes_common]
    lambda_mat_ref = lambda_mat_ref[genes_common,]

    D = data.frame(gene = genes_common, Y_obs, lambda_mat_ref, check.names = F) %>%
                left_join(gtf_transcript, by = 'gene')
    
    errs = lapply(colnames(lambda_mat_ref), function(ref) {

            lambdas_ref = D[,ref]

            D %>% 
            group_by(CHROM) %>%
            summarise(
                l_total = sum(dpois(x = Y_obs, lambda = get(ref) * d, log = TRUE)),
                mse = mean((log2((Y_obs/d + 1e-6)/get(ref)))^2),
                n = length(genes_common)
            ) %>%
            mutate(l_bar = l_total/n) %>%
            mutate(ref = ref)

        }) %>% 
        bind_rows() 
        # mutate(ref = colnames(lambda_mat_ref))
    
    if (debug) {
        return(errs)
    }

    return(errs$ref[which.min(-errs$l_bar)])
    
}

# choose ref based on likelihood
#' @export
choose_ref = function(Y_obs, lambdas_ref, d, gtf_transcript, debug = F) {

    genes_common = gtf_transcript$gene %>% 
        intersect(names(Y_obs)) %>%
        intersect(rownames(lambdas_ref)[rowSums(lambdas_ref == 0) == 0])
    
    Y_obs = Y_obs[genes_common]
    lambdas_ref = lambdas_ref[genes_common,]

    D = data.frame(gene = genes_common, Y_obs, lambdas_ref, check.names = F) %>%
                left_join(gtf_transcript, by = 'gene')
    
    errs = lapply(colnames(lambdas_ref), function(ref) {

            lambdas_ref = D[,ref]

            D %>% 
            summarise(
                l_total = sum(dpois(x = Y_obs, lambda = get(ref) * d, log = TRUE)),
                mse = mean((log2((Y_obs/d + 1e-6)/get(ref)))^2),
                n = length(genes_common)
            ) %>%
            mutate(l_bar = l_total/n) %>%
            mutate(ref = ref)

        }) %>% 
        bind_rows() %>%
        arrange(-l_bar)
    
    if (debug) {
        return(errs)
    }

    return(errs$ref[which.min(-errs$l_bar)])
    
}

# choose ref based on correlation
#' @export
choose_ref_cor = function(count_mat, lambdas_ref, gtf_transcript) {
    
    genes_annotated = gtf_transcript %>% 
        pull(gene) %>% 
        intersect(rownames(count_mat)) %>%
        intersect(rownames(lambdas_ref))

    count_mat = count_mat[genes_annotated,,drop=F]
    lambdas_ref = lambdas_ref[genes_annotated,,drop=F]
    
    # keep highly expressed genes in at least one of the references
    count_mat = count_mat[rowSums(lambdas_ref * 1e6 > 2) > 0,,drop=F]

    exp_mat = scale(count_mat, center=FALSE, scale=colSums(count_mat))
    
    cors = cor(log(exp_mat * 1e6 + 1), log(lambdas_ref * 1e6 + 1)[rownames(exp_mat),])
    best_refs = apply(cors, 1, function(x) {colnames(cors)[which.max(x)]})
    
    return(best_refs)
}

#' @export
fit_multi_ref_lnpois = function(Y_obs, lambdas_ref, d, gtf_transcript, min_lambda = 2e-6, verbose = FALSE) {

    if (length(dim(lambdas_ref)) == 1 | is.null(dim(lambdas_ref))) {
        return(list('w' = 1, 'lambdas_bar' = lambdas_ref))
    }

    # take the union of expressed genes across cell type
    genes_common = gtf_transcript$gene %>% 
        intersect(names(Y_obs)) %>%
        intersect(rownames(lambdas_ref)[rowMeans(lambdas_ref) > min_lambda])

    if (verbose) {
        log_info(glue('{length(genes_common)} genes common in reference and observation'))
    }

    Y_obs = Y_obs[genes_common]
    lambdas_ref = lambdas_ref[genes_common,,drop=F]

    n_ref = ncol(lambdas_ref)

    n = length(Y_obs)

    fit = optim(
        fn = function(params) {
            sig = params[1]
            w = params[2:length(params)]
            -sum(log(dpoilog(Y_obs, log(d * lambdas_ref %*% w), rep(sig, n))))
        },
        method = 'L-BFGS-B',
        par = c(1, rep(1/n_ref, n_ref)),
        lower = c(0.01, rep(1e-6, n_ref))
    )

    sig = fit$par[1]
    w = fit$par[2:length(fit$par)]
    mu = log(1/sum(w))
    w = w * exp(mu)
    w = set_names(w, colnames(lambdas_ref))

    lambdas_bar = lambdas_ref %*% w %>% {set_names(as.vector(.), rownames(.))}

    return(list('mu' = mu, 'sig' = sig, 'w' = w, 'lambdas_bar' = lambdas_bar))
}

#' @export
fit_multi_ref = function(Y_obs, lambdas_ref, d, gtf_transcript, min_lambda = 2e-6, verbose = FALSE) {

    if (length(dim(lambdas_ref)) == 1 | is.null(dim(lambdas_ref))) {
        return(list('w' = 1, 'lambdas_bar' = lambdas_ref))
    }

    # take the union of expressed genes across cell type
    genes_common = gtf_transcript$gene %>% 
        intersect(names(Y_obs)) %>%
        intersect(rownames(lambdas_ref)[rowMeans(lambdas_ref) > min_lambda])

    if (verbose) {
        log_info(glue('{length(genes_common)} genes common in reference and observation'))
    }

    Y_obs = Y_obs[genes_common]
    lambdas_ref = lambdas_ref[genes_common,,drop=F]

    n_ref = ncol(lambdas_ref)
    
    fit = optim(
        fn = function(params) {
            alpha = params[1]
            w = params[2:length(params)]
            -sum(dgpois(Y_obs, shape = alpha, rate = 1/(d * (lambdas_ref %*% w)), log = TRUE))
        },
        method = 'L-BFGS-B',
        par = c(1, rep(1/n_ref, n_ref)),
        lower = c(0.01, rep(1e-6, n_ref))
    )

    alpha = fit$par[1]
    w = fit$par[2:length(fit$par)]
    beta = 1/sum(w)
    w = w * beta
    w = set_names(w, colnames(lambdas_ref))

    lambdas_bar = lambdas_ref %*% w %>% {set_names(as.vector(.), rownames(.))}

    return(list('alpha' = alpha, 'beta' = beta, 'w' = w, 'lambdas_bar' = lambdas_bar))
}


Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

# for gamma poisson model
#' @export
process_exp_fc = function(count_mat_obs, lambdas_ref, gtf_transcript, verbose = TRUE) {
    
    genes_annotated = gtf_transcript$gene %>% 
        intersect(rownames(count_mat_obs)) %>%
        intersect(names(lambdas_ref))

    depth_obs = sum(count_mat_obs)
    
    count_mat_obs = count_mat_obs[genes_annotated,,drop=F]
    lambdas_ref = lambdas_ref[genes_annotated]

    lambdas_obs = rowSums(count_mat_obs)/depth_obs

    # filter for mutually expressed genes
    min_both = 2
    
    mut_expressed = ((lambdas_ref * 1e6 > min_both & lambdas_obs * 1e6 > min_both) |
        (lambdas_ref > mean(lambdas_ref[lambdas_ref != 0])) |
        (lambdas_obs > mean(lambdas_obs[lambdas_obs != 0]))) &
        (lambdas_ref > 0)
    
    count_mat_obs = count_mat_obs[mut_expressed,,drop=F]
    lambdas_ref = lambdas_ref[mut_expressed]

    if(verbose){message(paste0(nrow(count_mat_obs), ' genes remain after filtering'))}

    bulk_obs = count_mat_obs %>%
        rowSums() %>%
        data.frame() %>%
        setNames('Y_obs') %>%
        tibble::rownames_to_column('gene') %>%
        mutate(lambda_obs = (Y_obs/depth_obs)) %>%
        mutate(lambda_ref = lambdas_ref[gene]) %>%
        mutate(d_obs = depth_obs) %>%
        left_join(gtf_transcript, by = "gene") 
    
    # annotate using GTF
    bulk_obs = bulk_obs %>%
        mutate(gene = droplevels(factor(gene, gtf_transcript$gene))) %>%
        mutate(gene_index = as.integer(gene)) %>%
        arrange(gene) %>%
        mutate(CHROM = factor(CHROM)) %>%
        mutate(
            logFC = log2(lambda_obs) - log2(lambda_ref),
            lnFC = log(lambda_obs) - log(lambda_ref),
            logFC = ifelse(is.infinite(logFC), NA, logFC),
            lnFC = ifelse(is.infinite(lnFC), NA, lnFC)
        ) %>%
        group_by(CHROM) %>%
        mutate(
            phi_hat_roll = phi_hat_roll(Y_obs, lambda_ref, depth_obs, h = 50)
        ) %>%
        ungroup()
    
    return(list('bulk' = bulk_obs,
                'count_mat_obs' = count_mat_obs,
               'depth_obs' = depth_obs))
}

# match cells with best reference
#' @export
process_exp2 = function(count_mat_obs, lambdas_ref, gtf_transcript, window = 101, verbose = T) {

    count_mat_obs = as.matrix(count_mat_obs)

    lambdas_ref = as.matrix(lambdas_ref)
    
    genes_annotated = gtf_transcript$gene %>% 
        intersect(rownames(count_mat_obs)) %>%
        intersect(rownames(lambdas_ref))

    depth_obs = sum(count_mat_obs)

    count_mat_obs = count_mat_obs[genes_annotated,,drop=F]
    lambdas_ref = lambdas_ref[genes_annotated,]

    lambdas_obs = rowSums(count_mat_obs)/depth_obs

    # filter for mutually expressed genes
    min_both = 2

    mut_expressed = ((rowSums(lambdas_ref * 1e6 > min_both) > 0 & lambdas_obs * 1e6 > min_both) | 
        (lambdas_obs > mean(lambdas_obs[lambdas_obs != 0])))

    count_mat_obs = count_mat_obs[mut_expressed,,drop=F]
    lambdas_ref = lambdas_ref[mut_expressed,]
    
    if(verbose){message(nrow(count_mat_obs))}

    exp_mat = scale(count_mat_obs, center=FALSE, scale=colSums(count_mat_obs))

    if (verbose) {message('choosing best references by correlation')}
    cors = cor(log2(exp_mat * 1e6 + 1), log2(lambdas_ref * 1e6 + 1)[rownames(exp_mat),])
    best_refs = apply(cors, 1, function(x) {colnames(cors)[which.max(x)]})
    exp_mat_ref = lambdas_ref[,best_refs] %>% set_colnames(names(best_refs))
    
    exp_mat_norm = log2(exp_mat * 1e6 + 1) - log2(exp_mat_ref * 1e6 + 1)
    
    gexp.norm = exp_mat_norm %>% as.data.frame() %>%
        tibble::rownames_to_column('gene') %>%
        inner_join(gtf_transcript, by = "gene") %>%
        filter(!(CHROM == 6 & gene_start < 33480577 & gene_end > 28510120)) %>%
        mutate(gene = droplevels(factor(gene, gtf_transcript$gene))) %>%
        mutate(gene_index = as.integer(gene)) %>% 
        mutate(CHROM = as.integer(CHROM)) %>%
        arrange(CHROM, gene)

    gexp.norm.long = gexp.norm %>% 
        reshape2::melt(
            id.var = c('gene', 'gene_index', 'region', 'gene_start', 'gene_end', 'CHROM', 'gene_length'),
            variable.name = 'cell',
            value.name = 'exp') %>%
        ungroup()

    gexp.norm.long = gexp.norm.long %>%
        group_by(cell, CHROM) %>%
        mutate(exp_rollmean = caTools::runmean(exp, k = window, align = "center")) %>%
        ungroup()
    
    return(list('gexp.norm.long' = gexp.norm.long, 'gexp.norm' = gexp.norm, 'exp_mat' = exp_mat, 'best_refs' = best_refs))
}

# compare bulk vs single cells
#' @export
process_exp = function(count_mat_obs, lambdas_ref, gtf_transcript, window = 101, verbose = T) {

    count_mat_obs = as.matrix(count_mat_obs)
    
    genes_annotated = gtf_transcript$gene %>% 
        intersect(rownames(count_mat_obs)) %>%
        intersect(names(lambdas_ref))

    depth_obs = sum(count_mat_obs)

    count_mat_obs = count_mat_obs[genes_annotated,,drop=F]
    lambdas_ref = lambdas_ref[genes_annotated]

    lambdas_obs = rowSums(count_mat_obs)/depth_obs

    # filter for mutually expressed genes
    min_both = 2

    mut_expressed = ((lambdas_ref * 1e6 > min_both & lambdas_obs * 1e6 > min_both) |
        (lambdas_ref > mean(lambdas_ref[lambdas_ref != 0])) |
        (lambdas_obs > mean(lambdas_obs[lambdas_obs != 0]))) &
        (lambdas_ref > 0)

    count_mat_obs = count_mat_obs[mut_expressed,,drop=F]
    lambdas_ref = lambdas_ref[mut_expressed]
    
    if(verbose){message(nrow(count_mat_obs))}
    
    exp_mat = scale(count_mat_obs, center=FALSE, scale=colSums(count_mat_obs))

    exp_mat_norm = log2(exp_mat * 1e6 + 1) - log2(lambdas_ref * 1e6 + 1)

    # capping 
    cap = 3
    exp_mat_norm[exp_mat_norm > cap] = cap
    exp_mat_norm[exp_mat_norm < -cap] = -cap

    # centering by cell
    row_mean = apply(exp_mat_norm, 2, function(x) { mean(x, na.rm=TRUE) } )
    exp_mat_norm = t(apply(exp_mat_norm, 1, "-", row_mean))
    
    gexp.norm = exp_mat_norm %>% as.data.frame() %>%
        tibble::rownames_to_column('gene') %>%
        inner_join(gtf_transcript, by = "gene") %>%
        filter(!(CHROM == 6 & gene_start < 33480577 & gene_end > 28510120)) %>%
        mutate(gene = droplevels(factor(gene, gtf_transcript$gene))) %>%
        mutate(gene_index = as.integer(gene)) %>% 
        mutate(CHROM = as.integer(CHROM)) %>%
        arrange(CHROM, gene)

    gexp.norm.long = gexp.norm %>% 
        reshape2::melt(
            id.var = c('gene', 'gene_index', 'region', 'gene_start', 'gene_end', 'CHROM', 'gene_length'),
            variable.name = 'cell',
            value.name = 'exp')    

    gexp.norm.long = gexp.norm.long %>% 
        group_by(cell, CHROM) %>%
        mutate(exp_rollmean = caTools::runmean(exp, k = window, align = "center")) %>%
        ungroup()
    
    return(list('gexp.norm.long' = gexp.norm.long, 'gexp.norm' = gexp.norm, 'exp_mat' = exp_mat))
}

#' @export
get_bulk = function(count_mat, lambdas_ref, df_allele, gtf_transcript, genetic_map, min_depth = 0, lambda = 2, verbose = TRUE) {

    if (nrow(df_allele) == 0) {
        stop('empty allele dataframe - check cell names')
    }

    Y_obs = rowSums(count_mat)

    fit = fit_multi_ref(Y_obs, lambdas_ref, sum(Y_obs), gtf_transcript)

    if (verbose) {
        message('Fitted reference proportions:')
        message(paste0(paste0(names(fit$w), ':', signif(fit$w, 2)), collapse = ','))
    }
    
    gexp_bulk = process_exp_fc(
        count_mat,
        fit$lambdas_bar,
        gtf_transcript,
        verbose = verbose
    )$bulk %>%
    filter((logFC < 5 & logFC > -5) | Y_obs == 0) %>%
    mutate(w = paste0(paste0(names(fit$w), ':', signif(fit$w, 2)), collapse = ','))

    allele_bulk = get_allele_bulk(df_allele, genetic_map, lambda = lambda, min_depth = min_depth)
            
    bulk = combine_bulk(
        allele_bulk = allele_bulk,
        gexp_bulk = gexp_bulk
    )

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

# phase switch probablity as a function of genetic distance
switch_prob_cm = function(d, lambda = 1, min_p = 1e-10) {
    p = (1-exp(-2*lambda*d))/2
    p = pmax(p, min_p)
    p = ifelse(is.na(d), 0, p)
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

#' @export
get_allele_bulk = function(df_allele, genetic_map, lambda = 1, min_depth) {
    pseudobulk = df_allele %>%
        filter(GT %in% c('1|0', '0|1')) %>%
        group_by(snp_id, CHROM, POS, REF, ALT, GT, gene) %>%
        summarise(
            AD = sum(AD),
            DP = sum(DP),
            AR = AD/DP,
            .groups = 'drop'
        ) %>%
        arrange(CHROM, POS) %>%
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

#' @export
combine_bulk = function(allele_bulk, gexp_bulk) {
    
    Obs = allele_bulk %>% 
        full_join(
            gexp_bulk,
            by = c("CHROM", "gene")
        ) %>%
        mutate(
            snp_id = ifelse(is.na(snp_id), gene, snp_id),
            # levels will be missing if not expressed
            gene = factor(gene, levels(gexp_bulk$gene)),
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
    Obs = Obs %>% 
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
    
    return(Obs)
    
}


#' Utility function to make reference gene expression profiles
#' @param count_mat raw count matrices where rownames are genes and column names are cells
#' @param cell_annot dataframe with columns "cell" and "cell_type"
#' @return a matrix of reference gene expression profiles
#' @export
make_psbulk = function(count_mat, cell_annot, verbose = T) {

    cell_annot = cell_annot %>% mutate(cell_type = factor(cell_type))

    cells = intersect(colnames(count_mat), cell_annot$cell)
    
    count_mat = count_mat[,cells]

    cell_annot = cell_annot %>% filter(cell %in% cells)
    
    cell_dict = cell_annot %>% {setNames(.$cell_type, .$cell)} %>% droplevels

    cell_dict = cell_dict[cells]

    if (verbose) {
        message(table(cell_dict))
    }
   
    M = model.matrix(~ 0 + cell_dict) %>% magrittr::set_colnames(levels(cell_dict))
    count_mat_clust = count_mat %*% M
    exp_mat_clust = count_mat_clust %*% diag(1/colSums(count_mat_clust)) %>% magrittr::set_colnames(colnames(count_mat_clust))
    
    return(list('exp_mat' = as.matrix(exp_mat_clust), 'count_mat' = as.matrix(count_mat_clust)))
}

#' @export
fetch_results = function(out_dir, i = 2, max_cost = 150, verbose = F) {

    res = list()

    res[['mut_tree_file']] = glue('{out_dir}/mut_tree_{i}.gml')
    res[['joint_post']] = fread(glue('{out_dir}/joint_post_{i}.tsv'))
    res[['exp_post']] = fread(glue('{out_dir}/exp_post_{i}.tsv'))
    res[['allele_post']] = fread(glue('{out_dir}/allele_post_{i}.tsv'))
    # bulk0 = fread(glue('{out_dir}/bulk_subtrees_0.tsv'))
    res[['geno']] = fread(glue('{out_dir}/geno_{i}.tsv')) %>% tibble::column_to_rownames('V1')


    res[['tree_post']] = readRDS(glue('{out_dir}/tree_post_{i}.rds'))

    f = glue('{out_dir}/scistree_out_{i}.rds')

    if (file.exists(f)) {
        res[['scistree_out']] = readRDS(f)
    }
    
    f = glue('{out_dir}/segs_consensus_{i-1}.tsv')
    if (file.exists(f)) {
        segs_consensus = fread(f)
    } else {
        f = glue('{out_dir}/segs_consensus_{i}.tsv')
        segs_consensus = fread(f)
    }

    if ('cnv_states' %in% colnames(segs_consensus)) {
        segs_consensus = segs_consensus  %>% 
            mutate(cnv_states = str_split(cnv_states, '\\|'))
    }

    res[['segs_consensus']] = segs_consensus

    f = glue('{out_dir}/bulk_subtrees_{i}.tsv.gz')
    if (file.exists(f)) {
        res[['bulk_subtrees']] = fread(f)
    }

    res[['bulk_initial']] = fread(glue('{out_dir}/bulk_subtrees_0.tsv.gz'))
        # mutate(CHROM = factor(CHROM, unique(CHROM)))

    f = glue('{out_dir}/bulk_clones_{i}.tsv.gz')
    if (file.exists(f)) {
        res[['bulk_clones']] = fread(f)
    }

    f = glue('{out_dir}/mut_graph_{i}.rds')
    if (file.exists(f)) {
        res$G_m = readRDS(f)
    }

    # update tree
    f = glue('{out_dir}/tree_final_{i}.rds')
    if (file.exists(f)) {
        res$gtree = readRDS(f)
    }
    # } else {
        # res$gtree = mut_to_tree(res$tree_post$gtree, mut_nodes)
    # }

    f = glue('{out_dir}/clone_post_{i}.tsv')
    if (file.exists(f)) {
        res[['clone_post']] = fread(f)
    } else {
        res[['clone_post']] = cell_to_clone(res$gtree, res$exp_post, res$allele_post)
    }

    return(res)

}


########################### Analysis ############################

#' @export
analyze_bulk = function(
    bulk, t = 1e-5, gamma = 20, theta_min = 0.08, bal_cnv = TRUE, prior = NULL,
    exp_only = FALSE, allele_only = FALSE, run_hmm = TRUE, retest = TRUE, lambda = 1,
    phasing = TRUE, roll_phi = TRUE, verbose = TRUE, debug = FALSE, diploid_chroms = NULL,
    find_diploid = TRUE, classify_allele = FALSE
) {
    
    if (!is.numeric(t)) {
        stop('transition probability is not numeric')
    }

    # update transition probablity
    bulk = bulk %>% mutate(p_s = switch_prob_cm(inter_snp_cm, lambda = lambda))

    if (exp_only | allele_only) {
        bulk$diploid = TRUE
    } else if (!is.null(diploid_chroms)) {
        log_info(glue('Using diploid chromosomes given: {paste0(diploid_chroms, collapse = ",")}'))
        bulk = bulk %>% mutate(diploid = CHROM %in% diploid_chroms)
    } else if (find_diploid) {
        out = find_common_diploid(bulk, gamma = gamma, t = t)
        bulk = out$bulks
        bal_cnv = out$bamp
    } else if (!'diploid' %in% colnames(bulk)) {
        stop('Must define diploid region if not given')
    }

    # fit expression baseline
    fit = bulk %>%
        filter(!is.na(Y_obs)) %>%
        filter(logFC < 8 & logFC > -8) %>%
        filter(diploid) %>%
        {fit_lnpois(.$Y_obs, .$lambda_ref, unique(.$d_obs), approx = F)}
        
    bulk = bulk %>% mutate(mu = fit@coef[1], sig = fit@coef[2])

    if (run_hmm) {
        bulk = bulk %>% 
            group_by(CHROM) %>%
            mutate(state = 
                run_hmm_mv_inhom(
                    pAD = pAD,
                    DP = DP, 
                    p_s = p_s,
                    Y_obs = Y_obs, 
                    lambda_ref = lambda_ref, 
                    d_total = na.omit(unique(d_obs)),
                    phi_neu = 1,
                    phi_amp = 2^(0.25),
                    phi_del = 2^(-0.25),
                    phi_bamp = 2^(0.25),
                    phi_bdel = 2^(-0.25),
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
                    phasing = phasing,
                    exp_model = 'lnpois'
                )
            ) %>% 
            mutate(cnv_state = str_remove(state, '_down|_up')) %>%
            annot_segs %>%
            smooth_segs %>%
            annot_segs %>%
            ungroup() %>%
            group_by(CHROM) %>%
            mutate(seg_start = min(POS), seg_end = max(POS)) %>%
            ungroup()
    }

    # rolling theta estimates
    bulk = annot_theta_roll(bulk)
    
    # retest CNVs
    if (retest & (!exp_only)) {
        
        if (verbose) {
            message('Retesting CNVs..')
        }

        segs_post = retest_cnv(bulk, gamma = gamma, exp_model = 'lnpois', allele_only = allele_only)
        
        bulk = bulk %>% 
            select(-any_of(colnames(segs_post)[!colnames(segs_post) %in% c('seg', 'CHROM')])) %>%
            left_join(segs_post, by = c('seg', 'CHROM')) %>%
            mutate(
                cnv_state_post = tidyr::replace_na(cnv_state_post, 'neu'),
                cnv_state = tidyr::replace_na(cnv_state, 'neu')
            ) %>%
            mutate(state_post = ifelse(
                cnv_state_post %in% c('amp', 'del', 'loh') & (!cnv_state %in% c('bamp', 'bdel')),
                paste0(cnv_state_post, '_', str_extract(state, 'up_1|down_1|up_2|down_2|up|down|1_up|2_up|1_down|2_down')),
                cnv_state_post
            )) %>%
            mutate(state_post = str_remove(state_post, '_NA'))

        # note that this uses restest cnv states
        bulk = bulk %>% classify_alleles()

        bulk = annot_theta_roll(bulk)
        
    } else {
        bulk = bulk %>% mutate(state_post = state, cnv_state_post = cnv_state)
    }
    
    if (verbose) {
        message('Finishing..')
    }

    # annotate phi MLE for all segments
    bulk = bulk %>%
        group_by(seg) %>%
        mutate(
            approx_lik_exp(
                Y_obs[!is.na(Y_obs)], lambda_ref[!is.na(Y_obs)], unique(na.omit(d_obs)),
                mu = mu[!is.na(Y_obs)],
                sig = sig[!is.na(Y_obs)]
            )
        ) %>%
        ungroup()

    if (roll_phi) {
        # rolling phi estimates
        bulk = bulk %>% 
            select(-any_of('phi_mle_roll')) %>%
            left_join(
                bulk %>% 
                    group_by(CHROM) %>%
                    filter(!is.na(Y_obs)) %>%
                    filter(Y_obs > 0) %>%
                    mutate(
                        phi_mle_roll = phi_hat_lnpois_roll(Y_obs, lambda_ref, d_obs, mu, sig, h = 50)
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

# classify alleles using viterbi and forward-backward
#' @export
classify_alleles = function(Obs) {

    if (all(Obs$cnv_state_post %in% c('neu', 'bdel', 'bamp'))) {
        return(Obs)
    }
    
    allele_post = Obs %>%
        filter(!cnv_state_post %in% c('neu', 'bdel', 'bamp')) %>%
        filter(!is.na(AD)) %>%
        group_by(CHROM, seg) %>%
        filter(n() > 1) %>%
        mutate(
            p_up = forward_back_allele(get_allele_hmm(pAD, DP, p_s, theta = unique(theta_mle), gamma = 20)),
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

    Obs = Obs %>% 
            select(-any_of(colnames(allele_post)[!colnames(allele_post) %in% c('snp_id')])) %>%
            left_join(
                allele_post,
            by = c('snp_id')
        )

    return(Obs)
}

annot_theta_roll = function(Obs) {   

    Obs = Obs %>% 
        mutate(haplo_theta_min = case_when(
            str_detect(state, 'up') ~ 'major',
            str_detect(state, 'down') ~ 'minor',
            T ~ ifelse(pBAF > 0.5, 'major', 'minor')
        )) %>% 
        mutate(
            major_count = ifelse(haplo_theta_min == 'major', pAD, DP - pAD),
            minor_count = DP - major_count
        )

    Obs = Obs %>%
        select(-any_of('theta_hat_roll')) %>%
        left_join(
            Obs %>%
                group_by(CHROM) %>%
                filter(!is.na(pAD)) %>%
                mutate(theta_hat_roll = theta_hat_roll(major_count, minor_count, h = 100)) %>%
                select(theta_hat_roll, CHROM, snp_id),
            by = c('CHROM', 'snp_id')
        ) %>%
        group_by(CHROM) %>%
        mutate(theta_hat_roll = zoo::na.locf(theta_hat_roll, na.rm=FALSE)) %>%
        ungroup()

    return(Obs)
}

annot_segs = function(Obs) {

    Obs = Obs %>% 
            group_by(CHROM) %>%
            arrange(CHROM, snp_index) %>%
            mutate(boundary = c(0, cnv_state[2:length(cnv_state)] != cnv_state[1:(length(cnv_state)-1)])) %>%
            group_by(CHROM) %>%
            mutate(seg = paste0(CHROM, letters[cumsum(boundary)+1])) %>%
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

    return(Obs)
}

annot_haplo_segs = function(Obs) {

    Obs = Obs %>% 
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

    return(Obs)
}

smooth_segs = function(bulk, min_genes = 10) {
    bulk %>% group_by(seg) %>%
        mutate(
            cnv_state = ifelse(n_genes <= min_genes, NA, cnv_state)
        ) %>%
        ungroup() %>%
        group_by(CHROM) %>%
        mutate(cnv_state = zoo::na.locf(cnv_state, fromLast = F, na.rm=FALSE)) %>%
        mutate(cnv_state = zoo::na.locf(cnv_state, fromLast = T, na.rm=FALSE)) %>%
        ungroup()
}

#' @export
annot_consensus = function(bulk, segs_consensus) {

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
        set_names(c('marker_index', 'seg_index')) %>%
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
        select(-any_of(c('sample')))
    
    bulk = bulk %>% 
        select(-any_of(colnames(marker_seg)[!colnames(marker_seg) %in% c('snp_id')])) %>% 
        inner_join(marker_seg, by = c('snp_id'))  %>%
        mutate(seg = seg_cons) %>%
        mutate(seg = factor(seg, gtools::mixedsort(unique(seg)))) 

    return(bulk)
}

simes_p = function(p.vals, n_dim) {
    n_dim * min(sort(p.vals)/seq_along(p.vals))
}

# multi-sample generalization
#' @export
find_common_diploid = function(
    bulks, grouping = 'clique', gamma = 20, theta_min = 0.08, t = 1e-5, fc_min = 1.25, alpha = 1e-4, 
    bal_consensus = NULL, ncores = 5, debug = F, verbose = T) {

    if (!'sample' %in% colnames(bulks)) {
        bulks$sample = 1
    }
    
    # define balanced regions in each sample
    bulks = mclapply(
        bulks %>% split(.$sample),
        mc.cores = ncores,
        function(bulk) {
            
            bulk %>% 
                group_by(CHROM) %>%
                mutate(state = 
                    run_hmm_inhom2(
                        pAD = pAD,
                        DP = DP, 
                        p_s = p_s,
                        t = t,
                        theta_min = theta_min,
                        gamma = gamma
                    )
                ) %>% ungroup() %>%
                mutate(cnv_state = str_remove(state, '_down|_up')) %>%
                annot_segs()

        }) %>%
        bind_rows()

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

    bulks_bal = bulks %>% filter(cnv_state == 'neu') %>% filter(!is.na(lnFC))
    segs_bal = bulks_bal %>% count(seg) %>% filter(n > 50) %>% pull(seg) %>% as.character

    if (length(segs_bal) == 0) {
        msg = 'No balanced segments, using all segments as baseline'
        log_warn(msg)
        diploid_segs = bulks %>% pull(seg) %>% unique
        bamp = TRUE
        tests = data.frame()
        FC = data.frame()
        G = NULL
        test_dat = data.frame()
    } else if (length(segs_bal) == 1) {
        diploid_segs = segs_bal
        bamp = FALSE
        tests = data.frame()
        FC = data.frame()
        G = NULL
        test_dat = data.frame()
    } else {
        test_dat = bulks_bal %>%
            select(gene, seg, lnFC, sample) %>%
            reshape2::dcast(seg+gene ~ sample, value.var = 'lnFC') %>%
            na.omit() %>%
            mutate(seg = as.character(seg)) %>%
            select(-gene)
        
        # tests = t(combn(segs_bal, 2)) %>% 
        #     as.data.frame() %>%
        #     set_names(c('i', 'j')) %>%
        #     rowwise() %>%
        #     mutate(
        #         p = Hotelling::hotelling.test(
        #             .~seg,
        #             data = test_dat,
        #             pair = c(i, j)
        #         )$pval
        #     ) %>%
        #     ungroup() %>%
        #     mutate(p = p.adjust(p))

        # sime's p
        samples = unique(bulks_bal$sample)

        tests = lapply(
                samples,
                function(sample) {
                    t(combn(segs_bal, 2)) %>% 
                        as.data.frame() %>%
                        set_names(c('i', 'j')) %>%
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

        # E = tests %>% filter(q > alpha | delta_max < log(fc_min))
        E = tests %>% filter(q > alpha)

        G = igraph::graph_from_data_frame(d=E, vertices=V, directed=F)

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
            tibble::column_to_rownames('sample')

        }

        # choose diploid clique based on min FC in case there's a tie
        diploid_cluster = as.integer(names(sort(apply(FC[,Modes(apply(FC, 1, which.min)),drop=F], 2, min))))[1]

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
    
    return(list('bamp' = bamp, 'bulks' = bulks, 'diploid_segs' = diploid_segs, 'segs_consensus' = segs_consensus, 'G' = G, 'tests' = tests, 'test_dat' = test_dat, 'FC' = FC, 'cliques' = cliques))
}

get_segs_neu = function(bulks_all) {
    segs_neu = bulks_all %>% filter(cnv_state == "neu") %>%
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

# get all internal nodes, for cluster tree
get_internal_nodes = function(den, node) {

    membership = data.frame(
        cluster = dendextend::get_leaves_attr(den, attribute = 'label'),
        node = node
    )
    
    if (is.leaf(den)) {
        return(data.frame())
    }
                
    sub_dens = dendextend::get_subdendrograms(den, k = 2)
    
    membership_l = get_internal_nodes(
        sub_dens[[1]],
        paste0(node, '.', 1)
    )
    
    membership_r = get_internal_nodes(
        sub_dens[[2]],
        paste0(node, '.', 2)
    )
        
    return(rbind(membership, membership_l, membership_r))
}

phi_hat_seg = function(y_vec, lambda_vec, depth) {
    sum(y_vec)/(sum(lambda_vec) * depth)
}

phi_hat_roll = function(y_vec, lambda_vec, depth, h) {
    n = length(y_vec)
    sapply(
        1:n,
        function(c) {
            slice = max(c - h - 1, 1):min(c + h, n)
            phi_hat_seg(y_vec[slice], lambda_vec[slice], depth)
        }
    )
}

# homoskd
phi_hat_lnpois = function(Y_obs, lambda_ref, d, mu, sig) {
    logFC = log((Y_obs/d)/lambda_ref)
    exp(mean(logFC - mu))
}

# heterskd generalization
phi_hat_lnpois = function(Y_obs, lambda_ref, d, mu, sig) { 
        
    if (length(mu) == 1 & length(sig) == 1) {
        
        mu = rep(mu, length(Y_obs))
        sig = rep(sig, length(Y_obs))
    }
    
    logFC = log((Y_obs/d)/lambda_ref)
    weights = 1/sig^2
    exp(sum((logFC - mu) * weights)/sum(weights))
}

phi_hat_lnpois_roll = function(Y_obs, lambda_ref, d_obs, mu, sig, h) {
    n = length(Y_obs)
    
    if (length(mu) == 1 & length(sig) == 1) {
        mu = rep(mu, n)
        sig = rep(sig, n)
    }
    
    sapply(
        1:n,
        function(c) {
            slice = max(c - h - 1, 1):min(c + h, n)
            phi_hat_lnpois(Y_obs[slice], lambda_ref[slice], unique(d_obs), mu[slice], sig[slice])
        }
    )
}


theta_hat_seg = function(major_count, minor_count) {
    major_total = sum(major_count)
    minor_total = sum(minor_count)
    MAF = major_total/(major_total + minor_total)
    return(MAF - 0.5)
}

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

# vectorized version
phi_hat_seg_sc = function(y_mat, lambda_total, d_vec) {
    colSums(y_mat)/(lambda_total * d_vec)
}

phi_hat_roll_sc = function(y_mat, lambda_vec, d_vec, h) {
    m = nrow(y_mat)
    sapply(
        1:m,
        function(c) {
            slice = max(c - h - 1, 1):min(c + h, m)
            lambda_total = sum(lambda_vec[slice])
            phi_hat_seg_sc(y_mat[slice,], lambda_total, d_vec)
        }
    )
}

#' @export
calc_phi_mle = function(Y_obs, lambda_ref, d, alpha, beta, lower = 0.2, upper = 10) {
    
    if (length(Y_obs) == 0) {
        return(1)
    }
        
    start = max(min(1, upper), lower)
    
    res = stats4::mle(
        minuslogl = function(phi) {
            res = -l_gpois(Y_obs, lambda_ref, d, alpha, beta, phi)
            return(res)
        },
        start = start,
        lower = lower,
        upper = upper
    )
    
    return(res@coef)
}


calc_phi_mle_loglink = function(Y_obs, lambda_ref, d, alpha, beta) {
    
    if (length(Y_obs) == 0) {
        return(1)
    }
            
    res = stats4::mle(
        minuslogl = function(beta_0) {
            res = -l_gpois(Y_obs, lambda_ref, d, alpha, beta, exp(beta_0-1))
            return(res)
        },
        start = 0
    )
    
    return(res@coef)
}

phi_mle_roll = function(Y_obs, lambda_ref, alpha, beta, d_obs, h) {
    n = length(Y_obs)
    sapply(
        1:n,
        function(c) {
            slice = max(c - h - 1, 1):min(c + h, n)
            calc_phi_mle(Y_obs[slice], lambda_ref[slice], unique(d_obs), alpha[slice], beta[slice])
        }
    )
}

l_bbinom = function(AD, DP, alpha, beta) {
    sum(dbbinom(AD, DP, alpha, beta, log = TRUE))
}

#' @export
fit_bbinom = function(AD, DP) {

    fit = stats4::mle(
        minuslogl = function(alpha, beta) {
            -l_bbinom(AD, DP, alpha, beta)
        },
        start = c(5, 5),
        lower = c(0, 0)
    )

    alpha = fit@coef[1]
    beta = fit@coef[2]
    
    return(fit)
}

l_gpois = function(Y_obs, lambda_ref, d, alpha, beta, phi = 1) {
    sum(dgpois(Y_obs, shape = alpha, rate = beta/(phi * d * lambda_ref), log = TRUE))
}

#' @export
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

# Poisson LogNormal model
l_lnpois = function(Y_obs, lambda_ref, d, mu, sig, phi = 1) {
    if (any(sig <= 0)) {stop(glue('negative sigma. value: {sig}'))}
    if (length(sig) == 1) {sig = rep(sig, length(Y_obs))}
    sum(log(dpoilog(Y_obs, mu + log(phi * d * lambda_ref), sig)))
}

#' @export
fit_lnpois = function(Y_obs, lambda_ref, d, approx = F) {

    Y_obs = Y_obs[lambda_ref > 0]
    lambda_ref = lambda_ref[lambda_ref > 0]

    l_func = ifelse(approx, l_lnpois_mix, l_lnpois)
    
    fit = stats4::mle(
        minuslogl = function(mu, sig) {
            if (sig < 0) {stop('optim is trying negative sigma')}
            -l_func(Y_obs, lambda_ref, d, mu, sig)
        },
        start = c(0, 1),
        lower = c(-Inf, 0.01),
        control = list('trace' = FALSE)
    )

    return(fit)
}

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

binary_entropy = function(p) {
    H = -p*log2(p)-(1-p)*log2(1-p)
    H[is.na(H)] = 0
    return(H)
}

#' @export
retest_cnv = function(bulk, exp_model = 'lnpois', gamma = 20, allele_only = FALSE) {
    
    G = c('20' = 1/5, '10' = 1/5, '21' = 1/10, '31' = 1/10, '22' = 1/5, '00' = 1/5)
    
    theta_min = 0.065
    delta_phi_min = 0.15

    if (allele_only) {
        segs_post = bulk %>% 
            filter(cnv_state != 'neu') %>%
            group_by(CHROM, seg, cnv_state) %>%
            summarise(
                n_genes = length(na.omit(unique(gene))),
                n_snps = sum(!is.na(pAD)),
                seg_start = min(POS),
                seg_end = max(POS),
                theta_hat = theta_hat_seg(major_count[!is.na(major_count)], minor_count[!is.na(minor_count)]),
                approx_lik_ar(pAD[!is.na(pAD)], DP[!is.na(pAD)], p_s[!is.na(pAD)], gamma = unique(gamma), start = theta_hat),
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
            filter(cnv_state != 'neu') %>%
            group_by(CHROM, seg, cnv_state) %>%
            summarise(
                n_genes = length(na.omit(unique(gene))),
                n_snps = sum(!is.na(pAD)),
                seg_start = min(POS),
                seg_end = max(POS),
                theta_hat = theta_hat_seg(major_count[!is.na(major_count)], minor_count[!is.na(minor_count)]),
                approx_lik_ar(pAD[!is.na(pAD)], DP[!is.na(pAD)], p_s[!is.na(pAD)], gamma = unique(gamma), start = theta_hat),
                L_y_n = pnorm.range(0, theta_min, theta_mle, theta_sigma),
                L_y_d = pnorm.range(theta_min, 0.499, theta_mle, theta_sigma),
                L_y_a = pnorm.range(theta_min, 0.375, theta_mle, theta_sigma),
                approx_lik_exp(
                    Y_obs[!is.na(Y_obs)], lambda_ref[!is.na(Y_obs)], unique(na.omit(d_obs)),
                    alpha = alpha[!is.na(Y_obs)],
                    beta = beta[!is.na(Y_obs)],
                    mu = mu[!is.na(Y_obs)],
                    sig = sig[!is.na(Y_obs)],
                    model = exp_model
                ),
                L_x_n = pnorm.range(1 - delta_phi_min, 1 + delta_phi_min, phi_mle, phi_sigma),
                L_x_d = pnorm.range(0.1, 1 - delta_phi_min, phi_mle, phi_sigma),
                L_x_a = pnorm.range(1 + delta_phi_min, 3, phi_mle, phi_sigma),
                Z = sum(G['20'] * L_x_n * L_y_d,
                        G['10'] * L_x_d * L_y_d,
                        G['21'] * L_x_a * L_y_a,
                        G['31'] * L_x_a * L_y_a,
                        G['22'] * L_x_a * L_y_n, 
                        G['00'] * L_x_d * L_y_n),
                p_loh = (G['20'] * L_x_n * L_y_d)/Z,
                p_amp = ((G['31'] + G['21']) * L_x_a * L_y_a)/Z,
                p_del = (G['10'] * L_x_d * L_y_d)/Z,
                p_bamp = (G['22'] * L_x_a * L_y_n)/Z,
                p_bdel = (G['00'] * L_x_d * L_y_n)/Z,
                LLR_x = calc_exp_LLR(
                    Y_obs[!is.na(Y_obs)],
                    lambda_ref[!is.na(Y_obs)], 
                    unique(na.omit(d_obs)),
                    phi_mle,
                    alpha = alpha[!is.na(Y_obs)],
                    beta = beta[!is.na(Y_obs)],
                    mu = mu[!is.na(Y_obs)],
                    sig = sig[!is.na(Y_obs)],
                    model = exp_model
                ),
                LLR_y = calc_allele_LLR(pAD[!is.na(pAD)], DP[!is.na(pAD)], p_s[!is.na(pAD)], theta_mle, gamma = unique(gamma)),
                LLR = LLR_x + LLR_y,
                .groups = 'drop'
            ) %>%
            rowwise() %>%
            mutate(cnv_state_post = c('loh', 'amp', 'del', 'bamp', 'bdel')[
                which.max(c(p_loh, p_amp, p_del, p_bamp, p_bdel))
            ]) %>%
            ungroup()
    }

    return(segs_post)
}

approx_lik_ar = function(pAD, DP, p_s, lower = 0.001, upper = 0.499, start = 0.25, gamma = 20) {
    
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


approx_maxlik_ar = function(pAD, DP, p_s, lower = 0.001, upper = 0.499, start = 0.25, gamma = 20) {
    
    if (length(pAD) <= 10) {
        return(tibble('theta_mle' = 0, 'theta_sigma' = 0))
    }

    fit = optim(
        start, 
        function(theta) {-calc_allele_maxlik(pAD, DP, p_s, theta, gamma)},
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
    
    return(tibble('theta_mle' = mu, 'theta_sigma' = sigma, 'l_max' = -fit$value))
}

approx_lik_exp = function(Y_obs, lambda_ref, d, alpha = NULL, beta = NULL, mu = NULL, sig = NULL, model = 'lnpois', lower = 0.2, upper = 10, start = 1) {
    
    if (length(Y_obs) == 0) {
        return(tibble('phi_mle' = 1, 'phi_sigma' = 0))
    }
    
    start = max(min(1, upper), lower)

    if (model == 'gpois') {
        l = function(phi) {l_gpois(Y_obs, lambda_ref, d, alpha, beta, phi = phi)}
    } else {
        l = function(phi) {l_lnpois(Y_obs, lambda_ref, d, mu, sig, phi = phi)}
    }

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

approx_post_exp = function(Y_obs, lambda_ref, d, alpha, beta, sigma_0, lower = 0.2, upper = 10, start = 1) {
    
    if (length(Y_obs) == 0) {
        return(tibble('phi_map' = 1, 'phi_sigma' = 0))
    }
    
    start = max(min(1, upper), lower)
    
    fit = optim(
        start,
        function(phi) {
            - (dlnorm(phi, meanlog = 0, sdlog = sigma_0, log = TRUE) + l_gpois(Y_obs, lambda_ref, d, alpha, beta, phi = phi))
        },
        method = 'L-BFGS-B',
        lower = lower,
        upper = upper,
        hessian = TRUE
    )

    mu = fit$par
    sigma = sqrt(as.numeric(1/(fit$hessian)))
    
    return(tibble('phi_map' = mu, 'phi_sigma' = sigma))
}

pnorm.range = function(lower, upper, mu, sd) {
    if (sd == 0) {
        return(1)
    }
    pnorm(upper, mu, sd) - pnorm(lower, mu, sd)
}

# for a cell tree
get_internal_nodes_sc = function(den, node, labels) {

    membership = data.frame(
        cell = dendextend::get_leaves_attr(den, attribute = 'label'),
        node = node
    )
    
    is_leaf = length(unique(labels[dendextend::get_leaves_attr(den, attribute = 'label')])) == 1
    
    if (is_leaf) {
        return(data.frame())
    }
                
    sub_dens = dendextend::get_subdendrograms(den, k = 2)
    
    membership_l = get_internal_nodes_sc(
        sub_dens[[1]],
        paste0(node, '.', 1),
        labels
    )
    
    membership_r = get_internal_nodes_sc(
        sub_dens[[2]],
        paste0(node, '.', 2),
        labels
    )
        
    return(rbind(membership, membership_l, membership_r))
}

get_nodes_celltree = function(hc, clusters) {
        
    # internal nodes
    nodes = get_internal_nodes_sc(as.dendrogram(hc), '0', clusters)
    
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
        map(function(node){list(sample = unique(node$node), members = unique(node$cluster), cells = node$cell, size = length(node$cell))})
    
    return(nodes)
    
}


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

calc_theta_mle = function(pAD, DP, p_s, lower = 0.01, upper = 0.49, start = 0.25, gamma = 20) {
    
    if(length(pAD) == 0) {
        return(0)
    }
    
    res = optim(
        start, 
        function(theta) {-calc_allele_lik(pAD, DP, p_s, theta, gamma)},
        method = 'L-BFGS-B',
        lower = lower,
        upper = upper,
        control = list('factr' = 0.001/(.Machine$double.eps))
    )
    
    return(res$par)
}

calc_allele_LLR = function(pAD, DP, p_s, theta_mle, theta_0 = 0, gamma = 20) {
    if (length(pAD) <= 1) {
        return(0)
    }
    l_1 = calc_allele_lik(pAD, DP, p_s, theta = theta_mle, gamma = gamma) 
    l_0 = l_bbinom(pAD, DP, gamma*0.5, gamma*0.5)
    return(l_1 - l_0)
}

calc_exp_LLR = function(Y_obs, lambda_ref, d, phi_mle, alpha = NULL, beta = NULL, mu = NULL, sig = NULL, model = 'gpois') {
    if (length(Y_obs) == 0) {
        return(0)
    }
    if (model == 'gpois') {
        l_1 = l_gpois(Y_obs, lambda_ref, d, alpha, beta, phi = phi_mle)
        l_0 = l_gpois(Y_obs, lambda_ref, d, alpha, beta, phi = 1)
    } else {
        l_1 = l_lnpois(Y_obs, lambda_ref, d, mu, sig, phi = phi_mle)
        l_0 = l_lnpois(Y_obs, lambda_ref, d, mu, sig, phi = 1)
    }
    
    return(l_1 - l_0)
}

sc_exp_post = function(exp_sc) {

    eps = 0.15
    phi_max = 3.5

    prior_neu = (2*eps)^(-1) * (1/2)
    prior_amp = (phi_max - (1 + eps))^(-1) * (1/4)
    prior_del = (1 - eps)^(-1) * (1/4)
    
    depth_obs = sum(exp_sc$Y_obs)
    exp_sc = exp_sc %>% filter(lambda_ref > 0)

    fit = exp_sc %>% filter(cnv_state %in% c('neu', 'loh')) %>% {fit_gpois(.$Y_obs, .$lambda_ref, depth_obs)}
    
    exp_sc = exp_sc %>% mutate(alpha = fit@coef[1], beta = fit@coef[2])
        
    res = exp_sc %>% 
        group_by(seg) %>%
        filter(cnv_state != 'neu') %>%
        summarise(
            n = n(),
            alpha = unique(alpha),
            beta = unique(beta),
            cnv_state = unique(cnv_state),
            approx_lik_exp(Y_obs, lambda_ref, depth_obs, alpha, beta),
            L_del = pnorm.range(0, 1-eps, phi_mle, phi_sigma) * prior_del,
            L_neu = pnorm.range(1-eps, 1+eps, phi_mle, phi_sigma) * prior_neu,
            L_amp = pnorm.range(1+eps, phi_max, phi_mle, phi_sigma) * prior_amp,
            p_del = L_del/(L_del + L_neu + L_amp),
            p_amp = L_amp/(L_del + L_neu + L_amp),
            p_neu = L_neu/(L_del + L_neu + L_amp)
        )
        
    return(res)
}

rename_seg = function(seg) {
    chr = str_split(seg, '_')[[1]][1]
    id = as.integer(str_split(seg, '_')[[1]][2])
    paste0(chr, letters[id])
}

#' @export
robust.system <- function(cmd) {
  stderrFile = tempfile(pattern="R_robust.system_stderr", fileext=as.character(Sys.getpid()))
  stdoutFile = tempfile(pattern="R_robust.system_stdout", fileext=as.character(Sys.getpid()))

  retval = list()
  retval$exitStatus = system(paste0(cmd, " 2> ", shQuote(stderrFile), " > ", shQuote(stdoutFile)))
  retval$stdout = readLines(stdoutFile)
  retval$stderr = readLines(stderrFile)

  unlink(c(stdoutFile, stderrFile))
  return(retval)
}

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
get_bulk = function(count_mat, lambdas_ref, df_allele, gtf_transcript, genetic_map, min_depth = 0, lambda = 0.52, verbose = TRUE) {

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


# phase switch probablity as a function of genomic distance
switch_prob = function(x, pad = 0.01) {
    1-pmax(pmin(2.8 - 0.38 * log10(x), 1 - pad), 0.5 + pad)
}

switch_prob_bp = function(x, lambda = 1.0728, min_p = 0) {
    p = 0.5*(1 - exp(-lambda * x/1e6))
    pmax(p, min_p)
}
# phase switch probablity as a function of genetic distance
# switch_prob_cm = function(x, lambda = 2, min_p = 1e-10) {
#     p = 0.5*(1 - exp(-lambda * x))
#     pmax(p, min_p)
# }
switch_prob_cm = function(d, lambda = 2, min_p = 1e-10) {
    p = (1-exp(-2*lambda*d))/2
    pmax(p, min_p)
}

fit_switch_prob = function(y, d) {
    
    eta = function(d, lambda, min_p = 1e-10) {
        pmax(switch_prob(d, lambda), min_p)
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
        set_names(c('marker_index', 'map_index')) %>%
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
get_allele_bulk = function(df_allele, genetic_map, lambda = 0.52, min_depth) {
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
            p_s = switch_prob_cm(inter_snp_cm, lambda)
        ) %>%
        # mutate(
        #     inter_snp_dist = c(NA, POS[2:length(POS)] - POS[1:(length(POS)-1)]),
        #     p_s = switch_prob_bp(inter_snp_dist)
        # ) %>%
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
    
    count_mat = count_mat %>% extract(,cells)
    cell_annot = cell_annot %>% filter(cell %in% cells)
    
    cell_dict = cell_annot %>% {setNames(.$cell_type, .$cell)} %>% droplevels

    cell_dict = cell_dict[cells]

    if (verbose) {
        message(table(cell_dict))
    }
   
    M = model.matrix(~ 0 + cell_dict) %>% set_colnames(levels(cell_dict))
    count_mat_clust = count_mat %*% M
    exp_mat_clust = count_mat_clust %*% diag(1/colSums(count_mat_clust)) %>% set_colnames(colnames(count_mat_clust))
    
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
        res[['segs_consensus']] = fread(f)
    } else {
        f = glue('{out_dir}/segs_consensus_{i}.tsv')
        res[['segs_consensus']] = fread(f)
    }

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
    Obs, t = 1e-5, gamma = 20, theta_min = 0.08, bal_cnv = TRUE, prior = NULL,
    exp_only = FALSE, allele_only = FALSE, retest = TRUE, hskd = TRUE,
    phasing = TRUE, roll_phi = TRUE, verbose = TRUE, debug = FALSE, diploid_chroms = NULL,
    find_diploid = TRUE, classify_allele = FALSE
) {
    
    if (!is.numeric(t)) {
        stop('transition probability is not numeric')
    }

    if (exp_only) {
        Obs$diploid = TRUE
    } else if (!is.null(diploid_chroms)) {
        log_info(glue('Using diploid chromosomes given: {paste0(diploid_chroms, collapse = ",")}'))
        Obs = Obs %>% mutate(diploid = CHROM %in% diploid_chroms)
    } else if (find_diploid) {
        out = find_common_diploid(Obs, gamma = gamma, t = t)
        Obs = out$bulks
        bal_cnv = out$bamp
    } else if (!'diploid' %in% colnames(Obs)) {
        stop('Must define diploid region if not given')
    }

    fit = Obs %>%
        filter(!is.na(Y_obs)) %>%
        filter(logFC < 8 & logFC > -8) %>%
        filter(diploid) %>%
        {fit_lnpois(.$Y_obs, .$lambda_ref, unique(.$d_obs), approx = F)}
        
    Obs = Obs %>% mutate(mu = fit@coef[1], sig = fit@coef[2])
    
    Obs = Obs %>% 
        mutate(gamma = gamma) %>%
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

    # rolling theta estimates
    Obs = annot_theta_roll(Obs)
    
    # retest CNVs
    if (retest & (!exp_only)) {
        
        if (verbose) {
            message('Retesting CNVs..')
        }

        segs_post = retest_cnv(Obs, exp_model = 'lnpois', allele_only = allele_only)
        
        Obs = Obs %>% 
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

        Obs = Obs %>% classify_alleles()

        Obs = annot_theta_roll(Obs)
        
    } else {
        Obs = Obs %>% mutate(state_post = state, cnv_state_post = cnv_state)
    }
    
    if (verbose) {
        message('Finishing..')
    }

    if (roll_phi) {
        # rolling phi estimates
        Obs = Obs %>% 
            select(-any_of('phi_mle_roll')) %>%
            left_join(
                Obs %>% 
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
    
    return(Obs)
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
            # mpc = HiddenMarkov::Viterbi(get_allele_hmm(pAD, DP, p_s, theta = unique(theta_mle), gamma = 20))$y,
            # haplo_mpc = case_when(
            #     mpc == 1 & GT == '1|0' ~ 'major',
            #     mpc == 1 & GT == '0|1' ~ 'minor',
            #     mpc == 2 & GT == '1|0' ~ 'minor',
            #     mpc == 2 & GT == '0|1' ~ 'major'
            # ),
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
            mutate(seg = paste0(CHROM, '_', (cumsum(boundary)+1))) %>%
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
retest_cnv = function(bulk, exp_model = 'lnpois', allele_only = FALSE) {
    
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
        map(function(node){list(label = unique(node$node), members = unique(node$cluster), cells = node$cell, size = length(node$cell))})
    
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

########################### Visualization ############################

pal = RColorBrewer::brewer.pal(n = 8, 'Set1')


cnv_colors = c("neu" = "gray", "neu_up" = "gray", "neu_down" = "gray20",
        "del_up" = "royalblue", "del_down" = "darkblue", 
        "loh_up" = "darkgreen", "loh_down" = "olivedrab4",
        "amp_up" = "red", "amp_down" = "tomato3",
        "del_1_up" = "royalblue", "del_1_down" = "darkblue", 
        "loh_1_up" = "darkgreen", "loh_1_down" = "olivedrab4",
        "amp_1_up" = "red", "amp_1_down" = "tomato3",
        "del_2_up" = "royalblue", "del_2_down" = "darkblue", 
        "loh_2_up" = "darkgreen", "loh_2_down" = "olivedrab4",
        "amp_2_up" = "red", "amp_2_down" = "tomato3",
        "del_up_1" = "royalblue", "del_down_1" = "darkblue", 
        "loh_up_1" = "darkgreen", "loh_down_1" = "olivedrab4",
        "amp_up_1" = "red", "amp_down_1" = "tomato3",
        "del_up_2" = "royalblue", "del_down_2" = "darkblue", 
        "loh_up_2" = "darkgreen", "loh_down_2" = "olivedrab4",
        "amp_up_2" = "red", "amp_down_2" = "tomato3",
        "bamp" = "salmon", "bdel" = "skyblue",
        "amp" = "tomato3", "loh" = "olivedrab4", "del" = "royalblue", "neu2" = "gray30",
        "theta_up" = "darkgreen", "theta_down" = "olivedrab4",
        "theta_1_up" = "darkgreen", "theta_1_down" = "olivedrab4",
        "theta_2_up" = "darkgreen", "theta_2_down" = "olivedrab4",
        "theta_up_1" = "darkgreen", "theta_down_1" = "olivedrab4",
        "theta_up_2" = "darkgreen", "theta_down_2" = "olivedrab4",
        '0|1' = 'red', '1|0' = 'blue'
    )

#' @export
show_phasing = function(bulk, min_depth = 8, dot_size = 0.5, h = 50) {

    D = bulk %>% 
        filter(!is.na(pAD)) %>%
        group_by(CHROM) %>%
        mutate(pBAF = 1-pBAF, pAD = DP - pAD) %>%
        mutate(theta_ar_roll = theta_hat_roll(AD, DP-AD, h = h)) %>%
        mutate(theta_hat_roll = theta_hat_roll(pAD, DP-pAD, h = h)) %>%
        filter(DP >= min_depth) %>%
        mutate(snp_index = 1:n()) %>%
        mutate(dAR = 0.5+abs(AR-0.5)) %>%
        mutate(dAR_roll = caTools::runmean(abs(AR-0.5), align = 'center', k = 30))

    boundary = D %>% filter(boundary == 1) %>% pull(snp_index)
    
    p1 = D %>%
        mutate(state_post = 'neu') %>%
        ggplot(
            aes(x = snp_index, y = AR),
            na.rm=TRUE
        ) +
        geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'gray') +
        geom_point(
            aes(color = state_post),
            size = dot_size,
        ) +
        # geom_line(
        #     aes(x = snp_index, y = 0.5 + dAR_roll), color = 'red'
        # ) +
        # geom_line(
        #     aes(x = snp_index, y = 0.5 - dAR_roll), color = 'red'
        # ) +
        # geom_line(
        #     aes(x = snp_index, y = 0.5 + theta_ar_roll), color = 'red'
        # ) +
        theme_classic() +
        theme(
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
            strip.background = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank()
        ) +
        scale_color_manual(values = cnv_colors) +
        ylim(0,1) +
        facet_grid(.~CHROM, space = 'free_x', scale = 'free_x') +
        geom_vline(xintercept = boundary - 1, color = 'red', size = 0.5, linetype = 'dashed') +
        guides(color = 'none')

    p2 = D %>%
        mutate(state_post = 'neu') %>%
        ggplot(
            aes(x = snp_index, y = pBAF),
            na.rm=TRUE
        ) +
        geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'gray') +
        geom_point(
            aes(color = ifelse(theta_hat_roll > 0, 'loh_1_up', 'loh_1_down')),
            size = dot_size,
        ) +
        geom_line(
            aes(x = snp_index, y = 0.5 + theta_hat_roll), color = 'red'
        ) +
        theme_classic() +
        theme(
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
            strip.background = element_blank()
        ) +
        scale_color_manual(values = cnv_colors) +
        ylim(0,1) +
        facet_grid(.~CHROM, space = 'free_x', scale = 'free_x') +
        geom_vline(xintercept = boundary - 1, color = 'red', size = 0.5, linetype = 'dashed') +
        guides(color = 'none')

    (p1 / p2) + plot_layout(guides = 'auto')
}

#' @export
plot_psbulk = function(Obs, dot_size = 0.8, exp_limit = 2, min_depth = 10, theta_roll = FALSE, fc_correct = TRUE, allele_only = FALSE, phi_mle = FALSE, use_pos = FALSE, legend = TRUE) {

    if (!'state_post' %in% colnames(Obs)) {
        Obs = Obs %>% mutate(state_post = state)
    }

    if (use_pos) {
        marker = 'POS'
    } else {
        marker = 'snp_index'
    }

    # correct for baseline bias
    if (fc_correct & !allele_only) {
        Obs = Obs %>% mutate(logFC = logFC - mu)
    }

    D = Obs %>% 
        mutate(logFC = ifelse(logFC > exp_limit | logFC < -exp_limit, NA, logFC)) %>%
        mutate(pBAF = ifelse(DP >= min_depth, pBAF, NA)) %>%
        reshape2::melt(measure.vars = c('logFC', 'pBAF'))

    if (allele_only) {
        D = D %>% filter(variable == 'pBAF')
    }

    p = ggplot(
        D,
        aes(x = get(marker), y = value, color = state_post),
        na.rm=TRUE
    ) +
    geom_point(
        aes(shape = str_detect(state_post, '_2'), alpha = str_detect(state_post, '_2')),
        size = dot_size,
        # stroke = 0.3,
        # alpha = 0.5
    ) +
    scale_alpha_discrete(range = c(0.5, 1)) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 15)) +
    theme_classic() +
    theme(
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
        strip.background = element_blank(),
        axis.text.x = element_blank()
    ) +
    facet_grid(variable ~ CHROM, scale = 'free', space = 'free_x') +
    scale_color_manual(values = cnv_colors) +
    guides(color = guide_legend(title = "", override.aes = aes(size = 3)), fill = FALSE, alpha = FALSE, shape = FALSE) +
    xlab(marker) +
    ylab('')

    if (!legend) {
        p = p + guides(color = FALSE, fill = FALSE, alpha = FALSE, shape = FALSE)
    }

    if (phi_mle) {
        segs = Obs %>% 
            distinct(CHROM, seg, seg_start, seg_start_index, seg_end, seg_end_index, phi_mle) %>%
            mutate(variable = 'logFC') %>%
            filter(log2(phi_mle) < exp_limit)

        if (use_pos) {
            start = 'seg_start'
            end = 'seg_end'
        } else {
            start = 'seg_start_index'
            end = 'seg_end_index'
        }

        p = p + geom_segment(
            inherit.aes = FALSE,
            data = segs,
            aes(x = get(start), xend = get(end), y = log2(phi_mle), yend = log2(phi_mle)),
            color = 'darkred',
            size = 0.5
        ) +
        geom_hline(data = data.frame(variable = 'logFC'), aes(yintercept = 0), color = 'gray30', linetype = 'dashed')
    } else if (!allele_only) {
        p = p + geom_line(
            inherit.aes = FALSE,
            data = Obs %>% mutate(variable = 'logFC') %>% filter(log2(phi_mle_roll) < exp_limit),
            aes(x = get(marker), y = log2(phi_mle_roll), group = '1'),
            color = 'darkred',
            size = 0.35
        ) +
        geom_hline(data = data.frame(variable = 'logFC'), aes(yintercept = 0), color = 'gray30', linetype = 'dashed')
    }

    if (theta_roll) {
        p = p + geom_line(
            inherit.aes = FALSE,
            data = D %>% mutate(variable = 'pBAF'),
            aes(x = snp_index, y = 0.5 - theta_hat_roll, color = paste0(cnv_state_post, '_down')),
            # color = 'black',
            size = 0.35
        ) +
        geom_line(
            inherit.aes = FALSE,
            data = D %>% mutate(variable = 'pBAF'),
            aes(x = snp_index, y = 0.5 + theta_hat_roll, color = paste0(cnv_state_post, '_up')),
            # color = 'gray',
            size = 0.35
        )
    } 
    
    return(p)
}

#' @export
plot_bulks = function(bulk_all, min_depth = 8, fc_correct = TRUE, phi_mle = FALSE, allele_only = FALSE, use_pos = FALSE, ncol = 1, legend = TRUE, title = TRUE) {

    options(warn = -1)
    plot_list = bulk_all %>%
        split(.$sample) %>%
        lapply(
            function(bulk) {

                sample = unique(bulk$sample)
                n_cells = unique(bulk$n_cells)

                p = plot_psbulk(
                        bulk, 
                        min_depth = min_depth, fc_correct = fc_correct,
                        phi_mle = phi_mle, use_pos = use_pos, legend = legend,
                        allele_only = allele_only
                    ) + 
                    theme(
                        title = element_text(size = 8),
                        axis.text.x = element_blank(),
                        axis.title = element_blank()
                    )

                if (title) {
                    p = p + ggtitle(glue('{sample} (n={n_cells})'))
                }
                    

                return(p)
            }
        )
    options(warn = 0)

    panel = wrap_plots(plot_list, ncol = ncol, guides = 'collect')

    return(panel)
}

#' @export
plot_exp = function(gexp_bulk, exp_limit = 3) {
    ggplot(
        gexp_bulk,
        aes(x = gene_index, y = logFC)
    ) +
    theme_classic() +
    geom_point(size = 1, color = 'gray', alpha = 0.8, pch = 16) +
    theme(
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
        strip.background = element_blank()
    ) +
    facet_grid(~CHROM, space = 'free_x', scale = 'free_x') +
    geom_hline(yintercept = 0, color = 'gray30', linetype = 'dashed') +
    geom_line(
        inherit.aes = FALSE,
        aes(x = gene_index, y = log2(phi_hat_roll), group = '1'),
        color = 'darkred',
        size = 0.35
    ) +
    ylim(-exp_limit, exp_limit)
}

#' @export
plot_segs_post = function(segs_consensus) {
    segs_consensus %>% 
        filter(cnv_state != 'neu') %>%
        mutate(seg_label = paste0(seg_cons, '_', cnv_state_post)) %>%
        mutate(seg_label = factor(seg_label, unique(seg_label))) %>%
        reshape2::melt(measure.vars = c('p_loh', 'p_amp', 'p_del', 'p_bamp', 'p_bdel'), value.name = 'p') %>%
        ggplot(
            aes(x = seg_label, y = variable, fill = p, label = round(p, 2))
        ) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
        geom_tile() +
        geom_text(color = 'white')
}

# model diagnostics
#' @export
plot_exp_post = function(exp_post, jitter = TRUE) {
    if (!'annot' %in% colnames(exp_post)) {
        exp_post$annot = '0'
    }
    p = exp_post %>%
        filter(n > 20) %>%
        mutate(seg_label = paste0(seg, '(', cnv_state, ')')) %>%
        mutate(seg_label = factor(seg_label, gtools::mixedsort(unique(seg_label)))) %>%
        ggplot(
            aes(x = seg_label, y = log2(phi_mle), fill = cnv_state, color = p_cnv)
        ) +
        geom_violin(size = 0) +
        geom_hline(yintercept = 0, color = 'green', linetype = 'dashed') +
        geom_hline(yintercept = log2(1.5), color = 'red', linetype = 'dashed') +
        geom_hline(yintercept = -1, color = 'blue', linetype = 'dashed') +
        facet_grid(annot~cnv_state, scale = 'free', space = 'free') +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
        scale_fill_manual(values = cnv_colors)

    if (jitter) {
        p = p + geom_jitter(size = 0.1)
    }

    return(p)
}

#' @export
plot_clones = function(p_matrix, gtree, annot = TRUE, n_sample = 1e4, bar_ratio = 0.1, pal_clone = NULL, pal_annot = NULL) {

    if (is.null(pal_clone)) {
        pal_clone = c('gray', RColorBrewer::brewer.pal(n = 8, 'Set1'))
    }

    if (is.null(pal_annot)) {
        pal_annot = c('gray', RColorBrewer::brewer.pal(n = 8, 'Set1'))
    }

    p_matrix = p_matrix %>% 
        group_by(group) %>%
        filter(cell %in% sample(unique(cell), min(n_sample, length(unique(cell))))) %>%
        ungroup() %>%
        filter(cnv_state != 'neu')

    # ordering cells
    if (annot) {
        p_matrix = p_matrix %>% 
            arrange(group, annot) %>%
            mutate(cell = factor(cell, unique(cell)))
    } else {
        set.seed(0)
        p_matrix = p_matrix %>% 
            group_by(group) %>%
            sample_frac(1) %>%
            ungroup() %>%
            mutate(cell = factor(cell, unique(cell)))
    }
    
    # ordering cnvs
    cnv_order = gtree %>% 
            activate(nodes) %>%
            mutate(rank = dfs_rank(root = node_is_root())) %>%
            data.frame() %>%
            filter(!is.na(site)) %>%
            arrange(-rank) %>%
            pull(site) %>%
            map(function(x){rev(unlist(str_split(x, ',')))}) %>%
            unlist

    p_matrix = p_matrix %>%
        mutate(seg = factor(seg, cnv_order)) %>%
        arrange(seg) %>%
        mutate(seg_label = factor(seg_label, unique(seg_label)))  %>%
        mutate(group = factor(group))
    
    p = ggplot(
            p_matrix,
            aes(x = cell, y = seg_label, fill = logBF)
        ) +
        geom_tile(width=0.1, height=0.9) +
        theme_classic() +
        scale_y_discrete(expand = expansion(0)) +
        scale_x_discrete(expand = expansion(0)) +
        theme(
            panel.spacing = unit(0.1, 'mm'),
            panel.border = element_rect(size = 0.2, color = 'black', fill = NA),
            panel.background = element_rect(fill = 'white'),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            plot.margin = margin(3,0,0,0, unit = 'pt')
            # strip.text = element_text(angle = 0, size = 8, vjust = 0.5)
        ) +
        facet_grid(.~group, scale = 'free', space = 'free') +
        scale_fill_gradient2(low = pal[2], high = pal[1], midpoint = 0, limits = c(-5, 5), oob = scales::oob_squish) +
        xlab('') +
        ylab('') +
        guides(fill = guide_colorbar(barwidth = unit(3, 'mm'), barheight = unit(15, 'mm')))

    p_clones = ggplot(
            p_matrix %>% distinct(cell, group),
            aes(x = cell, y = 'clone', fill = group)
        ) +
        geom_tile(width=1, height=0.9) +
        theme_void() +
        scale_y_discrete(expand = expansion(0)) +
        scale_x_discrete(expand = expansion(0)) +
        theme(
            panel.spacing = unit(0.1, 'mm'),
            panel.border = element_rect(size = 0, color = 'black', fill = NA),
            panel.background = element_rect(fill = 'white'),
            strip.background = element_blank(),
            strip.text = element_text(angle = 0, size = 10, vjust = 0.5),
            axis.text.y = element_text(size = 8)
        ) +
        facet_grid(.~group, scale = 'free', space = 'free') +
        xlab('') +
        ylab('') + 
        scale_fill_manual(values = pal_clone) +
        guides(fill = 'none')

    if (annot) {

        p_annot = ggplot(
                p_matrix %>% distinct(cell, group, annot),
                aes(x = cell, y = '', fill = annot)
            ) +
            geom_tile(width=1, height=0.9) +
            theme_void() +
            scale_y_discrete(expand = expansion(0)) +
            scale_x_discrete(expand = expansion(0)) +
            theme(
                panel.spacing = unit(0.1, 'mm'),
                panel.border = element_rect(size = 0, color = 'black', fill = NA),
                panel.background = element_rect(fill = 'white'),
                strip.background = element_blank(),
                strip.text = element_blank(),
                axis.text.y = element_text(size = 8),
                plot.margin = margin(3,0,0,0, unit = 'pt')
            ) +
            facet_grid(.~group, scale = 'free', space = 'free') +
            xlab('') +
            ylab('') +
            scale_fill_manual(values = pal_annot) +
            guides(fill = guide_legend(keywidth = unit(2, 'mm'), keyheight = unit(2, 'mm'), title = ''))

        return((p_clones / p_annot / p) + plot_layout(height = c(bar_ratio, bar_ratio, 1), guides = 'collect'))
        
    } else {
        return((p_clones / p) + plot_layout(height = c(1,10)))
    }
    
}

#' @export
plot_mut_history = function(G_m, horizontal = TRUE, label = TRUE, pal_clone = NULL) {

    G_m = label_genotype(G_m)

    if (is.null(pal_clone)) {
        getPalette = colorRampPalette(RColorBrewer::brewer.pal(n = 5, 'Spectral'))
        pal_clone = c('gray', getPalette(length(V(G_m))))
    }

    G_df = G_m %>% as_tbl_graph() %>% mutate(clone = factor(clone))

    if (!label) {
        G_df = G_df %>% activate(edges) %>% mutate(to_label = '')
    }

    p = G_df %>% 
        ggraph(layout = 'tree') + 
        geom_edge_link(
            aes(label = str_trunc(to_label, 20, side = 'center')),
            vjust = -1,
            arrow = arrow(length = unit(3, "mm")),
            end_cap = circle(4, 'mm'),
            start_cap = circle(4, 'mm')
        ) + 
        geom_node_point(aes(color = clone), size = 10) +
        geom_node_text(aes(label = clone), size = 6) +
        theme_void() +
        scale_x_continuous(expand = expansion(0.2)) +
        scale_y_continuous(expand = expansion(0.2)) + 
        scale_color_manual(values = pal_clone) +
        guides(color = 'none')

    if (horizontal) {
        p = p + coord_flip() + scale_y_reverse(expand = expansion(0.2))
    }

    return(p)
}

getPalette = colorRampPalette(pal)

#' @export
plot_clone_panel = function(res, label = NULL, cell_annot = NULL, type = 'joint', ratio = 1, tvn = FALSE, tree = TRUE, p_min = 0.5, bar_ratio = 0.1, pal_clone = NULL, pal_annot = NULL) {

    if (is.null(pal_clone)) {
        n_clones = length(unique(res$clone_post$clone_opt))
        pal_clone = getPalette(max(V(res$G_m)-1, 8)) %>% c('gray', .) %>% setNames(1:n_clones)
    } 
    
    if (is.null(pal_annot) & !is.null(cell_annot)) {
        pal_annot = getPalette(length(unique(cell_annot$annot)))
    }
    

    if (type == 'joint') {
        p_matrix = res$joint_post
    } else if (type == 'allele') {
        p_matrix = res$allele_post
    } else {
        p_matrix = res$exp_post
    }

    if (!is.null(cell_annot)) {
        p_matrix = p_matrix %>% left_join(cell_annot, by = 'cell')
        annot = TRUE
    } else {
        annot = FALSE
    }

    if (!'p_opt' %in% colnames(res$clone_post)) {
        res$clone_post = res$clone_post %>% 
            rowwise() %>%
            mutate(p_opt = get(paste0('p_', clone_opt))) %>%
            ungroup()
    }

    if (tvn) {
        res$clone_post = res$clone_post %>%
            mutate(
                clone_opt = ifelse(clone_opt == 1, 'normal', 'tumor'),
                p_opt = ifelse(clone_opt == 'normal', p_1, 1-p_1)
            )
    }

    p_clones = p_matrix %>% 
        filter(seg %in% colnames(res$geno)) %>%
        inner_join(
            res$clone_post %>% filter(p_opt > p_min),
            by = 'cell'
        ) %>%
        mutate(group = clone_opt) %>%
        plot_clones(res$gtree, pal_clone = pal_clone, pal_annot = pal_annot, annot = annot, bar_ratio = bar_ratio)

    plot_title = plot_annotation(title = label, theme = theme(plot.title = element_text(hjust = 0.1)))

    if (tvn | (!tree)) {
        return(p_clones + plot_title)
    }

    p_mut = res$G_m %>% plot_mut_history(pal_clone = pal_clone) 

    (p_mut / p_clones) + plot_layout(heights = c(ratio, 1)) + plot_title
}

#' @export
tree_heatmap = function(joint_post, gtree, ratio = 1, limit = 5, cell_dict = NULL, cnv_order = NULL, label_mut = TRUE, cnv_type = TRUE, branch_width = 0.2, tip = T, tip_length = 0.5, pal_annot = NULL, pal_clone = NULL, layout = 'rect', tvn = FALSE, legend = T) {
    
    if (!'clone' %in% colnames(as.data.frame(activate(gtree, 'nodes')))) {
        gtree = gtree %>% activate(nodes) %>% mutate(clone = as.integer(as.factor(GT)))
    }

    gtree = mark_tumor_lineage(gtree)

    joint_post = joint_post %>% filter(cnv_state != 'neu')

    if (!'seg_label' %in% colnames(joint_post)) {
        joint_post = joint_post %>% mutate(seg_label = paste0(seg, '(', cnv_state, ')')) %>%
            mutate(seg_label = factor(seg_label, unique(seg_label)))
    }

    if (!'logBF' %in% colnames(joint_post)) {
        joint_post = joint_post %>% mutate(logBF = Z_cnv - Z_n)
    }

    if (tvn) {
        clone_dict = gtree %>%
            activate(nodes) %>%
            data.frame %>%
            mutate(compartment = factor(compartment)) %>%
            {setNames(.$compartment, .$name)}
    } else {
        clone_dict = gtree %>%
            activate(nodes) %>%
            data.frame %>%
            mutate(
                GT = ifelse(compartment == 'normal', '', GT),
                GT = factor(GT),
                clone = as.factor(clone)
            ) %>%
            {setNames(.$clone, .$name)}
    }

    getPalette = colorRampPalette(pal)

    if (is.null(pal_annot)) {
        pal_annot = getPalette(length(unique(cell_dict)))
    }

    if (is.null(pal_clone)) {
        pal_clone = getPalette(length(unique(clone_dict)))
    }

    OTU_dict = lapply(levels(clone_dict), function(x) names(clone_dict[clone_dict == x])) %>% setNames(levels(clone_dict))

    mut_nodes = gtree %>% activate(nodes) %>% filter(!is.na(site)) %>% data.frame() %>% select(name, site)

    gtree = gtree %>% activate(edges) %>% mutate(length = ifelse(leaf, pmax(length, tip_length), length))
    
    p_tree = gtree %>% 
        to_phylo() %>%
        groupOTU(
            OTU_dict,
            'clone'
        ) %>%
        ggtree(ladderize = T, size = branch_width, layout = layout) %<+%
        mut_nodes +
        layout_dendrogram() +
        # geom_rootedge(size = branch_width) +
        theme(
            plot.margin = margin(0,0,0,0),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line.y = element_line(size = 0.2),
            axis.ticks.y = element_line(size = 0.2),
            # axis.text.y = element_text(size = 5)
            axis.text.y = element_blank()
        ) +
        guides(color = F) 

    if (tip) {
        p_tree = p_tree + geom_tippoint(aes(color = clone), size=0, stroke = 0.2) +
            scale_color_manual(values = c('gray', pal_clone))
    }

    if (label_mut) {
        p_tree = p_tree + geom_point2(aes(subset = !is.na(site), x = branch), shape = 21, size = 1, fill = 'red') +
            geom_text2(
                aes(x = branch, label = str_trunc(site, 20, side = 'center')),
                size = 2, hjust = 0, vjust = -0.5, nudge_y = 1, color = 'darkred'
            )
    }
    
    if (legend) {
        p_tree = p_tree + 
            guides(color = guide_legend(keywidth = unit(3, 'mm'), override.aes = list(size = 2), keyheight = unit(1, 'mm'), title = NULL))
    }

    cell_order = p_tree$data %>% filter(isTip) %>% arrange(y) %>% pull(label)

    if (is.null(cnv_order)) {
        cnv_order = gtree %>% 
            activate(nodes) %>%
            mutate(rank = dfs_rank(root = node_is_root())) %>%
            data.frame() %>%
            filter(!is.na(site)) %>%
            arrange(-rank) %>%
            pull(site) %>%
            map(function(x){rev(unlist(str_split(x, ',')))}) %>%
            unlist
    }

    p_map = cell_heatmap(joint_post, cnv_order, cell_order, limit, cnv_type = cnv_type)

    p_clones = data.frame(
            cell = names(clone_dict),
            annot = unname(clone_dict)
        ) %>%
        mutate(cell = factor(cell, cell_order)) %>%
        annot_bar(transpose = F) +
        scale_fill_manual(values = c('gray', pal_clone))

    if (!is.null(cell_dict)) {
        p_annot = data.frame(
                cell = names(cell_dict),
                annot = unname(cell_dict)
            ) %>%
            mutate(cell = factor(cell, cell_order)) %>%
            annot_bar(transpose = F)
            # scale_fill_manual(values = pal_annot)

        panel = (p_tree / p_clones / p_annot / p_map) + plot_layout(heights = c(ratio,0.06,0.06,1), guides = 'collect')
    } else {
        panel = (p_tree / p_clones / p_map) + plot_layout(heights = c(ratio,0.1,1), guides = 'collect')
    }

    return(panel)
}

#' @export
plot_sc_joint = function(
        gtree, joint_post, segs_consensus, 
        cell_dict = NULL, size = 0.02, branch_width = 0.2, tip_length = 0.2, logBF_min = 1, logBF_max = 5, clone_bar = FALSE, clone_legend = TRUE, pal_clone = NULL
    ) {

    if (!'clone' %in% colnames(as.data.frame(activate(gtree, 'nodes')))) {
        gtree = gtree %>% activate(nodes) %>% mutate(clone = as.integer(as.factor(GT)))
    }
          
    gtree = mark_tumor_lineage(gtree)

    gtree = gtree %>% activate(edges) %>% mutate(length = ifelse(leaf, pmax(length, tip_length), length))

    # tree visualization
    p_tree = gtree %>% 
            to_phylo() %>%
            ggtree(ladderize = T, size = branch_width) +
            # geom_rootedge(size = branch_width) +
            theme(
                plot.margin = margin(0,1,0,0, unit = 'mm'),
                axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank(),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                # axis.text.y = element_text(size = 5)
                axis.text.y = element_blank()
            ) +
            guides(color = F)

    cell_order = p_tree$data %>% filter(isTip) %>% arrange(y) %>% pull(label)

    # cell heatmap
    D = joint_post %>% 
            inner_join(
                segs_consensus %>% select(seg = seg_cons, CHROM, seg_start, seg_end),
                by = c('seg', 'CHROM')
            ) %>%
            mutate(cell = factor(cell, cell_order)) %>%
            mutate(cell_index = as.integer(droplevels(cell))) 

    tumor_cells = gtree %>% 
        activate(nodes) %>% filter(leaf) %>%
        as.data.frame() %>% 
        filter(compartment == 'tumor') %>%
        pull(name)

    first_tumor_index = which(cell_order %in% tumor_cells)[1]

    p_segs = ggplot(
            D %>% mutate(logBF = pmax(pmin(logBF, logBF_max), logBF_min))
        ) +
        theme_classic() +
        geom_segment(
            aes(x = seg_start, xend = seg_end, y = cell_index, yend = cell_index, color = cnv_state, alpha = logBF),
            size = size
        ) +
        geom_segment(
            inherit.aes = F,
            aes(x = seg_start, xend = seg_end, y = 1, yend = 1),
            data = segs_consensus, size = 0, color = 'white', alpha = 0
        ) +
        geom_hline(yintercept = first_tumor_index, color = 'royalblue', size = 0.5, linetype = 'dashed') +
        theme(
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
            strip.background = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = margin(0,0,5,0, unit = 'mm'),
            axis.line = element_blank()
        ) +
        scale_x_continuous(expand = expansion(0)) +
        scale_y_continuous(expand = expansion(0)) +
        facet_grid(.~CHROM, space = 'free', scale = 'free') +
        scale_alpha_continuous(range = c(0,1)) +
        guides(
            alpha = 'none',
            color = guide_legend(override.aes = c('size' = 1))
        ) +
        scale_color_manual(
            values = c('amp' = 'darkred', 'del' = 'darkblue', 'bamp' = cnv_colors[['bamp']], 'loh' = 'darkgreen', 'bdel' = 'blue', 'neu' = 'white')
        )

    # clone annotation
    clone_dict = gtree %>%
        activate(nodes) %>%
        data.frame %>%
        mutate(
            GT = ifelse(compartment == 'normal', '', GT),
            GT = factor(GT),
            clone = as.factor(clone)
        ) %>%
        {setNames(.$clone, .$name)}

    if (is.null(pal_clone)) {
        getPalette = colorRampPalette(RColorBrewer::brewer.pal(n = 5, 'Spectral'))
        pal_clone = c('gray', getPalette(length(unique(clone_dict))))
    }

    p_clone = data.frame(
            cell = names(clone_dict),
            annot = unname(clone_dict)
        ) %>%
        mutate(cell = factor(cell, cell_order)) %>%
        annot_bar(transpose = T, legend = clone_legend) +
        scale_fill_manual(values =pal_clone)

    # external annotation
    if (!is.null(cell_dict)) {
        
        p_annot = data.frame(
                cell = names(cell_dict),
                annot = unname(cell_dict)
            ) %>%
            filter(cell %in% joint_post$cell) %>%
            mutate(cell = factor(cell, cell_order)) %>%
            annot_bar(transpose = T)

        if (clone_bar) {
            (p_tree | p_clone | p_annot | p_segs) + plot_layout(widths = c(1, 0.25, 0.25, 15), guides = 'collect')
        } else {
            (p_tree | p_annot | p_segs) + plot_layout(widths = c(1, 0.25, 15), guides = 'collect')
        }
    } else if (clone_bar) {
        (p_tree | p_clone | p_segs) + plot_layout(widths = c(1, 0.25, 15), guides = 'collect')
    } else {
        (p_tree | p_segs) + plot_layout(widths = c(1, 15), guides = 'collect')
    }
}

annot_bar = function(D, transpose = FALSE, legend = TRUE) {
    p = ggplot(
        D,
        aes(x = cell, y = '', fill = annot)
    ) +
    geom_tile(width=1, height=0.9) +
    theme_void() +
    scale_y_discrete(expand = expansion(0)) +
    scale_x_discrete(expand = expansion(0)) +
    theme(
        panel.spacing = unit(0.1, 'mm'),
        panel.border = element_rect(size = 0, color = 'black', fill = NA),
        panel.background = element_rect(fill = 'white'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        # axis.text = element_text(size = 8),
        axis.text = element_blank(),
        plot.margin = margin(0.5,0,0.5,0, unit = 'mm')
    ) 

    if (transpose) {
        p = p + coord_flip() +
            theme(plot.margin = margin(0,0.5,0,0.5, unit = 'mm'))
    }

    if (legend) {
        p = p + guides(fill = guide_legend(keywidth = unit(3, 'mm'), keyheight = unit(1, 'mm'), title = NULL))
    } else {
        p = p + guides(fill = 'none')
    }

    return(p)
}

#' @export
cell_heatmap = function(geno, cnv_order = NULL, cell_order = NULL, limit = 5, cnv_type = TRUE) {

    # geno = geno %>% mutate(logBF = Z_cnv - Z_n)

    if (is.null(cnv_order)) {
        cnv_order = unique(geno$seg)
    }

    if (is.null(cell_order)) {
        cell_order = unique(geno$cell)
    }

    geno = geno %>% 
        filter(cell %in% cell_order) %>%
        filter(cnv_state != 'neu') %>%
        mutate(seg = factor(seg, cnv_order)) %>%
        arrange(seg) %>%
        mutate(seg_label = factor(seg_label, unique(seg_label))) %>%
        mutate(cell = factor(cell, cell_order))

    pal = RColorBrewer::brewer.pal(n = 8, 'Set1')

    if (!cnv_type) {
        geno = geno %>% mutate(seg_label = seg)
    }

    p_map = ggplot(
            geno,
            aes(x = cell, y = seg_label, fill = logBF)
        ) +
        geom_tile(width=0.4, height=0.9) +
        theme_classic() +
        scale_y_discrete(expand = expansion(0)) +
        scale_x_discrete(expand = expansion(add = 0.5)) +
        theme(
            panel.spacing = unit(0.1, 'mm'),
            # panel.border = element_rect(size = 0.5, color = 'black', fill = NA),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.x = element_blank(),
            panel.border = element_blank(),
            panel.background = element_rect(fill = 'white'),
            strip.background = element_blank(),
            axis.text.x = element_blank(),
            strip.text = element_text(angle = 90, size = 8, vjust = 0.5),
            plot.margin = margin(0,0,0,0, unit = 'mm')
        ) +
        scale_fill_gradient2(low = pal[2], high = pal[1], midpoint = 0, limits = c(-limit, limit), oob = scales::oob_squish) +
        # xlab('') +
        theme(plot.title = element_blank()) +
        ylab('') +
        guides(fill = guide_colorbar(barwidth = unit(3, 'mm'), barheight = unit(15, 'mm')))

    return(p_map)
}

#' @export
plot_sc_roll = function(gexp.norm.long, hc, k, lim = 0.8, n_sample = 50) {

    cells = unique(gexp.norm.long$cell)
    
    cell_sample = sample(cells, min(n_sample, length(cells)), replace = FALSE)
        
    p_tree = ggtree(hc, size = 0.2)

    cell_order = p_tree$data %>% filter(isTip) %>% arrange(y) %>% pull(label)
    
    p_heatmap = gexp.norm.long %>%
        filter(cell %in% cell_sample) %>%
        mutate(cell = factor(cell, cell_order)) %>%
        mutate(cluster = cutree(hc, k = k)[as.character(cell)]) %>%
        arrange(cell) %>%
        mutate(cluster = factor(cluster, rev(unique(cluster)))) %>%
        ggplot(
            aes(x = gene_index, y = cell, fill = exp_rollmean)
        ) +
        geom_tile() +
        scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = 0, limits = c(-lim,lim), oob = scales::squish) +
        theme_void() +
        scale_x_discrete(expand = expansion(0)) +
        scale_y_discrete(expand = expansion(0)) +
        theme(
            axis.text.y = element_blank(),
            legend.position = 'top',
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'wheat4', fill = NA)
        ) +
        facet_grid(cluster~CHROM, scale = 'free', space = 'free')
    
    (p_tree | p_heatmap) + plot_layout(widths = c(1,10))
    # p_heatmap
    
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

#' @export
plot_sc_exp = function(exp_post, segs_consensus, size = 0.05, censor = 0) {
    
    # cell_order = exp_post %>% 
    #     filter(!cnv_state %in% c('neu', 'loh')) %>%
    #     left_join(cell_annot, by = 'cell') %>%
    #     group_by(cell_group) %>%
    #     do(
    #         reshape2::dcast(., cell ~ seg, value.var = 'phi_mle') %>%
    #         tibble::column_to_rownames('cell') %>%
    #         dist() %>%
    #         hclust %>%
    #         {.$labels[.$order]} %>%
    #         as.data.frame()
    #     ) %>%
    #     set_names(c('cell_group', 'cell'))

    exp_post = exp_post %>% filter(n > 15)

    exp_post = exp_post %>% 
            inner_join(
                segs_consensus %>% select(seg = seg_cons, CHROM, seg_start, seg_end),
                by = 'seg'
            ) %>%
            mutate(phi_mle = ifelse(phi_mle > 1-censor & phi_mle < 1+censor, 1, phi_mle))
            # mutate(cell = factor(cell, cell_order$cell))

    ggplot(
        exp_post,
        aes(x = seg_start, xend = seg_end, y = cell, yend = cell, color = phi_mle)
    ) +
    theme_classic() +
    geom_segment(size = size) +
    theme(
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
        strip.background = element_blank(),
        axis.text.y = element_blank()
    ) +
    scale_x_continuous(expand = expansion(0)) +
    facet_grid(group~CHROM, space = 'free', scale = 'free') +
    scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 1, limits = c(0.5, 2), oob = scales::oob_squish)
}

#' @export
plot_sc_allele = function(df_allele, bulk_subtrees, clone_post) {
    
    snp_seg = bulk_subtrees %>%
        filter(!is.na(pAD)) %>%
        mutate(haplo = case_when(
            str_detect(state, 'up') ~ 'major',
            str_detect(state, 'down') ~ 'minor',
            T ~ ifelse(pBAF > 0.5, 'major', 'minor')
        )) %>%
        select(snp_id, snp_index, sample, seg, haplo) %>%
        inner_join(
            segs_consensus,
            by = c('sample', 'seg')
        )

    snp_neu_haplo = bulk_subtrees %>% filter(sample == 1) %>%
        mutate(haplo = ifelse(pBAF > 0.5, 'major', 'minor')) %>% 
        filter(!is.na(haplo)) %>%
        {setNames(.$haplo, .$snp_id)}
    
    
    pal = RColorBrewer::brewer.pal(n = 8, 'Set1')

    p = df_allele %>% 
        left_join(
            snp_seg %>% select(snp_id, haplo, cnv_state),
            by = 'snp_id'
        ) %>% 
        mutate(haplo = ifelse(is.na(haplo), snp_neu_haplo[snp_id], haplo)) %>%
        mutate(cnv_state = ifelse(is.na(cnv_state), 'neu', cnv_state)) %>%
        mutate(pBAF = ifelse(GT == '1|0', AR, 1-AR)) %>%
        filter(!is.na(haplo)) %>%
        mutate(MAF = ifelse(haplo == 'major', pBAF, 1-pBAF)) %>%
        left_join(clone_post, by = 'cell') %>%
        arrange(clone) %>%
        arrange(CHROM, POS) %>%
        mutate(snp_index = as.integer(factor(snp_id, unique(snp_id)))) %>%
        ggplot(
            aes(x = snp_index, y = cell, color = MAF)
        ) +
        theme_classic() +
        geom_point(alpha = 0.5, pch = 16, size = 1) +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'white', fill = NA),
        ) +
        facet_grid(clone~CHROM, space = 'free', scale = 'free') +
        scale_x_discrete(expand = expansion(0)) +
        scale_color_gradient(low = pal[1], high = pal[2])
    
    return(p)
}

#' @export
clone_vs_annot = function(clone_post, cell_annot) {
    clone_post %>% 
    rename(clone = clone_opt) %>%
    filter(!is.na(clone)) %>%
    left_join(cell_annot, by = 'cell') %>%
    count(clone, cell_type) %>%
    arrange(cell_type) %>%
    mutate(clone = factor(clone, rev(unique(clone)))) %>%
    group_by(cell_type) %>%
    mutate(frac = n/sum(n)) %>%
    ggplot(
        aes(x = cell_type, y = clone, fill = frac, label = n)
    ) +
    theme_classic() +
    geom_tile() +
    geom_text() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    scale_fill_gradient(low = 'white', high = 'red') +
    xlab('')
}

#' @export
plot_markers = function(sample, count_mat, cell_annot, markers, clone_post, pal_annot = NULL) {

    if (is.null(pal_annot)) {
        pal_annot = getPalette(length(unique(cell_annot$annot)))
    }
    
    D = as.matrix(count_mat[,markers$gene]) %>%
        scale %>%
        reshape2::melt() %>%
        set_colnames(c('cell', 'gene', 'exp')) %>%
        inner_join(
            cell_annot, by = 'cell'
        ) %>%
        mutate(exp = ifelse(is.na(exp), 0, exp)) %>%
        inner_join(
            clone_post, by = 'cell'
        ) %>%
        left_join(markers, by = 'gene') %>%
        arrange(p_1) %>%
        mutate(cell = factor(cell, unique(cell)))

    p_markers = ggplot(
            D,
            aes(x = cell, y = gene, fill = exp)
        ) +
        geom_tile() +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 7),
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.2, fill = NA),
            strip.background.x = element_blank(),
            strip.text.x = element_blank(),
            strip.background.y = element_rect(size = 0, fill = NA),
            strip.text.y.left = element_text(size = 6, angle = 0),
        ) +
        ylab('marker') +
        facet_grid(marker_type ~ cell_group, space = 'free_y', scale = 'free', switch="y") +
        scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', limits = c(-1.5,1.5), oob = scales::oob_squish)

    p_annot = ggplot(
            D,
            aes(x = cell, y = 'annot', fill = annot)
        ) +
        geom_tile() +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 7),
            axis.ticks.x = element_blank(),
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.2, fill = NA),
            strip.background = element_rect(size = 0, fill = NA),
            strip.text = element_text(size = 6),
            axis.title.x = element_blank(),
            strip.text.x = element_blank()
        ) +
        ylab('') +
        facet_grid(~ cell_group, space = 'free_y', scale = 'free', switch="y") +
        scale_fill_manual(values = pal_annot) +
        guides(fill = guide_legend())

    p_cnv = ggplot(
            D %>% mutate(p_cnv = 1-p_1),
            aes(x = cell, y = 'cnv', fill = p_cnv)
        ) +
        geom_tile() +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 7),
            axis.ticks.x = element_blank(),
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.2, fill = NA),
            strip.background = element_rect(size = 0, fill = NA),
            strip.text = element_text(size = 6),
            axis.title.x = element_blank()
        ) +
        ylab('') +
        facet_grid(~cell_group, space = 'free_y', scale = 'free', switch="y") +
        scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0.5, limits = c(0,1), oob = scales::oob_squish) + 
        ggtitle(sample) +
        guides(fill = 'none')

    p_cnv/p_annot/p_markers + plot_layout(heights = c(0.5,0.5,10), guides = 'collect')
    
}
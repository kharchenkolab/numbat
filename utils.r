source('~/Armadillo/hmm.r')

########################### Data processing ############################

process_exp_fc = function(count_mat_obs, lambdas_ref, gtf_transcript, verbose = TRUE) {
    
    genes_annotated = gtf_transcript$gene %>% 
        intersect(rownames(count_mat)) %>%
        intersect(names(lambdas_ref))

    count_mat_obs = count_mat_obs[genes_annotated,]
    lambdas_ref = lambdas_ref[genes_annotated]

    depths = colSums(count_mat_obs)
    depth_obs = sum(depths)
    lambdas_obs = rowSums(count_mat_obs)/depth_obs

    # filter for mutually expressed genes
    min_both = 2
    
    mut_expressed = ((lambdas_ref * 1e6 > min_both & lambdas_obs * 1e6 > min_both) |
        (lambdas_ref > mean(lambdas_ref[lambdas_ref != 0])) |
        (lambdas_obs > mean(lambdas_obs[lambdas_obs != 0]))) &
        (lambdas_ref > 0 & lambdas_obs > 0)
    
    count_mat_obs = count_mat_obs[mut_expressed,]
    lambdas_ref = lambdas_ref[mut_expressed]

    if(verbose){display(nrow(count_mat_obs))}

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


preprocess_data = function(
    sample,
    vcf_pu,
    vcf_phased,
    AD,
    DP,
    barcodes,
    cell_annot,
    gtf_transcript
) {

    # transcript bed
    transcript_regions = gtf_transcript %>%
        pull(region) %>%
        unique %>% bedr::bedr.sort.region(verbose = F)

    # pileup VCF
    vcf_pu = vcf_pu %>%
        mutate(INFO = str_remove_all(INFO, '[:alpha:]|=')) %>%
        tidyr::separate(col = 'INFO', into = c('AD', 'DP', 'OTH'), sep = ';') %>%
        mutate_at(c('AD', 'DP', 'OTH'), as.integer) %>%
        mutate(snp_id = paste(CHROM, POS, REF, ALT, sep = '_'))

    # pileup count matrices
    DP = as.data.frame(summary(DP)) %>%
        mutate(
            cell = cell_barcodes[j],
            snp_id = vcf_pu$snp_id[i]
        ) %>%
        select(-i, -j) %>%
        rename(DP = x) %>%
        select(cell, snp_id, DP)

    AD = as.data.frame(summary(AD)) %>%
        mutate(
            cell = cell_barcodes[j],
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

    df = df %>% filter(DP_all > 1 & OTH_all == 0)

    # vcf has duplicated records sometimes
    df = df %>% distinct() 

    df = df %>% mutate(
        snp_index = as.integer(factor(snp_id, unique(snp_id))),
        cell_index = as.integer(factor(cell, sample(unique(cell))))
    )
    
    # phased VCF
    vcf_phased = vcf_phased %>% mutate(INFO = str_remove_all(INFO, '[:alpha:]|=')) %>%
        tidyr::separate(col = 'INFO', into = c('AD', 'DP', 'OTH'), sep = ';') %>%
        mutate_at(c('AD', 'DP', 'OTH'), as.integer)

    vcf_phased = vcf_phased %>% mutate(snp_id = paste(CHROM, POS, REF, ALT, sep = '_')) %>%
        mutate(GT = get(sample))

    vcf_phased = vcf_phased %>% mutate(region = paste0('chr', CHROM, ':', POS, '-', format(POS+1, scientific = F, trim = T)))

    # intersect with gene model
    vcf_phased_regions = vcf_phased %>%
        pull(region) %>%
        bedr::bedr.sort.region(verbose = F)

    overlap_transcript = bedr::bedr(
        input = list(a = vcf_phased_regions, b = transcript_regions), 
        method = "intersect", 
        params = "-loj -sorted",
        verbose = F
    ) %>% filter(V4 != '.')

    # annotate SNP by gene
    vcf_phased = vcf_phased %>% left_join(
        overlap_transcript %>%
        rename(CHROM = V4, gene_start = V5, gene_end = V6) %>%
        mutate(
            gene_start = as.integer(gene_start), 
            gene_end = as.integer(gene_end),
            CHROM = as.integer(str_remove(CHROM, 'chr')),
        ) %>%
        left_join(
            gtf_transcript %>% select(CHROM, gene_start, gene_end, gene),
            by = c('CHROM', 'gene_start', 'gene_end')
        ) %>%
        arrange(index, gene) %>%
        distinct(index, `.keep_all` = T),
        by = c('region' = 'index')
    )

    # add annotation to cell counts and pseudobulk
    df = df %>% left_join(vcf_phased %>% select(snp_id, gene, gene_start, gene_end, GT), by = 'snp_id') 
    df = df %>% mutate(CHROM = factor(CHROM, unique(CHROM)))

    df = df %>% 
        select(-any_of(c('cell_type', 'group'))) %>%
        left_join(
            cell_annot,
            by = "cell"
        ) 
    
    df_obs = df %>% filter(group == 'obs')
    
    df_ref = df %>% filter(group == 'ref')

    return(list('df_obs' = df_obs, 'df_ref' = df_ref))
}

get_bulk = function(count_mat, lambdas_ref, df, gtf_transcript, min_depth = 2) {
    
    gexp_bulk = process_exp_fc(
        count_mat,
        lambdas_ref,
        gtf_transcript,
        verbose = FALSE
        )$bulk %>%
        filter(logFC < 8 & logFC > -8)
            
    combine_bulk(
        df = df,
        gexp_bulk = gexp_bulk,
        min_depth = min_depth
    )
}

combine_bulk = function(df, gexp_bulk, min_depth = 2) {
    
    pseudobulk = df %>%
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
        group_by(CHROM) %>%
        filter(n() > 1) %>%
        mutate(inter_snp_dist = c(NA, POS[2:length(POS)] - POS[1:(length(POS)-1)])) %>%
        mutate(p_s = switch_prob(inter_snp_dist)) %>%
        ungroup()

    Obs = pseudobulk %>% 
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


########################### Analysis ############################

analyze_bulk_gpois = function(Obs, p_0, gamma = 16, bal_cnv = FALSE, prior = NULL, exp_only = FALSE, allele_only = FALSE, update_t = TRUE, verbose = TRUE) {
    
    # doesn't work with 0s in the ref
    Obs = Obs %>% filter(lambda_ref != 0 | is.na(gene)) 
        
    # define diploid regions
    Obs = Obs %>% group_by(CHROM) %>%
        mutate(state = 
            run_hmm_mv_inhom_gpois2(
                pAD = pAD,
                DP = DP, 
                p_s = p_s,
                Y_obs = Y_obs, 
                lambda_ref = lambda_ref, 
                d_total = 0,
                p_0 = p_0,
                allele_only = TRUE
            )
        ) %>% ungroup()
    
    diploid_regions = find_diploid_regions(Obs)
    
    Obs = Obs %>% mutate(cnv_state = ifelse(seg %in% diploid_regions, 'neu', 'aberrant'))
    
    if (verbose) {
        display(glue('diploid regions: {paste0(diploid_regions, collapse = ",")}'))
    }
    
    converge = FALSE
    
    i = 0
    max_iter = 10
    
    while (!(converge | i > max_iter)) {
        
        i = i + 1
        
        # update parameters
        fit = Obs %>%
            # exclude outliers
            filter(logFC < 8 & logFC > -8) %>%
            filter(!is.na(Y_obs) & cnv_state == 'neu') %>%
            {fit_gpois(.$Y_obs, .$lambda_ref, unique(.$d_obs))}
            
        alpha_hat = fit@coef[1]
        beta_hat = fit@coef[2]
        
        if (i > 1 & update_t) {
            p_0 = 1 - max(sum(Obs$boundary), 1)/nrow(Obs)
        }
        
        if (verbose) {
            display(glue('iteration {i}, alpha {signif(alpha_hat, 3)}, beta {signif(beta_hat, 3)}, t {signif(1-p_0, 3)}'))
        }
        
        state_old = Obs$cnv_state
        
        if (bal_cnv) {
            Obs = Obs %>% 
                group_by(CHROM) %>%
                mutate(state = 
                    run_hmm_mv_inhom_gpois2(
                        pAD = pAD,
                        DP = DP, 
                        p_s = p_s,
                        Y_obs = Y_obs, 
                        lambda_ref = lambda_ref, 
                        d_total = na.omit(unique(d_obs)),
                        phi_neu = 1,
                        phi_del = 2^(-0.25),
                        phi_amp = 2^(0.25),
                        phi_bamp = 2^(0.5),
                        phi_bdel = 2^(-0.5),
                        alpha = alpha_hat,
                        beta = beta_hat,
                        p_0 = p_0,
                        gamma = gamma,
                        prior = prior,
                        exp_only = exp_only,
                        allele_only = allele_only
                    )
                ) %>% ungroup()
            
        } else {
            Obs = Obs %>% 
                group_by(CHROM) %>%
                mutate(state = 
                    run_hmm_mv_inhom_gpois(
                        pAD = pAD,
                        DP = DP, 
                        p_s = p_s,
                        Y_obs = Y_obs, 
                        lambda_ref = lambda_ref, 
                        d_total = na.omit(unique(d_obs)),
                        phi_neu = 1,
                        phi_del = 2^(-0.25),
                        phi_amp = 2^(0.25),
                        alpha = alpha_hat,
                        beta = beta_hat,
                        p_0 = p_0,
                        gamma = gamma,
                        prior = prior,
                        exp_only = exp_only,
                        allele_only = allele_only
                    )
                ) %>% ungroup()
            
        }
        
        Obs = Obs %>% 
            mutate(cnv_state = str_remove(state, '_down|_up')) %>%
            group_by(CHROM) %>%
            arrange(CHROM, snp_index) %>%
            mutate(boundary = c(0, cnv_state[2:length(cnv_state)] != cnv_state[1:(length(cnv_state)-1)])) %>%
            group_by(CHROM) %>%
            mutate(seg = paste0(CHROM, '_', cumsum(boundary))) %>%
            arrange(CHROM) %>%
            mutate(seg = factor(seg, unique(seg))) %>%
            ungroup()
        
        converge = all(state_old == Obs$cnv_state)
        
    }
    
    if (verbose) {
        display('Finishing')
    }
        
    Obs = Obs %>% 
        select(-any_of('phi_mle_roll')) %>%
        left_join(
        Obs %>% 
            group_by(CHROM) %>%
            filter(!is.na(Y_obs)) %>%
            mutate(
                phi_mle_roll = phi_mle_roll(
                    Y_obs, lambda_ref, alpha_hat, beta_hat, d_obs, h = 50)
            ) %>%
            select(phi_mle_roll, CHROM, gene),
        by = c('CHROM', 'gene')
    )
    
    Obs = Obs %>% mutate(alpha = alpha_hat, beta = beta_hat)

    
    return(Obs)
}


find_diploid_regions = function(bulk, debug = F) {
    
    bulk_balanced = bulk %>% 
        filter(state == 'neu') %>% 
        filter(!is.na(lnFC)) %>%
        mutate(FC = exp(lnFC)) 

    segs = bulk_balanced %>% 
        group_by(seg) %>%
        summarise(
            n_genes = sum(!is.na(Y_obs))
        ) %>%
        ungroup() %>%
        filter(n_genes > 50)
    
    # pairwise K-S test
    tests = data.frame()
    
    options(warn = -1)
    for (i in segs$seg) {
        for (j in segs$seg) {
            tests = rbind(
                tests,
                data.frame(
                    i = i,
                    j = j,
                    p = ks.test(x = bulk_balanced$FC[bulk_balanced$seg == i], y = bulk_balanced$FC[bulk_balanced$seg == j])$p.value
                )
            )
        }       
    }
    options(warn = 0)
    
    # build graph
    V = segs

    E = tests %>% 
        filter(p > 0.01) %>%
        filter(i > j)

    G = igraph::graph_from_data_frame(d=E, vertices=V, directed=F)
    
    # find confident diploid clique
    cliques = igraph::largest_cliques(G)

    FC = sapply(cliques, function(c) {
        bulk_balanced %>% filter(seg %in% names(c)) %>% {phi_hat_seg(.$Y_obs, .$lambda_ref, unique(.$d_obs))}
    })

    diploid_clique = which.min(FC)
    
    diploid_segs = names(cliques[[diploid_clique]])
    
    if (debug) {
        options(repr.plot.width = 4, repr.plot.height = 4, repr.plot.res = 200)
        par(mar=c(0,0,0,0)+.1)
        plot(G, vertex.size = 30, vertex.label.cex = 1)
        display(FC)
    }
    
    return(diploid_segs)
    
}


# get all internal nodes 
get_internal_nodes = function(den, node, cuts, min_size = 400) {

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
        paste0(node, '.', 1),
        min_size
    )
    
    membership_r = get_internal_nodes(
        sub_dens[[2]],
        paste0(node, '.', 2),
        min_size
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


phi_mle = function(Y_obs, lambda_ref, d, alpha, beta, lower = 0.2, upper = 10) {
    
    res = stats4::mle(
        minuslogl = function(phi) {
            res = -sum(dgpois(Y_obs, shape = alpha, rate = beta/(phi * d * lambda_ref), log = TRUE))
            return(res)
        },
        start = 1,
        lower = lower,
        upper = upper
    )
    
    return(res)
}

phi_mle_roll = function(Y_obs, lambda_ref, alpha, beta, d_obs, h) {
    n = length(Y_obs)
    sapply(
        1:n,
        function(c) {
            slice = max(c - h - 1, 1):min(c + h, n)
            phi_mle(Y_obs[slice], lambda_ref[slice], unique(d_obs), alpha, beta)@coef
        }
    )
}

l_gpois = function(Y_obs, lambda_ref, d, alpha, beta, phi = 1) {
    sum(dgpois(Y_obs, shape = alpha, rate = beta/(phi * d * lambda_ref), log = TRUE))
}

fit_gpois = function(Y_obs, lambda_ref, d) {
    
    res = stats4::mle(
        minuslogl = function(alpha, beta) {
            -l_gpois(Y_obs, lambda_ref, d, alpha, beta)
        },
        start = c(1, 1),
        lower = c(0.01, 0.01),
        upper = c(5, 5)
    )
    
    return(res)
}

LLR_exp = function(Y_obs, lambda_ref, d, alpha, beta) {
    
    if (length(Y_obs) == 0) {
        return(0)
    }
    
    phi_mle = phi_mle(Y_obs, lambda_ref, d, alpha, beta)@coef
    l_1 = l_gpois(Y_obs, lambda_ref, d, alpha, beta, phi = phi_mle)
    l_0 = l_gpois(Y_obs, lambda_ref, d, alpha, beta, phi = 1)
    return(l_1 - l_0)
}


walk_tree = function(den, cut_prefix, cuts, segs_tested, gexp.norm.long, df, p_cutoff = 0.01, min_depth = 10) {
    
    # split into two clusters
    membership = dendextend::cutree(den, k = 2) %>%
        data.frame() %>%
        tibble::rownames_to_column('cell') %>%
        setNames(c('cell', 'cut')) %>%
        mutate(cut_prefix = cut_prefix)
        
    IRdisplay::display(cut_prefix)
    
    min_size = membership %>% count(cut) %>% pull(n) %>% min()
    
    accept = FALSE
    
    if (min_size > 400) {
        
        IRdisplay::display('Running HMMs..')
        Obs_cut = analyze_cut(df, gexp.norm.long, membership, min_depth)
        
        IRdisplay::display('Comparing clusters..')
        seg_diff = compare_cluster(Obs_cut) 
        
        if (nrow(seg_diff) > 0) {
            # accept the split if any significantly different segments
            accept = any(seg_diff$p_diff <= p_cutoff)
            
            segs_tested = rbind(
                segs_tested, 
                seg_diff %>% mutate(cut_prefix = cut_prefix))
        }
    } 
    
    if (!accept) {
        # this is a fringe node
        cuts = rbind(cuts, membership)
        return(list('cuts' = cuts, 'segs_tested' = segs_tested))
    } else {
        
        # keep splitting
        sub_dens = dendextend::get_subdendrograms(den, k = 2)
        
        res_1 = walk_tree(
            sub_dens[[1]],
            paste0(cut_prefix, '.', 1),
            cuts,
            segs_tested,
            gexp.norm.long,
            df,
            p_cutoff,
            min_depth
        )
        
        res_2 = walk_tree(
            sub_dens[[2]],
            paste0(cut_prefix, '.', 2),
            cuts,
            segs_tested, 
            gexp.norm.long,
            df,
            p_cutoff,
            min_depth
        )
        
        res = list(
            'cuts' = rbind(res_1$cuts, res_2$cuts), 
            'segs_tested' = rbind(res_1$segs_tested, res_2$segs_tested))
        
        return(res)
    }
}
    

calc_cluster_tree = function(count_mat, cell_annot) {
    
    exp_mat = t(count_mat/colSums(count_mat))
    
    exp_mat = exp_mat[,colSums(exp_mat) > 0]
        
    cell_dict = cell_annot %>% filter(group == 'obs') %>% {setNames(.$cell_type, .$cell)}

    exp_mat = exp_mat %>% data.frame() %>%
        tibble::rownames_to_column('cell') %>%
        mutate(cell_type = cell_dict[cell]) %>%
        filter(!is.na(cell_type)) %>%
        select(-cell)

    exp_mat_clust = aggregate(select(exp_mat, -cell_type), list('cell_type' = exp_mat$cell_type), mean) %>%
        tibble::column_to_rownames('cell_type')
    
    exp_mat_clust = log10(exp_mat_clust * 1e6 + 1)
    
    tree = hclust(as.dist(1-cor(t(exp_mat_clust))), method = "ward.D2")
    
    return(tree)
    
}

calc_alelle_lik = function(pAD, DP, p_s, theta) {
    
    # states
    states = c("theta_up", "theta_down")
    
    # transition matrices
    calc_trans_mat = function(p_s, p_0) {
        matrix(
            c(1 - p_s, p_s,
              p_s, 1 - p_s),
            ncol = 2,
            byrow = TRUE
        )
    }
    
    As = lapply(
        p_s,
        function(p_s) {calc_trans_mat(p_s, p_0)}
    )
    
    # intitial probabilities
    prior = c(0.5, 0.5)
    
    alpha_up = (0.5 + theta) * 16
    beta_up = (0.5 - theta) * 16
    alpha_down = beta_up
    beta_down = alpha_up
        
    hmm = HiddenMarkov::dthmm(
        x = pAD, 
        Pi = As, 
        delta = prior, 
        distn = "bbinom",
        pm = list(alpha=c(alpha_up, alpha_down), beta=c(beta_up, beta_down)),
        pn = list(size = DP),
        discrete = TRUE)

    class(hmm) = 'dthmm.inhom'
        
    solution = HiddenMarkov::Viterbi(hmm)
        
    return(solution$max.loglik)
}

calc_theta_mle = function(pAD, DP, p_s, lower = 0, upper = 0.49) {
    res = optim(
        0.25, 
        function(theta) {-calc_alelle_lik(pAD, DP, p_s, theta)},
        method = 'L-BFGS-B',
        lower = lower, upper = upper)
    return(list('theta_mle' = res$par, 'l' = -res$value))
}

LLR_allele = function(pAD, DP, p_s) {
    
    if (length(pAD) <= 1) {
        return(0)
    }
    
    l_1 = calc_theta_mle(pAD, DP, p_s)$l
    l_0 = calc_alelle_lik(pAD, DP, p_s, 0)
    return(l_1 - l_0)
}

# multi-state model
get_exp_post = function(exp_sc, pi) {
    
    depth_obs = sum(exp_sc$Y_obs)
        
    fit = exp_sc %>% filter(cnv_state == 'neu') %>%
        {fit_gpois(.$Y_obs, .$lambda_ref, depth_obs)}
    
    alpha = fit@coef[1]
    beta = fit@coef[2]
        
    res = exp_sc %>% 
        filter(cnv_state %in% c('del', 'amp')) %>%
        group_by(seg) %>%
        summarise(
            n = n(),
            cnv_state = unique(cnv_state),
            phi_mle = phi_mle(Y_obs, lambda_ref, depth_obs, alpha, beta, lower = 0.1, upper = 10)@coef,
            l11 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = 1),
            l10 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = 0.5),
            l21 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = 1.5),
            l31 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = 2)
        ) %>%
        mutate(alpha = alpha, beta = beta)
    
    return(res)
}

permute_phi_vec = function(lambda_mat_obs, lambda_ref, n_perm) {
    perm_sample = replicate(
        n = n_perm, {
            
            n_cells = ncol(lambda_mat_obs)
            n_genes = nrow(lambda_mat_obs)
            
            swap_index = sample(1:n_genes, size = n_genes/2, replace = FALSE)

            lambda_mat_ref = matrix(lambda_ref, ncol = n_cells, nrow = n_genes)
            
            lambda_mat_obs_perm = lambda_mat_obs
            lambda_mat_ref_perm = lambda_mat_ref
            lambda_mat_obs_perm[swap_index,] = lambda_mat_ref[swap_index,]
            lambda_mat_ref_perm[swap_index,] = lambda_mat_obs[swap_index,]
            
            phi_perm = colSums(lambda_mat_obs_perm)/colSums(lambda_mat_ref_perm)
    })
        
    return(perm_sample)
}


########################### Plotting ############################

cnv_colors = scale_color_manual(
    values = c("neu" = "gray", "del_up" = "royalblue", "del_down" = "darkblue", 
               "loh_up" = "darkgreen", "loh_down" = "olivedrab4",
               "amp_up" = "red", "amp_down" = "tomato3",
               "bamp" = "salmon", "bdel" = "skyblue",
              "amp" = "darkred", "loh" = "darkgreen", "del" = "darkblue", "neu2" = "gray30")
)


plot_psbulk = function(Obs, dot_size = 0.8, exp_limit = 2, min_depth = 10) {
    
    p = ggplot(
        Obs %>% 
            mutate(logFC = ifelse(logFC > exp_limit | logFC < -exp_limit, NA, logFC)) %>%
            mutate(pBAF = ifelse(DP >= min_depth, pBAF, NA)) %>%
            reshape2::melt(measure.vars = c('logFC', 'pBAF')),
        aes(x = snp_index, y = value, color = state),
        na.rm=TRUE
    ) +
    geom_point(size = dot_size, alpha = 0.5, pch = 16) +
    geom_hline(data = data.frame(variable = 'logFC'), aes(yintercept = 0), color = 'gray30', linetype = 'dashed') +
    geom_line(
        inherit.aes = FALSE,
        data = Obs %>% mutate(variable = 'logFC'),
        aes(x = snp_index, y = log2(phi_mle_roll), group = '1'),
        color = 'darkred'
    ) +
    theme_classic() +
    theme(
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
        strip.background = element_blank()
    ) +
    facet_grid(variable ~ CHROM, scale = 'free', space = 'free_x') +
    cnv_colors
    
    return(p)
}


########################### Main ############################

#' @param exp_mat scaled expression matrix
armadillo = function(count_mat, lambdas_ref, df, cell_annot, ncores, p_0 = 1e-6, verbose = TRUE) {
    
    #### 1. Build expression tree ####
    
    if (verbose) {
        display('Building expression tree ..')
    }
    
    tree = calc_cluster_tree(count_mat, cell_annot)
    
    # internal nodes
    nodes = get_internal_nodes(as.dendrogram(tree), '0', data.frame())
    
    # add terminal modes
    nodes = nodes %>% rbind(
        data.frame(
            node = cell_annot %>% filter(group != 'ref') %>% pull(cell_type) %>% unique
        ) %>%
        mutate(cluster = node)
    )
    
    # convert to list
    nodes = nodes %>%
        group_by(node) %>%
        summarise(
            members = list(cluster)
        ) %>%
        {setNames(.$members, .$node)} %>%
        lapply(function(members) {

            node = list(
                members = members,
                cells = cell_annot %>% filter(cell_type %in% members) %>% pull(cell)
            ) %>%
            inset('size', length(.$cells))

            return(node)
        })
    
    nodes = lapply(names(nodes), function(n) {nodes[[n]] %>% inset('label', n)}) %>% setNames(names(nodes))
    
    nodes = nodes %>% extract(map(., 'size') > 300)
    
    #### 2. Run pseudobulk HMM ####
    
    if (verbose) {
        display('Running HMMs ..')
    }
            
    bulk_all = mclapply(
            nodes,
            mc.cores = ncores,
            function(node) {

                bulk = get_bulk(
                    count_mat = count_mat[,node$cells],
                    df = df %>% filter(cell %in% node$cells),
                    lambdas_ref = lambdas_ref,
                    gtf_transcript = gtf_transcript
                ) %>%
                mutate(node = node$label) %>%
                analyze_bulk_gpois(p_0 = p_0, verbose = TRUE)

                return(bulk)
        }) %>%
        Reduce(rbind, .) %>%
        arrange(CHROM, POS) %>%
        mutate(snp_id = factor(snp_id, unique(snp_id))) %>%
        mutate(snp_index = as.integer(snp_id)) %>%
        arrange(node)
    
    #### 3. Find consensus CNVs ####
    
    if (verbose) {
        display('Finding consensus CNVs ..')
    }
        
    segs_all = bulk_all %>% 
        filter(cnv_state != 'neu') %>%
        split(.$node) %>%
        mclapply(
            mc.cores = ncores,
            function(bulk_node) {          
                bulk_node %>% 
                filter(cnv_state != 'neu') %>%
                group_by(CHROM, seg, cnv_state, node) %>%
                summarise(
                    n_genes = length(na.omit(unique(gene))),
                    n_snps = sum(!is.na(pAD)),
                    seg_start = min(POS),
                    seg_end = max(POS),
                    LLR_exp = LLR_exp(
                        Y_obs[!is.na(Y_obs)],
                        lambda_ref[!is.na(Y_obs)],
                        unique(na.omit(d_obs)),
                        unique(alpha_hat),
                        unique(beta_hat)
                    ),
                    LLR_allele = LLR_allele(
                        pAD[!is.na(pAD)],
                        DP[!is.na(pAD)],
                        p_s[!is.na(pAD)]
                    ),
                    LLR = LLR_exp + LLR_allele
                ) %>%
                ungroup()       
            }
        ) %>%
        bind_rows()
    
    # plot HMMs
    options(warn = -1)
    plot_list = bulk_all %>%
        split(.$node) %>%
        lapply(
            function(node_bulk) {
                node = unique(node_bulk$node)

                p = plot_psbulk(node_bulk) + 
                    theme(
                        title = element_text(size = 8),
                        axis.text.x = element_blank(),
                        axis.title = element_blank()
                    ) +
                    ggtitle(paste0(node, '(n=', nodes[[node]]$size, ')', ': ', paste0(nodes[[node]]$members, collapse = ', ')))

                return(p)
            }
        )
    options(warn = 0)
        
    segs_filtered = segs_all %>% filter(LLR_allele > 10)
    
    # resolve overlapping calls by graph
    V = segs_filtered %>% mutate(vertex = 1:n(), .before = 'CHROM')

    E = segs_filtered %>% {GenomicRanges::GRanges(
            seqnames = .$CHROM,
            IRanges::IRanges(start = .$seg_start,
                   end = .$seg_end)
        )} %>%
        GenomicRanges::findOverlaps(., .) %>%
        as.data.frame %>%
        setNames(c('from', 'to')) %>% 
        filter(from != to)

    G = igraph::graph_from_data_frame(d=E, vertices=V, directed=F)

    segs_filtered = segs_filtered %>% mutate(group = igraph::components(G)$membership)
    
    segs_consensus = segs_filtered %>% arrange(CHROM, group, -LLR) %>% distinct(group, `.keep_all` = TRUE) 
    
    #### 4. Per-cell CNV evaluations ####
    
    if (verbose) {
        display('Calculating per-cell CNV posteriors ..')
    }
    
    chrom_sizes = fread('~/ref/hg38.chrom.sizes.txt') %>% 
        set_names(c('CHROM', 'LEN')) %>%
        mutate(CHROM = str_remove(CHROM, 'chr')) %>%
        filter(CHROM %in% 1:22) %>%
        mutate(CHROM = as.integer(CHROM))

    segs_consensus = c(
            segs_consensus %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$seg_start,
                       end = .$seg_end)
            )},
            chrom_sizes %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = 1,
                       end = .$LEN)
            )}
        ) %>%
        GenomicRanges::disjoin() %>%
        as.data.frame() %>%
        select(CHROM = seqnames, seg_start = start, seg_end = end) %>%
        left_join(
            segs_consensus,
            by = c("CHROM", "seg_start", "seg_end")
        ) %>%
        mutate(cnv_state = tidyr::replace_na(cnv_state, 'neu')) %>%
        group_by(CHROM) %>%
        mutate(seg_cons = paste0(CHROM, '_', 1:n())) %>%
        ungroup() %>%
        mutate(CHROM = as.integer(CHROM))
    
    # gene expression posteriors
    gene_seg = GenomicRanges::findOverlaps(
            gtf_transcript %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$gene_start,
                       end = .$gene_end)
            )}, 
            segs_consensus %>% {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$seg_start,
                       end = .$seg_end)
            )}
        ) %>%
        as.data.frame() %>%
        set_names(c('gene_index', 'seg_index')) %>%
        left_join(
            gtf_transcript %>% mutate(gene_index = 1:n()),
            by = c('gene_index')
        ) %>%
        left_join(
            segs_consensus %>% mutate(seg_index = 1:n()),
            by = c('seg_index', 'CHROM')
        ) %>%
        distinct(gene, `.keep_all` = TRUE) 
    
    cells = colnames(count_mat)

    exp_sc = count_mat %>%
        as.data.frame() %>%
        tibble::rownames_to_column('gene') %>% 
        mutate(lambda_ref = lambdas_ref[gene]) %>%
        filter(lambda_ref > 0) %>%
        left_join(
            gene_seg %>% select(gene, seg = seg_cons, cnv_state),
            by = "gene"
        )
    
    priors = list(
        'amp' = c('21' = 1/4, '31' = 1/4, '11' = 1/2, '10' = 0),
        'del' = c('21' = 0, '31' = 0, '11' = 1/2, '10' = 1/2),
        'loh' = c('21' = 0, '31' = 0, '11' = 1/2, '10' = 1/2),
        'neu' = c('21' = 1/6, '31' = 1/6, '11' = 1/3, '10' = 1/3)
    )

    pi = lapply(priors, log)

    exp_post = mclapply(
            cells,
            mc.cores = ncores,
            function(cell) {
                exp_sc[,c('gene', 'lambda_ref', 'seg', 'cnv_state', cell)] %>%
                set_names(c('gene', 'lambda_ref', 'seg', 'cnv_state', 'Y_obs')) %>%
                get_exp_post(pi) %>%
                mutate(cell = cell)
            }
        ) %>%
        bind_rows()
    
    exp_post = exp_post %>% 
        rowwise() %>%
        mutate(
            l = matrixStats::logSumExp(c(l11 + pi[[cnv_state]]['11'], l10 + pi[[cnv_state]]['10'], l21 + pi[[cnv_state]]['21'], l31 + pi[[cnv_state]]['31'])),
            p_amp = exp(matrixStats::logSumExp(c(l21 + pi[[cnv_state]]['21'], l31 + pi[[cnv_state]]['31'])) - l),
            p_amp_mul = exp(l31 + pi[[cnv_state]]['31'] - l),
            p_amp_sin = exp(l21 + pi[[cnv_state]]['21'] - l),
            p_neu = exp(l11 + pi[[cnv_state]]['11'] - l),
            p_del = exp(l10 + pi[[cnv_state]]['10'] - l),
            l_amp = matrixStats::logSumExp(c(l21 + log(0.5), l31 + log(0.5))),
            l_del = l10,
            l_n = l11,
            l_t = ifelse(cnv_state == 'amp', l_amp, l_del)
        ) %>%
        ungroup()
    
    # alele posteriors
    snp_seg = bulk_all %>%
        filter(!is.na(pAD)) %>%
        mutate(haplo = case_when(
            str_detect(state, 'up') ~ 'major',
            str_detect(state, 'down') ~ 'minor',
            T ~ ifelse(pBAF > 0.5, 'major', 'minor')
        )) %>%
        select(snp_id, snp_index, node, seg, haplo) %>%
        inner_join(
            segs_consensus,
            by = c('node', 'seg')
        )
    
    allele_sc = df %>%
        mutate(pAD = ifelse(GT == '1|0', AD, DP - AD)) %>%
        select(-snp_index) %>% 
        inner_join(
            snp_seg %>% select(snp_id, snp_index, haplo, seg = seg_cons, cnv_state),
            by = c('snp_id')
        ) %>%
        mutate(
            major_count = ifelse(haplo == 'major', pAD, DP - pAD),
            minor_count = DP - major_count,
            MAF = major_count/DP
        ) %>%
        group_by(cell, CHROM) %>% 
        arrange(cell, CHROM, POS) %>%
        mutate(
            n_chrom_snp = n(),
            inter_snp_dist = ifelse(n_chrom_snp > 1, c(NA, POS[2:length(POS)] - POS[1:(length(POS)-1)]), NA)
        ) %>%
        ungroup() %>%
        filter(inter_snp_dist > 250 | is.na(inter_snp_dist))
    
    allele_post = allele_sc %>%
        group_by(cell, seg, cnv_state) %>%
        summarise(
            major = sum(major_count),
            minor = sum(minor_count),
            total = major + minor,
            MAF = major/total,
            .groups = 'drop'
        ) %>%
        rowwise() %>%
        mutate(
            l11 = dbinom(major, total, p = 0.5, log = TRUE),
            l10 = dbinom(major, total, p = 0.9, log = TRUE),
            l21 = dbinom(major, total, p = 0.66, log = TRUE),
            l31 = dbinom(major, total, p = 0.75, log = TRUE),
            l = matrixStats::logSumExp(c(l11 + pi[[cnv_state]]['11'], l10 + pi[[cnv_state]]['10'], l21 + pi[[cnv_state]]['21'], l31 + pi[[cnv_state]]['31'])),
            p_amp = exp(matrixStats::logSumExp(c(l21 + pi[[cnv_state]]['21'], l31 + pi[[cnv_state]]['31'])) - l),
            p_amp_mul = exp(l31 + pi[[cnv_state]]['31'] - l),
            p_amp_sin = exp(l21 + pi[[cnv_state]]['21'] - l),
            p_neu = exp(l11 + pi[[cnv_state]]['11'] - l),
            p_del = exp(l10 + pi[[cnv_state]]['10'] - l),
            p_imbalance = exp(matrixStats::logSumExp(c(l10 + pi[[cnv_state]]['10'], l21 + pi[[cnv_state]]['21'], l31 + pi[[cnv_state]]['31'])) - l),
            l_amp = matrixStats::logSumExp(c(l21 + log(0.5), l31 + log(0.5))),
            l_del = l10,
            l_n = l11,
            l_t = ifelse(cnv_state == 'amp', l_amp, l_del)
        ) %>%
        ungroup()
    
    # joint posteriors
    joint_post = exp_post %>%
        rename(l_n_exp = l_n, l_t_exp = l_t, l_exp = l) %>%
        full_join(
            allele_post %>% select(
                cell, seg, l_n_ar = l_n, l_t_ar = l_t, l_ar = l,
                n_snp = total
            ),
            c("cell", "seg")
        ) %>%
        mutate(cnv_state = tidyr::replace_na(cnv_state, 'loh')) %>%
        mutate_at(
            c('l_n_ar', 'l_t_ar', 'l_n_exp', 'l_t_exp', 'l_exp', 'l_ar'),
            function(x) tidyr::replace_na(x, 0)
        ) %>%
        rowwise() %>%
        mutate(
            l_n = l_n_exp + l_n_ar,
            l_t = l_t_exp + l_t_ar,
            l = matrixStats::logSumExp(c(l_n + log(0.5), l_t + log(0.5))),
            p_t = exp(l_t + log(0.5) - l),
            p_t_ar = exp(l_t_ar + log(0.5) - l_ar),
            p_t_exp = exp(l_t_exp + log(0.5) - l_exp)
        ) %>%
        ungroup() %>%
        left_join(cell_annot, by = 'cell')
    
    if (verbose) {
        display('All done!')
    }
    
    res = list(
        'nodes' = nodes,
        'segs_filtered' = segs_filtered,
        'segs_all' = segs_all,
        'segs_consensus' = segs_consensus,
        'bulk_all' = bulk_all,
        'tree' = tree,
        'plot_list' = plot_list,
        'graph' = G,
        'joint_post' = joint_post,
        'exp_post' = exp_post,
        'allele_post' = allele_post
    )
    
    return(res)
}
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

analyze_bulk_gpois = function(Obs, t, gamma = 20, bal_cnv = TRUE, prior = NULL, exp_only = FALSE, allele_only = FALSE, retest = TRUE, verbose = TRUE) {
    
    # doesn't work with 0s in the ref
    Obs = Obs %>% filter(lambda_ref != 0 | is.na(gene)) 
    
    x = find_diploid(Obs, gamma = gamma, verbose = verbose)
    
    Obs = x$bulk
    bal_cnv = x$bamp
    
    fit = Obs %>%
            filter(!is.na(Y_obs)) %>%
            filter(logFC < 8 & logFC > -8) %>%
            filter(diploid) %>%
            {fit_gpois(.$Y_obs, .$lambda_ref, unique(.$d_obs))}
            
    alpha_hat = fit@coef[1]
    beta_hat = fit@coef[2]
        
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
                phi_amp = 2^(0.25),
                phi_del = 2^(-0.25),
                phi_bamp = 2^(0.25),
                phi_bdel = 2^(-0.25),
                alpha = alpha_hat,
                beta = beta_hat,
                t = t,
                gamma = gamma,
                prior = prior,
                bal_cnv = bal_cnv,
                exp_only = exp_only,
                allele_only = allele_only
            )
        ) %>% 
        annot_segs %>%
        ungroup()
    
    Obs = Obs %>% mutate(alpha = alpha_hat, beta = beta_hat, gamma = gamma)
    
    if (retest) {
        
        if (verbose) {
            display('Retesting CNVs..')
        }

        segs_post = retest_cnv(Obs)

        Obs = Obs %>% 
            left_join(segs_post, by = c('seg', 'CHROM')) %>%
            mutate(cnv_state_post = tidyr::replace_na(cnv_state_post, 'neu')) %>%
            mutate(state_post = ifelse(
                cnv_state_post %in% c('amp', 'del', 'loh'),
                paste0(cnv_state_post, '_', str_extract(state, 'up|down')),
                cnv_state_post
            ))
        
    }
    
    if (verbose) {
        display('Finishing..')
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
    
    return(Obs)
}

annot_segs = function(Obs) {
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
}

find_diploid = function(bulk, gamma = 20, t = 1e-5, fc_min = 1.25, debug = F, verbose = T) {
    
    # define diploid regions
    bulk = bulk %>% group_by(CHROM) %>%
        mutate(state = 
            run_hmm_mv_inhom_gpois(
                pAD = pAD,
                DP = DP, 
                p_s = p_s,
                Y_obs = Y_obs, 
                lambda_ref = lambda_ref, 
                d_total = 0,
                t = t,
                allele_only = TRUE
            )
        ) %>% ungroup() %>%
        annot_segs()
    
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
    tests = t(combn(as.character(segs$seg), 2)) %>% 
        as.data.frame() %>%
        set_names(c('i', 'j')) %>%
        rowwise() %>%
        mutate(
            p = ks.test(x = bulk_balanced$FC[bulk_balanced$seg == i],
                        y = bulk_balanced$FC[bulk_balanced$seg == j])$p.value
        ) %>%
        ungroup()
    options(warn = 0)
    
    # build graph
    V = segs

    E = tests %>% filter(p > 0.01) 

    G = igraph::graph_from_data_frame(d=E, vertices=V, directed=F)
    
    # find confident diploid clique
    cliques = igraph::maximal.cliques(G)

    FC = sapply(cliques, function(c) {
        bulk_balanced %>% filter(seg %in% names(c)) %>% {phi_hat_seg(.$Y_obs, .$lambda_ref, unique(.$d_obs))}
    })
    
    # seperation needs to be clear (2 vs 4 copy)
    if (max(FC)/min(FC) > fc_min) {
        if (verbose) {display('quadruploid state found!')}
        diploid_segs = names(cliques[[which.min(FC)]])
        bamp = TRUE
    } else {
        diploid_segs = segs$seg
        bamp = FALSE
    }
    
    if (verbose) {
        display(glue('diploid regions: {paste0(gtools::mixedsort(diploid_segs), collapse = ",")}'))
    }
    
    bulk = bulk %>% mutate(diploid = seg %in% diploid_segs)
    
    if (debug) {
        return(list('G' = G, 'tests' = tests, 'segs' = segs, 'FC' = FC,
                    'bulk_balanced' = bulk_balanced, 'bulk' = bulk))
    }
    
    return(list('bulk' = bulk, 'bamp' = bamp))
    
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


calc_phi_mle = function(Y_obs, lambda_ref, d, alpha, beta, lower = 0.2, upper = 10) {
    
    if (length(Y_obs) == 0) {
        return(1)
    }
        
    start = max(min(1, upper), lower)
    
    res = stats4::mle(
        minuslogl = function(phi) {
            res = -sum(dgpois(Y_obs, shape = alpha, rate = beta/(phi * d * lambda_ref), log = TRUE))
            return(res)
        },
        start = start,
        lower = lower,
        upper = upper
    )
    
    return(res@coef)
}

phi_mle_roll = function(Y_obs, lambda_ref, alpha, beta, d_obs, h) {
    n = length(Y_obs)
    sapply(
        1:n,
        function(c) {
            slice = max(c - h - 1, 1):min(c + h, n)
            calc_phi_mle(Y_obs[slice], lambda_ref[slice], unique(d_obs), alpha, beta)
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
    
    phi_mle = calc_phi_mle(Y_obs, lambda_ref, d, alpha, beta)
    l_1 = l_gpois(Y_obs, lambda_ref, d, alpha, beta, phi = phi_mle)
    l_0 = l_gpois(Y_obs, lambda_ref, d, alpha, beta, phi = 1)
    return(l_1 - l_0)
}

approx_lik_exp = function(Y_obs, lambda_ref, d_obs, alpha, beta, lower, upper, step_size = 0.01) {
    
    if (length(Y_obs) == 0) {
        return(0)
    }
    
    phi_grid = seq(lower, upper, step_size)
        
    l_grid = sapply(phi_grid, function(phi) {
        l_gpois(Y_obs, lambda_ref, d_obs, alpha, beta, phi = phi)
    })
    
    l = matrixStats::logSumExp(l_grid)
    
    return(l)
}

approx_lik_ar = function(pAD, DP, p_s, lower, upper, step_size = 0.01, gamma = 16) {
    
    if (length(pAD) == 0) {
        return(0)
    }
    
    theta_grid = seq(lower, upper, step_size)
        
    l_grid = sapply(theta_grid, function(theta) {
        calc_alelle_lik(pAD, DP, p_s, theta = theta, gamma = gamma)
    })
    
    l = matrixStats::logSumExp(l_grid)
    
    return(l)
}

retest_cnv = function(bulk) {
    
    G = c('20' = 1/5, '10' = 1/5, '21' = 1/10, '31' = 1/10, '22' = 1/5, '00' = 1/5) %>% log
    
    theta_min = 0.065
    delta_phi_min = 0.15
    
    segs_post = bulk %>% 
        group_by(CHROM, seg) %>%
        filter(state != 'neu') %>%
        summarise(
            n_genes = length(na.omit(unique(gene))),
            n_snps = sum(!is.na(pAD)),
            seg_start = min(POS),
            seg_end = max(POS),
            phi_mle = calc_phi_mle(Y_obs[!is.na(Y_obs)], lambda_ref[!is.na(Y_obs)], unique(na.omit(d_obs)), unique(alpha), unique(beta)),
            theta_mle = calc_theta_mle(pAD[!is.na(pAD)], DP[!is.na(pAD)], p_s[!is.na(pAD)]),
            l_x_n = approx_lik_exp(Y_obs[!is.na(Y_obs)], lambda_ref[!is.na(Y_obs)], unique(na.omit(d_obs)), unique(alpha), unique(beta), lower = 1 - delta_phi_min, upper = 1 + delta_phi_min),
            l_x_a = approx_lik_exp(Y_obs[!is.na(Y_obs)], lambda_ref[!is.na(Y_obs)], unique(na.omit(d_obs)), unique(alpha), unique(beta), lower = 1 + delta_phi_min, upper = 5),
            l_x_d = approx_lik_exp(Y_obs[!is.na(Y_obs)], lambda_ref[!is.na(Y_obs)], unique(na.omit(d_obs)), unique(alpha), unique(beta), lower = 0.1, upper = 1 - delta_phi_min),
            l_y_n = approx_lik_ar(pAD[!is.na(pAD)], DP[!is.na(pAD)], p_s[!is.na(pAD)], gamma = unique(gamma), lower = 0, upper = theta_min),
            l_y_d = approx_lik_ar(pAD[!is.na(pAD)], DP[!is.na(pAD)], p_s[!is.na(pAD)], gamma = unique(gamma), lower = theta_min, upper = 0.49),
            l_y_a = approx_lik_ar(pAD[!is.na(pAD)], DP[!is.na(pAD)], p_s[!is.na(pAD)], gamma = unique(gamma), lower = theta_min, upper = 0.375),
            Z = matrixStats::logSumExp(
                c(G['20'] + l_x_n + l_y_d,
                 G['10'] + l_x_d + l_y_d,
                 G['21'] + l_x_a + l_y_a,
                 G['31'] + l_x_a + l_y_a,
                 G['22'] + l_x_a + l_y_n, 
                 G['00'] + l_x_d + l_y_n)),
            p_loh = exp(G['20'] + l_x_n + l_y_d - Z),
            p_amp = exp(matrixStats::logSumExp(c(G['21'] + l_x_a + l_y_a, G['31'] + l_x_a + l_y_a)) - Z),
            p_del = exp(G['10'] + l_x_d + l_y_d - Z),
            p_bamp = exp(G['22'] + l_x_a + l_y_n - Z),
            p_bdel = exp(G['00'] + l_x_d + l_y_n - Z),
            LLR_ar = l_y_d - l_y_n,
            LLR = Z - (l_x_n + l_y_n),
            .groups = 'drop'
        ) %>%
        rowwise() %>%
        mutate(cnv_state_post = c('loh', 'amp', 'del', 'bamp', 'bdel')[
            which.max(c(p_loh, p_amp, p_del, p_bamp, p_bdel))
        ]) %>%
        ungroup()

    return(segs_post)
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

calc_allele_lik = function(pAD, DP, p_s, theta, gamma = 20) {
            
    # states
    states = c("theta_up", "theta_down")
        
    # transition matrices
    calc_trans_mat = function(p_s) {
        matrix(
            c(1 - p_s, p_s,
              p_s, 1 - p_s),
            ncol = 2,
            byrow = TRUE
        )
    }
    
    As = lapply(
        p_s,
        function(p_s) {calc_trans_mat(p_s)}
    )
    
    # intitial probabilities
    prior = c(0.5, 0.5)
    
    alpha_up = (0.5 + theta) * gamma
    beta_up = (0.5 - theta) * gamma
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
    
    if(length(pAD) == 0) {
        return(0)
    }
    
    res = optim(
        0.25, 
        function(theta) {-calc_allele_lik(pAD, DP, p_s, theta)},
        method = 'L-BFGS-B',
        lower = lower, upper = upper)
    return(res$par)
}


# multi-state model
get_exp_likelihoods = function(exp_sc) {
    
    depth_obs = sum(exp_sc$Y_obs)
        
    fit = exp_sc %>% filter(cnv_state == 'neu') %>%
        {fit_gpois(.$Y_obs, .$lambda_ref, depth_obs)}
    
    alpha = fit@coef[1]
    beta = fit@coef[2]
        
    res = exp_sc %>% 
        filter(cnv_state %in% c('del', 'amp', 'bamp', 'bdel')) %>%
        group_by(seg) %>%
        summarise(
            n = n(),
            cnv_state = unique(cnv_state),
            phi_mle = calc_phi_mle(Y_obs, lambda_ref, depth_obs, alpha, beta, lower = 0.1, upper = 10),
            l11 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = 1),
            l20 = l11,
            l10 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = 0.5),
            l21 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = 1.5),
            l31 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = 2),
            l22 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = 2),
            l00 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = 0.25)
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
        aes(x = snp_index, y = value, color = state_post),
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
    cnv_colors +
    guides(color = guide_legend(""))
    
    return(p)
}
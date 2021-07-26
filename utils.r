########################### Data processing ############################

# take in a reference matrix, decide which one is best 
choose_ref = function(count_mat_obs, count_mat_ref, gtf_transcript, debug = F) {
    
    errs = apply(count_mat_ref, 2, function(lambdas_ref) {
            
            genes_common = gtf_transcript$gene %>% 
                intersect(rownames(count_mat_obs)) %>%
                intersect(names(lambdas_ref)[lambdas_ref > 0])

            Y_obs = rowSums(count_mat_obs)[genes_common]
            lambdas_ref = lambdas_ref[genes_common]

            tibble(
                l_total = sum(dpois(x = Y_obs, lambda = lambdas_ref * sum(Y_obs), log = TRUE)),
                mse = mean((log2((Y_obs/sum(Y_obs) + 1e-6)/lambdas_ref))^2),
                n = length(genes_common)
            ) %>%
            mutate(l_bar = l_total/n)

        }) %>% 
        bind_rows() %>%
        mutate(ref = colnames(count_mat_ref))
    
    if (debug) {
        return(errs)
    }

    return(errs$ref[which.min(-errs$l_bar)])
    
}

process_exp_fc = function(count_mat_obs, lambdas_ref, gtf_transcript, verbose = TRUE) {
    
    genes_annotated = gtf_transcript$gene %>% 
        intersect(rownames(count_mat_obs)) %>%
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

get_bulk = function(count_mat, lambdas_ref, df, gtf_transcript, min_depth = 2, verbose = FALSE) {

    
    
    gexp_bulk = process_exp_fc(
        count_mat,
        lambdas_ref,
        gtf_transcript,
        verbose = verbose
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
    
    if (!is.numeric(t)) {
        stop('transition probability is not numeric')
    }
    
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
        select(-any_of(c('gamma'))) %>%
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

    Obs = Obs %>%
        mutate(alpha = alpha_hat, beta = beta_hat, gamma = gamma) %>%
        mutate(haplo_post = case_when(
            str_detect(state, 'up') ~ 'major',
            str_detect(state, 'down') ~ 'minor',
            T ~ ifelse(pBAF > 0.5, 'major', 'minor')
        )) %>% 
        mutate(
            major_count = ifelse(haplo_post == 'major', pAD, DP - pAD),
            minor_count = DP - major_count
        )
    
    # rolling theta estimates
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
    

    if (retest) {
        
        if (verbose) {
            display('Retesting CNVs..')
        }

        segs_post = retest_cnv(Obs)
        
        Obs = Obs %>% 
            select(-any_of(colnames(segs_post)[!colnames(segs_post) %in% c('seg', 'CHROM')])) %>%
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
    
    # rolling phi estimates
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
    ) %>%
    group_by(CHROM) %>%
    mutate(phi_mle_roll = zoo::na.locf(phi_mle_roll, na.rm=FALSE)) %>%
    ungroup()

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

retest_cnv = function(bulk) {
    
    G = c('20' = 1/5, '10' = 1/5, '21' = 1/10, '31' = 1/10, '22' = 1/5, '00' = 1/5)
    
    theta_min = 0.065
    delta_phi_min = 0.15
    
    segs_post = bulk %>% 
        filter(cnv_state != 'neu') %>%
        group_by(CHROM, seg) %>%
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
            approx_lik_exp(Y_obs[!is.na(Y_obs)], lambda_ref[!is.na(Y_obs)], unique(na.omit(d_obs)), unique(alpha), unique(beta)),
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
            LLR_x = calc_exp_LLR(Y_obs[!is.na(Y_obs)], lambda_ref[!is.na(Y_obs)], unique(na.omit(d_obs)), unique(alpha), unique(beta), phi_mle),
            LLR_y = calc_allele_LLR(pAD[!is.na(pAD)], DP[!is.na(pAD)], p_s[!is.na(pAD)], theta_mle),
            LLR = LLR_x + LLR_y,
            .groups = 'drop'
        ) %>%
        rowwise() %>%
        mutate(cnv_state_post = c('loh', 'amp', 'del', 'bamp', 'bdel')[
            which.max(c(p_loh, p_amp, p_del, p_bamp, p_bdel))
        ]) %>%
        ungroup()

    return(segs_post)
}

approx_lik_ar = function(pAD, DP, p_s, lower = 0.001, upper = 0.499, start = 0.25, gamma = 20) {
    
    if (length(pAD) == 0) {
        return(tibble('theta_mle' = 0, 'theta_sigma' = 0))
    }

    fit = optim(
        start, 
        function(theta) {-calc_allele_lik(pAD, DP, p_s, theta, gamma)},
        method = 'L-BFGS-B',
        lower = lower,
        upper = upper,
        hessian = TRUE
#         control = list('factr' = 1e-4/(.Machine$double.eps))
    )
    
    mu = fit$par
    sigma = sqrt(as.numeric(1/(fit$hessian)))
    
    return(tibble('theta_mle' = mu, 'theta_sigma' = sigma))
}

approx_lik_exp = function(Y_obs, lambda_ref, d, alpha, beta, lower = 0.2, upper = 10, start = 1) {
    
    if (length(Y_obs) == 0) {
        return(tibble('phi_mle' = 1, 'phi_sigma' = 0))
    }
    
    start = max(min(1, upper), lower)
    
    fit = optim(
        start,
        function(phi) {
            -l_gpois(Y_obs, lambda_ref, d, alpha, beta, phi = phi)
        },
        method = 'L-BFGS-B',
        lower = lower,
        upper = upper,
        hessian = TRUE
    )

    mu = fit$par
    sigma = sqrt(as.numeric(1/(fit$hessian)))
    
    return(tibble('phi_mle' = mu, 'phi_sigma' = sigma))
}

pnorm.range = function(lower, upper, mu, sd) {
    if (sd == 0) {
        return(1)
    }
    pnorm(upper, mu, sd) - pnorm(lower, mu, sd)
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


calc_cluster_dist = function(count_mat, cell_annot) {

    cells = intersect(colnames(count_mat), cell_annot$cell)
    
    count_mat = count_mat %>% extract(,cells)
    cell_annot = cell_annot %>% filter(cell %in% cells)
    
    cell_dict = cell_annot %>% {setNames(.$cell_type, .$cell)} %>% droplevels()
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

calc_allele_LLR = function(pAD, DP, p_s, theta_mle, gamma = 20) {
    if (length(pAD) == 0) {
        return(0)
    }
    l_1 = calc_allele_lik(pAD, DP, p_s, theta = theta_mle) 
    l_0 = calc_allele_lik(pAD, DP, p_s, theta = 0)
    return(l_1 - l_0)
}

calc_exp_LLR = function(Y_obs, lambda_ref, depth_obs, alpha, beta, phi_mle) {
    if (length(Y_obs) == 0) {
        return(0)
    }
    l_1 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = phi_mle)
    l_0 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = 1)
    return(l_1 - l_0)
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
            l22 = l31,
            l00 = l10,
#             l11 = approx_lik_exp(Y_obs, lambda_ref, depth_obs, alpha, beta, lower = 0.75, upper = 1.25),
#             l20 = l11,
#             l10 = approx_lik_exp(Y_obs, lambda_ref, depth_obs, alpha, beta, lower = 0.1, upper = 0.75),
#             l21 = approx_lik_exp(Y_obs, lambda_ref, depth_obs, alpha, beta, lower = 1.25, upper = 1.75),
#             l31 = approx_lik_exp(Y_obs, lambda_ref, depth_obs, alpha, beta, lower = 1.75, upper = 3.5),
#             l22 = l31,
#             l00 = l10
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
              "amp" = "darkred", "loh" = "#34d834", "del" = "darkblue", "neu2" = "gray30")
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
        color = 'darkred',
        size = 0.35
    ) +
    theme_classic() +
    theme(
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
        strip.background = element_blank()
    ) +
    facet_grid(variable ~ CHROM, scale = 'free', space = 'free_x') +
    # geom_ribbon(
    #     inherit.aes = FALSE,
    #     data = Obs %>% filter(cnv_state_post != 'neu') %>% mutate(variable = 'pBAF'),
    #     aes(x = snp_index, ymin = 0.5 - theta_hat_roll, ymax = 0.5 + theta_hat_roll, fill = cnv_state_post),
    #     fill = NA,
    #     alpha = 0.2,
    #     color = 'white',
    #     size = 0.25
    # ) +
    cnv_colors +
    guides(color = guide_legend(""), fill = FALSE)
    
    return(p)
}

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

plot_cnv_post = function(segs_consensus) {
    segs_consensus %>% 
        filter(cnv_state != 'neu') %>%
        mutate(seg_label = paste0(seg, '_', cnv_state_post)) %>%
        mutate(seg_label = factor(seg_label, unique(seg_label))) %>%
        reshape2::melt(measure.vars = c('p_loh', 'p_amp', 'p_del', 'p_bamp', 'p_bdel'), value.name = 'p') %>%
        ggplot(
            aes(x = seg_label, y = variable, fill = p, label = round(p, 2))
        ) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
        geom_tile() +
        geom_text(color = 'white')
}


plot_cnv_cell = function(joint_post, n_sample = 100, rev_cnv = F, rev_cell = F, mut_clust = NULL) {
    
    cell_group_order = c('Adrenergic', 'SCP-like', 'Mesenchymal', 'Endothelial', 'Mast', 'Myofibroblasts', 'Immune')
    
    joint_post = joint_post %>% 
        mutate(p_cnv_y = exp(Z_cnv_y - Z_y)) %>%
        mutate(cell_group = case_when(
            cell_type %in% c('Th', 'NK', 'Tcyto', 'Monocytes', 'B', 'Treg', 'pDC', 'Plasma', 'mDC', 'Macrophages', 'ILC3') ~ 'Immune',
            T ~ as.character(cell_type)
        )) %>%
        mutate(cell_group = factor(cell_group, cell_group_order))
    
    tree_cnv = joint_post %>% 
        reshape2::dcast(seg ~ cell, value.var = 'p_cnv') %>%
        tibble::column_to_rownames('seg') %>%
        dist() %>%
        hclust
    
    if (rev_cnv) {tree_cnv = rev(tree_cnv)}
    
    cnv_order = tree_cnv %>% {.$labels[.$order]}
    
    cell_order = joint_post %>% 
            group_by(cell_group) %>%
            do(
                reshape2::dcast(., cell ~ seg, value.var = 'p_cnv') %>%
                tibble::column_to_rownames('cell') %>%
                dist() %>%
                hclust %>%
                {.$labels[.$order]} %>%
                as.data.frame()
            ) %>%
            set_names(c('cell_group', 'cell'))
    
    if (rev_cell) {cell_order = cell_order %>% arrange(-row_number())}
    
    cnv_clusters = dendextend::cutree(tree_cnv, 2) 
    
    joint_post_cell = joint_post %>% 
        mutate(cluster = cnv_clusters[seg]) %>%
        rbind(
            mutate(., cluster = 'all')
        ) %>%
        group_by(cell, cell_group, cluster) %>%
        summarise_at(
            vars(contains('Z_')),
            sum
        ) %>%
        rowwise() %>%
        mutate(
            Z = matrixStats::logSumExp(c(Z_cnv, Z_n)),
            Z_x = matrixStats::logSumExp(c(Z_cnv_x, Z_n_x)),
            Z_y = matrixStats::logSumExp(c(Z_cnv_y, Z_n_y)),
            p_cnv = exp(Z_cnv - Z),
            p_n = exp(Z_n - Z),
            p_cnv_x = exp(Z_cnv_x - Z_x),
            p_n_x = exp(Z_n_x - Z_x),
            p_cnv_y = exp(Z_cnv_y - Z_y),
            p_n_y = exp(Z_n_y - Z_y),
        ) %>%
        ungroup() %>%
        mutate(
            p_cnv = 1 - p.adjust(p_n, method = 'BH'),
            p_cnv_x = 1 - p.adjust(p_n_x, method = 'BH'),
            p_cnv_y = 1 - p.adjust(p_n_y, method = 'BH')
        )
    
    D = joint_post %>% 
        mutate(cluster = 'segment') %>%
        bind_rows(
            joint_post_cell %>% 
                reshape2::melt(measure.vars = c('p_cnv_x', 'p_cnv_y', 'p_cnv'), value.name = 'p_cnv') %>%
                mutate(seg_label = c('p_cnv' = 'joint', 'p_cnv_x' = 'expr', 'p_cnv_y' = 'allele')[as.character(variable)]) %>%
                filter(seg_label == 'joint' | cluster == 'all') %>%
                select(seg_label, cell, cell_group, p_cnv, cluster) %>%
                arrange(seg_label)
        ) %>%
        mutate(seg = factor(seg, cnv_order)) %>%
        arrange(seg) %>%
        mutate(seg_label = factor(seg_label, unique(seg_label))) %>%
        mutate(cell = factor(cell, cell_order$cell)) %>%
        group_by(cell_group) %>%
        filter(cell %in% sample(unique(cell), min(n_sample, length(unique(cell))))) %>%
        ungroup()
    
    p = ggplot(
            D %>% filter(cluster %in% c('segment')),
            aes(x = cell, y = seg_label, fill = p_cnv)
        ) +
        geom_tile() +
        theme_classic() +
        scale_y_discrete(expand = expansion(0)) +
        scale_x_discrete(expand = expansion(0)) +
        theme(
            panel.spacing = unit(0.1, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'white', fill = NA),
            panel.background = element_rect(fill = 'skyblue'),
            strip.background = element_blank(),
            axis.text.x = element_blank(),
            strip.text = element_text(angle = 90, size = 8, vjust = 0.5)
        ) +
        facet_grid(cluster~cell_group, scale = 'free', space = 'free') +
        scale_fill_gradient2(low = 'skyblue', high = 'red', mid = 'skyblue', midpoint = 0.5, limits = c(0,1)) +
        ylab('')
    
    # p = ggplot(
    #         D %>% filter(cluster %in% c('segment')),
    #         aes(x = cell, y = seg_label, fill = logBF)
    #     ) +
    #     geom_tile() +
    #     theme_classic() +
    #     scale_y_discrete(expand = expansion(0)) +
    #     scale_x_discrete(expand = expansion(0)) +
    #     theme(
    #         panel.spacing = unit(0.1, 'mm'),
    #         panel.border = element_rect(size = 0.5, color = 'white', fill = NA),
    #         panel.background = element_rect(fill = 'skyblue'),
    #         strip.background = element_blank(),
    #         axis.text.x = element_blank(),
    #         strip.text = element_text(angle = 90, size = 8, vjust = 0.5)
    #     ) +
    #     facet_grid(cluster~cell_group, scale = 'free', space = 'free') +
    #     scale_fill_gradient(low = 'skyblue', high = 'red', limits = c(0, 10), oob = scales::oob_squish) +
    #     ylab('')
    
    p_cnv_tree = tree_cnv %>%
        ggtree::ggtree(ladderize=F) +
        scale_x_discrete(expand = expansion(mult = 0)) +
        theme(plot.margin = margin(0,0,0,0)) 
    
    panel_1 = (p_cnv_tree | p) + plot_layout(widths = c(1,20))
    
    ## split plot
    D = joint_post %>%
        mutate(seg = factor(seg, rev(cnv_order))) %>%
        arrange(seg) %>%
        mutate(
            seg_label = factor(seg_label, unique(seg_label))
        ) %>%
        mutate(cell = factor(cell, cell_order$cell)) %>%
        group_by(cell_group) %>%
        filter(cell %in% sample(unique(cell), min(200, length(unique(cell))))) %>%
        ungroup()

    p = D %>% 
        reshape2::melt(measure.vars = c('p_cnv', 'p_cnv_x', 'p_cnv_y'), value.name = 'p_cnv') %>%
        mutate(source = c('p_cnv' = 'joint', 'p_cnv_x' = 'expr', 'p_cnv_y' = 'allele')[as.character(variable)]) %>%
        # filter(!(cnv_state %in% c('bamp') & source %in% c('allele', 'joint'))) %>%
        # filter(!(cnv_state %in% c('loh') & source %in% c('expr', 'joint'))) %>%
        filter(source != 'joint') %>%
        ggplot(
            aes(x = cell, y = source, fill = p_cnv)
        ) +
        geom_tile() +
        theme_classic() +
        scale_y_discrete(expand = expansion(0), position = "right") +
        scale_x_discrete(expand = expansion(0)) +
        facet_grid(seg_label~cell_group, scale = 'free', space = 'free_x', switch = 'y') +
        theme(
            panel.spacing = unit(0.1, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'white', fill = NA),
            panel.background = element_rect(fill = 'skyblue'),
            strip.background = element_blank(),
            axis.text.x = element_blank(),
            strip.text.x = element_text(angle = 90, size = 8),
            strip.text.y.left = element_text(angle = 0, size = 8),
            plot.margin = margin(0,0,0,0)
        ) +
        scale_fill_gradient2(low = 'skyblue', high = 'red', mid = 'skyblue', midpoint = 0.5, limits = c(0,1)) +
        ylab('')
    
    panel_2 = (p_cnv_tree | p) + plot_layout(widths = c(1,20))
    
    return(panel_1)
    
}
source('~/Armadillo/hmm.r')
require(LaplacesDemon)

annotate_segs = function(snps, segs_union) {
    
    snps = snps %>% select(-contains('seg_union'))
    
    snp = snps %>% 
        {GenomicRanges::GRanges(
            seqnames = .$CHROM,
            IRanges::IRanges(start = .$snp_index,
                   end = .$snp_index)
        )}
    
    seg = segs_union %>% {GenomicRanges::GRanges(
        seqnames = .$CHROM,
        IRanges::IRanges(start = .$start,
               end = .$end)
    )}
        
    overlaps = IRanges::mergeByOverlaps(snp, seg) %>%
        as.data.frame() %>%
        select(CHROM = snp.seqnames, snp_index = snp.start, seg_union_start = seg.start, seg_union_end = seg.end) %>%
        mutate(seg_union = paste0(CHROM, '_', seg_union_start, '_', seg_union_end))
    
    snps %>% left_join(overlaps, by = c("CHROM", "snp_index"))
    
}

compare_event = function(chr, start, end, state_1, state_2, Obs_cut, plot = FALSE) {
                
    D1 = Obs_cut %>% 
        filter(cut == 1) %>%
        filter(CHROM == chr & snp_index >= start & snp_index <= end)
    
    D2 = Obs_cut %>% 
        filter(cut == 2) %>%
        filter(CHROM == chr & snp_index >= start & snp_index <= end)
    
    if (state_1 == 'neu') {
        D1 = D1 %>% mutate(
            state = run_hmm_inhom(
                pAD, DP, p_s,
                p_0 = 1,
                prior = c(0.5, 0, 0.5)
        ))
    }
    
    if (state_2 == 'neu') {
        D2 = D2 %>% mutate(
            state = run_hmm_inhom(
                pAD, DP, p_s,
                p_0 = 1,
                prior = c(0.5, 0, 0.5)
        ))
    }
    
    D1 = D1 %>% mutate(pBAF_post = ifelse(str_detect(state, 'down'), 1-pBAF, pBAF))

    D2 = D2 %>% mutate(pBAF_post = ifelse(str_detect(state, 'down'), 1-pBAF, pBAF))

    p = t.test(D1$pBAF_post, D2$pBAF_post)$p.value
        
    if (plot) {
        
        p1 = ggplot(
            D1,
            aes(x = pBAF_post)
        ) +
        geom_density()

        p2 = ggplot(
            D1,
            aes(x = snp_index, y = pBAF, color = state)
        ) +
        geom_point()

        p3 = ggplot(
            D2,
            aes(x = snp_index, y = pBAF, color = state)
        ) +
        geom_point()

        p4 = ggplot(
            D2,
            aes(x = pBAF_post)
        ) +
        geom_density()
        
        print(p1 | p2 | p3 | p4)
    }
    
    return(p)
    
}

compare_cluster = function(Obs_cut, min_snps = 30, plot = FALSE) {
    
    Obs_cut = Obs_cut %>% filter(!is.na(pBAF)) %>%
        group_by(CHROM) %>%
        arrange(POS) %>%
        mutate(snp_index = as.integer(factor(snp_id, unique(snp_id)))) %>%
        ungroup()
    
    # coerce cut label into factor
    Obs_cut = Obs_cut %>% mutate(cut = as.integer(as.factor(cut)))

    segs_union = Obs_cut %>% 
        group_by(CHROM, seg, cut, cnv_state) %>%
        summarise(
            n_markers = n(),
            seg_start = min(snp_index),
            seg_end = max(snp_index),
            .groups = 'drop'
        ) %>%
        do(
            GenomicRanges::disjoin(
                GenomicRanges::GRanges(
                seqnames = .$CHROM,
                ranges = IRanges::IRanges(start = .$seg_start, end = .$seg_end),
                cnv_state = .$cnv_state
            )) %>% as.data.frame
        ) %>%
        select(-strand) %>%
        rename(CHROM = seqnames)

    Obs_cut = Obs_cut %>% group_by(cut) %>%
        do(annotate_segs(., segs_union)) %>%
        ungroup() %>%
        arrange(CHROM) %>%
        mutate(seg_union = factor(seg_union, unique(seg_union))) 

    seg_diff = Obs_cut %>% 
        filter(!is.na(seg_union)) %>%
        group_by(CHROM, seg_union, seg_union_start, seg_union_end) %>%
        summarise(
            state_1 = unique(cnv_state[cut == 1]),
            state_2 = unique(cnv_state[cut == 2]),
            min_markers = min(sum(cut == 1), sum(cut == 2)),
            .groups = 'drop'
        ) %>%
        filter(
            ((state_1 == 'neu' & state_2 != 'neu') | (state_1 != 'neu' & state_2 == 'neu')) &
            min_markers >= 30
        )
    
    if (nrow(seg_diff) == 0) {
        return(seg_diff)
    }
    
#     IRdisplay::display(seg_diff)
    
    seg_diff = seg_diff %>% 
        rowwise() %>%
        mutate(
            p_diff = compare_event(CHROM, seg_union_start, seg_union_end, state_1, state_2, Obs_cut, plot)
        ) %>%
        ungroup()
    
    return(seg_diff)
}

analyze_cut = function(df, gexp.norm.long, membership, min_depth = 10, exp_delta = 0.1, k = 20, p_trans = 1e-4, exp_model = 'normal') {
    
    gexp.norm.long = gexp.norm.long %>% 
        select(-contains('cut')) %>%
        left_join(
            membership,
            by = "cell"
        )

    df = df %>% 
        select(-contains('cut')) %>%
        left_join(
            membership,
            by = "cell"
        )
        
    pseudobulk_cut = df %>%
        filter(!is.na(cut)) %>%
        filter(GT %in% c('1|0', '0|1')) %>%
        group_by(snp_id, cut, CHROM, POS, REF, ALT, GT, gene) %>%
        summarise(
            AD = sum(AD),
            DP = sum(DP),
            AR = AD/DP,
            .groups = 'drop'
        ) %>%
        arrange(CHROM, POS) %>%
        mutate(snp_index = as.integer(factor(snp_id, unique(snp_id)))) %>%
        ungroup()

    gexp_cut = gexp.norm.long %>% 
        filter(!is.na(cut)) %>%
        group_by(gene, CHROM, gene_start, gene_end, cut) %>% 
        summarise(
            exp = mean(exp),
            `.groups` = 'drop'
        ) %>%
        group_by(CHROM, cut) %>%
        arrange(gene) %>%
        mutate(exp_smooth = caTools::runmean(exp, k = k, align = 'center')) %>%
        arrange(CHROM) %>%
        mutate(CHROM = factor(CHROM, unique(CHROM)))

    Obs_cut = pseudobulk_cut %>% 
        filter(DP >= min_depth) %>%
        mutate(pBAF = ifelse(GT == '1|0', AR, 1-AR)) %>%
        mutate(pAD = ifelse(GT == '1|0', AD, DP - AD)) %>%
        mutate(CHROM = factor(CHROM, unique(CHROM))) %>%
        group_by(cut, CHROM) %>%
        arrange(cut, CHROM, POS) %>%
        mutate(inter_snp_dist = c(NA, POS[2:length(POS)] - POS[1:(length(POS)-1)])) %>%
        mutate(p_s = switch_prob(inter_snp_dist)) %>%
        ungroup() %>%
        select(-any_of('exp')) %>%
        full_join(
            gexp_cut,
            by = c("cut", "CHROM", "gene")
        ) %>%
        mutate(
            snp_id = ifelse(is.na(snp_id), gene, snp_id),
            # levels will be missing if not expressed
            gene = factor(gene, levels(gexp_cut$gene)),
            POS = ifelse(is.na(POS), gene_start, POS),
            # phase switch is forbidden if not heteroSNP
            p_s = ifelse(is.na(p_s), 0, p_s)
        ) %>%
        arrange(CHROM, POS) %>%
        group_by(CHROM) %>%
        mutate(snp_index = as.integer(factor(snp_id, unique(snp_id)))) %>%
        ungroup()
    
    # get rid of duplicate gene expression values
    Obs_cut = Obs_cut %>% 
        group_by(cut, CHROM, gene) %>%
        mutate(
            exp = ifelse(
                !is.na(gene) & n() > 1,
                c(unique(exp), rep(NA, n()-1)), exp
            ),
            exp_smooth = ifelse(
            !is.na(gene) & n() > 1,
            c(unique(exp_smooth), rep(NA, n()-1)), exp_smooth
        )) %>%
        ungroup() 

    if (exp_model == 'normal') {
        # gene expression model parameters
        sigma = 0.08
        mu_neu = 0.008
        mu_del = mu_neu - exp_delta
        mu_gain = mu_neu + exp_delta

        Obs_cut = Obs_cut %>% group_by(CHROM, cut) %>%
            mutate(state = run_hmm_mv_inhom(
                pAD, DP, exp, p_s,
                sigma = sigma,
                mu_neu = mu_neu,
                mu_del = mu_del,
                mu_gain = mu_gain,
                p_0 = 1 - p_trans
            ))
    } else if (exp_model == 'laplace') {
        Obs_cut = Obs_cut %>% group_by(CHROM, cut) %>% 
            mutate(state = run_hmm_mv_inhom_laplace(
                pAD, DP, exp, p_s,
                scale = 0.03,
                mu_neu = 0,
                mu_del = 0,
                mu_amp = 0,
                kappa_neu = 1,
                kappa_del = 1.43,
                kappa_amp = 0.7,
                p_0 = 1 - p_trans
            ))
    }
    
    
    Obs_cut = Obs_cut %>% 
        mutate(cnv_state = str_remove(state, '_down|_up')) %>%
        group_by(CHROM, cut) %>%
        arrange(snp_index) %>%
        mutate(boundary = c(0, cnv_state[2:length(cnv_state)] != cnv_state[1:(length(cnv_state)-1)])) %>%
        group_by(CHROM, cut) %>%
        mutate(seg = paste0(CHROM, '_', cumsum(boundary))) %>%
        arrange(CHROM) %>%
        mutate(seg = factor(seg, unique(seg))) %>%
        ungroup()
    
    return(Obs_cut)
    
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

analyze_node = function(df, gexp_bulk, min_depth = 8, p_trans = 1e-4, exp_only = FALSE) {
    
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
        group_by(CHROM) %>%
        arrange(CHROM, POS) %>%
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
        ungroup()
    
    Obs = Obs %>% 
        mutate(
            lambda_obs = Y_obs/d_obs,
            logFC = log2(lambda_obs/lambda_ref),
            lnFC = log(lambda_obs/lambda_ref)
        )
    
    if (exp_only) {
        Obs = Obs %>% mutate(pAD = NA, p_s = 0)
    }
    
    sigma_hat = Obs %>% pull(logFC) %>% sd(na.rm = TRUE)
    
    Obs = Obs %>% 
        mutate(sd = lambda_ref^(-1/4)/10) %>%
#         mutate(
#             precision_j = (n() * sigma_hat^(-2))*sqrt(lambda_ref)/sum(sqrt(lambda_ref)),
#             sd_j = precision_j ^ (-1/2)
#         ) %>%
        group_by(CHROM) %>%
        mutate(state = run_hmm_mv_inhom_norm(
            pAD = pAD,
            DP,
            p_s = p_s,
            logFC, 
            sd, 
            mu_neu = 0,
            mu_del = -0.3,
            mu_amp = 0.3,
            p_0 = 1-p_trans
        ))
    
    Obs = Obs %>% 
        mutate(cnv_state = str_remove(state, '_down|_up')) %>%
        group_by(CHROM) %>%
        arrange(snp_index) %>%
        mutate(boundary = c(0, cnv_state[2:length(cnv_state)] != cnv_state[1:(length(cnv_state)-1)])) %>%
        group_by(CHROM) %>%
        mutate(seg = paste0(CHROM, '_', cumsum(boundary))) %>%
        arrange(CHROM) %>%
        mutate(seg = factor(seg, unique(seg))) %>%
        ungroup()
    
    
    return(Obs)
    
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

mu_mle = function(x_vec, lambda_vec) {
    weighted.mean(x = x_vec, w = sqrt(lambda_vec), na.rm = T)
}

mu_mle_roll = function(x_vec, lambda_vec, h) {
    n = length(x_vec)
    sapply(
        1:n,
        function(c) {
            slice = max(c - h - 1, 1):min(c + h, n)
            mu_mle(x_vec[slice], lambda_vec[slice])
        }
    )
}

process_exp_fc2 = function(count_mat, depths, cell_annot, gtf_transcript) {
    
    ref_cells = cell_annot %>% filter(group == 'ref') %>% pull(cell)
    obs_cells = cell_annot %>% filter(group == 'obs') %>% pull(cell)
        
    names(depths) = str_remove(names(depths), paste0(sample, '_'))
    colnames(count_mat) = str_remove(colnames(count_mat), paste0(sample, '_'))
    
    count_mat = count_mat[,colnames(count_mat) %in% cell_annot$cell]
    depths = depths[names(depths) %in% cell_annot$cell]
    
    genes_annotated = intersect(gtf_transcript$gene, rownames(count_mat))
    count_mat = count_mat[genes_annotated,]
    count_mat = count_mat[rowSums(count_mat) > 0,]

    depth_ref = sum(depths[ref_cells])
    depth_obs = sum(depths[obs_cells])

    # lambda matrices
#     lambda_mat = sapply(
#         1:length(depths),
#         function (i) {
#             count_mat[,i]/depths[i]
#         }
#     )

#     colnames(lambda_mat) = colnames(count_mat)

#     lambda_mat_ref = lambda_mat[,colnames(lambda_mat) %in% ref_cells]
#     lambda_mat_obs = lambda_mat[,colnames(lambda_mat) %in% obs_cells]

    # count matrices
    count_mat_ref = count_mat[,colnames(count_mat) %in% ref_cells]
    count_mat_obs = count_mat[,colnames(count_mat) %in% obs_cells, drop = FALSE]
    
    lambdas_ref = rowSums(count_mat_ref)/depth_ref
    lambdas_obs = rowSums(count_mat_obs)/depth_obs
    
    # filter for mutually expressed genes
#     min_both = 0
    
#     mut_expressed = (rowSums(count_mat_ref)/depth_ref * 1e6 > min_both & rowSums(count_mat_obs)/depth_obs * 1e6 > min_both) |
#         (rowMeans(lambda_mat_ref) > mean(lambda_mat_ref[lambda_mat_ref != 0])) |
#         (rowMeans(lambda_mat_obs) > mean(lambda_mat_obs[lambda_mat_obs != 0]))
    
#     mut_expressed = rowSums(count_mat) > 0

#     count_mat_ref = count_mat_ref[mut_expressed,]
#     count_mat_obs = count_mat_obs[mut_expressed,]

    display(nrow(count_mat))
    
    bulk_ref = count_mat_ref %>%
        rowSums() %>%
        data.frame() %>%
        setNames('Y_ref') %>%
        tibble::rownames_to_column('gene') %>%
        mutate(lambda_ref = (Y_ref/depth_ref)) %>%
        left_join(gtf_transcript, by = "gene") %>%
        mutate(d_obs = depth_obs)

    bulk_obs = count_mat_obs %>%
        rowSums() %>%
        data.frame() %>%
        setNames('Y_obs') %>%
        tibble::rownames_to_column('gene') %>%
        mutate(lambda_obs = (Y_obs/depth_obs)) %>%
        left_join(gtf_transcript, by = "gene") %>%
        mutate(d_ref = depth_ref)
    
    bulk_all = bulk_obs %>% left_join(bulk_ref)

    # annotate using GTF
    bulk_all = bulk_all %>%
        mutate(gene = droplevels(factor(gene, gtf_transcript$gene))) %>%
        mutate(gene_index = as.integer(gene)) %>%
        arrange(gene) %>%
        mutate(CHROM = factor(CHROM)) %>%
        group_by(CHROM) %>%
        mutate(
            logFC = log2(lambda_obs) - log2(lambda_ref),
            logFC = ifelse(is.infinite(logFC), NA, logFC),
            logFC_pad = log2(lambda_obs * 1e6 + 1) - log2(lambda_ref * 1e6 + 1),
            mu_mle_roll = mu_mle_roll(logFC, lambda_ref, h = 50),
            phi_hat_roll = phi_hat_roll(Y_obs, lambda_ref, depth_obs, h = 50)
        )%>%
        ungroup()
    
    sigma_hat = sd(bulk_all$logFC)
    
    bulk_all = bulk_all %>%
        mutate(sd = lambda_ref^(-1/4)/10) %>%
        mutate(
            precision_j = (n() * sigma_hat^(-2))*sqrt(lambda_ref)/sum(sqrt(lambda_ref)),
            sd_j = precision_j ^ (-1/2)
        )
    
    return(bulk_all)
}


process_exp_fc = function(count_mat, depths, cell_annot, gtf_transcript) {
    
    ref_cells = cell_annot %>% filter(group == 'ref') %>% pull(cell)
    obs_cells = cell_annot %>% filter(group == 'obs') %>% pull(cell)
        
    names(depths) = str_remove(names(depths), paste0(sample, '_'))
    colnames(count_mat) = str_remove(colnames(count_mat), paste0(sample, '_'))
    
    count_mat = count_mat[,colnames(count_mat) %in% cell_annot$cell]
    depths = depths[names(depths) %in% cell_annot$cell]
    
    genes_annotated = intersect(gtf_transcript$gene, rownames(count_mat))
    count_mat = count_mat[genes_annotated,]

    depth_ref = sum(depths[ref_cells])
    depth_obs = sum(depths[obs_cells])

    # lambda matrices
#     lambda_mat = sapply(
#         1:length(depths),
#         function (i) {
#             count_mat[,i]/depths[i]
#         }
#     )

#     colnames(lambda_mat) = colnames(count_mat)

#     lambda_mat_ref = lambda_mat[,colnames(lambda_mat) %in% ref_cells]
#     lambda_mat_obs = lambda_mat[,colnames(lambda_mat) %in% obs_cells]

    # count matrices
    count_mat_ref = count_mat[,colnames(count_mat) %in% ref_cells]
    count_mat_obs = count_mat[,colnames(count_mat) %in% obs_cells]

    lambdas_ref = rowSums(count_mat_ref)/depth_ref
    lambdas_obs = rowSums(count_mat_obs)/depth_obs

    # filter for mutually expressed genes
    min_both = 2
    
    mut_expressed = (rowSums(count_mat_ref)/depth_ref * 1e6 > min_both & rowSums(count_mat_obs)/depth_obs * 1e6 > min_both) |
        (lambdas_ref > mean(lambdas_ref[lambdas_ref != 0])) |
        (lambdas_obs > mean(lambdas_obs[lambdas_obs != 0]))
    
    count_mat_ref = count_mat_ref[mut_expressed,]
    count_mat_obs = count_mat_obs[mut_expressed,]

    display(nrow(count_mat_obs))
    
    bulk_ref = count_mat_ref %>%
        rowSums() %>%
        data.frame() %>%
        setNames('Y_ref') %>%
        tibble::rownames_to_column('gene') %>%
        mutate(lambda_ref = (Y_ref/depth_ref)) %>%
        left_join(gtf_transcript, by = "gene") %>%
        mutate(d_obs = depth_obs)

    bulk_obs = count_mat_obs %>%
        rowSums() %>%
        data.frame() %>%
        setNames('Y_obs') %>%
        tibble::rownames_to_column('gene') %>%
        mutate(lambda_obs = (Y_obs/depth_obs)) %>%
        left_join(gtf_transcript, by = "gene") %>%
        mutate(d_ref = depth_ref)

    bulk_all = bulk_obs %>% left_join(bulk_ref)

    # annotate using GTF
    bulk_all = bulk_all %>%
        mutate(gene = droplevels(factor(gene, gtf_transcript$gene))) %>%
        mutate(gene_index = as.integer(gene)) %>%
        arrange(gene) %>%
        mutate(CHROM = factor(CHROM)) %>%
        mutate(
            logFC = log2(lambda_obs) - log2(lambda_ref),
            logFC = ifelse(is.infinite(logFC), NA, logFC),
            logFC_pad = log2(lambda_obs * 1e6 + 1) - log2(lambda_ref * 1e6 + 1)
        ) %>%
        group_by(CHROM) %>%
#         filter(logFC < 5 & logFC > -5) %>%
        mutate(
            mu_mle_roll = mu_mle_roll(logFC, lambda_ref, h = 50),
            phi_hat_roll = phi_hat_roll(Y_obs, lambda_ref, depth_obs, h = 50)
        ) %>%
        ungroup()
    
    return(bulk_all)
}


cnv_colors = scale_color_manual(
    values = c("neu" = "gray", "del_up" = "royalblue", "del_down" = "darkblue", 
               "loh_up" = "darkgreen", "loh_down" = "olivedrab4",
               "amp_up" = "red", "amp_down" = "tomato3",
              "amp" = "darkred", "loh" = "darkgreen", "del" = "darkblue", "neu2" = "gray30")
)

plot_hmm = function(Obs_cut, dot_size = 0.8, exp_limit = 1, rollmean = FALSE, k = 51) {
    
    segs_cut = Obs_cut %>% group_by(CHROM, cut, seg, cnv_state) %>%
        filter(state != 'neu') %>%
        mutate(haplo = case_when(
            str_detect(state, 'up') ~ 'major',
            str_detect(state, 'down') ~ 'minor',
#             state == 'neu' & pBAF >= 0.5 ~ 'major',
#             T ~ 'minor'
        )) %>%
        summarise(
            exp = mean(exp, na.rm = TRUE),
            AF_major = median(pBAF[haplo == 'major']),
            AF_minor = median(pBAF[haplo == 'minor']),
            seg_start = min(snp_index),
            seg_end = max(snp_index),
            .groups = 'drop'
        ) %>%
        mutate(exp = ifelse(exp > exp_limit | exp < -exp_limit, NA, exp))
    
    exp_roll = Obs_cut %>% group_by(CHROM, cut) %>%
        mutate(
            exp = caTools::runmean(exp, k = k, align = 'center')) %>%
        mutate(exp = ifelse(exp > exp_limit | exp < -exp_limit, NA, exp)) %>%
        mutate(cnv_state = ifelse(cnv_state == 'neu', 'neu2', cnv_state))

    p = ggplot(
        Obs_cut %>% 
            mutate(exp = ifelse(exp > exp_limit | exp < -exp_limit, NA, exp)) %>%
            reshape2::melt(measure.vars = c('exp', 'pBAF')),
        aes(x = snp_index, y = value, color = state)
    ) +
    geom_point(size = dot_size, alpha = 0.5, pch = 16) +
#     geom_segment(
#         inherit.aes = FALSE,
#         data = segs_cut %>% mutate(variable = 'pBAF'),
#         aes(y = AF_major, yend = AF_major, x = seg_start, xend = seg_end, color = cnv_state),
#         lineend = 'round',
#         size = 0.8
#     ) +
#     geom_segment(
#         inherit.aes = FALSE,
#         data = segs_cut %>% mutate(variable = 'pBAF'),
#         aes(y = AF_minor, yend = AF_minor, x = seg_start, xend = seg_end, color = cnv_state),
#         lineend = 'round',
#         size = 0.8
#     ) +
#     geom_segment(
#         inherit.aes = FALSE,
#         data = segs_cut %>% mutate(variable = 'exp'),
#         aes(y = exp, yend = exp, x = seg_start, xend = seg_end, color = cnv_state),
#         lineend = 'round',
#         size = 0.8
#     ) +
    geom_line(
        inherit.aes = FALSE,
        data = exp_roll %>% mutate(variable = 'exp'),
        aes(x = snp_index, y = exp, color = cnv_state, group = '1')
    ) +
    theme_classic() +
    theme(
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA)
    ) +
    facet_grid(cut + variable ~ CHROM, scale = 'free', space = 'free_x') +
    cnv_colors
    
    return(p)
}


plot_psbulk = function(Obs, dot_size = 0.8, exp_limit = 2, min_depth = 10) {
    
    p = ggplot(
        Obs %>% 
            mutate(logFC = ifelse(logFC > exp_limit | logFC < -exp_limit, NA, logFC)) %>%
            mutate(pBAF = ifelse(DP >= min_depth, pBAF, NA)) %>%
            reshape2::melt(measure.vars = c('logFC', 'pBAF')),
        aes(x = snp_index, y = value, color = state)
    ) +
    geom_point(size = dot_size, alpha = 0.5, pch = 16) +
    geom_hline(data = data.frame(variable = 'logFC'), aes(yintercept = 0), color = 'gray30', linetype = 'dashed') +
    geom_line(
        inherit.aes = FALSE,
        data = Obs %>% mutate(variable = 'logFC'),
        aes(x = snp_index, y = mu_mle_roll, group = '1'),
        color = 'darkred'
    ) +
    theme_classic() +
    theme(
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA)
    ) +
    facet_grid(variable ~ CHROM, scale = 'free', space = 'free_x') +
    cnv_colors
    
    return(p)
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

process_exp = function(exp_mat, cell_annot, controls, gtf_transcript, minMeanBoth = 0, window = 51) {
    
    # Normal reference
    ref_cells = cell_annot %>% filter(group == 'ref') %>% pull(cell)
    obs_cells = cell_annot %>% filter(group == 'obs') %>% pull(cell)

    ref_mat = exp_mat[rownames(exp_mat) %in% ref_cells,]
    obs_mat = exp_mat[rownames(exp_mat) %in% obs_cells,]

    obs_mat = t(obs_mat)
    ref_mat = t(ref_mat)

    gexp.ref.init = ref_mat
    gexp.sc.init = obs_mat

    # insersect gene list
    vi <- intersect(rownames(gexp.sc.init), rownames(gexp.ref.init))
    
    gexp.sc <- gexp.sc.init[vi,]
    gexp.ref <- gexp.ref.init[vi,,drop=FALSE]

    # compute thresholds
    minMeanTest=mean(gexp.sc.init[gexp.sc.init!=0])
    minMeanRef=mean(gexp.ref.init[gexp.ref.init!=0])

    # Filter for highly expressed genes
    vi <- (rowMeans(gexp.sc) > minMeanBoth & rowMeans(gexp.ref) > minMeanBoth) |
        rowMeans(gexp.sc) > minMeanTest | rowMeans(gexp.ref) > minMeanRef

    gexp.sc <- gexp.sc[vi,]
    gexp.ref <- gexp.ref[vi,,drop=FALSE]

    cat(paste0(sum(vi), " genes passed filtering ... \n"))

    # log-transform
    gexp.sc = log2(gexp.sc + 1)
    gexp.ref = log2(gexp.ref + 1)

    # scaling for library size
#     gexp.sc <- scale(gexp.sc)
#     gexp.ref <- scale(gexp.ref)

    # substract reference background
    refmean <- rowMeans(gexp.ref)
    gexp.norm <- gexp.sc - refmean

    gexp.norm = gexp.norm %>% as.data.frame() %>%
        tibble::rownames_to_column('gene') %>%
#         setNames(str_remove(colnames(.), paste0(sample, '_'))) %>%
        inner_join(gtf_transcript, by = "gene") %>%
        mutate(gene = droplevels(factor(gene, gtf_transcript$gene))) %>%
        mutate(gene_index = as.integer(gene)) %>% 
        mutate(CHROM = as.integer(CHROM)) %>%
        arrange(CHROM, gene)

    gexp.norm.long = gexp.norm %>% 
        reshape2::melt(
            id.var = c('gene', 'gene_index', 'region', 'gene_start', 'gene_end', 'CHROM', 'length'),
            variable.name = 'cell',
            value.name = 'exp')    

    gexp.norm.long = gexp.norm.long %>% 
        group_by(cell, CHROM) %>%
#         arrange(gene) %>%
        mutate(exp_rollmean = caTools::runmean(exp, k = window, align = "center")) %>%
        ungroup()

    gexp.norm.long = gexp.norm.long %>% 
        left_join(
            cell_annot,
            by = 'cell'
        )
    
    return(list('gexp.norm' = gexp.norm, 'gexp.norm.long' = gexp.norm.long, 'gexp.ref' = gexp.ref))
}

calc_cluster_tree = function(exp_mat, cell_annot) {
    
    cell_dict = cell_annot %>% {setNames(.$cell_type, .$cell)}

    exp_mat = exp_mat %>% data.frame() %>%
        tibble::rownames_to_column('cell') %>%
        mutate(cell_type = cell_dict[cell]) %>%
        filter(!is.na(cell_type)) %>%
        select(-cell)

    exp_mat_clust = aggregate(select(exp_mat, -cell_type), list('cell_type' = exp_mat$cell_type), mean) %>%
        tibble::column_to_rownames('cell_type')
    
    hc_exp_clust = hclust(as.dist(1-cor(t(exp_mat_clust))), method = "ward.D2")
    
    return(hc_exp_clust)
    
}

analyze_node_iter = function(Obs, p_0, prior = NULL) {
    
    Obs = Obs %>% 
        mutate(
            lnFC = log(lambda_obs/lambda_ref),
            state = 'neu'
        )
    
    converge = FALSE
    
    i = 0
    
    while (!converge) {
        
        i = i + 1
        
        mu_hat = mean(Obs$lnFC[Obs$state == 'neu'], na.rm = TRUE)
        sigma_hat = sd(Obs$lnFC[Obs$state == 'neu'], na.rm = TRUE)
        
        display(glue('iteration {i}, mu_hat {signif(mu_hat, 3)}, sigma_hat {signif(sigma_hat, 3)}'))
        
        state_old = Obs$state
        
        Obs = Obs %>% 
            group_by(CHROM) %>%
            mutate(
                state = run_hmm_mv_inhom_poilog(
                    pAD, DP, p_s, Y_obs, lambda_ref, d_obs,
                    phi_neu = exp(mu_hat),
                    phi_del = exp(qnorm(0.4, mean = mu_hat, sd = sigma_hat)),
                    phi_amp = exp(qnorm(0.6, mean = mu_hat, sd = sigma_hat)),
                    sigma = sigma_hat,
                    p_0, 
                    prior)
        )
        
        converge = all(state_old == Obs$state)
        
    }
    
    return(Obs)
}
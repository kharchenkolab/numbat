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

compare_cluster = function(Obs_cut, plot = FALSE) {

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

analyze_cut = function(df, gexp.norm.long, membership, min_depth = 10) {
    
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
        group_by(snp_id, cut) %>%
        summarise(
            AD = sum(AD),
            DP = sum(DP),
            AR = AD/DP,
            .groups = 'drop'
        ) %>%
        left_join(
            vcf_phased %>% select(CHROM, POS, REF, ALT, snp_id, GT, gene, gene_start, gene_end),
            by = "snp_id"
        ) %>%
        arrange(CHROM, POS) %>%
        mutate(snp_index = as.integer(factor(snp_id, unique(snp_id)))) %>%
        ungroup()

    gexp_cut = gexp.norm.long %>% 
        filter(!is.na(cut)) %>%
        group_by(gene, CHROM, start, end, cut) %>% 
        summarise(
            exp = mean(exp),
            `.groups` = 'drop'
        ) %>%
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
            gene = factor(gene, levels(gexp_all$gene)),
            POS = ifelse(is.na(POS), start, POS),
            # phase switch is forbidden if not heteroSNP
            p_s = ifelse(is.na(p_s), 0, p_s)
        ) %>%
        arrange(CHROM, POS) %>%
        group_by(CHROM) %>%
        mutate(snp_index = 1:n()) %>%
        ungroup()
    
    # get rid of duplicate gene expression values
    Obs_cut = Obs_cut %>% 
        group_by(cut, CHROM, gene) %>%
        mutate(exp = ifelse(
            !is.na(gene) & n() > 1,
            c(unique(exp), rep(NA, n()-1)), exp
        )) %>%
        ungroup() 

    # gene expression model parameters
    sigma = 1
    mu_neu = 0
    mu_del = mu_neu - 0.1
    mu_gain = mu_neu + 0.1

    Obs_cut = Obs_cut %>% group_by(CHROM, cut) %>%
        mutate(state = run_hmm_mv_inhom(
            pAD, DP, exp, p_s,
            sigma = sigma,
            mu_neu = mu_neu,
            mu_del = mu_del,
            mu_gain = mu_gain,
            p_0 = 1-1e-4
        ))
    
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


cnv_colors = scale_color_manual(
    values = c("neu" = "gray", "del_up" = "royalblue", "del_down" = "darkblue", 
               "loh_up" = "darkgreen", "loh_down" = "olivedrab4",
               "amp_up" = "red", "amp_down" = "tomato3",
              "amp" = "darkred", "loh" = "darkgreen", "del" = "darkblue")
)

plot_hmm = function(Obs_cut) {
    
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
        mutate(exp = ifelse(exp > 1 | exp < -1, NA, exp))

    p = ggplot(
        Obs_cut %>% 
            mutate(exp = ifelse(exp > 1 | exp < -1, NA, exp)) %>%
            reshape2::melt(measure.vars = c('exp', 'pBAF')),
        aes(x = snp_index, y = value, color = state)
    ) +
    geom_point(size = 0.5, alpha = 0.5, pch = 16) +
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
    theme_classic() +
    theme(
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA)
    ) +
    facet_grid(cut + variable ~ CHROM, scale = 'free', space = 'free_x') +
    cnv_colors
    
    return(p)
}


walk_tree = function(den, cut_prefix, cuts, segs_tested, gexp.norm.long, df, p_cutoff = 0.01) {
    
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
        Obs_cut = analyze_cut(df, gexp.norm.long, membership)
        
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
            p_cutoff)
        
        res_2 = walk_tree(
            sub_dens[[2]],
            paste0(cut_prefix, '.', 2),
            cuts,
            segs_tested, 
            gexp.norm.long,
            df,
            p_cutoff)
        
        res = list(
            'cuts' = rbind(res_1$cuts, res_2$cuts), 
            'segs_tested' = rbind(res_1$segs_tested, res_2$segs_tested))
        
        return(res)
    }
}
    
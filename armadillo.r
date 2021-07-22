source('~/Armadillo/hmm.r')
source('~/Armadillo/utils.r')
########################### Main ############################

#' @param exp_mat scaled expression matrix
armadillo = function(count_mat, lambdas_ref, df, cell_annot, ncores, t = 1e-5, verbose = TRUE) {
    
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
                analyze_bulk_gpois(t = t, verbose = TRUE)

                return(bulk)
        }) %>%
        bind_rows() %>%
        arrange(CHROM, POS) %>%
        mutate(snp_id = factor(snp_id, unique(snp_id))) %>%
        mutate(snp_index = as.integer(snp_id)) %>%
        arrange(node)
    
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
    
    #### 3. Find consensus CNVs ####
    
    if (verbose) {
        display('Finding consensus CNVs ..')
    }
    
    segs_all = bulk_all %>% 
        filter(state != 'neu') %>%
        distinct(node, CHROM, seg, cnv_state, cnv_state_post, seg_start, seg_end,
                 theta_mle, phi_mle, p_loh, p_del, p_amp, p_bamp, p_bdel, LLR, LLR_ar)
    
    segs_filtered = segs_all %>% filter(!(LLR_ar < 10 & cnv_state %in% c('del', 'amp') & (phi_mle > 2 | phi_mle < 0.5)))
        
    segs_consensus = resolve_cnvs(segs_filtered)
    
    #### 4. Per-cell CNV evaluations ####
    
    if (verbose) {
        display('Calculating per-cell CNV posteriors ..')
    }
    
    exp_post = get_exp_post(
        segs_consensus %>%
            mutate(cnv_state = ifelse(cnv_state == 'neu', cnv_state, cnv_state_post)),
        count_mat,
        lambdas_ref,
        ncores = ncores
    )
    
    allele_post = get_allele_post(bulk_all, segs_consensus, df)
    
    joint_post = get_joint_post(exp_post, allele_post, segs_consensus)
    
    if (verbose) {
        display('All done!')
    }
    
    res = list(
        'tree' = tree,
        'nodes' = nodes,
        'bulk_all' = bulk_all,
        'plot_list' = plot_list,
        'segs_all' = segs_all,
        'segs_filtered' = segs_filtered,
        'graph' = G,
        'segs_consensus' = segs_consensus,
        'exp_post' = exp_post,
        'allele_post' = allele_post,
        'joint_post' = joint_post
    )
    
    return(res)
}

resolve_cnvs = function(segs_all) {
            
    # resolve overlapping calls by graph
    V = segs_all %>% mutate(vertex = 1:n(), .before = 'CHROM')

    E = segs_all %>% {GenomicRanges::GRanges(
            seqnames = .$CHROM,
            IRanges::IRanges(start = .$seg_start,
                   end = .$seg_end)
        )} %>%
        GenomicRanges::findOverlaps(., .) %>%
        as.data.frame %>%
        setNames(c('from', 'to')) %>% 
        filter(from != to)

    G = igraph::graph_from_data_frame(d=E, vertices=V, directed=F)

    segs_all = segs_all %>% mutate(group = igraph::components(G)$membership)

    segs_consensus = segs_all %>% arrange(CHROM, group, -LLR) %>% distinct(group, `.keep_all` = TRUE) 
    
    segs_consensus = fill_neu_segs(segs_consensus)
    
    return(segs_consensus)
}



# retrieve neutral segments
fill_neu_segs = function(segs_consensus) {
    
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
        mutate(CHROM = as.factor(CHROM))
    
    return(segs_consensus)
}

get_exp_post = function(segs_consensus, count_mat, lambdas_ref, ncores = 30) {

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
        mutate(CHROM = as.factor(CHROM)) %>%
        left_join(
            segs_consensus %>% mutate(seg_index = 1:n()),
            by = c('seg_index', 'CHROM')
        ) %>%
        distinct(gene, `.keep_all` = TRUE) 

    exp_sc = count_mat %>%
        as.data.frame() %>%
        tibble::rownames_to_column('gene') %>% 
        mutate(lambda_ref = lambdas_ref[gene]) %>%
        filter(lambda_ref > 0) %>%
        left_join(
            gene_seg %>% select(gene, seg = seg_cons, cnv_state),
            by = "gene"
        )
    
    cells = colnames(count_mat)
    
    exp_post = mclapply(
        cells,
        mc.cores = ncores,
        function(cell) {
            exp_sc[,c('gene', 'lambda_ref', 'seg', 'cnv_state', cell)] %>%
            set_names(c('gene', 'lambda_ref', 'seg', 'cnv_state', 'Y_obs')) %>%
            get_exp_likelihoods() %>%
            mutate(cell = cell)
        }
    ) %>% 
    bind_rows() %>%
    mutate(seg = factor(seg, gtools::mixedsort(unique(seg))))
    
    exp_post = exp_post %>% 
        rowwise() %>%
        left_join(
            segs_consensus %>% select(seg = seg_cons, p_loh, p_amp, p_del, p_bamp, p_bdel),
            by = 'seg'
        ) %>%
        rowwise() %>%
        mutate(
            Z = matrixStats::logSumExp(
                c(l11 + log(1/2),
                  l20 + log(p_loh/2),
                  l10 + log(p_del/2),
                  l21 + log(p_amp/4),
                  l31 + log(p_amp/4),
                  l22 + log(p_bamp/2),
                  l00 + log(p_bdel/2))
            ),
            p_amp = exp(matrixStats::logSumExp(c(l21 + log(p_amp/4), l31 + log(p_amp/4))) - Z),
            p_amp_sin = exp(l21 + log(p_amp/4) - Z),
            p_amp_mul = exp(l31 + log(p_amp/4) - Z),
            p_neu = exp(l11 + log(1/2) - Z),
            p_del = exp(l10 + log(p_del/2) - Z),
            p_loh = exp(l20 + log(p_loh/2) - Z),
            p_bamp = exp(l22 + log(p_bamp/2) - Z),
            p_bdel = exp(l22 + log(p_bdel/2) - Z),
            p_cnv = p_amp + p_del + p_loh + p_bamp + p_bdel
        ) %>%
        ungroup()
    
    return(exp_post)
}

get_allele_post = function(bulk_all, segs_consensus, df) {
    
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
            segs_consensus %>% mutate(cnv_state = cnv_state_post),
            by = c('node', 'seg')
        )

    allele_sc = df %>%
        mutate(pAD = ifelse(GT == '1|0', AD, DP - AD)) %>%
        select(-snp_index) %>% 
        inner_join(
            snp_seg %>% select(snp_id, snp_index, haplo, seg = seg_cons, cnv_state),
            by = c('snp_id')
        ) %>%
        filter(cnv_state %in% c('amp', 'del', 'loh')) %>%
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
        left_join(
            segs_consensus %>% select(seg = seg_cons, p_loh, p_amp, p_del, p_bamp, p_bdel),
            by = 'seg'
        ) %>%
        rowwise() %>%
        mutate(
            l11 = dbinom(major, total, p = 0.5, log = TRUE),
            l10 = dbinom(major, total, p = 0.9, log = TRUE),
            l20 = dbinom(major, total, p = 0.9, log = TRUE),
            l21 = dbinom(major, total, p = 0.66, log = TRUE),
            l31 = dbinom(major, total, p = 0.75, log = TRUE),
            l22 = l11,
            l00 = l11,
            Z = matrixStats::logSumExp(
                c(l11 + log(1/2),
                  l20 + log(p_loh/2),
                  l10 + log(p_del/2),
                  l21 + log(p_amp/4),
                  l31 + log(p_amp/4),
                  l22 + log(p_bamp/2),
                  l00 + log(p_bdel/2)
                 )
            ),
            p_amp = exp(matrixStats::logSumExp(c(l21 + log(p_amp/4), l31 + log(p_amp/4))) - Z),
            p_amp_sin = exp(l21 + log(p_amp/4) - Z),
            p_amp_mul = exp(l31 + log(p_amp/4) - Z),
            p_neu = exp(l11 + log(1/2) - Z),
            p_del = exp(l10 + log(p_del/2) - Z),
            p_loh = exp(l20 + log(p_loh/2) - Z),
            p_bamp = exp(l22 + log(p_bamp/2) - Z),
            p_bdel = exp(l22 + log(p_bdel/2) - Z),
            p_cnv = p_amp + p_del + p_loh + p_bamp + p_bdel
        ) %>%
        ungroup()
}

get_joint_post = function(exp_post, allele_post, segs_consensus) {
    joint_post = exp_post %>%
        select(
            cell, seg, cnv_state, l11_x = l11, l20_x = l20, l10_x = l10, l21_x = l21, l31_x = l31, l22_x = l22, l00_x = l00,
            Z_x = Z
        ) %>%
        full_join(
            allele_post %>% select(
                cell, seg, l11_y = l11, l20_y = l20, l10_y = l10, l21_y = l21, l31_y = l31, l22_y = l22, l00_y = l00,
                n_snp = total,
                Z_y = Z
            ),
            c("cell", "seg")
        ) %>%
        mutate(cnv_state = tidyr::replace_na(cnv_state, 'loh')) %>%
        mutate_at(
                vars(matches("_x|_y")),
                function(x) tidyr::replace_na(x, 0)
            ) %>%
        left_join(
            segs_consensus %>% select(seg = seg_cons, CHROM, p_loh, p_amp, p_del, p_bamp, p_bdel),
            by = 'seg'
        ) %>%
        rowwise() %>%
        mutate(
            Z = matrixStats::logSumExp(
                c(l11_x + l11_y + log(1/2),
                  l20_x + l20_y + log(p_loh/2),
                  l10_x + l10_y + log(p_del/2),
                  l21_x + l21_y + log(p_amp/4),
                  l31_x + l31_y + log(p_amp/4),
                  l22_x + l22_y + log(p_bamp/2),
                  l00_x + l00_y + log(p_bdel/2))
            ),
            Z_cnv_x = matrixStats::logSumExp(
                c(l20_x + log(p_loh/2),
                  l10_x + log(p_del/2),
                  l21_x + log(p_amp/4),
                  l31_x + log(p_amp/4),
                  l22_x + log(p_bamp/2),
                  l00_x + log(p_bdel/2))
            ),
            Z_cnv_y = matrixStats::logSumExp(
                c(l20_y + log(p_loh/2),
                  l10_y + log(p_del/2),
                  l21_y + log(p_amp/4),
                  l31_y + log(p_amp/4),
                  l22_y + log(p_bamp/2),
                  l00_y + log(p_bdel/2))
            ),
            Z_cnv = matrixStats::logSumExp(
                c(l20_x + l20_y + log(p_loh/2),
                  l10_x + l10_y + log(p_del/2),
                  l21_x + l21_y + log(p_amp/4),
                  l31_x + l31_y + log(p_amp/4),
                  l22_x + l22_y + log(p_bamp/2),
                  l00_x + l00_y + log(p_bdel/2))
            ),
            Z_n = l11_x + l11_y + log(1/2),
            Z_n_x = l11_x + log(1/2),
            Z_n_y = l11_y + log(1/2),
            p_cnv = exp(Z_cnv - Z),
            p_cnv_x = exp(Z_cnv_x - Z_x),
            p_cnv_y = exp(Z_cnv_y - Z_y)
        ) %>%
        ungroup()

    joint_post = joint_post %>% 
        mutate(seg = factor(seg, gtools::mixedsort(unique(seg)))) %>%
        mutate(seg_label = paste0(seg, '(', cnv_state, ')')) %>%
        mutate(seg_label = factor(seg_label, unique(seg_label)))
    
    return(joint_post)
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
        
#     if (is.null(mut_clust)) {
#         mut_clust = joint_post %>% 
#             pull(seg) %>% unique
#     }
    
    cnv_clusters = dendextend::cutree(tree_cnv, 2) 
#     cnv_clusters = lapply(unique(cnv_clusters), function(c) {names(cnv_clusters[cnv_clusters == c])})
    
#     joint_post_cell = joint_post %>% 
#         group_by(cell, cell_group) %>%
#         summarise_at(
#             vars(contains('Z_')),
#             sum
#         ) %>%
#         mutate(
#             Z = matrixStats::logSumExp(c(Z_cnv + log(0.5), Z_n + log(0.5))),
#             Z_x = matrixStats::logSumExp(c(Z_cnv_x + log(0.5), Z_n_x + log(0.5))),
#             Z_y = matrixStats::logSumExp(c(Z_cnv_y + log(0.5), Z_n_y + log(0.5))),
#             p_cnv = exp(Z_cnv + log(0.5) - Z),
#             p_cnv_x = exp(Z_cnv_x + log(0.5) - Z_x),
#             p_cnv_y = exp(Z_cnv_y + log(0.5) - Z_y)
#         ) %>%
#         ungroup() %>%
#         mutate(panel = 'all')
        
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
            Z = matrixStats::logSumExp(c(Z_cnv + log(0.5), Z_n + log(0.5))),
            Z_x = matrixStats::logSumExp(c(Z_cnv_x + log(0.5), Z_n_x + log(0.5))),
            Z_y = matrixStats::logSumExp(c(Z_cnv_y + log(0.5), Z_n_y + log(0.5))),
            p_cnv = exp(Z_cnv + log(0.5) - Z),
            p_cnv_x = exp(Z_cnv_x + log(0.5) - Z_x),
            p_cnv_y = exp(Z_cnv_y + log(0.5) - Z_y)
        ) %>%
        ungroup()
    
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
#         mutate(cluster = factor(cluster, c('all', '1', '2', 'segment'))) %>%
        mutate(seg = factor(seg, cnv_order)) %>%
        arrange(seg) %>%
        mutate(seg_label = factor(seg_label, unique(seg_label))) %>%
        mutate(cell = factor(cell, cell_order$cell)) %>%
        group_by(cell_group) %>%
        filter(cell %in% sample(unique(cell), min(n_sample, length(unique(cell))))) %>%
        ungroup()
    
    p = ggplot(
            D,
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
    
    p_cnv_tree = tree_cnv %>%
        ggtree::ggtree(ladderize=F) +
        scale_x_discrete(expand = expansion(mult = 0)) +
        theme(plot.margin = margin(0,0,0,0)) 
    
    panel_1 = p
    
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
#         filter(!(cnv_state %in% c('bamp') & source %in% c('allele', 'joint'))) %>%
#         filter(!(cnv_state %in% c('loh') & source %in% c('expr', 'joint'))) %>%
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
    
    return(list(panel_1, panel_2))
    
}
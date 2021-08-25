source('/home/tenggao/Armadillo/hmm.r')
source('/home/tenggao/Armadillo/utils.r')
source('/home/tenggao/Armadillo/graphs.r')


#' @param exp_mat scaled expression matrix
armadillo = function(count_mat, lambdas_ref, df, cell_annot, ncores, t = 1e-5, gbuild = 'hg38', verbose = TRUE) {
    
    res = list()
    #### 1. Build expression tree ####
    
    if (verbose) {
        display('Building expression tree ..')
    }
    
    dist_mat = calc_cluster_dist(
        count_mat, 
        cell_annot %>% filter(group == 'obs')
    )

    tree = hclust(as.dist(dist_mat), method = "ward.D2")
    
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
    
    nodes = nodes %>% extract(purrr::map(., 'size') > 300)
    
    res[['tree']] = tree
    res[['nodes']] = nodes
    
    #### 2. Run pseudobulk HMM ####
    
    if (verbose) {
        display('Running HMMs ..')
    }

    bulk_all = mclapply(
            nodes,
            mc.cores = ncores,
            function(node) {
                get_bulk(
                    count_mat = count_mat[,node$cells],
                    df = df %>% filter(cell %in% node$cells),
                    lambdas_ref = lambdas_ref,
                    gtf_transcript = gtf_transcript
                ) %>%
                mutate(node = node$label) %>%
                analyze_bulk_gpois(t = t, verbose = TRUE)
        }) %>%
        bind_rows() %>%
        arrange(CHROM, POS) %>%
        mutate(snp_id = factor(snp_id, unique(snp_id))) %>%
        mutate(snp_index = as.integer(snp_id)) %>%
        group_by(seg, node) %>%
        mutate(
            seg_start_index = min(snp_index),
            seg_end_index = max(snp_index)
        ) %>%
        ungroup() %>%
        mutate(sample = node) %>%
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

    res[['bulk_all']] = bulk_all
    res[['plot_list']] = plot_list
    
    #### 3. Find consensus CNVs ####
    
    if (verbose) {
        display('Finding consensus CNVs ..')
    }
    
    segs_all = bulk_all %>% 
        filter(state != 'neu') %>%
        distinct(node, CHROM, seg, cnv_state, cnv_state_post, seg_start, seg_end,
                 theta_mle, theta_sigma, phi_mle, phi_sigma, p_loh, p_del, p_amp, p_bamp, p_bdel, LLR, LLR_y, n_genes, n_snps)
    
    segs_filtered = segs_all %>% filter(!(LLR_y < 20 & cnv_state %in% c('del', 'amp'))) %>% filter(n_genes >= 20)
        
    segs_consensus = segs_filtered %>% mutate(sample = node) %>% resolve_cnvs() %>% fill_neu_segs(gbuild = gbuild)

    res[['segs_all']] = segs_all
    res[['segs_filtered']] = segs_filtered
    res[['segs_consensus']] = segs_consensus

    #### 4. Per-cell CNV evaluations ####
    
    if (verbose) {
        display('Calculating per-cell CNV posteriors ..')
    }

    exp_post = get_exp_post(
        segs_consensus %>% mutate(cnv_state = ifelse(cnv_state == 'neu', cnv_state, cnv_state_post)),
        count_mat,
        lambdas_ref,
        gtf_transcript = gtf_transcript,
        ncores = ncores)
    
    
    allele_post = get_allele_post(
        bulk_all,
        segs_consensus %>% mutate(cnv_state = ifelse(cnv_state == 'neu', cnv_state, cnv_state_post)),
        df)
    
    joint_post = get_joint_post(
        exp_post,
        allele_post,
        segs_consensus
    )

    joint_post = joint_post %>% left_join(cell_annot)
    
    if (verbose) {
        display('All done!')
    }
    
    res[['exp_post']] = exp_post
    res[['allele_post']] = allele_post
    res[['joint_post']] = joint_post
    
    return(res)
}

armadillo_subclone = function(count_mat, lambdas_ref, df, gtf_transcript, out_dir = './', t = 1e-5, sample_size = 400, min_cells = 200, min_depth = 0, ncores = 30, gbuild = 'hg38', verbose = TRUE) {
    
    dir.create(out_dir, showWarnings = FALSE)

    res = list()

    if (verbose) {
        display('Running initial HMM ..')
    }

    # Initialize bulk HMM 
    bulk_subtrees = get_bulk(
            count_mat = count_mat,
            df = df,
            lambdas_ref = lambdas_ref,
            gtf_transcript = gtf_transcript,
            min_depth = min_depth
        ) %>%
        analyze_bulk_gpois(t = t) %>%
        mutate(sample = 0)

    fwrite(bulk_subtrees, glue('{out_dir}/bulk_subtrees_0.tsv'), sep = '\t')

    # resolve CNVs
    segs_consensus = get_segs_consensus(bulk_subtrees, gbuild = gbuild)

    res[['0']] = list(bulk_subtrees, segs_consensus)

    for (i in 1:2) {

        if (verbose) {
            display(glue('Iteration {i}'))
        }

        ######## Evaluate CNV per cell ########

        if (verbose) {
            display('Evaluating CNV per cell ..')
        }

        exp_post = get_exp_post(
            segs_consensus %>% mutate(cnv_state = ifelse(cnv_state == 'neu', cnv_state, cnv_state_post)),
            count_mat,
            lambdas_ref,
            gtf_transcript = gtf_transcript,
            ncores = ncores
        )

        allele_post = get_allele_post(
            bulk_subtrees,
            segs_consensus %>% mutate(cnv_state = ifelse(cnv_state == 'neu', cnv_state, cnv_state_post)),
            df
        )

        joint_post = get_joint_post(
            exp_post,
            allele_post,
            segs_consensus)

        fwrite(exp_post, glue('{out_dir}/exp_post_{i}.tsv'), sep = '\t')
        fwrite(allele_post, glue('{out_dir}/allele_post_{i}.tsv'), sep = '\t')
        fwrite(joint_post, glue('{out_dir}/joint_post_{i}.tsv'), sep = '\t')

        ######## Build phylogeny ########

        if (verbose) {
            display('Building phylogeny ..')
        }

        cell_sample = colnames(count_mat) %>% sample(min(sample_size, length(.)))

        p_min = 1e-10

        geno = joint_post %>%
            group_by(seg) %>%
            filter(mean(p_cnv > 0.95) > 0.05) %>%
            ungroup() %>%
            filter(cell %in% cell_sample) %>%
            filter(cnv_state != 'neu') %>%
            mutate(p_n = 1 - p_cnv) %>%
            mutate(p_n = pmax(pmin(p_n, 1-p_min), p_min)) %>%
            reshape2::dcast(seg ~ cell, value.var = 'p_n', fill = 0.5) %>%
            tibble::column_to_rownames('seg')

        fwrite(geno, glue('{out_dir}/geno_{i}.tsv'), row.names = T, sep = '\t')

        n_cnvs = geno %>% nrow
        n_cells = ncol(geno)
        in_file = glue('{out_dir}/geno_{i}.txt')

        fwrite(geno, in_file, row.names = T, sep = ' ', quote = F)

        cmd = glue('sed -i "1s;^;HAPLOID {n_cnvs} {n_cells};" {in_file}')
        system(cmd)

        # run 
        out_file = glue('{out_dir}/out_{i}.txt')
        mut_tree_file = glue('{out_dir}/mut_tree_{i}.gml')
        cmd = glue('scistree -k 10 -v -e -t 0.9 -o {mut_tree_file} {in_file} > {out_file}')
        system(cmd, wait = T)

        # read output
        scistree_out = parse_scistree(out_file, geno, joint_post)
        tree_post = get_tree_post(scistree_out$MLtree, geno)

        saveRDS(scistree_out, glue('{out_dir}/scistree_out_{i}.rds'))
        saveRDS(tree_post, glue('{out_dir}/tree_post_{i}.rds'))

        # simplify mutational history
        G_m = mut_tree_file %>% 
            read_mut_tree(tree_post$mut_nodes) %>%
            simplify_history(tree_post$l_matrix, max_cost = 100) %>%
            label_genotype()

        mut_nodes = G_m %>% as_data_frame('vertices') %>% select(name = node, site = label)

        # update tree
        gtree = mut_to_tree(tree_post$gtree, mut_nodes)

        # map cells to the phylogeny
        clone_post = cell_to_clone(gtree, exp_post, allele_post)

        clones = clone_post %>% split(.$clone) %>%
            map(function(c){list(clone = unique(c$clone), GT = unique(c$GT), cells = c$cell, size = length(c$cell))})

        clones = keep(clones, function(x) x$size > min_cells)

        subtrees = lapply(1:vcount(G_m), function(v) {
            G_m %>%
            as_tbl_graph %>% 
            mutate(rank = dfs_rank(root = v)) %>%
            filter(!is.na(rank)) %>%
            data.frame() %>%
            inner_join(clone_post, by = c('GT')) %>%
            {list(subtree = v, members = unique(.$GT), cells = .$cell, size = length(.$cell))}
        })

        subtrees = keep(subtrees, function(x) x$size > min_cells)

        saveRDS(clones, glue('{out_dir}/clones_{i}.rds'))
        saveRDS(subtrees, glue('{out_dir}/subtrees_{i}.rds'))

        ######## Run HMMs ########
        
        if (verbose) {
            display('Running HMMs ..')
        }

        bulk_clones = mclapply(
                clones,
                mc.cores = max(length(clones), 1),
                function(clone) {
                    get_bulk(
                        count_mat = count_mat[,clone$cells],
                        df = df %>% filter(cell %in% clone$cells),
                        lambdas_ref = lambdas_ref,
                        gtf_transcript = gtf_transcript,
                        min_depth = min_depth
                    ) %>%
                    mutate(
                        sample = clone$clone,
                        n_cells = clone$size
                    ) %>%
                    analyze_bulk_gpois(t = t)
            }) %>%
            bind_rows() %>%
            arrange(CHROM, POS) %>%
            mutate(snp_id = factor(snp_id, unique(snp_id))) %>%
            mutate(snp_index = as.integer(snp_id)) %>%
            arrange(sample)
        
        bulk_subtrees = mclapply(
                subtrees,
                mc.cores = max(length(subtrees), 1),
                function(subtree) {
                    get_bulk(
                        count_mat = count_mat[,subtree$cells],
                        df = df %>% filter(cell %in% subtree$cells),
                        lambdas_ref = lambdas_ref,
                        gtf_transcript = gtf_transcript,
                        min_depth = min_depth
                    ) %>%
                    mutate(
                        sample = subtree$subtree,
                        n_cells = subtree$size
                    ) %>%
                    analyze_bulk_gpois(t = t)
            }) %>%
            bind_rows() %>%
            arrange(CHROM, POS) %>%
            mutate(snp_id = factor(snp_id, unique(snp_id))) %>%
            mutate(snp_index = as.integer(snp_id)) %>%
            group_by(seg, sample) %>%
            mutate(
                seg_start_index = min(snp_index),
                seg_end_index = max(snp_index)
            ) %>%
            ungroup() %>%
            arrange(sample) 

        fwrite(bulk_clones, glue('{out_dir}/bulk_clones_{i}.tsv'), sep = '\t')
        fwrite(bulk_subtrees, glue('{out_dir}/bulk_subtrees_{i}.tsv'), sep = '\t')

        ######## Find consensus CNVs ########

        segs_consensus = get_segs_consensus(bulk_subtrees, gbuild = gbuild)

        res[[as.character(i)]] = list(cell_sample, exp_post, allele_post, joint_post, tree_post, G_m, gtree, subtrees, bulk_subtrees, bulk_clones, segs_consensus)

    }

    return(res)
}

get_segs_consensus = function(bulk_all, gbuild = 'hg38') {

    segs_all = bulk_all %>% 
        filter(state != 'neu') %>%
        distinct(sample, CHROM, seg, cnv_state, cnv_state_post, seg_start, seg_end, seg_start_index, seg_end_index,
                theta_mle, theta_sigma, phi_mle, phi_sigma, p_loh, p_del, p_amp, p_bamp, p_bdel, LLR, LLR_y, n_genes, n_snps)
    
    segs_filtered = segs_all %>% filter(!(LLR_y < 20 & cnv_state %in% c('del', 'amp'))) %>% filter(n_genes >= 20)
    
    segs_consensus = segs_filtered %>% resolve_cnvs() %>% fill_neu_segs(gbuild = gbuild)

    return(segs_consensus)

}

cell_to_clone = function(gtree, exp_post, allele_post) {
    
    clones = gtree %>% data.frame() %>% 
        group_by(GT) %>%
        summarise(
            clone_size = n()
        ) %>%
        mutate(prior_clone = clone_size/sum(clone_size)) %>%
        mutate(seg = GT) %>%
        tidyr::separate_rows(seg, sep = ',') %>%
        mutate(I = 1) %>%
        tidyr::complete(seg, tidyr::nesting(GT, prior_clone, clone_size), fill = list('I' = 0)) %>%
        mutate(clone = as.integer(factor(GT))) 

    clone_post = inner_join(
            exp_post %>%
                filter(cnv_state != 'neu') %>%
                inner_join(clones, by = c('seg' = 'seg')) %>%
                mutate(l_clone = ifelse(I == 1, Z_cnv, Z_n)) %>%
                group_by(cell, clone, GT, prior_clone) %>%
                summarise(
                    l_clone_x = sum(l_clone),
                    .groups = 'drop'
                ),
            allele_post %>%
                filter(cnv_state != 'neu') %>%
                inner_join(clones, by = c('seg' = 'seg')) %>%
                mutate(l_clone = ifelse(I == 1, Z_cnv, Z_n)) %>%
                group_by(cell, clone, GT, prior_clone) %>%
                summarise(
                    l_clone_y = sum(l_clone),
                    .groups = 'drop'
                ),
            by = c("cell", "clone", "GT", "prior_clone")
        ) %>%
        group_by(cell) %>%
        mutate(
            Z_clone = log(prior_clone) + l_clone_x + l_clone_y,
            post_clone = exp(Z_clone - matrixStats::logSumExp(Z_clone))
        ) %>%
        summarise(
            clone = clone[which.max(post_clone)],
            GT = GT[which.max(post_clone)]
        )
    
    return(clone_post)
    
}



# resolve overlapping calls by graph reduction
resolve_cnvs = function(segs_all, debug = FALSE) {
            
    V = segs_all %>% mutate(vertex = 1:n(), .before = 'CHROM')

    E = segs_all %>% {GenomicRanges::GRanges(
            seqnames = .$CHROM,
            IRanges::IRanges(start = .$seg_start_index,
                end = .$seg_end_index)
        )} %>%
        GenomicRanges::findOverlaps(., .) %>%
        as.data.frame %>%
        setNames(c('from', 'to')) %>% 
        filter(from != to) %>%
        rowwise() %>%
        mutate(vp = paste0(sort(c(from, to)), collapse = ',')) %>%
        ungroup() %>%
        distinct(vp, .keep_all = T)

    # cut some edges with weak overlaps
    E = E %>% 
        left_join(
            V %>% select(from = vertex, start_x = seg_start_index, end_x = seg_end_index),
            by = 'from'
        ) %>%
        left_join(
            V %>% select(to = vertex, start_y = seg_start_index, end_y = seg_end_index),
            by = 'to'
        ) %>%
        mutate(
            len_x = end_x - start_x,
            len_y = end_y - start_y,
            len_overlap = pmin(end_x, end_y) - pmax(start_x, start_y),
            frac_overlap_x = len_overlap/len_x,
            frac_overlap_y = len_overlap/len_y
        ) %>%
        filter(!(frac_overlap_x < 0.5 & frac_overlap_y < 0.5))

    G = igraph::graph_from_data_frame(d=E, vertices=V, directed=F)

    segs_all = segs_all %>% mutate(component = igraph::components(G)$membership)

    segs_consensus = segs_all %>% group_by(component, sample) %>%
        mutate(LLR_sample = max(LLR)) %>%
        arrange(CHROM, component, -LLR_sample) %>%
        group_by(component) %>%
        filter(sample == sample[which.max(LLR_sample)])

    # segs_consensus = segs_all %>% arrange(CHROM, group, -LLR) %>% distinct(group, `.keep_all` = TRUE) 

    segs_consensus = segs_consensus %>% arrange(CHROM, seg_start) %>%
        mutate(CHROM = factor(CHROM, 1:22))
    
    if (debug) {
        return(list('G' = G, 'segs_consensus' = segs_consensus))
    }
    
    return(segs_consensus)
}


# # retrieve neutral segments
# fill_neu_segs = function(segs_consensus, gbuild = 'hg38') {
    
#     chrom_sizes = fread(glue('~/ref/{gbuild}.chrom.sizes.txt')) %>% 
#             set_names(c('CHROM', 'LEN')) %>%
#             mutate(CHROM = str_remove(CHROM, 'chr')) %>%
#             filter(CHROM %in% 1:22) %>%
#             mutate(CHROM = as.factor(as.integer(CHROM)))

#     out_of_bound = segs_consensus %>% left_join(chrom_sizes, by = "CHROM") %>% 
#         filter(seg_end > LEN) %>% pull(seg)

#     if (length(out_of_bound) > 0) {
#         warning(glue('Segment end exceeds genome length: {out_of_bound}'))
#         chrom_sizes = chrom_sizes %>%
#             left_join(
#                 segs_filtered %>% group_by(CHROM) %>%
#                     summarise(seg_end = max(seg_end)),
#                 by = "CHROM"
#             ) %>%
#             mutate(LEN = ifelse(is.na(seg_end), LEN, pmax(LEN, seg_end)))
#     }

#     segs_consensus = c(
#             segs_consensus %>% {GenomicRanges::GRanges(
#                 seqnames = .$CHROM,
#                 IRanges::IRanges(start = .$seg_start,
#                        end = .$seg_end)
#             )},
#             chrom_sizes %>% {GenomicRanges::GRanges(
#                 seqnames = .$CHROM,
#                 IRanges::IRanges(start = 1,
#                        end = .$LEN)
#             )}
#         ) %>%
#         GenomicRanges::disjoin() %>%
#         as.data.frame() %>%
#         select(CHROM = seqnames, seg_start = start, seg_end = end) %>%
#         left_join(
#             segs_consensus,
#             by = c("CHROM", "seg_start", "seg_end")
#         ) %>%
#         mutate(cnv_state = tidyr::replace_na(cnv_state, 'neu')) %>%
#         group_by(CHROM) %>%
#         mutate(seg_cons = paste0(CHROM, '_', 1:n())) %>%
#         ungroup() %>%
#         mutate(CHROM = as.factor(CHROM))
    
#     return(segs_consensus)
# }

# retrieve neutral segments
fill_neu_segs = function(segs_consensus, gbuild = 'hg38') {
    
    chrom_sizes = fread(glue('~/ref/{gbuild}.chrom.sizes.txt')) %>% 
            set_names(c('CHROM', 'LEN')) %>%
            mutate(CHROM = str_remove(CHROM, 'chr')) %>%
            filter(CHROM %in% 1:22) %>%
            mutate(CHROM = as.factor(as.integer(CHROM)))

    out_of_bound = segs_consensus %>% left_join(chrom_sizes, by = "CHROM") %>% 
        filter(seg_end > LEN) %>% pull(seg)

    if (length(out_of_bound) > 0) {
        warning(glue('Segment end exceeds genome length: {out_of_bound}'))
        chrom_sizes = chrom_sizes %>%
            left_join(
                segs_consensus %>% group_by(CHROM) %>%
                    summarise(seg_end = max(seg_end)),
                by = "CHROM"
            ) %>%
            mutate(LEN = ifelse(is.na(seg_end), LEN, pmax(LEN, seg_end)))
    }

    gaps = GenomicRanges::setdiff(
        chrom_sizes %>% {GenomicRanges::GRanges(
            seqnames = .$CHROM,
            IRanges::IRanges(start = 1,
                   end = .$LEN)
        )},
        segs_consensus %>% 
            group_by(CHROM) %>%
            summarise(seg_start = min(seg_start), seg_end = max(seg_end)) %>%
            ungroup() %>%
            {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$seg_start,
                       end = .$seg_end)
            )},
        ) %>%
        as.data.frame() %>%
        select(CHROM = seqnames, seg_start = start, seg_end = end)

    segs_consensus = segs_consensus %>%
        bind_rows(gaps) %>% 
        mutate(cnv_state = tidyr::replace_na(cnv_state, 'neu')) %>%
        arrange(CHROM, seg_start) %>%
        group_by(CHROM) %>%
        mutate(seg_cons = paste0(CHROM, '_', 1:n())) %>%
        ungroup() %>%
        mutate(CHROM = as.factor(CHROM)) 
    
    return(segs_consensus)
}


# multi-state model
get_exp_likelihoods = function(exp_sc, alpha = NULL, beta = NULL) {
    
    depth_obs = sum(exp_sc$Y_obs)
    exp_sc = exp_sc %>% filter(lambda_ref > 0)
    
    if (is.null(alpha) & is.null(beta)) {

        # if (hskd) {

        #     exp_sc = exp_sc %>% mutate(exp_bin = as.factor(ntile(lambda_ref, 4)))

        #     fits = exp_sc %>% 
        #         filter(cnv_state == 'neu') %>%
        #         group_by(exp_bin) %>%
        #         do({
        #             display(depth_obs)
        #             coef = fit_gpois(.$Y_obs, .$lambda_ref, depth_obs)@coef
        #             alpha = coef[1]
        #             beta = coef[2]
        #             data.frame(
        #                 alpha = alpha,
        #                 beta = beta,
        #                 mean = alpha/beta,
        #                 var = alpha/beta^2
        #             )
        #         })

        #     exp_sc = exp_sc %>%
        #         select(-any_of(c('alpha', 'beta'))) %>%
        #         left_join(fits, by = 'exp_bin')

        fit = exp_sc %>% filter(cnv_state %in% c('neu', 'loh')) %>% {fit_gpois(.$Y_obs, .$lambda_ref, depth_obs)}
        
        exp_sc = exp_sc %>% mutate(alpha = fit@coef[1], beta = fit@coef[2])
    }
        
    res = exp_sc %>% 
        group_by(seg) %>%
        summarise(
            n = n(),
            alpha = unique(alpha),
            beta = unique(beta),
            cnv_state = unique(cnv_state),
            phi_mle = calc_phi_mle(Y_obs, lambda_ref, depth_obs, alpha, beta, lower = 0.1, upper = 10),
            l11 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = 1),
            l20 = l11,
            l10 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = 0.5),
            l21 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = 1.5),
            l31 = l_gpois(Y_obs, lambda_ref, depth_obs, alpha, beta, phi = 2),
            l22 = l31,
            l00 = l10
        )
        
        
    return(res)
}


get_exp_post = function(segs_consensus, count_mat, gtf_transcript, lambdas_ref = NULL, lambdas_fit = NULL, alpha = NULL, beta = NULL, ncores = 30, debug = F) {

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
        inner_join(
            gene_seg %>% select(CHROM, gene, seg = seg_cons, seg_start, seg_end, gene_start, cnv_state),
            by = "gene"
        ) %>%
        arrange(CHROM, gene_start) %>%
        mutate(gene_index = 1:n()) %>%
        group_by(seg) %>%
        mutate(
            seg_start_index = min(gene_index),
            seg_end_index = max(gene_index),
            n_genes = n()
        ) %>%
        ungroup()
    
    cells = colnames(count_mat)
    
    # should write a try catch to handle failures in some cells..
    results = mclapply(
        cells,
        mc.cores = ncores,
        function(cell) {

            exp_sc = exp_sc[,c('gene', 'seg', 'CHROM', 'cnv_state', 'seg_start', 'seg_end', cell)] %>%
                rename(Y_obs = ncol(.))
            
            if (is.null(lambdas_fit)) {
                lambdas_fit = exp_sc %>% {fit_multi_ref(setNames(.$Y_obs, .$gene), lambdas_ref, sum(.$Y_obs), gtf_transcript)}
            }

            exp_sc %>%
                mutate(
                    lambda_ref = lambdas_fit$lambdas_bar[gene],
                    lambda_obs = Y_obs/sum(Y_obs),
                    logFC = log2(lambda_obs/lambda_ref)
                ) %>%
                filter(logFC < 5) %>%
                get_exp_likelihoods(alpha, beta) %>%
                mutate(cell = cell, w = paste0(signif(lambdas_fit$w, 2), collapse = ','))
        }
    )

    exp_post = results[map(results, class) != 'try-error'] %>% bind_rows()
    
    exp_post = exp_post %>% 
        mutate(seg = factor(seg, gtools::mixedsort(unique(seg)))) %>%
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
            Z_cnv = matrixStats::logSumExp(
                c(l20 + log(p_loh/2),
                l10 + log(p_del/2),
                l21 + log(p_amp/4),
                l31 + log(p_amp/4),
                l22 + log(p_bamp/2),
                l00 + log(p_bdel/2))
            ),
            Z_n = l11 + log(1/2),
            logBF = Z_cnv - Z_n,
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

    if ((!'sample' %in% colnames(bulk_all)) | (!'sample' %in% colnames(segs_consensus))) {
        bulk_all['sample'] = '0'
        segs_consensus['sample'] = '0'
        warning('Sample column missing')
    }
    
    # allele posteriors
    snp_seg = bulk_all %>%
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

    allele_sc = df %>%
        mutate(pAD = ifelse(GT == '1|0', AD, DP - AD)) %>%
        select(-snp_index) %>% 
        inner_join(
            snp_seg %>% select(snp_id, snp_index, haplo, seg = seg_cons, cnv_state),
            by = c('snp_id')
        ) %>%
        filter(!cnv_state %in% c('neu', 'bamp', 'bdel')) %>%
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
            Z_cnv = matrixStats::logSumExp(
                c(l20 + log(p_loh/2),
                  l10 + log(p_del/2),
                  l21 + log(p_amp/4),
                  l31 + log(p_amp/4),
                  l22 + log(p_bamp/2),
                  l00 + log(p_bdel/2)
                 )
            ),
            Z_n = l11 + log(1/2),
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

    cells_common = intersect(exp_post$cell, allele_post$cell)
    exp_post = exp_post %>% filter(cell %in% cells_common)
    allele_post = allele_post %>% filter(cell %in% cells_common)

    joint_post = exp_post %>%
        filter(cnv_state != 'neu') %>%
        select(
            cell, seg, cnv_state, l11_x = l11, l20_x = l20,
            l10_x = l10, l21_x = l21, l31_x = l31, l22_x = l22, l00_x = l00,
            Z_x = Z, Z_cnv_x = Z_cnv, Z_n_x = Z_n
        ) %>%
        full_join(
            allele_post %>% select(
                cell, seg, l11_y = l11, l20_y = l20, l10_y = l10, l21_y = l21, l31_y = l31, l22_y = l22, l00_y = l00,
                n_snp = total,
                Z_y = Z, Z_cnv_y = Z_cnv, Z_n_y = Z_n
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
            Z_cnv = matrixStats::logSumExp(
                c(l20_x + l20_y + log(p_loh/2),
                  l10_x + l10_y + log(p_del/2),
                  l21_x + l21_y + log(p_amp/4),
                  l31_x + l31_y + log(p_amp/4),
                  l22_x + l22_y + log(p_bamp/2),
                  l00_x + l00_y + log(p_bdel/2))
            ),
            Z_n = l11_x + l11_y + log(1/2),
            p_cnv = exp(Z_cnv - Z),
            p_n = exp(Z_n - Z),
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

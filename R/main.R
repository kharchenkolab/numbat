#' @import logger
#' @import dplyr
#' @importFrom data.table fread fwrite as.data.table
#' @import stringr
#' @import glue
#' @importFrom parallel mclapply
#' @import tidygraph
#' @import ggplot2
#' @import ggtree
#' @import ggraph
#' @importFrom igraph vcount ecount E V V<- E<-
#' @import patchwork
#' @importFrom extraDistr dgpois
#' @useDynLib Numbat

#' @description Run workflow to decompose tumor subclones
#' @param count_mat raw count matrices where rownames are genes and column names are cells
#' @param lambdas_ref either a named vector with gene names as names and normalized expression as values, or a matrix where rownames are genes and columns are pseudobulk names
#' @param df_allele dataframe of allele counts per cell, produced by preprocess_allele
#' @param gtf_transcript gtf dataframe of transcripts 
#' @param genetic_map a genetic map dataframe
#' @param out_dir output directory
#' @param gamma dispersion parameter for the Beta-Binomial allele model
#' @param t transition probability
#' @param init_method initialization method, either "smooth" (clustering based on smoothed expression) or "bulk" (all cell pseudobulk)
#' @param init_k number of clusters in the initial clustering
#' @param sample_size subsampling size
#' @param min_cells minimum number of cells to run HMM on
#' @param max_cost maximum cost to simplify the phylogeny
#' @param max_iter maximum number of iterations to run the phyologeny optimization
#' @param min_dpeth minimum allele depth 
#' @param common_diploid whether to find common diploid regions in a group of peusdobulks
#' @param ncores number of threads to use
#' @param max_entropy entorpy threshold to filter CNVs by
#' @param min_LLR minimum LLR to filter CNVs by 
#' @param eps convergence threshold for ML tree search
#' @param multi_allelic whether to call multi-allelic CNVs
#' @param use_loh whether to include LOH regions in the expression baseline
#' @param skip_nj whether to skip NJ tree construction and only use UPGMA
#' @param exclude_normal whether to exclude normal cells in the tree construction
#' @param diploid_chroms known diploid chromosomes
#' @return a status code
#' @export
run_numbat = function(
        count_mat, lambdas_ref, df_allele, gtf_transcript, genetic_map, 
        out_dir = './', max_iter = 2, t = 1e-5, gamma = 20, min_LLR = 50, alpha = 1e-4, eps = 1e-5, max_entropy = 0.6, 
        init_method = 'smooth', init_k = 3, sample_size = 1e5, min_cells = 10, max_cost = ncol(count_mat) * 0.3, 
        min_depth = 0, common_diploid = TRUE, ncores = 30, exp_model = 'lnpois', 
        verbose = TRUE, diploid_chroms = NULL, use_loh = NULL,
        exclude_normal = FALSE, skip_nj = FALSE, multi_allelic = TRUE,
        plot = TRUE
    ) {
    
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    logfile = glue('{out_dir}/log.txt')
    if (file.exists(logfile)) {file.remove(logfile)}
    log_appender(appender_file(logfile))

    log_info(paste(
        'Running under parameters:',
        glue('t = {t}'), 
        glue('alpha = {alpha}'),
        glue('gamma = {gamma}'),
        glue('min_cells = {min_cells}'), 
        glue('init_k = {init_k}'),
        glue('sample_size = {sample_size}'),
        glue('max_cost = {max_cost}'),
        glue('max_iter = {max_iter}'),
        glue('min_depth = {min_depth}'),
        glue('use_loh = {ifelse(is.null(use_loh), "auto", use_loh)}'),
        glue('multi_allelic = {multi_allelic}'),
        glue('min_LLR = {min_LLR}'),
        glue('max_entropy = {max_entropy}'),
        glue('skip_nj = {skip_nj}'),
        glue('exclude_normal = {exclude_normal}'),
        glue('diploid_chroms = {paste0(diploid_chroms, collapse = ",")}'),
        glue('ncores = {ncores}'),
        glue('common_diploid = {common_diploid}'),
        'Input metrics:',
        glue('{ncol(count_mat)} cells'),
        sep = "\n"
    ))

    ######## Initialization ########
    i = 0

    if (init_method == 'bulk') {

        if (verbose) {log_info('Initializing using all-cell bulk ..')}   
        bulk_subtrees = get_bulk(
                count_mat = count_mat,
                df_allele = df_allele,
                lambdas_ref = lambdas_ref,
                gtf_transcript = gtf_transcript,
                genetic_map = genetic_map,
                min_depth = min_depth
            ) %>%
            analyze_bulk(t = t, gamma = gamma) %>%
            mutate(sample = 0)

    } else if (init_method == 'smooth') {

        if (verbose) {log_info('Approximating initial clusters using smoothed expression ..')}

        clust = exp_hclust(
            count_mat,
            lambdas_ref = lambdas_ref,
            gtf_transcript,
            k = init_k,
            ncores = ncores
        )

        fwrite(
            as.data.frame(clust$gexp_roll_wide) %>% tibble::rownames_to_column('cell'),
            glue('{out_dir}/gexp_roll_wide.tsv.gz'),
            sep = '\t',
            nThread = min(4, ncores)
        )
        saveRDS(clust$hc, glue('{out_dir}/hc.rds'))
        saveRDS(clust$nodes, glue('{out_dir}/hc_nodes.rds'))

        nodes = purrr::keep(clust$nodes, function(x) x$size > min_cells)

        bulk_subtrees = make_group_bulks(
                groups = nodes,
                count_mat = count_mat,
                df_allele = df_allele, 
                lambdas_ref = lambdas_ref,
                gtf_transcript = gtf_transcript,
                genetic_map = genetic_map,
                min_depth = min_depth,
                ncores = ncores
            )

        bulk_subtrees = bulk_subtrees %>% 
            run_group_hmms(
                t = t,
                gamma = gamma,
                alpha = alpha,
                common_diploid = common_diploid,
                diploid_chroms = diploid_chroms,
                exp_model = exp_model,
                ncores = ncores,
                verbose = verbose
            )

    } else {
        stop('init_method can be raw, bulk, or smooth')
    }

    fwrite(bulk_subtrees, glue('{out_dir}/bulk_subtrees_{i}.tsv.gz'), sep = '\t')

    # resolve CNVs
    segs_consensus = get_segs_consensus(bulk_subtrees)

    if (multi_allelic) {
        # only test the leaf nodes
        out = test_multi_allelic(
            bulk_subtrees %>% filter(sample %in% 1:init_k),
            segs_consensus,
            use_loh = use_loh,
            diploid_chroms = diploid_chroms)

        segs_consensus = out$segs_consensus

        fwrite(out$bulks, glue('{out_dir}/bulk_clones_retest_{i}.tsv.gz'), sep = '\t')

    }

    fwrite(segs_consensus, glue('{out_dir}/segs_consensus_0.tsv'), sep = '\t')

    if (plot) {

        log_info('Making initial visualization ..')

        p = plot_bulks(bulk_subtrees)

        ggsave(glue('{out_dir}/bulk_subtrees_{i}.png'), p, width = 14, height = 2*length(unique(bulk_subtrees$sample)), dpi = 200)
        
    }

    normal_cells = c()

    ######## Begin iterations ########
    for (i in 1:max_iter) {
        
        log_info('Iteration {i}')

        ######## Evaluate CNV per cell ########

        log_info('Evaluating CNV per cell ..')

        exp_post_res = get_exp_post(
            segs_consensus %>% mutate(cnv_state = ifelse(cnv_state == 'neu', cnv_state, cnv_state_post)),
            count_mat,
            lambdas_ref,
            use_loh = use_loh,
            gtf_transcript = gtf_transcript,
            ncores = ncores)

        log_info('Done with computing expression posterior ..')

        exp_post = exp_post_res$exp_post
        exp_sc = exp_post_res$exp_sc

        fwrite(exp_sc, glue('{out_dir}/exp_sc_{i}.tsv.gz'), sep = '\t')
        
        allele_post = get_allele_post(
            bulk_subtrees,
            segs_consensus %>% mutate(cnv_state = ifelse(cnv_state == 'neu', cnv_state, cnv_state_post)),
            df_allele
        )

        log_info('Done with computing allele posterior ..')

        joint_post = get_joint_post(
            exp_post,
            allele_post,
            segs_consensus)

        log_info('Done with computing joint posterior ..')

        joint_post = joint_post %>%
            group_by(seg) %>%
            mutate(
                avg_entropy = mean(binary_entropy(p_cnv), na.rm = TRUE)
            ) %>%
            ungroup()

        fwrite(exp_post, glue('{out_dir}/exp_post_{i}.tsv'), sep = '\t')
        fwrite(allele_post, glue('{out_dir}/allele_post_{i}.tsv'), sep = '\t')
        fwrite(joint_post, glue('{out_dir}/joint_post_{i}.tsv'), sep = '\t')

        if (multi_allelic) {
            log_info('Expanding allelic states..')
            # expand multi-allelic segments into multiple CNV states
            joint_post = expand_states(joint_post, segs_consensus)
            exp_post = expand_states(exp_post, segs_consensus)
            allele_post = expand_states(allele_post, segs_consensus)

            fwrite(exp_post, glue('{out_dir}/exp_post_expand_{i}.tsv'), sep = '\t')
            fwrite(allele_post, glue('{out_dir}/allele_post_expand_{i}.tsv'), sep = '\t')
            fwrite(joint_post, glue('{out_dir}/joint_post_expand_{i}.tsv'), sep = '\t')
        }

        ######## Build phylogeny ########

        log_info('Building phylogeny ..')

        cell_sample = colnames(count_mat) %>% sample(min(sample_size, length(.)), replace = FALSE)
        
        if (exclude_normal) {
            cell_sample = cell_sample[!cell_sample %in% normal_cells]
        }

        # filter CNVs
        joint_post_filtered = joint_post %>%
            filter(cnv_state != 'neu') %>%
            filter(cell %in% cell_sample) %>%
            filter(avg_entropy < max_entropy & LLR > min_LLR)

        if (nrow(joint_post_filtered) == 0) {
            msg = 'No CNVs remain after filtering! Check threshold'
            log_error(msg)
            stop(msg)
        } else {
            n_cnv = length(unique(joint_post_filtered$seg))
            log_info(glue('Using {n_cnv} CNVs to construct phylogeny'))
        }

        # construct genotype probability matrix
        p_min = 1e-10

        P = joint_post_filtered %>%
            mutate(p_n = 1 - p_cnv) %>%
            mutate(p_n = pmax(pmin(p_n, 1-p_min), p_min)) %>%
            reshape2::dcast(cell ~ seg, value.var = 'p_n', fill = 0.5) %>%
            tibble::column_to_rownames('cell') %>%
            as.matrix

        fwrite(as.data.frame(P), glue('{out_dir}/geno_{i}.tsv'), row.names = T, sep = '\t')

        # contruct initial tree
        dist_mat = parallelDist::parDist(rbind(P, 'outgroup' = 1), threads = ncores)

        treeUPGMA = phangorn::upgma(dist_mat) %>%
            ape::root(outgroup = 'outgroup') %>%
            ape::drop.tip('outgroup') %>%
            reorder(order = 'postorder')

        saveRDS(treeUPGMA, glue('{out_dir}/treeUPGMA_{i}.rds'))

        UPGMA_score = score_tree(treeUPGMA, as.matrix(P))$l_tree
        tree_init = treeUPGMA
        # note that dist_mat gets modified by NJ
        if(!skip_nj){
            treeNJ = phangorn::NJ(dist_mat) %>%
                ape::root(outgroup = 'outgroup') %>%
                ape::drop.tip('outgroup') %>%
                reorder(order = 'postorder')
            NJ_score = score_tree(treeNJ, as.matrix(P))$l_tree

            if (UPGMA_score > NJ_score) {
                log_info('Using UPGMA tree as seed..')
            } else {
                tree_init = treeNJ
                log_info('Using NJ tree as seed..')
            }
            saveRDS(treeNJ, glue('{out_dir}/treeNJ_{i}.rds'))
        } else{
            log_info('Only computing UPGMA..')
            log_info('Using UPGMA tree as seed..')
        }
        
        # maximum likelihood tree search with NNI
        tree_list = perform_nni(tree_init, P, ncores = ncores, eps = eps)
        saveRDS(tree_list, glue('{out_dir}/tree_list_{i}.rds'))

        tree_post = get_tree_post(tree_list[[length(tree_list)]], P)
        saveRDS(tree_post, glue('{out_dir}/tree_post_{i}.rds'))

        # simplify mutational history
        G_m = get_mut_tree(tree_post$gtree, tree_post$mut_nodes)  %>% 
            simplify_history(tree_post$l_matrix, max_cost = max_cost) %>% 
            label_genotype()

        mut_nodes = G_m %>% igraph::as_data_frame('vertices') %>% 
            select(name = node, site = label, clone = clone, GT = GT)

        # update tree
        gtree = mut_to_tree(tree_post$gtree, mut_nodes)
        gtree = mark_tumor_lineage(gtree)

        saveRDS(gtree, glue('{out_dir}/tree_final_{i}.rds'))
        saveRDS(G_m, glue('{out_dir}/mut_graph_{i}.rds'))

        clone_post = cell_to_clone(gtree, exp_post, allele_post)

        fwrite(clone_post, glue('{out_dir}/clone_post_{i}.tsv'), sep = '\t')

        normal_cells = clone_post %>% filter(GT_opt == '') %>% pull(cell)

        log_info('Found {length(normal_cells)} normal cells..')

        # form cell groupings
        clones = clone_post %>% split(.$clone_opt) %>%
            purrr::map(function(c){list(sample = unique(c$clone_opt), members = unique(c$GT_opt), cells = c$cell, size = length(c$cell))})

        saveRDS(clones, glue('{out_dir}/clones_{i}.rds'))
        clones = purrr::keep(clones, function(x) x$size > min_cells)

        clone_to_node = setNames(V(G_m)$id, V(G_m)$clone)

        subtrees = lapply(1:vcount(G_m), function(c) {
            G_m %>%
            as_tbl_graph %>% 
            mutate(rank = dfs_rank(root = clone_to_node[[as.character(c)]])) %>%
            filter(!is.na(rank)) %>%
            data.frame() %>%
            inner_join(clone_post, by = c('GT' = 'GT_opt')) %>%
            {list(sample = c, members = unique(.$GT), clones = unique(.$clone), cells = .$cell, size = length(.$cell))}
        })

        saveRDS(subtrees, glue('{out_dir}/subtrees_{i}.rds'))
        subtrees = purrr::keep(subtrees, function(x) x$size > min_cells)

        ######## Run HMMs ########

        bulk_clones = make_group_bulks(
                groups = clones,
                count_mat = count_mat,
                df_allele = df_allele, 
                lambdas_ref = lambdas_ref,
                gtf_transcript = gtf_transcript,
                genetic_map = genetic_map,
                min_depth = min_depth,
                ncores = ncores)

        bulk_clones = bulk_clones %>% 
            run_group_hmms(
                t = t,
                gamma = gamma,
                alpha = alpha,
                exp_model = exp_model,
                common_diploid = common_diploid,
                diploid_chroms = diploid_chroms,
                verbose = verbose,
                ncores = ncores
            )

        fwrite(bulk_clones, glue('{out_dir}/bulk_clones_{i}.tsv.gz'), sep = '\t')

        bulk_subtrees = make_group_bulks(
                groups = subtrees,
                count_mat = count_mat,
                df_allele = df_allele, 
                lambdas_ref = lambdas_ref,
                gtf_transcript = gtf_transcript,
                genetic_map = genetic_map,
                min_depth = min_depth,
                ncores = ncores)

        bulk_subtrees = bulk_subtrees %>%
            run_group_hmms(
                t = t,
                gamma = gamma,
                alpha = alpha,
                exp_model = exp_model,
                common_diploid = common_diploid,
                diploid_chroms = diploid_chroms,
                verbose = verbose
            )
        
        fwrite(bulk_subtrees, glue('{out_dir}/bulk_subtrees_{i}.tsv.gz'), sep = '\t')

        ######## plotting ########
        if (plot) {

            log_info('Making plots..')

            panel = plot_sc_joint(
                    gtree,
                    joint_post,
                    segs_consensus,
                    tip_length = 0.2,
                    branch_width = 0.2,
                    size = 0.15,
                    clone_bar = T,
                    logBF_min = 2,
                    logBF_max = 5
                )
            
            ggsave(glue('{out_dir}/panel_{i}.png'), panel, width = 8, height = 3.5, dpi = 200)

            p = plot_bulks(bulk_subtrees)

            ggsave(glue('{out_dir}/bulk_subtrees_{i}.png'), p, width = 14, height = 2*length(unique(bulk_subtrees$sample)), dpi = 200)

            p = plot_bulks(bulk_clones)

            ggsave(glue('{out_dir}/bulk_clones_{i}.png'), p, width = 14, height = 2*length(unique(bulk_clones$sample)), dpi = 200)

        }

        ######## Find consensus CNVs ########

        segs_consensus = get_segs_consensus(bulk_subtrees)

        if (multi_allelic) {

            out = test_multi_allelic(
                bulk_clones, 
                segs_consensus,
                use_loh = use_loh,
                diploid_chroms = diploid_chroms)

            segs_consensus = out$segs_consensus

            fwrite(out$bulks, glue('{out_dir}/bulk_clones_retest_{i}.tsv.gz'), sep = '\t')

        }

        fwrite(segs_consensus, glue('{out_dir}/segs_consensus_{i}.tsv'), sep = '\t')
    }
    # reflect change
    return('Success')
}

#' Run smoothed expression-based hclust
#' @keywords internal 
exp_hclust = function(count_mat_obs, lambdas_ref, gtf_transcript, k = 5, ncores = 10, verbose = T) {

    Y_obs = rowSums(count_mat_obs)

    fit = fit_multi_ref(Y_obs, lambdas_ref, sum(Y_obs), gtf_transcript)

    gexp_norm_long = process_exp_sc(
        count_mat_obs,
        fit$lambdas_bar,
        gtf_transcript,
        verbose = verbose
    )
    
    gexp_roll_wide = gexp_norm_long %>%
        reshape2::dcast(cell ~ gene, value.var = 'exp_rollmean') %>%
        tibble::column_to_rownames('cell') %>%
        as.matrix

    dist_mat = parallelDist::parDist(gexp_roll_wide, threads = ncores)

    log_info('running hclust...')
    hc = hclust(dist_mat, method = "ward.D2")

    cell_annot = data.frame(
        cell = colnames(count_mat_obs)
        ) %>%
        mutate(cluster = cutree(hc, k = k)[cell]) %>%
        mutate(group = 'obs')

    nodes = get_nodes_celltree(hc, cutree(hc, k = k))

    return(list('cell_annot' = cell_annot, 'nodes' = nodes, 'gexp_roll_wide' = gexp_roll_wide, 'hc' = hc, 'fit' = fit))
}

#' Make a group of pseudobulks
#' @keywords internal 
make_group_bulks = function(groups, count_mat, df_allele, lambdas_ref, gtf_transcript, genetic_map, min_depth = 0, ncores = NULL) {
    
    if (length(groups) == 0) {
        return(data.frame())
    }

    ncores = ifelse(is.null(ncores), length(groups), ncores)

    results = mclapply(
            groups,
            mc.cores = ncores,
            function(g) {
                get_bulk(
                    count_mat = count_mat[,g$cells],
                    df_allele = df_allele %>% filter(cell %in% g$cells),
                    lambdas_ref = lambdas_ref,
                    gtf_transcript = gtf_transcript,
                    genetic_map = genetic_map,
                    min_depth = min_depth
                ) %>%
                mutate(
                    n_cells = g$size,
                    members = paste0(g$members, collapse = ';'),
                    sample = g$sample
                )
        })

    bad = sapply(results, inherits, what = "try-error")

    if (any(bad)) {
        log_error(glue('job {paste(which(bad), collapse = ",")} failed'))
        log_error(results[bad][[1]])
    }

    # unify marker index
    bulks = results %>% 
        bind_rows() %>%
        arrange(CHROM, POS) %>%
        mutate(snp_id = factor(snp_id, unique(snp_id))) %>%
        mutate(snp_index = as.integer(snp_id)) %>%
        arrange(sample)

    return(bulks)
}

#' Run mutitple HMMs 
#' @export
run_group_hmms = function(
    bulks, t = 1e-4, gamma = 20, theta_min = 0.08,
    exp_model = 'lnpois', alpha = 1e-4,
    common_diploid = TRUE, diploid_chroms = NULL, allele_only = FALSE, retest = TRUE, run_hmm = TRUE,
    ncores = NULL, verbose = FALSE, debug = FALSE
) {

    if (nrow(bulks) == 0) {
        return(data.frame())
    }

    n_groups = length(unique(bulks$sample))

    if (verbose) {
        log_info('Running HMMs on {n_groups} cell groups..')
    }

    ncores = ifelse(is.null(ncores), n_groups, ncores)

    # find common diploid region
    if (!run_hmm) {
        find_diploid = FALSE
    } else if (common_diploid & is.null(diploid_chroms)) {
        diploid_out = find_common_diploid(bulks, gamma = gamma, alpha = alpha, ncores = ncores)
        bulks = diploid_out$bulks
        find_diploid = FALSE
    } else {
        find_diploid = TRUE
    }

    results = mclapply(
        bulks %>% split(.$sample),
        mc.cores = ncores,
        function(bulk) {
            bulk %>% analyze_bulk(
                t = t,
                gamma = gamma, 
                find_diploid = find_diploid, 
                run_hmm = run_hmm,
                allele_only = allele_only, 
                diploid_chroms = diploid_chroms,
                retest = retest, 
                verbose = verbose)
    })
                
    bad = sapply(results, inherits, what = "try-error")

    if (any(bad)) {
        log_error(glue('job {paste(which(bad), collapse = ",")} failed'))
        log_error(results[bad][[1]])
        message(results[bad][[1]])
    }

    bulks = results %>% bind_rows() %>%
        group_by(seg, sample) %>%
        mutate(
            seg_start_index = min(snp_index),
            seg_end_index = max(snp_index)
        ) %>%
        ungroup()

    return(bulks)
}

#' @description Extract consensus CNV segments
#' @param bulks pseudobulks dataframe
#' @param LLR_min LLR threshold to filter CNVs 
#' @param min_overlap minimum overlap fraction to determine count two events as as overlapping
#' @return consensus segments dataframe
#' @export
get_segs_consensus = function(bulks, LLR_min = 20, min_overlap = 0.45) {

    if (!'sample' %in% colnames(bulks)) {
        bulks$sample = 1
    }

    info_cols = c('sample', 'CHROM', 'seg', 'cnv_state', 'cnv_state_post',
            'seg_start', 'seg_end', 'seg_start_index', 'seg_end_index',
            'theta_mle', 'theta_sigma', 'phi_mle', 'phi_sigma', 
            'p_loh', 'p_del', 'p_amp', 'p_bamp', 'p_bdel',
            'LLR', 'LLR_y', 'LLR_x', 'n_genes', 'n_snps')

    segs_all = bulks %>% 
        group_by(sample, seg, CHROM) %>%
        mutate(seg_start = min(POS), seg_end = max(POS)) %>%
        ungroup() %>%
        select(any_of(info_cols)) %>%
        distinct()

    segs_filtered = segs_all %>% 
        filter(cnv_state != 'neu') %>%
        filter(LLR_y > LLR_min | LLR_x > 100 | theta_mle > 0.1 | cnv_state %in% c('bamp', 'bdel')) %>% 
        filter(n_genes >= 20)
    
    if(dim(segs_filtered)[[1]] == 0){
        stop('No CNV remains after filtering')
    }

    # reduce to unique intervals
    segs_neu = segs_all %>%
        filter(cnv_state == 'neu') %>%
        arrange(CHROM) %>%
        {GenomicRanges::GRanges(
            seqnames = .$CHROM,
            IRanges::IRanges(start = .$seg_start,
                end = .$seg_end)
        )} %>%
        GenomicRanges::reduce() %>%
        as.data.frame() %>%
        select(CHROM = seqnames, seg_start = start, seg_end = end) %>%
        mutate(seg_length = seg_end - seg_start)
    
    segs_consensus = segs_filtered %>% 
        resolve_cnvs(
            min_overlap = min_overlap
        ) %>% 
        fill_neu_segs(segs_neu) %>%
        mutate(cnv_state_post = ifelse(cnv_state == 'neu', cnv_state, cnv_state_post))

    return(segs_consensus)

}

#' @description Fill neutral regions into consensus segments
#' @param segs_consensus a dataframe containing info of all CNV segments from multiple samples
#' @param segs_neu neutral segments
#' @return collections of neutral and aberrant segments with no gaps
#' @keywords internal
fill_neu_segs = function(segs_consensus, segs_neu) {
    
    # take complement of consensus aberrant segs
    gaps = GenomicRanges::setdiff(
        segs_neu %>% {GenomicRanges::GRanges(
            seqnames = .$CHROM,
            IRanges::IRanges(start = .$seg_start,
                   end = .$seg_end)
        )},
        segs_consensus %>% 
            group_by(seg, CHROM) %>%
            summarise(
                seg_start = min(seg_start),
                seg_end = max(seg_end),
                .groups = 'drop'
            ) %>%
            ungroup() %>%
            {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$seg_start,
                       end = .$seg_end)
            )},
        ) %>%
        suppressWarnings() %>%
        as.data.frame() %>%
        select(CHROM = seqnames, seg_start = start, seg_end = end) %>%
        mutate(seg_length = seg_end - seg_start) %>%
        filter(seg_length > 0)

    segs_consensus = segs_consensus %>%
        mutate(seg_length = seg_end - seg_start) %>%
        bind_rows(gaps) %>% 
        mutate(cnv_state = tidyr::replace_na(cnv_state, 'neu')) %>%
        arrange(CHROM, seg_start) %>%
        group_by(CHROM) %>%
        mutate(seg_cons = paste0(CHROM, letters[1:n()])) %>%
        ungroup() %>%
        mutate(CHROM = factor(CHROM, 1:22)) %>%
        arrange(CHROM)
    
    return(segs_consensus)
}

#' @description Map cells to the phylogeny (or genotypes) based on CNV posteriors
#' @param gtree a cell lineage tree as a tidygraph object
#' @param exp_post expression posteriors
#' @param allele_post allele posteriors
#' @return clone posteriors
#' @export
cell_to_clone = function(gtree, exp_post, allele_post) {

    clones = gtree %>%
        activate(nodes) %>%
        data.frame() %>% 
        # mutate(GT = ifelse(compartment == 'normal', '', GT)) %>%
        group_by(GT, clone, compartment) %>%
        summarise(
            clone_size = sum(leaf),
            .groups = 'drop'
        )

    clone_segs = clones %>%
        mutate(
            prior_clone = ifelse(GT == '', 0.5, 0.5/(length(unique(GT)) - 1))
        ) %>%
        mutate(seg = GT) %>%
        tidyr::separate_rows(seg, sep = ',') %>%
        mutate(I = 1) %>%
        tidyr::complete(seg, tidyr::nesting(GT, clone, compartment, prior_clone, clone_size), fill = list('I' = 0))

    clone_post = inner_join(
            exp_post %>%
                filter(cnv_state != 'neu') %>%
                inner_join(clone_segs, by = c('seg' = 'seg')) %>%
                mutate(l_clone = ifelse(I == 1, Z_cnv, Z_n)) %>%
                group_by(cell, clone, GT, prior_clone) %>%
                summarise(
                    l_clone_x = sum(l_clone),
                    .groups = 'drop'
                ),
            allele_post %>%
                filter(cnv_state != 'neu') %>%
                inner_join(clone_segs, by = c('seg' = 'seg')) %>%
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
            Z_clone_x = log(prior_clone) + l_clone_x,
            Z_clone_y = log(prior_clone) + l_clone_y,
            p = exp(Z_clone - matrixStats::logSumExp(Z_clone)),
            p_x = exp(Z_clone_x - matrixStats::logSumExp(Z_clone_x)),
            p_y = exp(Z_clone_y - matrixStats::logSumExp(Z_clone_y))
        ) %>%
        mutate(
            clone_opt = clone[which.max(p)],
            GT_opt = GT[clone_opt],
            p_opt = p[which.max(p)]
        ) %>%
        as.data.table() %>%
        data.table::dcast(cell + clone_opt + GT_opt + p_opt ~ clone, value.var = c('p', 'p_x', 'p_y')) %>%
        as.data.frame()

    tumor_clones = clones %>% 
        filter(compartment == 'tumor') %>%
        pull(clone)

    clone_post['p_cnv'] = clone_post[paste0('p_', tumor_clones)] %>% rowSums
    clone_post['p_cnv_x'] = clone_post[paste0('p_x_', tumor_clones)] %>% rowSums
    clone_post['p_cnv_y'] = clone_post[paste0('p_y_', tumor_clones)] %>% rowSums

    clone_post = clone_post %>% mutate(compartment_opt = ifelse(p_cnv > 0.5, 'tumor', 'normal'))
    
    return(clone_post)
    
}

#' @param segs_all a dataframe containing info of all CNV segments from multiple samples
#' @return consensus CNV segments
#' @export
resolve_cnvs = function(segs_all, min_overlap = 0.5, debug = FALSE) {
            
    V = segs_all %>% ungroup() %>% mutate(vertex = 1:n(), .before = 1)

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
        filter(!(frac_overlap_x < min_overlap & frac_overlap_y < min_overlap))

    G = igraph::graph_from_data_frame(d=E, vertices=V, directed=F)

    segs_all = segs_all %>% mutate(component = igraph::components(G)$membership)

    segs_consensus = segs_all %>% group_by(component, sample) %>%
        mutate(LLR_sample = max(LLR)) %>%
        arrange(CHROM, component, -LLR_sample) %>%
        group_by(component) %>%
        filter(sample == sample[which.max(LLR_sample)])

    segs_consensus = segs_consensus %>% arrange(CHROM, seg_start) %>%
        mutate(CHROM = factor(CHROM, 1:22))
    
    if (debug) {
        return(list('G' = G, 'segs_consensus' = segs_consensus))
    }
    
    return(segs_consensus)
}

#' @description get the single cell expression likelihoods
#' @param exp_sc single-cell expression dataframe
#' @keywords internal
get_exp_likelihoods_lnpois = function(exp_sc, diploid_chroms = NULL, use_loh = FALSE, depth_obs = NULL, mu = NULL, sigma = NULL) {

    exp_sc = exp_sc %>% filter(lambda_ref > 0)
    
    if (is.null(depth_obs)){
        depth_obs = sum(exp_sc$Y_obs)
    }

    if (use_loh) {
        ref_states = c('neu', 'loh')
    } else {
        ref_states = c('neu')
    }

    if (is.null(mu) | is.null(sigma)) {

        if (!is.null(diploid_chroms)) {
            exp_sc_diploid = exp_sc %>% filter(CHROM %in% diploid_chroms)
        } else {
            exp_sc_diploid = exp_sc %>% filter(cnv_state %in% ref_states)
        }

        fit = exp_sc_diploid %>% {fit_lnpois(.$Y_obs, .$lambda_ref, depth_obs)}
        mu = fit@coef[1]
        sigma = fit@coef[2]
    }

    res = exp_sc %>% 
        filter(cnv_state != 'neu') %>%
        group_by(seg, cnv_state) %>%
        summarise(
            n = n(),
            phi_mle = calc_phi_mle_lnpois(Y_obs, lambda_ref, depth_obs, mu, sigma, lower = 0.1, upper = 10),
            l11 = l_lnpois(Y_obs, lambda_ref, depth_obs, mu, sigma, phi = 1),
            l20 = l11,
            l10 = l_lnpois(Y_obs, lambda_ref, depth_obs, mu, sigma, phi = 0.5),
            l21 = l_lnpois(Y_obs, lambda_ref, depth_obs, mu, sigma, phi = 1.5),
            l31 = l_lnpois(Y_obs, lambda_ref, depth_obs, mu, sigma, phi = 2),
            l22 = l31,
            l32 = l_lnpois(Y_obs, lambda_ref, depth_obs, mu, sigma, phi = 2.5),
            l00 = l_lnpois(Y_obs, lambda_ref, depth_obs, mu, sigma, phi = 0.25),
            mu = mu,
            sigma = sigma,
            .groups = 'drop'
        )
        
    return(res)
}

#' @description get the single cell expression dataframe
#' @param segs_consensus consensus segments
#' @param count_mat gene expression count matrix
#' @param gtf_transcript transcript dataframe
#' @return single cell expression dataframe
#' @keywords internal
get_exp_sc = function(segs_consensus, count_mat, gtf_transcript) {

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
        setNames(c('gene_index', 'seg_index')) %>%
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

    return(exp_sc)
}


#' @description compute single-cell expression posteriors
#' @param segs_consensus consensus segments
#' @param count_mat gene expression count matrix
#' @param gtf_transcript transcript dataframe
#' @param lambdas_ref reference expression matrix
#' @return expression posterior datatframe
#' @keywords internal
get_exp_post = function(segs_consensus, count_mat, gtf_transcript, lambdas_ref = NULL, diploid_chroms = NULL, use_loh = NULL, ncores = 30, verbose = TRUE, debug = F) {

    exp_sc = get_exp_sc(segs_consensus, count_mat, gtf_transcript) 

    if (is.null(use_loh)) {
        if (mean(exp_sc$cnv_state == 'neu') < 0.05) {
            use_loh = TRUE
            log_info('less than 5% genes are in neutral region - including LOH in baseline')
        } else {
            use_loh = FALSE
        }
    } else if (use_loh) {
        log_info('Including LOH in baseline as specified')
    }
    
    cells = colnames(count_mat)

    if ((!is.matrix(lambdas_ref))) {
        lambdas_ref = as.matrix(lambdas_ref) %>% set_colnames('ref')
        best_refs = setNames(rep('ref', length(cells)), cells)
    } else {
        best_refs = choose_ref_cor(count_mat, lambdas_ref, gtf_transcript)
    }

    results = mclapply(
        cells,
        mc.cores = ncores,
        function(cell) {
   
            ref = best_refs[cell]

            exp_sc = exp_sc[,c('gene', 'seg', 'CHROM', 'cnv_state', 'seg_start', 'seg_end', cell)] %>%
                rename(Y_obs = ncol(.))

            exp_sc %>%
                mutate(
                    lambda_ref = lambdas_ref[, ref][gene],
                    lambda_obs = Y_obs/sum(Y_obs),
                    logFC = log2(lambda_obs/lambda_ref)
                ) %>%
                get_exp_likelihoods_lnpois(
                    use_loh = use_loh,
                    diploid_chroms = diploid_chroms
                ) %>%
                mutate(cell = cell, ref = ref)

        }
    )

    bad = sapply(results, inherits, what = "try-error")

    if (any(bad)) {
        if (verbose) {log_warn(glue('{sum(bad)} jobs failed'))}
        log_warn(results[bad][1])
        log_warn(cells[bad][1])
    } else {
        log_info('All cells succeeded')
    }
    
    exp_post = results[!bad] %>%
        bind_rows() %>%
        mutate(seg = factor(seg, gtools::mixedsort(unique(seg)))) %>%
        rowwise() %>%
        left_join(
            segs_consensus %>% select(
                CHROM, seg = seg_cons, 
                prior_loh = p_loh, prior_amp = p_amp, prior_del = p_del, prior_bamp = p_bamp, prior_bdel = p_bdel
            ),
            by = 'seg'
        ) %>%
        # if the opposite state has a very small prior, and phi is in the opposite direction, then CNV posterior can still be high which is miselading
        mutate_at(
            vars(contains('prior')),
            function(x) {ifelse(x < 0.05, 0, x)}
        ) %>%
        compute_posterior() %>%
        mutate(seg_label = paste0(seg, '(', cnv_state, ')')) %>%
        mutate(seg_label = factor(seg_label, unique(seg_label)))
    
    return(list('exp_post' = exp_post, 'exp_sc' = exp_sc, 'best_refs' = best_refs))
}

#' @description do bayesian averaging to get posteriors
#' @keywords internal
compute_posterior = function(sc_post) {
    sc_post %>% 
    rowwise() %>%
    mutate(
        Z_amp = matrixStats::logSumExp(c(l21 + log(prior_amp/4), l31 + log(prior_amp/4))),
        Z_loh = l20 + log(prior_loh/2),
        Z_del = l10 + log(prior_del/2),
        Z_bamp = l22 + log(prior_bamp/2),
        Z_bdel = l00 + log(prior_bdel/2),
        Z_n = l11 + log(1/2),
        Z = matrixStats::logSumExp(
            c(Z_n, Z_loh, Z_del, Z_amp, Z_bamp, Z_bdel)
        ),
        Z_cnv = matrixStats::logSumExp(
            c(Z_loh, Z_del, Z_amp, Z_bamp, Z_bdel)
        ),
        p_amp = exp(Z_amp - Z),
        p_neu = exp(Z_n - Z),
        p_del = exp(Z_del - Z),
        p_loh = exp(Z_loh - Z),
        p_bamp = exp(Z_bamp - Z),
        p_bdel = exp(Z_bdel - Z),
        logBF = Z_cnv - Z_n,
        p_cnv = exp(Z_cnv - Z),
        p_n = exp(Z_n - Z)
    ) %>%
    ungroup()
}

#' get CNV allele posteriors
#' @keywords internal
get_allele_post = function(bulk_all, segs_consensus, df_allele, naive = FALSE) {

    if ((!'sample' %in% colnames(bulk_all)) | (!'sample' %in% colnames(segs_consensus))) {
        bulk_all['sample'] = '0'
        segs_consensus['sample'] = '0'
        # warning('Sample column missing')
    }

    if (naive) {
        bulk_all = bulk_all %>% mutate(haplo_post = ifelse(AR >= 0.5, 'major', 'minor'))
    }
    
    # allele posteriors
    snp_seg = bulk_all %>%
        filter(!is.na(pAD)) %>%
        select(snp_id, snp_index, sample, seg, haplo_post) %>%
        inner_join(
            segs_consensus,
            by = c('sample', 'seg')
        )

    allele_sc = df_allele %>%
        mutate(pAD = ifelse(GT == '1|0', AD, DP - AD)) %>%
        select(-snp_index) %>% 
        inner_join(
            snp_seg %>% select(snp_id, snp_index, haplo_post, seg = seg_cons, cnv_state),
            by = c('snp_id')
        ) %>%
        filter(!cnv_state %in% c('neu', 'bamp', 'bdel')) %>%
        mutate(
            major_count = ifelse(haplo_post == 'major', AD, DP - AD),
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
        group_by(cell, CHROM, seg, cnv_state) %>%
        summarise(
            major = sum(major_count),
            minor = sum(minor_count),
            total = major + minor,
            MAF = major/total,
            .groups = 'drop'
        ) %>%
        left_join(
            segs_consensus %>% select(seg = seg_cons, 
                prior_loh = p_loh, prior_amp = p_amp, prior_del = p_del, prior_bamp = p_bamp, prior_bdel = p_bdel
            ),
            by = 'seg'
        ) %>%
        rowwise() %>%
        mutate(
            l11 = dbinom(major, total, p = 0.5, log = TRUE),
            l10 = dbinom(major, total, p = 0.9, log = TRUE),
            l01 = dbinom(major, total, p = 0.1, log = TRUE),
            l20 = dbinom(major, total, p = 0.9, log = TRUE),
            l02 = dbinom(major, total, p = 0.1, log = TRUE),
            l21 = dbinom(major, total, p = 0.66, log = TRUE),
            l12 = dbinom(major, total, p = 0.33, log = TRUE),
            l31 = dbinom(major, total, p = 0.75, log = TRUE),
            l13 = dbinom(major, total, p = 0.25, log = TRUE),
            l32 = dbinom(major, total, p = 0.6, log = TRUE),
            l22 = l11,
            l00 = l11
        ) %>%
        ungroup() %>%
        compute_posterior() %>%
        mutate(seg_label = paste0(seg, '(', cnv_state, ')')) %>%
        mutate(seg_label = factor(seg_label, unique(seg_label)))

        return(allele_post)
}

#' get joint posteriors
#' @keywords internal
get_joint_post = function(exp_post, allele_post, segs_consensus) {
    
    joint_post = full_join(
            exp_post %>%
                filter(cnv_state != 'neu') %>%
                select(
                    cell, CHROM, seg, cnv_state, 
                    l11_x = l11, l20_x = l20, l10_x = l10, l21_x = l21, l31_x = l31, l22_x = l22, l00_x = l00,
                    Z_x = Z, Z_cnv_x = Z_cnv, Z_n_x = Z_n, logBF_x = logBF
                ),
            allele_post %>% 
                select(
                    cell, seg, 
                    l11_y = l11, l20_y = l20, l10_y = l10, l21_y = l21, l31_y = l31, l22_y = l22, l00_y = l00,
                    Z_y = Z, Z_cnv_y = Z_cnv, Z_n_y = Z_n, logBF_y = logBF,
                    n_snp = total, MAF, major, total
            ),
            c("cell", "seg")
        ) %>%
        mutate(cnv_state = tidyr::replace_na(cnv_state, 'loh')) %>%
        mutate_at(
            vars(matches("_x|_y")),
            function(x) tidyr::replace_na(x, 0)
        ) %>%
        left_join(
            segs_consensus %>% select(
                seg = seg_cons,
                any_of(c('n_genes', 'n_snps', 'prior_loh' = 'p_loh', 'prior_amp' = 'p_amp', 'prior_del' = 'p_del', 'prior_bamp' = 'p_bamp', 'prior_bdel' = 'p_bdel', 'LLR', 'LLR_x', 'LLR_y'))
            ),
            by = 'seg'
        ) %>%
        mutate(
            l11 = l11_x + l11_y,
            l20 = l20_x + l20_y,
            l10 = l10_x + l10_y,
            l21 = l21_x + l21_y,
            l31 = l31_x + l31_y,
            l22 = l22_x + l22_y,
            l00 = l00_x + l00_y,
        ) %>%
        compute_posterior() %>%
        mutate(
            p_cnv_x = 1/(1+exp(-logBF_x)),
            p_cnv_y = 1/(1+exp(-logBF_y))
        ) %>%
        rowwise() %>%
        mutate(
            cnv_state_mle = c('neu', 'loh', 'del', 'amp', 'amp', 'bamp')[which.max(c(l11, l20, l10, l21, l31, l22))],
            cnv_state_map = c('neu', 'loh', 'del', 'amp', 'bamp')[which.max(c(p_neu, p_loh, p_del, p_amp, p_bamp))],
        ) %>%
        ungroup()

    joint_post = joint_post %>% 
        mutate(seg = factor(seg, gtools::mixedsort(unique(seg)))) %>%
        mutate(seg_label = paste0(seg, '(', cnv_state, ')')) %>%
        mutate(seg_label = factor(seg_label, unique(seg_label)))
    
    return(joint_post)
}

#' test for multi-allelic CNVs
#' @keywords internal 
test_multi_allelic = function(bulks, segs_consensus, use_loh = FALSE, LLR_min = 100, p_min = 0.995, diploid_chroms = NULL) {

    log_info('Testing for multi-allelic CNVs ..')

    # default
    if (is.null(use_loh)) {
        length_neu = segs_consensus %>% filter(cnv_state == 'neu') %>% pull(seg_length) %>% sum
        if (length_neu > 1.5e8) {
            use_loh = TRUE
            log_info('less than 5% of genome is in neutral region - including LOH in baseline')
        } else {
            use_loh = FALSE
        }
    }

    if (use_loh) {
        ref_states = c('neu', 'loh')
    } else {
        ref_states = c('neu')
    }
    
    bulks = bulks %>% annot_consensus(segs_consensus)
    
    if (!is.null(diploid_chroms)) {
        bulks = bulks %>% mutate(diploid = CHROM %in% diploid_chroms)
    } else {
        bulks = bulks %>% mutate(diploid = cnv_state %in% ref_states)
    }
    
    bulks = bulks %>% 
        run_group_hmms(run_hmm = F) %>%
        mutate(state_post = ifelse(LLR < LLR_min | is.na(LLR), 'neu', state_post))
    
    segs_multi = bulks %>% 
        distinct(sample, CHROM, seg_cons, LLR, p_amp, p_del, p_loh, p_bamp, cnv_state_post) %>%
        rowwise() %>%
        mutate(p_max = max(c(p_amp, p_del, p_loh, p_bamp))) %>%
        filter(LLR > LLR_min & p_max > p_min) %>%
        group_by(seg_cons) %>%
        summarise(
            cnv_states = list(sort(unique(cnv_state_post))),
            n_states = length(unlist(cnv_states))
        ) %>%
        filter(n_states > 1)

    segs = segs_multi$seg_cons

    log_info(glue('{length(segs)} multi-allelic CNVs found: {paste(segs, collapse = ",")}'))

    if (length(segs) > 0) {
        segs_consensus = segs_consensus %>%
            left_join(
                segs_multi,
                by = 'seg_cons'
            ) %>%
            rowwise() %>%
            mutate(
                cnv_states = ifelse(is.null(cnv_states), list(cnv_state_post), list(cnv_states)),
                n_states = sum(cnv_states != 'neu')
            ) %>%
            mutate(
                p_del = ifelse(n_states > 1, ifelse('del' %in% cnv_states, 0.5, 0), p_del),
                p_amp = ifelse(n_states > 1, ifelse('amp' %in% cnv_states, 0.5, 0), p_amp),
                p_loh = ifelse(n_states > 1, ifelse('loh' %in% cnv_states, 0.5, 0), p_loh),
                p_bamp = ifelse(n_states > 1, ifelse('bamp' %in% cnv_states, 0.5, 0), p_bamp),
                p_bdel = ifelse(n_states > 1, ifelse('bdel' %in% cnv_states, 0.5, 0), p_bdel)
            ) %>%
            ungroup()
    } else {
        segs_consensus = segs_consensus %>% 
            mutate(
                n_states = ifelse(cnv_state == 'neu', 1, 0),
                cnv_states = cnv_state
            )
    }
    
    return(list('segs_consensus' = segs_consensus, 'bulks' = bulks, 'segs_multi' = segs_multi))
}

#' expand multi-allelic CNVs into separate entries in the single-cell posterior dataframe
#' @keywords internal 
expand_states = function(sc_post, segs_consensus) {

    if (any(segs_consensus$n_states > 1)) {
        sc_post = sc_post %>%
            left_join(
                segs_consensus %>% select(seg = seg_cons, cnv_states, n_states),
                by = 'seg'
            ) %>%
            reshape2::melt(
                measure.vars = c('p_amp', 'p_loh', 'p_del', 'p_bamp', 'p_bdel'),
                variable.name = 'cnv_state_expand',
                value.name = 'p_cnv_expand'
            ) %>%
            mutate(
                cnv_state_expand = c('p_amp' = 'amp', 'p_loh' = 'loh', 'p_del' = 'del', 'p_bamp' = 'bamp', 'p_bdel' = 'bdel')[cnv_state_expand]
            ) %>%
            rowwise() %>%
            filter(cnv_state_expand %in% cnv_states) %>%
            ungroup() %>%
            mutate(
                seg = ifelse(n_states > 1, paste0(seg, '_', cnv_state_expand), seg),
                p_cnv = ifelse(n_states > 1, p_cnv_expand, p_cnv),
                cnv_state = ifelse(n_states > 1, cnv_state_expand, cnv_state)
            ) %>%
            mutate(
                Z_cnv = case_when(
                    n_states <= 1 ~ Z_cnv,
                    cnv_state == 'loh' ~ Z_loh,
                    cnv_state == 'del' ~ Z_del,
                    cnv_state == 'bamp' ~ Z_bamp,
                    cnv_state == 'amp' ~ Z_amp,
                    cnv_state == 'bdel' ~ Z_bdel
                )
            ) %>%
            ungroup()
    } else {
        log_info('No multi-allelic CNVs, skipping ..')
    }

    return(sc_post)
}


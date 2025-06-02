#' @import logger
#' @import dplyr
#' @import Matrix
#' @importFrom data.table fread fwrite as.data.table
#' @import stringr
#' @import glue
#' @importFrom parallel mclapply
#' @import tidygraph
#' @import ggplot2
#' @import ggraph
#' @importFrom scistreer perform_nni get_mut_graph score_tree annotate_tree mut_to_tree ladderize to_phylo
#' @importFrom ggtree %<+%
#' @importFrom methods is as
#' @importFrom igraph vcount ecount E V V<- E<-
#' @import patchwork
#' @importFrom grDevices colorRampPalette
#' @importFrom stats as.dendrogram as.dist cor cutree dbinom dnbinom dnorm dpois end hclust integrate model.matrix na.omit optim p.adjust pnorm reorder rnorm setNames start t.test as.ts complete.cases is.leaf na.contiguous
#' @importFrom hahmmr logSumExp dpoilog dbbinom l_lnpois l_bbinom fit_lnpois_cpp likelihood_allele forward_back_allele run_joint_hmm_s15 run_allele_hmm_s5
#' @import tibble
#' @importFrom utils combn
#' @useDynLib numbat, .registration=TRUE
NULL

#' Run workflow to decompose tumor subclones
#'
#' @param count_mat dgCMatrix Raw count matrices where rownames are genes and column names are cells
#' @param lambdas_ref matrix Either a named vector with gene names as names and normalized expression as values, or a matrix where rownames are genes and columns are pseudobulk names
#' @param df_allele dataframe Allele counts per cell, produced by preprocess_allele
#' @param genome character Genome version (hg38, hg19, or mm10) 
#' @param out_dir string Output directory
#' @param gamma numeric Dispersion parameter for the Beta-Binomial allele model
#' @param t numeric Transition probability
#' @param init_k integer Number of clusters in the initial clustering
#' @param min_cells integer Minimum number of cells to run HMM on
#' @param min_genes integer Minimum number of genes to call a segment
#' @param max_cost numeric Likelihood threshold to collapse internal branches
#' @param n_cut integer Number of cuts on the phylogeny to define subclones
#' @param tau numeric Factor to determine max_cost as a function of the number of cells (0-1)
#' @param nu numeric Phase switch rate
#' @param alpha numeric P value cutoff for diploid finding
#' @param min_overlap numeric Minimum CNV overlap threshold
#' @param max_iter integer Maximum number of iterations to run the phyologeny optimization
#' @param max_nni integer Maximum number of iterations to run NNI in the ML phylogeny inference
#' @param min_depth integer Minimum allele depth 
#' @param common_diploid logical Whether to find common diploid regions in a group of peusdobulks
#' @param ncores integer Number of threads to use 
#' @param ncores_nni integer Number of threads to use for NNI
#' @param max_entropy numeric Entropy threshold to filter CNVs
#' @param min_LLR numeric Minimum LLR to filter CNVs
#' @param eps numeric Convergence threshold for ML tree search
#' @param multi_allelic logical Whether to call multi-allelic CNVs
#' @param p_multi numeric P value cutoff for calling multi-allelic CNVs
#' @param use_loh logical Whether to include LOH regions in the expression baseline
#' @param skip_nj logical Whether to skip NJ tree construction and only use UPGMA
#' @param diploid_chroms vector Known diploid chromosomes
#' @param segs_loh dataframe Segments of clonal LOH to be excluded
#' @param call_clonal_loh logical Whether to call segments with clonal LOH
#' @param segs_consensus_fix dataframe Pre-determined segmentation of consensus CNVs
#' @param check_convergence logical Whether to terminate iterations based on consensus CNV convergence 
#' @param random_init logical Whether to initiate phylogney using a random tree (internal use only)
#' @param exclude_neu logical Whether to exclude neutral segments from CNV retesting (internal use only)
#' @param plot logical Whether to plot results
#' @param verbose logical Verbosity
#' @return a status code
#' @export
run_numbat = function(
        count_mat, lambdas_ref, df_allele,gtf=NULL,genome = 'hg38', 
        out_dir = tempdir(), max_iter = 2, max_nni = 100, t = 1e-5, gamma = 20, min_LLR = 5,
        alpha = 1e-4, eps = 1e-5, max_entropy = 0.5, init_k = 3, min_cells = 50, tau = 0.3, nu = 1,
        max_cost = ncol(count_mat) * tau, n_cut = 0, min_depth = 0, common_diploid = TRUE, min_overlap = 0.45, 
        ncores = 1, ncores_nni = ncores, random_init = FALSE, segs_loh = NULL, call_clonal_loh = FALSE, 
        verbose = TRUE, diploid_chroms = NULL, segs_consensus_fix = NULL, use_loh = NULL, min_genes = 10,
        skip_nj = FALSE, multi_allelic = TRUE, p_multi = 1-alpha, 
        plot = TRUE, check_convergence = FALSE, exclude_neu = TRUE
    ) {

    ######### Setup output folder #########
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    logfile = glue('{out_dir}/log.txt')
    if (file.exists(logfile)) {file.remove(logfile)}
    log_appender(appender_file(logfile))

    ######### Basic checks #########
    if(is.null(gtf)){
       if (genome == 'hg38') {
        gtf = gtf_hg38
    } else if (genome == 'hg19') {
        gtf = gtf_hg19
    } else if (genome == 'mm10') {
        gtf = gtf_mm10
    } else {
        stop('Genome version must be hg38, hg19, or mm10')
    } 
    }
    

    count_mat = check_matrix(count_mat)
    df_allele = annotate_genes(df_allele, gtf)
    df_allele = check_allele_df(df_allele)
    lambdas_ref = check_exp_ref(lambdas_ref)

    # filter for annotated genes
    genes_annotated = unique(gtf$gene) %>% 
        intersect(rownames(count_mat)) %>%
        intersect(rownames(lambdas_ref))

    count_mat = count_mat[genes_annotated,,drop=FALSE]
    lambdas_ref = lambdas_ref[genes_annotated,,drop=FALSE]

    zero_cov = names(which(colSums(count_mat) == 0))
    if (length(zero_cov) > 0) {
        log_message(glue('Filtering out {length(zero_cov)} cells with 0 coverage'))
        count_mat = count_mat[,!colnames(count_mat) %in% zero_cov]
        df_allele = df_allele %>% filter(!cell %in% zero_cov)
    }

    # only keep cells that have a transcriptome
    df_allele = df_allele %>% filter(cell %in% colnames(count_mat))
    if (nrow(df_allele) == 0){
        stop('No matching cell names between count_mat and df_allele')
    }

    if ((!is.null(segs_loh)) & (!is.null(segs_consensus_fix))) {
        stop('Cannot specify both segs_loh and segs_consensus_fix')
    }

    # check provided consensus CNVs
    segs_consensus_fix = check_segs_fix(segs_consensus_fix)

    # check clonal LOH
    if (!is.null(segs_loh)) {
        if (call_clonal_loh) {
            stop('Cannot specify both segs_loh and call_clonal_loh')
        }
        segs_loh = check_segs_loh(segs_loh)
    }

    ######### Log parameters #########
    log_message(paste('\n',
        glue('numbat version: ', as.character(utils::packageVersion("numbat"))),
        glue('scistreer version: ', as.character(utils::packageVersion("scistreer"))),
        glue('hahmmr version: ', as.character(utils::packageVersion("hahmmr"))),
        'Running under parameters:',
        glue('t = {t}'), 
        glue('alpha = {alpha}'),
        glue('gamma = {gamma}'),
        glue('min_cells = {min_cells}'), 
        glue('init_k = {init_k}'),
        glue('max_cost = {max_cost}'),
        glue('n_cut = {n_cut}'),
        glue('max_iter = {max_iter}'),
        glue('max_nni = {max_nni}'),
        glue('min_depth = {min_depth}'),
        glue('use_loh = {ifelse(is.null(use_loh), "auto", use_loh)}'),
        glue('segs_loh = {ifelse(is.null(segs_loh), "None", "Given")}'),
        glue('call_clonal_loh = {call_clonal_loh}'),
        glue('segs_consensus_fix = {ifelse(is.null(segs_consensus_fix), "None", "Given")}'),
        glue('multi_allelic = {multi_allelic}'),
        glue('min_LLR = {min_LLR}'),
        glue('min_overlap = {min_overlap}'),
        glue('max_entropy = {max_entropy}'),
        glue('skip_nj = {skip_nj}'),
        glue('diploid_chroms = {ifelse(is.null(diploid_chroms), "None", "Given")}'),
        glue('ncores = {ncores}'),
        glue('ncores_nni = {ncores_nni}'),
        glue('common_diploid = {common_diploid}'),
        glue('tau = {tau}'),
        glue('check_convergence = {check_convergence}'),
        glue('plot = {plot}'),
        glue('genome = {genome}'),
        'Input metrics:',
        glue('{ncol(count_mat)} cells'),
        sep = "\n"
    ), verbose = verbose)

    log_mem()

    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    data.table::setDTthreads(1)

    ######## Initialization ########

    if (call_clonal_loh) {
        
        log_message('Calling segments with clonal LOH')

        bulk = get_bulk(
            count_mat = count_mat,
            lambdas_ref = lambdas_ref,
            df_allele = df_allele,
            gtf = gtf,
            min_depth = min_depth,
            nu = nu
        )

        segs_loh = bulk %>% detect_clonal_loh(t = t, min_depth = min_depth)

        if (!is.null(segs_loh)) {
            fwrite(segs_loh, glue('{out_dir}/segs_loh.tsv'), sep = '\t')
        } else {
            log_message('No segments with clonal LOH detected')
        }
    }

    sc_refs = choose_ref_cor(count_mat, lambdas_ref, gtf)
    saveRDS(sc_refs, glue('{out_dir}/sc_refs.rds'))

    if (random_init) {

        log_message('Initializing with random tree...')
        n_cells = ncol(count_mat)
        dist_mat = matrix(rnorm(n_cells^2),nrow=n_cells)
        colnames(dist_mat) = colnames(count_mat)
        hc = hclust(as.dist(dist_mat), method = "ward.D2")

        saveRDS(hc, glue('{out_dir}/hc.rds'))

        # extract cell groupings
        subtrees = get_nodes_celltree(hc, cutree(hc, k = init_k))

    } else if (init_k == 1) {

        log_message('Initializing with all-cell pseudobulk ..', verbose = verbose)
        log_mem()

        subtrees = list(list('cells' = colnames(count_mat), 'size' = ncol(count_mat), 'sample' = 1))

    } else {

        log_message('Approximating initial clusters using smoothed expression ..', verbose = verbose)
        log_mem()

        clust = exp_hclust(
            count_mat = count_mat,
            lambdas_ref = lambdas_ref,
            gtf = gtf,
            sc_refs = sc_refs,
            ncores = ncores
        )

        hc = clust$hc

        fwrite(
            as.data.frame(clust$gexp_roll_wide) %>% tibble::rownames_to_column('cell'),
            glue('{out_dir}/gexp_roll_wide.tsv.gz'),
            sep = '\t',
            nThread = min(4, ncores)
        )

        saveRDS(hc, glue('{out_dir}/hc.rds'))

        # extract cell groupings
        subtrees = get_nodes_celltree(hc, cutree(hc, k = init_k))

        if (plot) {

            p = plot_exp_roll(
                gexp_roll_wide = clust$gexp_roll_wide,
                hc = hc,
                k = init_k,
                gtf = gtf,
                n_sample = 1e4
            )
        
            ggsave(glue('{out_dir}/exp_roll_clust.png'), p, width = 8, height = 4, dpi = 200)

        }

    }

    clones = purrr::keep(subtrees, function(x) x$sample %in% 1:init_k)

    normal_cells = c()
    segs_consensus_old = data.frame()

    ######## Begin iterations ########
    for (i in 1:max_iter) {

        log_message(glue('Iteration {i}'), verbose = verbose)
        log_mem()

        subtrees = purrr::keep(subtrees, function(x) x$size > min_cells)

        bulk_subtrees = make_group_bulks(
                groups = subtrees,
                count_mat = count_mat,
                df_allele = df_allele, 
                lambdas_ref = lambdas_ref,
                gtf = gtf,
                min_depth = min_depth,
                nu = nu,
                segs_loh = segs_loh,
                ncores = ncores)

        # diagnostics
        if (i == 1) {
            bulk_subtrees %>% filter(sample == 0) %>% check_contam()
            bulk_subtrees %>% filter(sample == 0) %>% check_exp_noise()
        }

        if (is.null(segs_consensus_fix)) {

            bulk_subtrees = bulk_subtrees %>%
                run_group_hmms(
                    t = t,
                    gamma = gamma,
                    alpha = alpha,
                    nu = nu,
                    min_genes = min_genes,
                    common_diploid = common_diploid,
                    diploid_chroms = diploid_chroms,
                    ncores = ncores,
                    verbose = verbose)

            fwrite(bulk_subtrees, glue('{out_dir}/bulk_subtrees_{i}.tsv.gz'), sep = '\t')
            
            if (plot) {
                p = plot_bulks(bulk_subtrees, min_LLR = min_LLR, use_pos = TRUE, genome = genome)
                ggsave(
                    glue('{out_dir}/bulk_subtrees_{i}.png'), p, 
                    width = 13, height = 2*length(unique(bulk_subtrees$sample)), dpi = 250
                )
            }

            # define consensus CNVs
            segs_consensus = bulk_subtrees %>% 
                get_segs_consensus(min_LLR = min_LLR, min_overlap = min_overlap, retest = TRUE)
            
            # check termination
            if (all(segs_consensus$cnv_state_post == 'neu')) {
                msg = 'No CNV remains after filtering by LLR in pseudobulks. Consider reducing min_LLR.'
                log_message(msg)
                return(msg)
            }

            # retest all segments on subtrees
            bulk_subtrees = retest_bulks(
                    bulk_subtrees,
                    segs_consensus, 
                    diploid_chroms = diploid_chroms, 
                    gamma = gamma,
                    min_LLR = min_LLR,
                    ncores = ncores
                )

            fwrite(bulk_subtrees, glue('{out_dir}/bulk_subtrees_retest_{i}.tsv.gz'), sep = '\t')
            
            # find consensus CNVs again
            segs_consensus = bulk_subtrees %>%
                get_segs_consensus(min_LLR = min_LLR, min_overlap = min_overlap, retest = FALSE)

            # check termination again
            if (all(segs_consensus$cnv_state_post == 'neu')) {
                msg = 'No CNV remains after filtering by LLR in pseudobulks. Consider reducing min_LLR.'
                log_message(msg)
                return(msg)
            }

        } else {

            log_message('Using fixed consensus CNVs')
            segs_consensus = segs_consensus_fix

            bulk_subtrees = bulk_subtrees %>% 
                annot_consensus(segs_consensus) %>%
                annot_theta_mle() %>%
                classify_alleles()

        }

        # retest on clones
        clones = purrr::keep(clones, function(x) x$size > min_cells)

        if (length(clones) == 0) {
            msg = 'No clones remain after filtering by size. Consider reducing min_cells.'
            log_message(msg)
            return(msg)
        }

        bulk_clones = make_group_bulks(
                groups = clones,
                count_mat = count_mat,
                df_allele = df_allele, 
                lambdas_ref = lambdas_ref,
                gtf = gtf,
                min_depth = min_depth,
                nu = nu,
                segs_loh = segs_loh,
                ncores = ncores)

        bulk_clones = bulk_clones %>% 
            run_group_hmms(
                t = t,
                gamma = gamma,
                alpha = alpha,
                nu = nu,
                min_genes = min_genes,
                common_diploid = common_diploid,
                diploid_chroms = diploid_chroms,
                ncores = ncores,
                verbose = verbose,
                retest = FALSE)

        bulk_clones = retest_bulks(
            bulk_clones,
            segs_consensus,
            gamma = gamma,
            use_loh = use_loh,
            min_LLR = min_LLR,
            diploid_chroms = diploid_chroms,
            ncores = ncores)
        
        fwrite(bulk_clones, glue('{out_dir}/bulk_clones_{i}.tsv.gz'), sep = '\t')

        if (plot) {
            p = plot_bulks(bulk_clones, min_LLR = min_LLR, use_pos = TRUE, genome = genome)
            ggsave(
                glue('{out_dir}/bulk_clones_{i}.png'), p, 
                width = 13, height = 2*length(unique(bulk_clones$sample)), dpi = 250
            )
        }

        # test for multi-allelic CNVs
        if (multi_allelic) {
            segs_consensus = test_multi_allelic(bulk_clones, segs_consensus, min_LLR = min_LLR, p_min = p_multi)
        }

        fwrite(segs_consensus, glue('{out_dir}/segs_consensus_{i}.tsv'), sep = '\t')
        
        ######## Evaluate CNV per cell ########
        log_message('Evaluating CNV per cell ..', verbose = verbose)
        log_mem()

        exp_post = get_exp_post(
            segs_consensus %>% mutate(cnv_state = ifelse(cnv_state == 'neu', cnv_state, cnv_state_post)),
            count_mat,
            lambdas_ref,
            use_loh = use_loh,
            segs_loh = segs_loh,
            gtf = gtf,
            sc_refs = sc_refs,
            ncores = ncores)

        haplotype_post = get_haplotype_post(
            bulk_subtrees,
            segs_consensus %>% mutate(cnv_state = ifelse(cnv_state == 'neu', cnv_state, cnv_state_post))
        )

        allele_post = get_allele_post(
            df_allele, 
            haplotype_post, 
            segs_consensus %>% mutate(cnv_state = ifelse(cnv_state == 'neu', cnv_state, cnv_state_post))
        )

        joint_post = get_joint_post(
            exp_post,
            allele_post,
            segs_consensus)

        joint_post = joint_post %>%
            group_by(seg) %>%
            mutate(
                avg_entropy = mean(binary_entropy(p_cnv), na.rm = TRUE)
            ) %>%
            ungroup()

        if (multi_allelic) {
            log_message('Expanding allelic states..', verbose = verbose)
            exp_post = expand_states(exp_post, segs_consensus)
            allele_post = expand_states(allele_post, segs_consensus)
            joint_post = expand_states(joint_post, segs_consensus)
        }

        fwrite(exp_post, glue('{out_dir}/exp_post_{i}.tsv'), sep = '\t')
        fwrite(allele_post, glue('{out_dir}/allele_post_{i}.tsv'), sep = '\t')
        fwrite(joint_post, glue('{out_dir}/joint_post_{i}.tsv'), sep = '\t')

        ######## Build phylogeny ########
        log_message('Building phylogeny ..', verbose = verbose)
        log_mem()

        # filter CNVs
        joint_post_filtered = joint_post %>%
            filter(cnv_state != 'neu') %>%
            filter(avg_entropy < max_entropy & LLR > min_LLR)

        if (nrow(joint_post_filtered) == 0) {
            msg = 'No CNV remains after filtering by entropy in single cells. Consider increasing max_entropy.'
            log_message(msg)
            return(msg)
        } else {
            n_cnv = length(unique(joint_post_filtered$seg))
            log_message(glue('Using {n_cnv} CNVs to construct phylogeny'), verbose = verbose)
        }

        # construct genotype probability matrix
        p_min = 1e-10

        P = joint_post_filtered %>%
            mutate(p_cnv = pmax(pmin(p_cnv, 1-p_min), p_min)) %>%
            as.data.table %>%
            data.table::dcast(cell ~ seg, value.var = 'p_cnv', fill = 0.5) %>%
            tibble::column_to_rownames('cell') %>%
            as.matrix

        fwrite(
            as.data.frame(P) %>% tibble::rownames_to_column('cell'),
            glue('{out_dir}/geno_{i}.tsv'),
            sep = '\t'
        )

        # contruct initial tree
        dist_mat = parallelDist::parDist(rbind(P, 'outgroup' = 0), threads = ncores)

        treeUPGMA = upgma(dist_mat) %>%
            ape::root(outgroup = 'outgroup') %>%
            ape::drop.tip('outgroup') %>%
            reorder(order = 'postorder')

        saveRDS(treeUPGMA, glue('{out_dir}/treeUPGMA_{i}.rds'))

        UPGMA_score = score_tree(treeUPGMA, as.matrix(P))$l_tree
        tree_init = treeUPGMA

        # note that dist_mat gets modified by NJ
        if(!skip_nj){
            treeNJ = ape::nj(dist_mat) %>%
                ape::root(outgroup = 'outgroup') %>%
                ape::drop.tip('outgroup') %>%
                reorder(order = 'postorder')
            NJ_score = score_tree(treeNJ, as.matrix(P))$l_tree

            if (UPGMA_score > NJ_score) {
                log_message('Using UPGMA tree as seed..', verbose = verbose)
            } else {
                tree_init = treeNJ
                log_message('Using NJ tree as seed..', verbose = verbose)
            }
            saveRDS(treeNJ, glue('{out_dir}/treeNJ_{i}.rds'))
        } else {
            log_message('Only computing UPGMA..', verbose = verbose)
            log_message('Using UPGMA tree as seed..', verbose = verbose)
        }
        log_mem()

        # maximum likelihood tree search with NNI
        tree_list = perform_nni(tree_init, P, ncores = ncores_nni, eps = eps, max_iter = max_nni)
        saveRDS(tree_list, glue('{out_dir}/tree_list_{i}.rds'))
        treeML = tree_list[[length(tree_list)]]
        saveRDS(treeML, glue('{out_dir}/treeML_{i}.rds'))

        gtree = get_gtree(treeML, P, n_cut = n_cut, max_cost = max_cost)
        saveRDS(gtree, glue('{out_dir}/tree_final_{i}.rds'))
        G_m = get_mut_graph(gtree) %>% label_genotype()
        saveRDS(G_m, glue('{out_dir}/mut_graph_{i}.rds'))

        clone_post = get_clone_post(gtree, exp_post, allele_post)
        fwrite(clone_post, glue('{out_dir}/clone_post_{i}.tsv'), sep = '\t')

        normal_cells = clone_post %>% filter(p_cnv < 0.5) %>% pull(cell)

        log_message(glue('Found {length(normal_cells)} normal cells..'), verbose = verbose)

        if (plot) {

            panel = plot_phylo_heatmap(
                gtree,
                joint_post,
                segs_consensus,
                clone_post,
                tip_length = 0.2,
                branch_width = 0.2,
                line_width = 0.1,
                clone_bar = TRUE
            )
        
            ggsave(glue('{out_dir}/panel_{i}.png'), panel, width = 7.5, height = 3.75, dpi = 250)

        }

        # form cell groupings using the obtained phylogeny
        clone_to_node = setNames(V(G_m)$id, V(G_m)$clone)

        subtrees = lapply(1:vcount(G_m), function(c) {

            nodes = na.omit(igraph::dfs(G_m, root = c, unreachable = F)$order)
            
            G_m %>% 
            igraph::as_data_frame('vertices') %>%
            filter(id %in% nodes) %>%
            inner_join(clone_post, by = c('GT' = 'GT_opt')) %>%
            {list(sample = c, members = unique(.$GT), clones = unique(.$clone), cells = .$cell, size = length(.$cell))}
        })

        saveRDS(subtrees, glue('{out_dir}/subtrees_{i}.rds'))

        clones = clone_post %>% split(.$clone_opt) %>%
            purrr::map(function(c){list(sample = unique(c$clone_opt), members = unique(c$GT_opt), cells = c$cell, size = length(c$cell))})

        saveRDS(clones, glue('{out_dir}/clones_{i}.rds'))

        #### check convergence ####
        if (check_convergence) {

            converge = segs_equal(segs_consensus_old, segs_consensus)

            if (converge) {
                log_message('Convergence reached')
                break
            } else {
                segs_consensus_old = segs_consensus
            }

        }
    }

    # Output final subclone bulk profiles - logic can be simplified
    bulk_clones = make_group_bulks(
        groups = clones,
        count_mat = count_mat,
        df_allele = df_allele, 
        lambdas_ref = lambdas_ref,
        gtf = gtf,
        min_depth = min_depth,
        nu = nu,
        segs_loh = segs_loh,
        ncores = ncores)

    bulk_clones = bulk_clones %>% 
        run_group_hmms(
            t = t,
            gamma = gamma,
            alpha = alpha,
            nu = nu,
            min_genes = min_genes,
            common_diploid = FALSE,
            diploid_chroms = diploid_chroms,
            ncores = ncores,
            verbose = verbose,
            retest = FALSE)

    bulk_clones = retest_bulks(
        bulk_clones,
        segs_consensus,
        gamma = gamma,
        use_loh = use_loh,
        min_LLR = min_LLR,
        diploid_chroms = diploid_chroms,
        ncores = ncores)
    
    fwrite(bulk_clones, glue('{out_dir}/bulk_clones_final.tsv.gz'), sep = '\t')

    if (plot) {
        p = plot_bulks(bulk_clones, min_LLR = min_LLR, use_pos = TRUE, genome = genome)
        ggsave(
            glue('{out_dir}/bulk_clones_final.png'), p, 
            width = 13, height = 2*length(unique(bulk_clones$sample)), dpi = 250
        )
    }
    
    log_message('All done!')

    return('Success')
}

segs_equal = function(segs_1, segs_2) {
    
    cols = c('CHROM', 'seg', 'seg_start', 'seg_end', 'cnv_state_post')
    
    equal = isTRUE(all.equal(
        segs_1 %>% select(any_of(cols)), 
        segs_2 %>% select(any_of(cols))
    ))
    
    return(equal)
    
}

subtrees_equal = function(subtrees_1, subtrees_2) {
    isTRUE(all.equal(subtrees_1, subtrees_2))
}

#' Run smoothed expression-based hclust
#' @param count_mat dgCMatrix Gene counts
#' @param lambdas_ref matrix Reference expression profiles
#' @param gtf dataframe Transcript GTF
#' @param sc_refs named list Reference choices for single cells
#' @param window integer Sliding window size
#' @param ncores integer Number of cores
#' @param verbose logical Verbosity
#' @keywords internal 
exp_hclust = function(count_mat, lambdas_ref, gtf, sc_refs = NULL, window = 101, ncores = 1, verbose = TRUE) {

    count_mat = check_matrix(count_mat)
    
    if (is.null(sc_refs)) {
        sc_refs = choose_ref_cor(count_mat, lambdas_ref, gtf)
    }

    lambdas_bar = get_lambdas_bar(lambdas_ref, sc_refs, verbose = FALSE)

    gexp_roll_wide = smooth_expression(
        count_mat,
        lambdas_bar,
        gtf,
        window = window,
        verbose = verbose
    ) %>% t

    dist_mat = parallelDist::parDist(gexp_roll_wide, threads = ncores)

    if (sum(is.na(dist_mat)) > 0) {
        log_warn('NAs in distance matrix, filling with 0s. Consider filtering out cells with low coverage.')
        dist_mat[is.na(dist_mat)] = 0
    }

    log_message('running hclust...')
    hc = hclust(dist_mat, method = "ward.D2")

    return(list('gexp_roll_wide' = gexp_roll_wide, 'hc' = hc))
}

#' Make a group of pseudobulks
#' @param groups list Contains fields named "sample", "cells", "size", "members"
#' @param count_mat dgCMatrix Gene counts
#' @param df_allele dataframe Alelle counts
#' @param lambdas_ref matrix Reference expression profiles
#' @param gtf dataframe Transcript GTF
#' @param min_depth integer Minimum allele depth to include
#' @param segs_loh dataframe Segments with clonal LOH to be excluded
#' @param ncores integer Number of cores
#' @return dataframe Pseudobulk profiles
#' @keywords internal 
make_group_bulks = function(groups, count_mat, df_allele, lambdas_ref, gtf, min_depth = 0, nu = 1, segs_loh = NULL, ncores = NULL) {
    
    if (length(groups) == 0) {
        return(data.frame())
    }

    ncores = ifelse(is.null(ncores), length(groups), ncores)

    results = mclapply(
            groups,
            mc.cores = ncores,
            function(g) {
                get_bulk(
                    count_mat = count_mat,
                    df_allele = df_allele,
                    subset = g$cells,
                    lambdas_ref = lambdas_ref,
                    gtf = gtf,
                    min_depth = min_depth,
                    nu = nu,
                    segs_loh = segs_loh
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

#' Run multiple HMMs 
#' @param bulks dataframe Pseudobulk profiles
#' @param gamma numeric Dispersion parameter for the Beta-Binomial allele model
#' @param t numeric Transition probability
#' @param alpha numeric P value cut-off to determine segment clusters in find_diploid
#' @param common_diploid logical Whether to find common diploid regions between pseudobulks
#' @param diploid_chroms character vector Known diploid chromosomes to use as baseline 
#' @param retest logcial Whether to retest CNVs
#' @param run_hmm logical Whether to run HMM segments or just retest
#' @param ncores integer Number of cores
#' @param allele_only logical Whether only use allele data to run HMM
#' @keywords internal
run_group_hmms = function(
    bulks, t = 1e-4, gamma = 20, alpha = 1e-4, min_genes = 10, nu = 1,
    common_diploid = TRUE, diploid_chroms = NULL, allele_only = FALSE, retest = TRUE, run_hmm = TRUE,
    exclude_neu = TRUE, ncores = 1, verbose = FALSE, debug = FALSE
) {

    # drop samples with no allele data
    bulks = bulks %>% group_by(sample) %>% filter(sum(!is.na(DP)) > 0) %>% ungroup()

    if (nrow(bulks) == 0) {
        return(data.frame())
    }

    n_groups = length(unique(bulks$sample))

    if (verbose) {
        log_message(glue('Running HMMs on {n_groups} cell groups..'))
    }

    # find common diploid region
    if (!run_hmm) {
        find_diploid = FALSE
    } else if (common_diploid & is.null(diploid_chroms)) {
        bulks = find_common_diploid(bulks, gamma = gamma, alpha = alpha, ncores = ncores)
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
                nu = nu,
                find_diploid = find_diploid, 
                run_hmm = run_hmm,
                allele_only = allele_only, 
                diploid_chroms = diploid_chroms,
                min_genes = min_genes,
                retest = retest, 
                verbose = verbose,
                exclude_neu = exclude_neu
            )
    })
                
    bad = sapply(results, inherits, what = "try-error")

    if (any(bad)) {
        log_error(results[bad][[1]])
        stop(results[bad][[1]])
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

#' Extract consensus CNV segments
#' @param bulks dataframe Pseudobulks
#' @param min_LLR numeric LLR threshold to filter CNVs 
#' @param min_overlap numeric Minimum overlap fraction to determine count two events as as overlapping
#' @return dataframe Consensus segments
#' @keywords internal
get_segs_consensus = function(bulks, min_LLR = 5, min_overlap = 0.45, retest = TRUE) {

    if (!'sample' %in% colnames(bulks)) {
        bulks$sample = 1
    }

    info_cols = c('sample', 'CHROM', 'seg', 'cnv_state', 'cnv_state_post',
            'seg_start', 'seg_end', 'seg_start_index', 'seg_end_index',
            'theta_mle', 'theta_sigma', 'phi_mle', 'phi_sigma', 
            'p_loh', 'p_del', 'p_amp', 'p_bamp', 'p_bdel',
            'LLR', 'LLR_y', 'LLR_x', 'n_genes', 'n_snps')

    # all possible aberrant segments
    segs_all = bulks %>% 
        group_by(sample, seg, CHROM) %>%
        mutate(seg_start = min(POS), seg_end = max(POS)) %>%
        filter(seg_start != seg_end) %>%
        ungroup() %>%
        select(any_of(info_cols)) %>%
        distinct() %>%
        mutate(cnv_state = ifelse(LLR < min_LLR | is.na(LLR), 'neu', cnv_state))

    # confident aberrant segments
    segs_star = segs_all %>% 
        filter(cnv_state != 'neu') %>%
        resolve_cnvs(
            min_overlap = min_overlap
        )

    if (retest) {

        # union of all aberrant segs
        segs_cnv = segs_all %>% 
            filter(cnv_state != 'neu') %>%
            arrange(CHROM) %>%
            {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$seg_start,
                    end = .$seg_end)
            )} %>%
            GenomicRanges::reduce() %>%
            as.data.frame() %>%
            select(CHROM = seqnames, seg_start = start, seg_end = end)

        # segs to be retested
        segs_retest = GenomicRanges::setdiff(
            segs_cnv %>% 
                {GenomicRanges::GRanges(
                    seqnames = .$CHROM,
                    IRanges::IRanges(start = .$seg_start,
                        end = .$seg_end)
            )},
            segs_star %>% 
                {GenomicRanges::GRanges(
                    seqnames = .$CHROM,
                    IRanges::IRanges(start = .$seg_start,
                        end = .$seg_end)
                )},
            ) %>%
            suppressWarnings() %>%
            as.data.frame() %>%
            select(CHROM = seqnames, seg_start = start, seg_end = end) %>%
            filter(seg_end - seg_start > 0) %>%
            mutate(cnv_state = 'retest', cnv_state_post = 'retest')

    } else {
        segs_retest = data.frame()
    }

    # union of neutral segments
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
    
    if (all(segs_all$cnv_state == 'neu')){
        segs_neu = segs_neu %>% 
            arrange(CHROM) %>%
            group_by(CHROM) %>%
            mutate(
                seg = paste0(CHROM, generate_postfix(1:n())),
                cnv_state = 'neu', 
                cnv_state_post = 'neu'
            ) %>% 
            ungroup()
        return(segs_neu)
    }
    
    segs_consensus = bind_rows(segs_star, segs_retest) %>%
        fill_neu_segs(segs_neu) %>%
        mutate(cnv_state_post = ifelse(cnv_state == 'neu', cnv_state, cnv_state_post))

    return(segs_consensus)

}

#' Fill neutral regions into consensus segments
#' @param segs_consensus dataframe CNV segments from multiple samples
#' @param segs_neu dataframe Neutral segments
#' @return dataframe Collections of neutral and aberrant segments with no gaps
#' @keywords internal
fill_neu_segs = function(segs_consensus, segs_neu) {
    
    # take complement of consensus aberrant segs
    gaps = GenomicRanges::setdiff(
        segs_neu %>% 
            {GenomicRanges::GRanges(
                seqnames = .$CHROM,
                IRanges::IRanges(start = .$seg_start,
                    end = .$seg_end)
        )},
        segs_consensus %>% 
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
        mutate(seg_cons = paste0(CHROM, generate_postfix(1:n()))) %>%
        ungroup() %>%
        mutate(CHROM = factor(CHROM, 1:22)) %>%
        arrange(CHROM)
    
    return(segs_consensus)
}

#' Map cells to the phylogeny (or genotypes) based on CNV posteriors
#' @param gtree tbl_graph A cell lineage tree
#' @param exp_post dataframe Expression posteriors
#' @param allele_post dataframe Allele posteriors
#' @return dataframe Clone posteriors
#' @keywords internal
get_clone_post = function(gtree, exp_post, allele_post) {

    clones = gtree %>%
        activate(nodes) %>%
        data.frame() %>% 
        # mutate(GT = ifelse(compartment == 'normal', '', GT)) %>%
        group_by(GT, clone, compartment) %>%
        summarise(
            clone_size = sum(leaf),
            .groups = 'drop'
        )
    
    # add the normal genotype if not in the tree
    if (min(clones$clone) > 1) {
        clones = clones %>% 
            add_row(.before = 1, GT = '', clone = 1, compartment = 'normal', clone_size = 0)
    }

    clone_segs = clones %>%
        mutate(
            prior_clone = ifelse(GT == '', 0.5, 0.5/(length(unique(GT)) - 1))
        ) %>%
        mutate(seg = GT) %>%
        tidyr::separate_rows(seg, sep = ',') %>%
        mutate(I = 1) %>%
        tidyr::complete(
            seg,
            tidyr::nesting(GT, clone, compartment, prior_clone, clone_size),
            fill = list('I' = 0)
        ) %>% 
        filter(seg != '')

    clone_post = full_join(
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
        mutate(
            l_clone_y = ifelse(is.na(l_clone_y), 0, l_clone_y),
            l_clone_x = ifelse(is.na(l_clone_x), 0, l_clone_x)
        ) %>%
        group_by(cell) %>%
        mutate(
            Z_clone = log(prior_clone) + l_clone_x + l_clone_y,
            Z_clone_x = log(prior_clone) + l_clone_x,
            Z_clone_y = log(prior_clone) + l_clone_y,
            p = exp(Z_clone - logSumExp(Z_clone)),
            p_x = exp(Z_clone_x - logSumExp(Z_clone_x)),
            p_y = exp(Z_clone_y - logSumExp(Z_clone_y))
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

#' Get unique CNVs from set of segments
#' @param segs_all dataframe CNV segments from multiple samples
#' @param min_overlap numeric scalar Minimum overlap fraction to determine count two events as as overlapping
#' @return dataframe Consensus CNV segments
#' @keywords internal 
resolve_cnvs = function(segs_all, min_overlap = 0.5, debug = FALSE) {

    if (nrow(segs_all) == 0) {
        return(segs_all)
    }
            
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
        distinct(vp, .keep_all = TRUE)

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

    G = igraph::graph_from_data_frame(d=E, vertices=V, directed=FALSE)

    segs_all = segs_all %>% mutate(component = igraph::components(G)$membership)

    segs_consensus = segs_all %>% group_by(component, sample) %>%
        mutate(LLR_sample = max(LLR_x + LLR_y)) %>%
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

#' get the single cell expression likelihoods
#' @param exp_counts dataframe Single-cell expression counts (CHROM, seg, cnv_state, gene, Y_obs, lambda_ref)
#' @param diploid_chroms character vector Known diploid chromosomes
#' @param use_loh logical Whether to include CNLOH regions in baseline
#' @return dataframe Single-cell CNV likelihood scores
#' @keywords internal
get_exp_likelihoods = function(exp_counts, diploid_chroms = NULL, use_loh = FALSE, depth_obs = NULL, mu = NULL, sigma = NULL) {
    
    exp_counts = exp_counts %>% filter(!is.na(Y_obs)) %>% filter(lambda_ref > 0)
    
    if (is.null(depth_obs)){
        depth_obs = sum(exp_counts$Y_obs)
    }

    if (use_loh) {
        ref_states = c('neu', 'loh')
    } else {
        ref_states = c('neu')
    }

    if (is.null(mu) | is.null(sigma)) {

        if (!is.null(diploid_chroms)) {
            exp_counts_diploid = exp_counts %>% filter(!loh) %>% filter(CHROM %in% diploid_chroms)
        } else {
            exp_counts_diploid = exp_counts %>% filter(!loh) %>% filter(cnv_state %in% ref_states)
        }

        fit = exp_counts_diploid %>% {fit_lnpois_cpp(.$Y_obs, .$lambda_ref, depth_obs)}
        mu = fit[1]
        sigma = fit[2]
    }

    res = exp_counts %>% 
        filter(cnv_state != 'neu') %>%
        group_by(CHROM, seg, cnv_state) %>%
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
            mu = UQ(mu),
            sigma = UQ(sigma),
            .groups = 'drop'
        )
        
    return(res)
}

#' get the single cell expression dataframe
#' @param segs_consensus dataframe Consensus segments
#' @param count_mat dgCMatrix gene expression count matrix
#' @param gtf dataframe Transcript gtf
#' @return dataframe single cell expression counts annotated with segments
#' @keywords internal
get_exp_sc = function(segs_consensus, count_mat, gtf, segs_loh = NULL) {

    gene_seg = GenomicRanges::findOverlaps(
            gtf %>% {GenomicRanges::GRanges(
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
            gtf %>% mutate(gene_index = 1:n()),
            by = c('gene_index')
        ) %>%
        mutate(CHROM = as.factor(CHROM)) %>%
        left_join(
            segs_consensus %>% mutate(seg_index = 1:n()),
            by = c('seg_index', 'CHROM')
        ) %>%
        distinct(gene, `.keep_all` = TRUE) 

    exp_sc = count_mat %>%
        as.matrix() %>%
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

    exp_sc = exclude_loh(exp_sc, segs_loh)

    return(exp_sc)
}

# exclude genes in clonal LOH regions
exclude_loh = function(exp_sc, segs_loh = NULL) {

    if (is.null(segs_loh)) {
        exp_sc = exp_sc %>% mutate(loh = FALSE)
    } else {
        log_message('Excluding clonal LOH regions .. ')
        
        genes_loh = GenomicRanges::findOverlaps(
                exp_sc %>% {GenomicRanges::GRanges(
                    seqnames = .$CHROM,
                    IRanges::IRanges(start = .$gene_start,
                        end = .$gene_start)
                )}, 
                segs_loh %>% {GenomicRanges::GRanges(
                    seqnames = .$CHROM,
                    IRanges::IRanges(start = .$seg_start,
                        end = .$seg_end)
                )}
            ) %>%
            as.data.frame() %>% 
            setNames(c('gene_index', 'seg_index')) %>%
            left_join(
                exp_sc %>% mutate(gene_index = 1:n()),
                by = c('gene_index')
            ) %>%
            pull(gene)
        
        exp_sc = exp_sc %>% mutate(loh = gene %in% genes_loh) 
    }
    
    return(exp_sc)
    
}

#' compute single-cell expression posteriors
#' @param segs_consensus dataframe Consensus segments
#' @param count_mat dgCMatrix gene expression count matrix
#' @param gtf dataframe transcript gtf
#' @param lambdas_ref matrix Reference expression profiles
#' @return dataframe Expression posteriors
#' @keywords internal
get_exp_post = function(segs_consensus, count_mat, gtf, lambdas_ref, sc_refs = NULL, diploid_chroms = NULL, use_loh = NULL, segs_loh = NULL, ncores = 30, verbose = TRUE, debug = FALSE) {

    exp_sc = get_exp_sc(segs_consensus, count_mat, gtf, segs_loh)

    if (is.null(use_loh)) {
        if (mean(exp_sc$cnv_state == 'neu' & (!exp_sc$loh)) < 0.05) {
            use_loh = TRUE
            log_message('less than 5% genes are in neutral region - including LOH in baseline')
        } else {
            use_loh = FALSE
        }
    } else if (use_loh) {
        log_message('Including LOH in baseline as specified')
    }
    
    if (is.null(sc_refs)) {
        sc_refs = choose_ref_cor(count_mat, lambdas_ref, gtf)
    }

    cells = names(sc_refs)

    results = mclapply(
        cells,
        mc.cores = ncores,
        function(cell) {
   
            ref = sc_refs[cell]

            exp_sc = exp_sc[,c('gene', 'seg', 'CHROM', 'cnv_state', 'loh', 'seg_start', 'seg_end', cell)] %>%
                rename(Y_obs = ncol(.))

            exp_sc %>%
                mutate(
                    lambda_ref = lambdas_ref[, ref][gene],
                    lambda_obs = Y_obs/sum(Y_obs),
                    logFC = log2(lambda_obs/lambda_ref)
                ) %>%
                get_exp_likelihoods(
                    use_loh = use_loh,
                    diploid_chroms = diploid_chroms
                ) %>%
                mutate(cell = cell, ref = ref)

        }
    )

    bad = sapply(results, inherits, what = "try-error")

    if (any(bad)) {
        if (verbose) {log_warn(glue('{sum(bad)} cell(s) failed'))}
        log_warn(results[bad][1])
        log_warn(cells[bad][1])
    } else {
        log_message('All cells succeeded')
    }
    
    exp_post = results[!bad] %>%
        bind_rows() %>%
        mutate(seg = factor(seg, mixedsort(unique(seg)))) %>%
        left_join(
            segs_consensus %>% select(
                CHROM, 
                seg = seg_cons, 
                seg_start,
                seg_end,
                prior_loh = p_loh, prior_amp = p_amp, prior_del = p_del, prior_bamp = p_bamp, prior_bdel = p_bdel
            ),
            by = c('seg', 'CHROM')
        ) %>%
        # if the opposite state has a very small prior, and phi is in the opposite direction, then CNV posterior can still be high which is miselading
        mutate_at(
            vars(contains('prior')),
            function(x) {ifelse(x < 0.05, 0, x)}
        ) %>%
        compute_posterior() %>%
        mutate(seg_label = paste0(seg, '(', cnv_state, ')')) %>%
        mutate(seg_label = factor(seg_label, unique(seg_label)))
    
    return(exp_post)
}

#' Do bayesian averaging to get posteriors
#' @param PL dataframe Likelihoods and priors
#' @return dataframe Posteriors
#' @keywords internal
compute_posterior = function(PL) {
    PL %>% 
    rowwise() %>%
    mutate(
        Z_amp = logSumExp(c(l21 + log(prior_amp/4), l31 + log(prior_amp/4))),
        Z_loh = l20 + log(prior_loh/2),
        Z_del = l10 + log(prior_del/2),
        Z_bamp = l22 + log(prior_bamp/2),
        Z_bdel = l00 + log(prior_bdel/2),
        Z_n = l11 + log(1/2),
        Z = logSumExp(
            c(Z_n, Z_loh, Z_del, Z_amp, Z_bamp, Z_bdel)
        ),
        Z_cnv = logSumExp(
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


#' Get phased haplotypes
#' @param bulks dataframe Subtree pseudobulk profiles
#' @param segs_consensus dataframe Consensus CNV segments
#' @param naive logical Whether to use naive haplotype classification
#' @return dataframe Posterior haplotypes
#' @keywords internal
get_haplotype_post = function(bulks, segs_consensus, naive = FALSE) {
    
    # add sample column if only one sample
    if ((!'sample' %in% colnames(bulks)) | (!'sample' %in% colnames(segs_consensus))) {
        bulks['sample'] = '0'
        segs_consensus['sample'] = '0'
    }

    # nothing to test
    if (all(segs_consensus$cnv_state_post == 'neu')) {
        stop('No CNVs')
    }
    
    if (naive) {
        bulks = bulks %>% mutate(haplo_post = ifelse(AR >= 0.5, 'major', 'minor'))
    }
    
    haplotypes = bulks %>%
        filter(!is.na(pAD)) %>%
        select(CHROM, seg, snp_id, sample, haplo_post) %>%
        inner_join(
            segs_consensus,
            by = c('sample', 'CHROM', 'seg')
        ) %>%
        select(CHROM, seg = seg_cons, cnv_state, snp_id, haplo_post)
    
    return(haplotypes)
}

#' get CNV allele posteriors
#' @param df_allele dataframe Allele counts
#' @param segs_consensus dataframe Consensus CNV segments
#' @param haplotypes dataframe Haplotype classification
#' @return dataframe Allele posteriors
#' @keywords internal
get_allele_post = function(df_allele, haplotypes, segs_consensus) {
    
    allele_counts = df_allele %>%
        mutate(pAD = ifelse(GT == '1|0', AD, DP - AD)) %>%
        inner_join(
            haplotypes %>% select(CHROM, seg, cnv_state, snp_id, haplo_post),
            by = c('CHROM', 'snp_id')
        )  %>%
        filter(cnv_state != 'neu') %>%
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
    
    allele_post = allele_counts %>%
        group_by(cell, CHROM, seg, cnv_state) %>%
        summarise(
            major = sum(major_count),
            minor = sum(minor_count),
            total = major + minor,
            MAF = major/total,
            .groups = 'drop'
        ) %>%
        left_join(
            segs_consensus %>% select(
                seg = seg_cons, seg_start, seg_end,
                prior_loh = p_loh, prior_amp = p_amp, prior_del = p_del, prior_bamp = p_bamp, prior_bdel = p_bdel
            ),
            by = 'seg'
        ) %>%
        rowwise() %>%
        mutate(
            l11 = dbinom(major, total, prob = 0.5, log = TRUE),
            l10 = dbinom(major, total, prob = 0.9, log = TRUE),
            l01 = dbinom(major, total, prob = 0.1, log = TRUE),
            l20 = dbinom(major, total, prob = 0.9, log = TRUE),
            l02 = dbinom(major, total, prob = 0.1, log = TRUE),
            l21 = dbinom(major, total, prob = 0.66, log = TRUE),
            l12 = dbinom(major, total, prob = 0.33, log = TRUE),
            l31 = dbinom(major, total, prob = 0.75, log = TRUE),
            l13 = dbinom(major, total, prob = 0.25, log = TRUE),
            l32 = dbinom(major, total, prob = 0.6, log = TRUE),
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
#' @param exp_post dataframe Expression single-cell CNV posteriors
#' @param allele_post dataframe Allele single-cell CNV posteriors
#' @param segs_consensus dataframe Consensus CNV segments
#' @return dataframe Joint single-cell CNV posteriors
#' @keywords internal
get_joint_post = function(exp_post, allele_post, segs_consensus) {
    
    joint_post = full_join(
            exp_post %>%
                filter(cnv_state != 'neu') %>%
                select(
                    any_of(c('cell', 'CHROM', 'seg', 'cnv_state', 
                    'l11_x' = 'l11', 'l20_x' = 'l20', 'l10_x' = 'l10', 'l21_x' = 'l21', 'l31_x' = 'l31', 'l22_x' = 'l22', 'l00_x' = 'l00',
                    'Z_x' = 'Z', 'Z_cnv_x' = 'Z_cnv', 'Z_n_x' = 'Z_n', 'logBF_x' = 'logBF'))
                ),
            allele_post %>% 
                select(
                    any_of(c('cell', 'CHROM', 'seg', 'cnv_state',
                    'l11_y' = 'l11', 'l20_y' = 'l20', 'l10_y' = 'l10', 'l21_y' = 'l21', 'l31_y' = 'l31', 'l22_y' = 'l22', 'l00_y' = 'l00',
                    'Z_y' = 'Z', 'Z_cnv_y' = 'Z_cnv', 'Z_n_y' = 'Z_n', 'logBF_y' = 'logBF',
                    'n_snp' = 'total', 'MAF', 'major', 'total'))
                ),
            c("cell", "CHROM", "seg", "cnv_state")
        ) %>%
        mutate_at(
            vars(matches("_x|_y")),
            function(x) tidyr::replace_na(x, 0)
        ) %>%
        left_join(
            segs_consensus %>% select(
                seg = seg_cons,
                seg_start,
                seg_end,
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
        mutate(seg = factor(seg, mixedsort(unique(seg)))) %>%
        mutate(seg_label = paste0(seg, '(', cnv_state, ')')) %>%
        mutate(seg_label = factor(seg_label, unique(seg_label)))
    
    return(joint_post)
}

#' retest consensus segments on pseudobulks
#' @param bulks dataframe Pseudobulk profiles
#' @param segs_consensus dataframe Consensus segments
#' @param use_loh logical Whether to use loh in the baseline
#' @param diploid_chroms vector User-provided diploid chromosomes
#' @return dataframe Retested pseudobulks 
#' @keywords internal 
retest_bulks = function(bulks, segs_consensus = NULL,
    t = 1e-5, min_genes = 10, gamma = 20, nu = 1,
    use_loh = FALSE, diploid_chroms = NULL, ncores = 1, exclude_neu = TRUE, min_LLR = 5) {

    if (is.null(segs_consensus)) {
        segs_consensus = get_segs_consensus(bulks)
    }

    # default
    if (is.null(use_loh)) {
        length_neu = segs_consensus %>% filter(cnv_state == 'neu') %>% pull(seg_length) %>% sum
        if (length_neu < 1.5e8) {
            use_loh = TRUE
            log_message('less than 5% of genome is in neutral region - including LOH in baseline')
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
    
    # retest CNVs
    bulks = bulks %>% 
        run_group_hmms(
            t = t, 
            gamma = gamma, 
            nu = nu,
            min_genes = min_genes,
            run_hmm = FALSE,
            exclude_neu = exclude_neu,
            ncores = ncores
        ) %>%
        mutate(
            LLR = ifelse(is.na(LLR), 0, LLR)
        ) %>%
        mutate(cnv_state_post = ifelse(LLR < min_LLR, 'neu', cnv_state_post)) %>%
        mutate(cnv_state = cnv_state_post)

    return(bulks)
}

#' test for multi-allelic CNVs
#' @param bulks dataframe Pseudobulk profiles
#' @param segs_consensus dataframe Consensus segments
#' @param min_LLR numeric CNV LLR threshold to filter events
#' @param p_min numeric Probability threshold to call multi-allelic events
#' @return dataframe Consensus segments annotated with multi-allelic events
#' @keywords internal 
test_multi_allelic = function(bulks, segs_consensus, min_LLR = 5, p_min = 0.999) {

    log_message('Testing for multi-allelic CNVs ..')
    
    segs_multi = bulks %>% 
        distinct(sample, CHROM, seg_cons, LLR, p_amp, p_del, p_bdel, p_loh, p_bamp, cnv_state_post) %>%
        rowwise() %>%
        mutate(p_max = max(c(p_amp, p_del, p_bdel, p_loh, p_bamp))) %>%
        filter(LLR > min_LLR & p_max > p_min) %>%
        group_by(seg_cons) %>%
        summarise(
            cnv_states = list(sort(unique(cnv_state_post))),
            n_states = length(unlist(cnv_states))
        ) %>%
        filter(n_states > 1)

    segs = segs_multi$seg_cons

    log_message(glue('{length(segs)} multi-allelic CNVs found: {paste(segs, collapse = ",")}'))

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
                n_states = ifelse(cnv_state == 'neu', 0, 1),
                cnv_states = cnv_state
            )
    }

    segs_consensus = segs_consensus %>%
        mutate(cnv_states = unlist(purrr::map(cnv_states, function(x){paste0(x, collapse = ',')})))
    
    return(segs_consensus)
}

#' expand multi-allelic CNVs into separate entries in the single-cell posterior dataframe
#' @param sc_post dataframe Single-cell posteriors
#' @param segs_consensus dataframe Consensus segments
#' @return dataframe Single-cell posteriors with multi-allelic CNVs split into different entries
#' @keywords internal 
expand_states = function(sc_post, segs_consensus) {

    segs_multi = segs_consensus %>% filter(n_states > 1) %>%
        select(seg = seg_cons, cnv_states, n_states) %>%
        tidyr::separate_longer_delim(cnv_states, ',') %>%
        rename(cnv_state = cnv_states)

    if (any(segs_consensus$n_states > 1)) {

        sc_post_multi = sc_post %>% 
            select(-cnv_state) %>%
            inner_join(
                segs_multi,
                by = 'seg',
                relationship = "many-to-many"
            ) %>%
            mutate(
                seg = paste0(seg, '_', cnv_state), 
            ) %>%
            rowwise() %>%
            mutate(
                p_cnv = get(glue('p_{cnv_state}')),
                p_n = 1 - p_cnv,
                Z_cnv = get(glue('Z_{cnv_state}'))
            ) %>%
            ungroup()

        # note that *_x and *_y columns are not updated .. to fix
        sc_post = sc_post %>% filter(!seg %in% segs_multi$seg) %>%
            mutate(n_states = 1) %>%
            bind_rows(sc_post_multi) %>%
            arrange(cell, CHROM, seg) %>%
            mutate(seg_label = paste0(seg, '(', cnv_state, ')')) %>%
            mutate(seg_label = factor(seg_label, unique(seg_label)))

    } else {
        log_message('No multi-allelic CNVs, skipping ..')
    }

    return(sc_post)
}


## gtools is orphaned
## https://github.com/cran/gtools/blob/master/R/mixedsort.R
#' @keywords internal 
mixedsort <- function(x,
                      decreasing = FALSE,
                      na.last = TRUE,
                      blank.last = FALSE,
                      numeric.type = c("decimal", "roman"),
                      roman.case = c("upper", "lower", "both"),
                      scientific = TRUE) {
  ord <- mixedorder(x,
    decreasing = decreasing,
    na.last = na.last,
    blank.last = blank.last,
    numeric.type = numeric.type,
    roman.case = roman.case,
    scientific = scientific
  )
  x[ord]
}

#' @keywords internal 
mixedorder <- function(x,
                       decreasing = FALSE,
                       na.last = TRUE,
                       blank.last = FALSE,
                       numeric.type = c("decimal", "roman"),
                       roman.case = c("upper", "lower", "both"),
                       scientific = TRUE) {
  # - Split each each character string into an vector of strings and
  #   numbers
  # - Separately rank numbers and strings
  # - Combine orders so that strings follow numbers

  numeric.type <- match.arg(numeric.type)
  roman.case <- match.arg(roman.case)

  if (length(x) < 1) {
    return(NULL)
  } else if (length(x) == 1) {
    return(1)
  }

  if (!is.character(x)) {
    return(order(x, decreasing = decreasing, na.last = na.last))
  }

  delim <- "\\$\\@\\$"

  if (numeric.type == "decimal") {
    if (scientific) {
      regex <- "((?:(?i)(?:[-+]?)(?:(?=[.]?[0123456789])(?:[0123456789]*)(?:(?:[.])(?:[0123456789]{0,}))?)(?:(?:[eE])(?:(?:[-+]?)(?:[0123456789]+))|)))"
    } # uses PERL syntax
    else {
      regex <- "((?:(?i)(?:[-+]?)(?:(?=[.]?[0123456789])(?:[0123456789]*)(?:(?:[.])(?:[0123456789]{0,}))?)))"
    } # uses PERL syntax

    numeric <- function(x) as.numeric(x)
  }
  else if (numeric.type == "roman") {
    regex <- switch(roman.case,
      "both"  = "([IVXCLDMivxcldm]+)",
      "upper" = "([IVXCLDM]+)",
      "lower" = "([ivxcldm]+)"
    )
    numeric <- function(x) roman2int(x)
  }
  else {
    stop("Unknown value for numeric.type: ", numeric.type)
  }

  nonnumeric <- function(x) {
    ifelse(is.na(numeric(x)), toupper(x), NA)
  }

  x <- as.character(x)

  which.nas <- which(is.na(x))
  which.blanks <- which(x == "")

  ####
  # - Convert each character string into an vector containing single
  #   character and  numeric values.
  ####

  # find and mark numbers in the form of +1.23e+45.67
  delimited <- gsub(regex,
    paste(delim, "\\1", delim, sep = ""),
    x,
    perl = TRUE
  )

  # separate out numbers
  step1 <- strsplit(delimited, delim)

  # remove empty elements
  step1 <- lapply(step1, function(x) x[x > ""])

  # create numeric version of data
  suppressWarnings(step1.numeric <- lapply(step1, numeric))

  # create non-numeric version of data
  suppressWarnings(step1.character <- lapply(step1, nonnumeric))

  # now transpose so that 1st vector contains 1st element from each
  # original string
  maxelem <- max(sapply(step1, length))

  step1.numeric.t <- lapply(
    1:maxelem,
    function(i) {
      sapply(
        step1.numeric,
        function(x) x[i]
      )
    }
  )

  step1.character.t <- lapply(
    1:maxelem,
    function(i) {
      sapply(
        step1.character,
        function(x) x[i]
      )
    }
  )

  # now order them
  rank.numeric <- sapply(step1.numeric.t, rank)
  rank.character <- sapply(
    step1.character.t,
    function(x) as.numeric(factor(x))
  )

  # and merge
  rank.numeric[!is.na(rank.character)] <- 0 # mask off string values

  rank.character <- t(
    t(rank.character) +
      apply(matrix(rank.numeric), 2, max, na.rm = TRUE)
  )

  rank.overall <- ifelse(is.na(rank.character), rank.numeric, rank.character)

  order.frame <- as.data.frame(rank.overall)
  if (length(which.nas) > 0) {
    if (is.na(na.last)) {
      order.frame[which.nas, ] <- NA
    } else if (na.last) {
      order.frame[which.nas, ] <- Inf
    } else {
      order.frame[which.nas, ] <- -Inf
    }
  }

  if (length(which.blanks) > 0) {
    if (is.na(blank.last)) {
      order.frame[which.blanks, ] <- NA
    } else if (blank.last) {
      order.frame[which.blanks, ] <- 1e99
    } else {
      order.frame[which.blanks, ] <- -1e99
    }
  }

  order.frame <- as.list(order.frame)
  order.frame$decreasing <- decreasing
  order.frame$na.last <- NA

  retval <- do.call("order", order.frame)

  return(retval)
}


#' @keywords internal
roman2int <- function(roman){
    roman <- trimws(toupper(as.character(roman)))

    roman2int.inner <- function(roman){
        results <- roman2int_internal(letters = as.character(roman), nchar = as.integer(nchar(roman)))
        return(results)
    }

    tryIt <- function(x){
            retval <- try(roman2int.inner(x), silent=TRUE)
            if(is.numeric(retval)){
                retval
            }
            else{
                NA
            }
        }

    retval <- sapply(roman, tryIt)

    retval
}


## for tidyverse (magrittr & dplyr) functions 
if (getRversion() >= "2.15.1"){
  utils::globalVariables(c(".", "AD", "AD_all", "ALT", "AR", "CHROM", "DP", "DP_all", "FILTER", "FP",
    "GT", "ID", "LLR", "LLR_sample","LLR_x", "LLR_y", "L_x_a", "L_x_d", "L_x_n", 
    "L_y_a", "L_y_d", "L_y_n", "MAF", "OTH", "OTH_all", "POS",
    "QUAL", "REF", "TP", "UQ", "Y_obs", "Z", "Z_amp", "Z_bamp", "Z_bdel", "Z_clone", "Z_clone_x",
    "Z_clone_y", "Z_cnv", "Z_del", "Z_loh", "Z_n", "acen_hg19", "acen_hg38", "annot", "avg_entropy",
    "branch", "cM", "cell", "cell_index", "chrom_sizes_hg19", "chrom_sizes_hg38", "clone",
    "clone_opt", "clone_size", "cluster", "cnv_state", "cnv_state_expand", "cnv_state_map",
    "cnv_state_post", "cnv_states", "compartment", "component", "cost", "d_obs", "diploid",
    "down", "edges", "end_x", "end_y", "exp_rollmean", "expected_colnames", "extract",
    "frac_overlap_x", "frac_overlap_y", "from", "from_label", "from_node", "gaps_hg19",
    "gaps_hg38", "gene", "gene_end", "gene_index", "gene_length", "gene_snps", "gene_start",
    "group", "groupOTU", "gtf_hg19", "gtf_hg38", "gtf_mm10", "haplo", "haplo_naive", "haplo_post",
    "haplo_theta_min", "het", "hom_alt", "i", "inter_snp_cm", "inter_snp_dist", "isTip",
    "is_desc", "j", "keep", "l", "l00", "l00_x", "l00_y", "l10", "l10_x", "l10_y", "l11", "l11_x",
    "l11_y", "l20", "l20_x", "l20_y", "l21", "l21_x", "l21_y", "l22", "l22_x", "l22_y", "l31",
    "l31_x", "l31_y", "l_clone", "l_clone_x", "l_clone_y", "label", "lambda_obs", "lambda_ref", "last_mut", "leaf",
    "len_overlap", "len_x", "len_y", "lnFC", "lnFC_i", "lnFC_j", "lnFC_max_i", "lnFC_max_j",
    "logBF", "logBF_x", "logBF_y", "logFC", "loh", "major", "major_count", "marker_index", 
    "minor", "minor_count", "mu", "mut_burden", "n_chrom_snp", "n_genes", "n_mut", "n_sibling", 
    "n_snps", "n_states", "name", "node", "node_phylo", "nodes", "p", "pAD", "pBAF", "p_0",
    "p_1", "p_amp", "p_bamp", "p_bdel", "p_cnv", "p_cnv_expand", "p_del", "p_loh", "p_max", "p_n", "p_neu", "p_s", "p_up",
    "phi_mle", "phi_mle_roll", "phi_sigma", "pnorm.range", "potential_missing_columns",
    "precision", "prior_amp", "prior_bamp", "prior_bdel", "prior_clone", "prior_del",
    "prior_loh", "s", "seg", "seg_cons", "seg_end", "seg_end_index", "seg_label", "seg_length",
    "seg_start", "seg_start_index", "segs_consensus", "seqnames", "set_colnames", "sig",
    "site", "size", "snp_id", "snp_index", "snp_rate", "start_x", "start_y", "state", "state_post",
    "superclone", "theta_hat", "theta_level", "theta_mle", "theta_sigma", "to", "to_label",
    "to_node", "total", "value", "variable", "vcf_meta", "vertex", "vp", "width", "write.vcf", "x", "y"))
}








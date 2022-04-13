##### R implementation of ScisTree perfect phylogeny #####

#' score a tree based on maximum likelihood
#' @param tree phylo object
#' @param P genotype probability matrix
#' @param get_l_matrix whether to compute the whole likelihood matrix
#' @return list Likelihood scores of a tree
#' @export
score_tree = function(tree, P, get_l_matrix = FALSE) {
 
    tree = reorder(tree, order = 'postorder')
        
    n = nrow(P)
    m = ncol(P)

    logQ = matrix(nrow = tree$Nnode * 2 + 1, ncol = m)

    logP_0 = log(P)
    logP_1 = log(1-P)
    
    node_order = c(tree$edge[,2], n+1)
    node_order = node_order[node_order > n]
    
    logQ[1:n,] = logP_1 - logP_0

    children_dict = allChildrenCPP(tree$edge)

    logQ = CgetQ(logQ, children_dict, node_order)

    if (get_l_matrix) {
        l_matrix = sweep(logQ, 2, colSums(logP_0), FUN = '+')
        l_tree = sum(apply(l_matrix, 2, max))
    } else {
        l_matrix = NULL
        l_tree = sum(apply(logQ, 2, max)) + sum(logP_0)
    }
    
    return(list('l_tree' = l_tree, 'logQ' = logQ, 'l_matrix' = l_matrix))
    
}

#' maximum likelihood tree search via NNI
#' @param tree_init iintial tree as phylo object
#' @param P genotype probability matrix
#' @param max_iter maximum number of iterations
#' @param eps tolerance threshold in likelihood difference for stopping
#' @param verbose verbosity
#' @param ncores number of cores to use
#' @return a list of trees corresponding to the rearrangement steps
#' @export
perform_nni = function(tree_init, P, max_iter = 100, eps = 0.01, ncores = 1, verbose = TRUE) {

    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    P = as.matrix(P)
    
    converge = FALSE

    i = 1
    max_current = score_tree(tree_init, P)$l_tree
    tree_current = tree_init
    tree_list = list()
    tree_list[[1]] = tree_current

    while (!converge & i <= max_iter) {
        
        i = i + 1
        
        ptm = proc.time()

        neighbours = nni(tree_current, ncores = ncores)

        scores = mclapply(
                mc.cores = ncores,
                neighbours,
                function(tree) {
                    score_tree(tree, P)$l_tree
                }
            ) %>%
            unlist()

        if (max(scores) > max_current + eps) {
            tree_list[[i]] = tree_current = neighbours[[which.max(scores)]]
            tree_list[[i]]$likelihood = max_current = max(scores)
            converge = FALSE
        } else {
            converge = TRUE
        }
        
        runtime = proc.time() - ptm
        
        if (verbose) {
            msg = glue('Iter {i} {max_current} {signif(unname(runtime[3]),2)}s')
            message(msg)
            log_info(msg)
        }
    }

    class(tree_list) = 'multiPhylo'

    return(tree_list)
}

# from phangorn
#' Perform the NNI at a specific branch
#' @param tree phylo Single-cell phylogenetic tree
#' @param n integer Branch ID
#' @keywords internal
nnin <- function(tree, n) {
  attr(tree, "order") <- NULL
  tree1 <- tree
  tree2 <- tree
  edge <- matrix(tree$edge, ncol = 2)
  parent <- edge[, 1]
  child <- tree$edge[, 2]
  k <- min(parent) - 1
  ind <- which(child > k)[n]
  if (is.na(ind)) return(NULL)
  p1 <- parent[ind]
  p2 <- child[ind]
  ind1 <- which(parent == p1)
  ind1 <- ind1[ind1 != ind][1]
  ind2 <- which(parent == p2)
  e1 <- child[ind1]
  e2 <- child[ind2[1]]
  e3 <- child[ind2[2]]
  tree1$edge[ind1, 2] <- e2
  tree1$edge[ind2[1], 2] <- e1
  tree2$edge[ind1, 2] <- e3
  tree2$edge[ind2[2], 2] <- e1
  if (!is.null(tree$edge.length)) {
    tree1$edge.length[c(ind1, ind2[1])] <- tree$edge.length[c(ind2[1], ind1)]
    tree2$edge.length[c(ind1, ind2[2])] <- tree$edge.length[c(ind2[2], ind1)]
  }
  tree1 <- reorder(tree1, "postorder")
  tree2 <- reorder(tree2, "postorder")

  tree1$tip.label <- tree2$tip.label <- NULL
  
  result <- list(tree1, tree2)
  
  result
}

# from phangorn
#' Generate tree neighbourhood
#' @param tree phylo Single-cell phylogenetic tree
#' @param ncores integer Number of cores to use
#' @keywords internal
nni <- function(tree, ncores = 1) {
  tip.label <- tree$tip.label
  attr(tree, "order") <- NULL
  k <- min(tree$edge[, 1]) - 1
  n <- sum(tree$edge[, 2] > k)
  result <- vector("list", 2 * n)
  l <- 1
  
  result = mclapply(
          mc.cores = ncores,
          seq(1,n),
          function(i) {
               nnin(tree, i)
          }
     ) %>% Reduce('c', .)

  attr(result, "TipLabel") <- tip.label
  class(result) <- "multiPhylo"

  return(result)
}

#' from phangorn
#' @rdname upgma
#' @export
"upgma" <- function(D, method = "average", ...) {
  DD <- as.dist(D)
  hc <- hclust(DD, method = method, ...)
  result <- ape::as.phylo(hc)
  result <- reorder(result, "postorder")
  result
}

#' mark the tumor lineage of a phylogeny
#' @param gtree a tidygraph tree
#' @return a tidygraph tree
#' @export
mark_tumor_lineage = function(gtree) {
  
    mut_nodes = gtree %>%
        activate(nodes) %>%
        filter(!is.na(site)) %>%
        as.data.frame() %>%
        pull(id)

    mut_burdens = lapply(
        mut_nodes,
        function(node) {
            gtree %>%
            activate(nodes) %>%
            mutate(
                mut_burden = ifelse(GT == '', 0, str_count(GT, ',') + 1)
            ) %>%
            ungroup() %>%
            mutate(seq = bfs_rank(root = node)) %>%
            data.frame %>%
            filter(leaf & !is.na(seq)) %>%
            pull(mut_burden) %>%
            sum
        }
    )

    tumor_root = mut_nodes[which.max(mut_burdens)]
    
    gtree = gtree %>%
        activate(nodes) %>%
        mutate(
            seq = bfs_rank(root = tumor_root),
            compartment = ifelse(is.na(seq), 'normal', 'tumor'),
            is_tumor_root = tumor_root == id
        )

    compartment_dict = gtree %>%
        activate(nodes) %>%
        as.data.frame() %>%
        {setNames(.$compartment, .$id)} 

    gtree = gtree %>%
        activate(edges) %>%
        mutate(compartment = compartment_dict[to]) 
    
    return(gtree)
    
}

#' Convert a single-cell phylogeny with mutation placements into a mutation graph
#' @param gtree tbl_graph The single-cell phylogeny
#' @param mut_nodes dataframe Mutation placements
#' @return igraph Mutation graph
#' @keywords internal
get_mut_tree = function(gtree, mut_nodes) {

    G = gtree %>%
        activate(nodes) %>%
        arrange(last_mut) %>%
        convert(to_contracted, last_mut) %>%
        mutate(label = last_mut, id = 1:n()) %>%
        as.igraph

    G = label_edges(G)

    V(G)$node = G %>%
        igraph::as_data_frame('vertices') %>%
        left_join(
            mut_nodes %>% dplyr::rename(node = name),
            by = c('label' = 'site')
        ) %>%
        pull(node)

    G = G %>% transfer_links()

    return(G)

}

#' Compute site branch likelihood (not used)
#' @param node integer Node id
#' @param site character Mutation site name
#' @param gtree tbl_graph The single-cell phylogeny
#' @param geno matrix Genotype probability matrix
#' @return numeric Likelihood of assigning the mutation to this node
#' @keywords internal
l_s_v = function(node, site, gtree, geno) {
  
    gtree %>%
    activate(nodes) %>%
    mutate(seq = bfs_rank(root = node)) %>%
    data.frame %>%
    filter(leaf) %>%
    mutate(
        p_0 = unlist(geno[site, name]),
        p_1 = 1 - p_0
    )  %>%
    mutate(is_desc = !is.na(seq)) %>%
    mutate(l = ifelse(is_desc, log(p_1), log(p_0))) %>%
    pull(l) %>%
    sum
}

#' Find maximum lilkelihood assignment of mutations on a tree
#' @param tree phylo Single-cell phylogenetic tree
#' @param P matrix Genotype probability matrix
#' @return list Mutation 
#' @keywords internal
get_tree_post = function(tree, P) {
   
    sites = colnames(P)
    n = nrow(P)
    tree_stats = score_tree(tree, P, get_l_matrix = TRUE)

    l_matrix = as.data.frame(tree_stats$l_matrix)

    colnames(l_matrix) = sites
    rownames(l_matrix) = c(tree$tip.label, paste0('Node', 1:tree$Nnode))

    mut_nodes = data.frame(
            site = sites,
            node_phylo = apply(l_matrix, 2, which.max),
            l = apply(l_matrix, 2, max)
        ) %>%
        mutate(
            name = ifelse(node_phylo <= n, tree$tip.label[node_phylo], paste0('Node', node_phylo - n))
        ) %>%
        group_by(name) %>%
        summarise(
            site = paste0(sort(site), collapse = ','),
            n_mut = n(),
            l = sum(l),
            .groups = 'drop'
        )

    gtree = tree %>%
        ape::ladderize() %>%
        as_tbl_graph() %>%
        mutate(
            leaf = node_is_leaf(),
            root = node_is_root(),
            depth = bfs_dist(root = 1),
            id = row_number()
        )

    # leaf annotation for edges
    gtree = gtree %>%
        activate(edges) %>%
        select(-any_of(c('leaf'))) %>%
        left_join(
            gtree %>%
                activate(nodes) %>%
                data.frame() %>%
                select(id, leaf),
            by = c('to' = 'id')
        )

    # annotate the tree
    gtree = mut_to_tree(gtree, mut_nodes)
    gtree = mark_tumor_lineage(gtree)

    return(list('mut_nodes' = mut_nodes, 'gtree' = gtree, 'l_matrix' = l_matrix))
}

#' transfer mutation assignment onto a single-cell phylogeny
#' @param gtree tbl_graph The single-cell phylogeny
#' @param mut_nodes dataframe Mutation placements
#' @return tbl_graph A single-cell phylogeny with mutation placements
#' @keywords internal
mut_to_tree = function(gtree, mut_nodes) {
   
    # transfer mutation to tree
    gtree = gtree %>%
        activate(nodes) %>%
        select(-any_of(c('n_mut', 'l', 'site', 'clone'))) %>%
        left_join(
            mut_nodes %>%
                mutate(n_mut = unlist(purrr::map(str_split(site, ','), length))) %>%
                select(name, n_mut, site),
            by = 'name'
        ) %>%
        mutate(n_mut = ifelse(is.na(n_mut), 0, n_mut))

    # get branch length
    gtree = gtree %>% 
        activate(edges) %>%
        select(-any_of(c('length'))) %>%
        left_join(
            gtree %>%
                activate(nodes) %>%
                data.frame() %>%
                select(id, length = n_mut),
            by = c('to' = 'id')
        ) %>%
        mutate(length = ifelse(leaf, pmax(length, 0.2), length))

    # label genotype on nodes
    node_to_mut = gtree %>% activate(nodes) %>% data.frame() %>% {setNames(.$site, .$id)}

    gtree = gtree %>%
        activate(nodes) %>%
        mutate(GT = unlist(
            map_bfs(node_is_root(),
            .f = function(path, ...) { paste0(na.omit(node_to_mut[path$node]), collapse = ',') })
            ),
            last_mut = unlist(
                map_bfs(node_is_root(),
                .f = function(path, ...) { 
                    past_muts = na.omit(node_to_mut[path$node])
                    if (length(past_muts) > 0) {
                        return(past_muts[length(past_muts)])
                    } else {
                        return('')
                    }
                })
            )
        ) %>%
        mutate(GT = ifelse(GT == '' & !is.na(site), site, GT))

    # preserve the clone ids
    if ('GT' %in% colnames(mut_nodes)) {
        gtree = gtree %>% activate(nodes) %>%
            left_join(
                mut_nodes %>% select(GT, clone),
                by = 'GT'
            )
    }
    
    return(gtree)
}

#' Convert the phylogeny from tidygraph to igraph object
#' @param gtree tbl_graph The single-cell phylogeny
#' @return phylo The single-cell phylogeny
#' @keywords internal
to_phylo = function(gtree) {
    
    phytree = gtree %>% ape::as.phylo()
    phytree$edge.length = gtree %>% activate(edges) %>% data.frame() %>% pull(length)
    
    n_mut_root = gtree %>% activate(nodes) %>% filter(node_is_root()) %>% pull(n_mut)
    phytree$root.edge = n_mut_root
    
    return(phytree)
}

#' Annotate the direct upstream or downstream mutations on the edges
#' @param G igraph Mutation graph
#' @return igraph Mutation graph 
#' @keywords internal
label_edges = function(G) {
    
    edge_df = G %>% igraph::as_data_frame('edges') %>%
        left_join(
            G %>% igraph::as_data_frame('vertices') %>% select(from_label = label, id),
            by = c('from' = 'id')
        ) %>%
        left_join(
            G %>% igraph::as_data_frame('vertices') %>% select(to_label = label, id),
            by = c('to' = 'id')
        ) %>%
        mutate(label = paste0(from_label, '->', to_label))
    
    E(G)$label = edge_df$label
    E(G)$from_label = edge_df$from_label
    E(G)$to_label = edge_df$to_label
    
    return(G)
}

#' Annotate the direct upstream or downstream node on the edges
#' @param G igraph Mutation graph
#' @return igraph Mutation graph 
#' @keywords internal
transfer_links = function(G) {
    
    edge_df = G %>% igraph::as_data_frame('edges') %>%
            left_join(
                G %>% igraph::as_data_frame('vertices') %>% select(from_node = node, id),
                by = c('from' = 'id')
            ) %>%
            left_join(
                G %>% igraph::as_data_frame('vertices') %>% select(to_node = node, id),
                by = c('to' = 'id')
            )

    E(G)$from_node = edge_df$from_node
    E(G)$to_node = edge_df$to_node
    
    return(G)
}

#' Label the genotypes on a mutation graph
#' @param G igraph Mutation graph
#' @return igraph Mutation graph
#' @keywords internal
label_genotype = function(G) {

    id_to_label = igraph::as_data_frame(G, 'vertices') %>% {setNames(.$label, .$id)}

    # for some reason, the output from all_simple_path is out of order if supplied directly
    # V(G)$GT = igraph::all_simple_paths(G, from = 1) %>% 
    V(G)$GT = lapply(
            2:length(V(G)), 
            function(v) {dplyr::first(igraph::all_simple_paths(G, from = 1, to = v), default = NULL)}
        ) %>%
        purrr::map(as.character) %>%
        purrr::map(function(x) {
            muts = id_to_label[x]
            muts = muts[muts != '']
            paste0(muts, collapse = ',')
        }) %>%
        c(id_to_label[[1]],.) %>%
        as.character

    visit_order = setNames(1:length(V(G)), as.numeric(igraph::dfs(G, root = 1)$order))
    V(G)$clone = visit_order[as.character(as.numeric(V(G)))]
    
    return(G)
}

#' Merge adjacent set of nodes
#' @param G igraph Mutation graph
#' @param vset vector Set of adjacent vertices to merge
#' @return igraph Mutation graph
#' @keywords internal
contract_nodes = function(G, vset, node_tar = NULL, debug = F) {
    
    vset = unlist(vset)
    
    if (length(vset) == 1) {
        return(G)
    }
    
    # reorder the nodes according to graph
    vorder = V(G)$label[igraph::dfs(G, root = 1)$order]
    vset = vorder[vorder %in% vset]
    
    vset_ids = V(G)[label %in% vset]
    
    ids_new = 1:vcount(G)
    
    # the indices before do not change
    ids_new[vset_ids] = min(vset_ids)
    # indices after might need to be reset
    if (max(vset_ids) != vcount(G)) {
        ids_new[(max(vset_ids)+1):length(ids_new)] = ids_new[(max(vset_ids)+1):length(ids_new)] - length(vset_ids) + 1
    }

    G = G %>% igraph::contract(
        ids_new,
        vertex.attr.comb = list(label = function(x){paste0(sort(x), collapse = ',')}, node = "first", "ignore")
    )
    
    if (!is.null(node_tar)) {
        V(G)[min(vset_ids)]$node = node_tar
    }

    V(G)$id = 1:vcount(G)

    G = igraph::simplify(G)

    if (debug) {
        return(G)
    }
    
    G = label_edges(G)
    
    return(G)
    
}

#' Simplify the mutational history based on likelihood evidence
#' @param G igraph Mutation graph 
#' @param l_matrix matrix Mutation placement likelihood matrix (node by mutation)
#' @return igraph Mutation graph
#' @export
simplify_history = function(G, l_matrix, max_cost = 150, verbose = T) {

    # moves = data.frame()

    for (i in 1:ecount(G)) {
    
        move_opt = get_move_opt(G, l_matrix)

        if (move_opt$cost < max_cost) {

            if (move_opt$direction == 'up') {
                G = G %>% contract_nodes(c(move_opt$from_label, move_opt$to_label), move_opt$from_node) %>% transfer_links()
                msg = glue('opt_move:{move_opt$to_label}->{move_opt$from_label}, cost={signif(move_opt$cost,3)}')
            } else {
                G = G %>% contract_nodes(c(move_opt$from_label, move_opt$to_label), move_opt$to_node) %>% transfer_links()
                msg = glue('opt_move:{move_opt$from_label}->{move_opt$to_label}, cost={signif(move_opt$cost,3)}')
            }

            # moves = moves %>% rbind(move_opt %>% mutate(i = i))

            log_info(msg)
            # if (verbose) {display(msg)}
        } else {
            break()
        }
    }
    
    return(G)
}

#' Get the cost of a mutation reassignment
#' @param muts character Mutations dlimited by comma
#' @param node_ori character Name of the "from" node
#' @param node_tar character Name of the "to" node
#' @return numeric Likelihood cost of the mutation reassignment
#' @keywords internal
get_move_cost = function(muts, node_ori, node_tar, l_matrix) {

    if (muts == '') {
        return(Inf)
    }

    if (str_detect(muts, ',')) {
        muts = unlist(str_split(muts, ','))
    }

    sum(l_matrix[node_ori, muts] - l_matrix[node_tar, muts])
}

#' Get the least costly mutation reassignment 
#' @param G igraph Mutation graph
#' @param l_matrix matrix Likelihood matrix of mutation placements
#' @return numeric Lieklihood cost of performing the mutation move
#' @keywords internal
get_move_opt = function(G, l_matrix) {
    
    move_opt = G %>% igraph::as_data_frame('edges') %>%
        group_by(from) %>%
        mutate(n_sibling = n()) %>%
        ungroup() %>%
        rowwise() %>%
        mutate(
            up = get_move_cost(to_label, to_node, from_node, l_matrix),
            down = get_move_cost(from_label, from_node, to_node, l_matrix)
        ) %>%
        ungroup() %>%
        # prevent a down move if branching. Technically it's fine but graph has to be modified correctly
        mutate(down = ifelse(n_sibling > 1, Inf, down)) %>%
        reshape2::melt(measure.vars = c('up', 'down'), variable.name = 'direction', value.name = 'cost') %>%
        arrange(cost) %>%
        head(1)

    return(move_opt)
}

#' Annotate superclones on a tree
#' @export
annot_superclones = function(gtree, geno, p_min = 0.95, precision_cutoff = 0.9) {
    
    anchors = score_mut(gtree, geno, p_min = p_min) %>%
        filter(precision > precision_cutoff) %>% 
        pull(site)
    
    mut_nodes = gtree %>%
        activate(nodes) %>%
        filter(!is.na(site)) %>%
        as.data.frame() %>%
        select(name, site) %>%
        tidyr::separate_rows(site, sep = ',') %>%
        filter(site %in% anchors) %>%
        group_by(name) %>%
        summarise(
            site = paste0(site, collapse = ','),
            .groups = 'drop'
        )
    
    superclone_dict = gtree %>% 
        mut_to_tree(mut_nodes) %>%
        activate(nodes) %>%
        as.data.frame() %>%
        mutate(clone = as.integer(as.factor(GT))) %>%
        {setNames(.$clone, .$name)}
    
    gtree = gtree %>% 
        activate(nodes) %>%
        mutate(superclone = superclone_dict[name])
    
    return(gtree)
}

#' Score mutations on goodness of fit on the tree
#' @keywords internal
score_mut = function(gtree, geno, p_min = 0.95) {

    mut_scores = gtree %>%
        activate(nodes) %>%
        filter(!is.na(site)) %>%
        as.data.frame() %>%
        distinct(site, id, name) %>%
        rename(node = id) %>%
        tidyr::separate_rows(site, sep = ',') %>%
        group_by(site, node, name) %>%
        summarise(
            score_mut_helper(gtree, geno, site, node, p_min),
            .groups = 'drop'
        ) %>%
        arrange(-precision)
    
    return(mut_scores)
    
}

#' Helper for mutation scoring
#' @keywords internal
score_mut_helper = function(gtree, geno, s, v, p_min = 0.95) {
    gtree %>%
        activate(nodes) %>%
        mutate(seq = bfs_rank(root = v)) %>%
        data.frame %>%
        filter(leaf) %>%
        mutate(
            p_0 = unlist(geno[name, s]),
            p_1 = 1 - p_0
        ) %>%
        mutate(
            is_desc = !is.na(seq),
            p = ifelse(is_desc, p_1, p_0)
        ) %>%
        summarise(
            p = mean(p),
            TP = sum(is_desc & p_1 > p_min),
            FP = sum((!is_desc) & p_1 > p_min),
            precision = TP/(TP + FP)
        ) %>%
        tibble()
}
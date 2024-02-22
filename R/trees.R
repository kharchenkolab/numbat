# from phangorn
#' UPGMA and WPGMA clustering
#'
#' @param D A distance matrix.
#' @param method The agglomeration method to be used. This should be (an
#' unambiguous abbreviation of) one of "ward", "single", "complete", "average",
#' "mcquitty", "median" or "centroid". The default is "average".
#' @param \dots Further arguments passed to or from other methods.
upgma <- function(D, method = "average", ...) {
  DD <- as.dist(D)
  hc <- hclust(DD, method = method, ...)
  result <- ape::as.phylo(hc)
  result <- reorder(result, "postorder")
  result
}

#' Mark the tumor lineage of a phylogeny
#' @param gtree tbl_graph Single-cell phylogeny
#' @return tbl_graph Phylogeny annotated with tumor versus normal compartment
#' @keywords internal
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
            filter(leaf & seq > 0) %>%
            pull(mut_burden) %>%
            sum
        }
    )

    tumor_root = mut_nodes[which.max(mut_burdens)]
    
    gtree = gtree %>%
        activate(nodes) %>%
        mutate(
            seq = bfs_rank(root = tumor_root),
            compartment = ifelse(seq > 0, 'tumor', 'normal'),
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

    gtree = annotate_tree(tree, P)

    return(list('gtree' = gtree, 'l_matrix' = l_matrix))
}

#' Get a tidygraph tree with simplified mutational history. 
#' @description Specify either max_cost or n_cut. 
#' max_cost works similarly as h and n_cut works similarly as k in stats::cutree.
#' The top-level normal diploid clone is always included.
#' @param tree phylo Single-cell phylogenetic tree
#' @param P matrix Genotype probability matrix
#' @param max_cost numeric Likelihood threshold to collapse internal branches
#' @param n_cut integer Number of cuts on the phylogeny to define subclones
#' @return tbl_graph Phylogeny annotated with branch lengths and mutation events
#' @export
get_gtree = function(tree, P, n_cut = 0, max_cost = 0) {
    
    # calculate L matrix and create tidygraph object
    tree_post = get_tree_post(tree, P)

    # simplify mutational history
    G_m = get_mut_graph(tree_post$gtree)  %>% 
        simplify_history(tree_post$l_matrix, max_cost = max_cost, n_cut = n_cut) %>% 
        label_genotype()

    mut_nodes = G_m %>% igraph::as_data_frame('vertices') %>% 
        select(name = node, site = label, clone = clone, GT = GT)

    # update tree
    gtree = mut_to_tree(tree_post$gtree, mut_nodes)
    gtree = mark_tumor_lineage(gtree)
    
    return(gtree)
    
}

#' Annotate the direct upstream or downstream node on the edges
#'
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
#'
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

#' Annotate the direct upstream or downstream mutations on the edges
#'
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

#' Merge adjacent set of nodes
#'
#' @param G igraph Mutation graph
#' @param vset vector Set of adjacent vertices to merge
#' @return igraph Mutation graph
#' @keywords internal
contract_nodes = function(G, vset, node_tar = NULL, debug = FALSE) {
    
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
#'
#' @param G igraph Mutation graph 
#' @param l_matrix matrix Mutation placement likelihood matrix (node by mutation)
#' @return igraph Mutation graph
#' @keywords internal
simplify_history = function(G, l_matrix, max_cost = 150, n_cut = 0, verbose = TRUE) {

    if (n_cut > 0) {
        max_cost = Inf
    }

    for (i in 1:ecount(G)) {
    
        move_opt = get_move_opt(G, l_matrix)

        if (move_opt$cost < max_cost & ecount(G) > n_cut) {

            if (move_opt$direction == 'up') {
                G = G %>% contract_nodes(c(move_opt$from_label, move_opt$to_label), move_opt$from_node) %>% transfer_links()
                msg = glue('opt_move:{move_opt$to_label}->{move_opt$from_label}, cost={signif(move_opt$cost,3)}')
            } else {
                G = G %>% contract_nodes(c(move_opt$from_label, move_opt$to_label), move_opt$to_node) %>% transfer_links()
                msg = glue('opt_move:{move_opt$from_label}->{move_opt$to_label}, cost={signif(move_opt$cost,3)}')
            }

            log_info(msg)
            
        } else {
            break()
        }
    }
    
    return(G)
}

#' Get the cost of a mutation reassignment
#'
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
#'
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
        as.data.table %>%
        data.table::melt(measure.vars = c('up', 'down'), variable.name = 'direction', value.name = 'cost') %>%
        arrange(cost) %>%
        head(1)

    return(move_opt)
}

#' Get ordered tips from a tree
#' @keywords internal
get_ordered_tips = function(tree) {
    is_tip <- tree$edge[,2] <= length(tree$tip.label)
    ordered_tips <- tree$edge[is_tip, 2]
    tree$tip.label[ordered_tips]
}

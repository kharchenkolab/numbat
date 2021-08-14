
label_edges = function(G) {
    
    edge_df = G %>% as_data_frame('edges') %>%
        left_join(
            G %>% as_data_frame('vertices') %>% select(from_label = label, id),
            by = c('from' = 'id')
        ) %>%
        left_join(
            G %>% as_data_frame('vertices') %>% select(to_label = label, id),
            by = c('to' = 'id')
        ) %>%
        mutate(label = paste0(from_label, '->', to_label))
    
    E(G)$label = edge_df$label
    E(G)$from_label = edge_df$from_label
    E(G)$to_label = edge_df$to_label
    
    return(G)
}

transfer_links = function(G) {
    
    edge_df = G %>% as_data_frame('edges') %>%
            left_join(
                G %>% as_data_frame('vertices') %>% select(from_node = node, id),
                by = c('from' = 'id')
            ) %>%
            left_join(
                G %>% as_data_frame('vertices') %>% select(to_node = node, id),
                by = c('to' = 'id')
            )

    E(G)$from_node = edge_df$from_node
    E(G)$to_node = edge_df$to_node
    
    return(G)
}

# merge adjacent set of nodes
contract_nodes = function(G, vset, node_tar = NULL, debug = F) {
    
    vset = unlist(vset)
    
    if (length(vset) == 1) {
        return(G)
    }
    
    # reorder the nodes according to graph
    vorder = V(G)$label[dfs(G, root = 1)$order]
    vset = vorder[vorder %in% vset]
    
    vset_ids = V(G)[label %in% vset]
    
    ids_new = 1:vcount(G)
    
    # the indices before do not change
    ids_new[vset_ids] = min(vset_ids)
    # indices after might need to be reset
    if (max(vset_ids) != vcount(G)) {
        ids_new[(max(vset_ids)+1):length(ids_new)] = ids_new[(max(vset_ids)+1):length(ids_new)] - length(vset_ids) + 1
    }
    
    if (debug) {
        print(vset)
        print(vset_ids)
        print(ids_new)
    }

    G = G %>% contract(
        ids_new,
        vertex.attr.comb = list(label = function(x){paste0(sort(x), collapse = ',')}, node = "first", "ignore")
    )
    
    if (!is.null(node_tar)) {
        V(G)[min(vset_ids)]$node = node_tar
    }

    V(G)$id = 1:vcount(G)

    G = simplify(G)
    
    G = label_edges(G)
    
    return(G)
    
}

get_move_cost = function(muts, node_ori, node_tar, l_matrix) {
    if (str_detect(muts, ',')) {
        muts = unlist(str_split(muts, ','))
    }
    sum(l_matrix[node_ori, muts] - l_matrix[node_tar, muts])
}

get_move_opt = function(G, l_matrix, verbose = T) {
    
    move_opt = G %>% as_data_frame('edges') %>%
        rowwise() %>%
        mutate(
            cost = get_move_cost(to_label, to_node, from_node, l_matrix)
        ) %>%
        ungroup() %>%
        arrange(cost) %>%
        head(1)
    
    if (verbose) {
        display(glue('opt_move:{move_opt$to_label}->{move_opt$from_label}, cost={signif(move_opt$cost,3)}'))
    }
    
    return(move_opt)
}

plot_mut_tree = function(G) {
    
    par(mar = c(0, 0, 0, 0))

    plot(G, 
         edge.label = NA,
         edge.label.cex = 0.5,
         vertex.label.cex = 0.5,
         edge.label.dist = 1,
         vertex.label.degree = -pi/2,
         vertex.label.dist = 1, 
         layout = layout_as_tree(G) %>% {-.[,2:1]},
         vertex.size = 2,
         edge.arrow.size = 0.25,
         asp = 0.25
    )
}
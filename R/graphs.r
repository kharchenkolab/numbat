##### R implementation of ScisTree perfect phylogeny #####

#' @description score a tree based on maximum likelihood
#' @param tree phylo object
#' @param P genotype probability matrix
#' @export
score_tree = function(tree, P) {
    
    tree = reorder(tree, order = 'postorder')
        
    n = nrow(P)
    m = ncol(P)

    logQ = matrix(nrow = tree$Nnode * 2 + 1, ncol = m)

    logP_0 = log(P)
    logP_1 = log(1-P)
    
    node_order = c(tree$edge[,2], n+1)
    node_order = node_order[node_order > n]
    
    logQ[1:n,] = logP_1 - logP_0

    children_dict = phangorn:::allChildren(tree)

    # for (i in node_order) {
    #     children = children_dict[[i]]
    #     logQ[i,] = logQ[children[1],] + logQ[children[2],]
    # }

    logQ = CgetQ(logQ, children_dict, node_order)
    
    l_matrix = sweep(logQ, 2, colSums(logP_0), FUN = '+')
    
    l_tree = sum(apply(l_matrix, 2, max))
    
    return(list('l_tree' = l_tree, 'logQ' = logQ, 'l_matrix' = l_matrix))
    
}

#' @description maximum likelihood tree search via NNI
#' @param tree_init iintial tree as phylo object
#' @param P genotype probability matrix
#' @export
perform_nni = function(tree_init, P, max_iter = 100, ncores = 20, eps = 0.01, verbose = TRUE) {

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

        neighbours = nni(tree_current)

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
            msg = glue('Iter {i} {max_current} {signif(unname(runtime[3]),2)}')
            message(msg)
            log_info(msg)
        }
    }

    class(tree_list) = 'multiPhylo'

    return(tree_list)
}

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

parse_scistree = function(out_file, geno, joint_post) {
    
    MLtree = fread(out_file, skip = 'Constructed single cell phylogeny', fill=TRUE, sep = ':', nrow = 1) %>%
        colnames %>% .[2] %>% paste0(';') %>%
        treeio::read.tree(text = .) %>%
        ape::compute.brlen(method = 1)

    NJtree = grep('Neighbor joining tree from noisy genotypes', readLines(out_file), value = T) %>%
        str_remove('Neighbor joining tree from noisy genotypes: ') %>%
        paste0(';') %>%
        treeio::read.tree(text = .) 

    NJtree$tip.label = colnames(geno)[as.integer(NJtree$tip.label)]

    cnv_order = grep('Mutation tree', readLines(out_file), value = T) %>% str_remove_all('\\(|\\)|\\^|Mutation tree:|,| ') %>%
        str_split('#') %>% unlist %>% str_subset('') 
    
    options(warn = -1)
    
    Gopt = fread(out_file, skip = 'Imputed genotypes', header = F, sep = ' ', nrows = nrow(geno)) %>%
        select(-V1) %>%
        mutate(V2 = as.integer(purrr::map(str_split(V2, '\t'),2))) %>%
        set_colnames(colnames(geno)) %>%
        mutate(seg = rownames(geno)) %>% 
        reshape2::melt(id.var = 'seg', variable.name = 'cell', value.name = 'p_cnv') %>%
        mutate(logBF = ifelse(p_cnv == 0, -5, 5)) %>%
        left_join(
            joint_post %>% distinct(seg, seg_label, cnv_state),
            by = 'seg'
        )
    options(warn = 0)
    
    return(list('MLtree' = MLtree, 'cnv_order' = cnv_order, 'G' = Gopt, 'NJtree' = NJtree))
}

label_mut_tree = function(G, mut_assign) {

    # fix the root node
    if (all(V(G)$label != ' ')) {
        
        V(G)$id = V(G)$id + 1
        
        G = G %>% add_vertices(1, attr = list('label' = ' ', 'id' = 1)) %>%
            add_edges(c(length(V(G))+1,1)) 
    }

    G = G %>% as_tbl_graph() %>% arrange(id != 1) %>%
        mutate(id = 1:length(V(G))) %>%
        mutate(label = ifelse(label == ' ', '', label)) 

    G = label_edges(G)

    for (i in 1:nrow(mut_assign)) {
        node_set = unlist(str_split(mut_assign[i,]$site, ','))
        G = G %>% contract_nodes(node_set)
    }

    V(G)$node = G %>%
        igraph::as_data_frame('vertices') %>%
        left_join(
            mut_assign %>% rename(node = name),
            by = c('label' = 'site')
        ) %>%
        pull(node)

    G = G %>% transfer_links()

    return(G)
}

get_mut_tree = function(gtree, mut_assign) {

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
            mut_assign %>% rename(node = name),
            by = c('label' = 'site')
        ) %>%
        pull(node)

    G = G %>% transfer_links()

    return(G)

}

# compute site branch likelihood
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

get_tree_post = function(tree, P) {
    
    sites = colnames(P)
    n = nrow(P)
    tree_stats = score_tree(tree, P)

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



mut_to_tree = function(gtree, mut_nodes) {
    
    gtree = gtree %>%
        activate(nodes) %>%
        select(-any_of(c('n_mut', 'l', 'site'))) %>%
        left_join(
            mut_nodes %>%
                mutate(n_mut = unlist(purrr::map(str_split(site, ','), length))) %>%
                select(name, n_mut, site),
            by = 'name'
        ) %>%
        mutate(n_mut = ifelse(is.na(n_mut), 0, n_mut))

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
    
    return(gtree)
}

to_phylo = function(gtree) {
    
    phytree = gtree %>% ape::as.phylo()
    phytree$edge.length = gtree %>% activate(edges) %>% data.frame() %>% pull(length)
    
    n_mut_root = gtree %>% activate(nodes) %>% filter(node_is_root()) %>% pull(n_mut)
    phytree$root.edge = -n_mut_root
    
    return(phytree)
}

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

label_genotype = function(G) {
    id_to_label = igraph::as_data_frame(G, 'vertices') %>% {setNames(.$label, .$id)}

    # for some reason, the output from all_simple_path is out of order if supplied directly
    # V(G)$GT = igraph::all_simple_paths(G, from = 1) %>% 
    V(G)$GT = lapply(
            2:length(V(G)), 
            function(v) {first(igraph::all_simple_paths(G, from = 1, to = v), default = NULL)}
        ) %>%
        purrr::map(as.character) %>%
        purrr::map(function(x) {
            muts = id_to_label[x]
            muts = muts[muts != '']
            paste0(muts, collapse = ',')
        }) %>%
        c(id_to_label[[1]],.) %>%
        as.character
    
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

    G = G %>% contract(
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

get_move_cost = function(muts, node_ori, node_tar, l_matrix) {

    if (muts == '') {
        return(Inf)
    }

    if (str_detect(muts, ',')) {
        muts = unlist(str_split(muts, ','))
    }

    sum(l_matrix[node_ori, muts] - l_matrix[node_tar, muts])
}

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


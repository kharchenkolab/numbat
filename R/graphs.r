library(igraph)
library(tidygraph)
library(ggtree)

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
        mutate(V2 = as.integer(map(str_split(V2, '\t'),2))) %>%
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

read_mut_tree = function(tree_file, mut_assign) {

    G = igraph::read_graph(tree_file, format = 'gml')
    V(G)$id = 1:length(igraph::V(G))
    V(G)$label[V(G)$label == ' '] = ''

    G = label_edges(G)

    for (i in 1:nrow(mut_assign)) {
        node_set = unlist(str_split(mut_assign[i,]$site, ','))
        G = G %>% contract_nodes(node_set)
    }

    V(G)$node = G %>%
        as_data_frame('vertices') %>%
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

get_tree_post = function(MLtree, geno) {
    
    gtree = MLtree %>%
        ape::ladderize() %>%
        as_tbl_graph() %>%
        mutate(
            leaf = node_is_leaf(),
            root = node_is_root(),
            depth = bfs_dist(root = 1),
            id = row_number()
        )

    sites = rownames(geno)

    nodes = mclapply(
            mc.cores = length(sites),
            sites,
            function(site) {
                gtree %>%
                data.frame() %>%
                rowwise() %>%
                mutate(l = l_s_v(id, site, gtree, geno)) %>%
                ungroup() %>%
                mutate(site = site)
            }
        ) %>%
        bind_rows()
    
    nodes = nodes %>% group_by(site) %>%
        mutate(
            p = exp(l - matrixStats::logSumExp(l)),
            mle = p == max(p)
        ) %>%
        ungroup()

    l_matrix = nodes %>% 
        reshape2::dcast(name ~ site, value.var = 'l') %>%
        tibble::column_to_rownames('name')

    # mutation assignments
    mut_nodes = nodes %>% 
        filter(mle) %>%
        group_by(name) %>%
        summarise(
            site = paste0(sort(site), collapse = ','),
            n_mut = n(),
            l = sum(l),
            .groups = 'drop'
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

    # annotate the tree with mutation assignment
    gtree = mut_to_tree(gtree, mut_nodes)

    return(list('mut_nodes' = mut_nodes, 'gtree' = gtree, 'nodes' = nodes, 'l_matrix' = l_matrix))
}



mut_to_tree = function(gtree, mut_nodes) {
    
    gtree = gtree %>%
        activate(nodes) %>%
        select(-any_of(c('n_mut', 'l', 'site'))) %>%
        left_join(
            mut_nodes %>%
                mutate(n_mut = unlist(map(str_split(site, ','), length))) %>%
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
        ))
    
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

label_genotype = function(G) {
    id_to_label = igraph::as_data_frame(G, 'vertices') %>% {setNames(.$label, .$id)}

    V(G)$GT = igraph::all_simple_paths(G, from = 1) %>% 
        map(as.character) %>%
        map(function(x) {
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

simplify_history = function(G, l_matrix, max_cost = 150, verbose = T) {

    moves = data.frame()

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

            moves = moves %>% rbind(move_opt %>% mutate(i = i))

            if (verbose) {
                log_info(msg)
                display(msg)
            }
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
    
    move_opt = G %>% as_data_frame('edges') %>%
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
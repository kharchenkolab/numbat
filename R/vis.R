########################### Visualizations ############################
# default color palette
pal = RColorBrewer::brewer.pal(n = 8, 'Set1')
getPalette = colorRampPalette(pal)

#' @export
cnv_colors = c("neu" = "gray", 
        "neu_up" = "darkgray", "neu_down" = "gray",
        "del_up" = "royalblue", "del_down" = "darkblue", 
        "loh_up" = "darkgreen", "loh_down" = "olivedrab4",
        "amp_up" = "red", "amp_down" = "tomato3",
        "del_1_up" = "royalblue", "del_1_down" = "darkblue", 
        "loh_1_up" = "darkgreen", "loh_1_down" = "olivedrab4",
        "amp_1_up" = "red", "amp_1_down" = "tomato3",
        "del_2_up" = "royalblue", "del_2_down" = "darkblue", 
        "loh_2_up" = "darkgreen", "loh_2_down" = "olivedrab4",
        "amp_2_up" = "red", "amp_2_down" = "tomato3",
        "del_up_1" = "royalblue", "del_down_1" = "darkblue", 
        "loh_up_1" = "darkgreen", "loh_down_1" = "olivedrab4",
        "amp_up_1" = "red", "amp_down_1" = "tomato3",
        "del_up_2" = "royalblue", "del_down_2" = "darkblue", 
        "loh_up_2" = "darkgreen", "loh_down_2" = "olivedrab4",
        "amp_up_2" = "red", "amp_down_2" = "tomato3",
        "bamp" = "salmon", "bdel" = "skyblue",
        "amp" = "tomato3", "loh" = "olivedrab4", "del" = "royalblue",
        "theta_up" = "darkgreen", "theta_down" = "olivedrab4",
        "theta_1_up" = "darkgreen", "theta_1_down" = "olivedrab4",
        "theta_2_up" = "darkgreen", "theta_2_down" = "olivedrab4",
        "theta_up_1" = "darkgreen", "theta_down_1" = "olivedrab4",
        "theta_up_2" = "darkgreen", "theta_down_2" = "olivedrab4",
        '0|1' = 'red', '1|0' = 'blue',
        'major' = '#66C2A5', 'minor' = '#FC8D62'
    )

#' @export
cnv_labels = names(cnv_colors) %>%
    str_remove_all('_') %>% 
    str_to_upper() %>%
    str_replace('UP', '(major)') %>%
    str_replace('DOWN', '(minor)') %>%
    str_replace('LOH', 'CNLoH') %>%
    setNames(names(cnv_colors))

#' plot a pseudobulk HMM profile
#' @param bulk dataframe Pseudobulk profile
#' @param use_pos logical Use marker position instead of index as x coordinate
#' @param allele_only logical Only plot alleles
#' @param min_LLR numeric LLR threshold for event filtering
#' @param min_depth numeric Minimum coverage depth for a SNP to be plotted
#' @param exp_limit numeric Expression logFC axis limit
#' @param phi_mle logical Whether to plot estimates of segmental expression fold change
#' @param theta_roll logical Whether to plot rolling estimates of allele imbalance
#' @param dot_size numeric Size of marker dots
#' @param dot_alpha numeric Transparency of the marker dots
#' @param legend logical Whether to show legend
#' @return ggplot Plot of pseudobulk HMM profile
#' @export
plot_psbulk = function(
    bulk, use_pos = FALSE, allele_only = FALSE, min_LLR = 20, min_depth = 8, exp_limit = 2, 
    phi_mle = TRUE, theta_roll = FALSE, dot_size = 0.8, dot_alpha = 0.5, legend = FALSE
    ) {

    if (!all(c('state_post', 'cnv_state_post') %in% colnames(bulk))) {
        bulk = bulk %>%
            mutate(
                state_post = state,
                cnv_state_post = cnv_state
            )
    }

    # filter events by LLR
    if (min_LLR != 0) {
        bulk = bulk %>% mutate(
            LLR = ifelse(is.na(LLR), 0, LLR),
            cnv_state_post = ifelse(LLR < min_LLR, 'neu', cnv_state_post),
            state_post = ifelse(LLR < min_LLR, 'neu', state_post)
        )
    }

    # mark clonal LOH
    if ('loh' %in% colnames(bulk)) {
        bulk = bulk %>% mutate(state_post = ifelse(loh, 'del', state_post))
    }

    if (use_pos) {
        marker = 'POS'
        marker_label = 'Genomic position'
    } else {
        marker = 'snp_index'
        marker_label = 'SNP index'
    }

    # fix retest states 
    bulk = bulk %>% 
        mutate(
            theta_level = ifelse(str_detect(state_post, '_2'), 2, 1),
            state_post = ifelse(
                cnv_state_post %in% c('amp', 'loh', 'del'),
                ifelse(p_up > 0.5, paste0(cnv_state_post, '_', theta_level, '_', 'up'), paste0(cnv_state_post, '_', theta_level, '_', 'down')),
                state_post
        ))

    # correct for baseline bias
    if (!allele_only) {
        bulk = bulk %>% mutate(logFC = logFC - mu)
    }

    D = bulk %>% 
        mutate(logFC = ifelse(logFC > exp_limit | logFC < -exp_limit, NA, logFC)) %>%
        mutate(pBAF = ifelse(DP >= min_depth, pBAF, NA)) %>%
        mutate(pHF = pBAF) %>%
        reshape2::melt(measure.vars = c('logFC', 'pHF'))

    if (allele_only) {
        D = D %>% filter(variable == 'pHF')
    }

    p = ggplot(
            D,
            aes(x = get(marker), y = value, color = state_post),
            na.rm=TRUE
        ) +
        geom_point(
            aes(
                shape = str_detect(state_post, '_2'),
                alpha = str_detect(state_post, '_2')
            ),
            size = dot_size,
            na.rm = TRUE
        ) +
        geom_hline(
            data = data.frame(y = c(0,1), variable = 'pHF'),
            aes(yintercept = y),
            size = 0
        ) +
        geom_hline(
            data = data.frame(y = c(-exp_limit, exp_limit), variable = 'logFC'),
            aes(yintercept = y),
            size = 0
        ) +
        scale_alpha_discrete(range = c(dot_alpha, 1)) +
        scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 15)) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, 'mm'),
            panel.spacing.y = unit(1, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
            strip.background = element_blank(),
            axis.text.x = element_blank()
        ) +
        facet_grid(variable ~ CHROM, scale = 'free', space = 'free_x') +
        # scale_x_continuous(expand = expansion(add = 5)) +
        scale_color_manual(
            values = cnv_colors,
            limits = force,
            labels = cnv_labels,
            na.translate = F
        ) +
        guides(color = guide_legend(title = "", override.aes = aes(size = 3)), fill = 'none', alpha = 'none', shape = 'none') +
        xlab(marker) +
        ylab('')

    if (!legend) {
        p = p + guides(color = 'none', fill = 'none', alpha = 'none', shape = 'none')
    }

    if (phi_mle) {
        segs = bulk %>% 
            distinct(CHROM, seg, seg_start, seg_start_index, seg_end, seg_end_index, phi_mle) %>%
            mutate(variable = 'logFC') %>%
            filter(log2(phi_mle) < exp_limit)

        if (use_pos) {
            start = 'seg_start'
            end = 'seg_end'
        } else {
            start = 'seg_start_index'
            end = 'seg_end_index'
        }

        p = p + geom_segment(
            inherit.aes = FALSE,
            data = segs,
            aes(x = get(start), xend = get(end), y = log2(phi_mle), yend = log2(phi_mle)),
            color = 'darkred',
            size = 0.5
        ) +
        geom_hline(data = data.frame(variable = 'logFC'), aes(yintercept = 0), color = 'gray30', linetype = 'dashed')
    } else if (!allele_only) {
        p = p + geom_line(
            inherit.aes = FALSE,
            data = bulk %>% mutate(variable = 'logFC') %>% filter(log2(phi_mle_roll) < exp_limit),
            aes(x = get(marker), y = log2(phi_mle_roll), group = '1'),
            color = 'darkred',
            size = 0.35
        ) +
        geom_hline(data = data.frame(variable = 'logFC'), aes(yintercept = 0), color = 'gray30', linetype = 'dashed')
    }

    if (theta_roll) {
        p = p + 
            geom_line(
                inherit.aes = FALSE,
                data = D %>% mutate(variable = 'pHF'),
                aes(x = snp_index, y = 0.5 - theta_hat_roll, color = paste0(cnv_state_post, '_down')),
                # color = 'black',
                size = 0.35
            ) +
            geom_line(
                inherit.aes = FALSE,
                data = D %>% mutate(variable = 'pHF'),
                aes(x = snp_index, y = 0.5 + theta_hat_roll, color = paste0(cnv_state_post, '_up')),
                # color = 'gray',
                size = 0.35
            )
    } 

    p = p + xlab(marker_label)
    
    return(p)
}

#' plot a group of pseudobulks HMM profile
#' @param bulks dataframe Pseudobulk profiles annotated with "sample" column
#' @param ncol integer Number of columns
#' @param title logical Whether to add titles to individual plots
#' @param ... additional parameters passed to plot_psbulk()
#' @return a ggplot object
#' @export
plot_bulks = function(
    bulks, ..., ncol = 1, title = TRUE
    ) {

    if (!'sample' %in% colnames(bulks)) {
        bulks$sample = 1
    }

    options(warn = -1)
    plot_list = bulks %>%
        split(.$sample) %>%
        lapply(
            function(bulk) {

                sample = unique(bulk$sample)
                n_cells = unique(bulk$n_cells)

                p = plot_psbulk(
                        bulk, ...
                    ) + 
                    theme(
                        title = element_text(size = 8),
                        axis.text.x = element_blank(),
                        axis.title = element_blank()
                    )

                if (title) {
                    if (is.null(n_cells)) {
                        title_text = sample
                    } else {
                        title_text = glue('{sample} (n={n_cells})')
                    }
                    p = p + ggtitle(title_text)
                }
                    
                return(p)
            }
        )
    options(warn = 0)

    panel = wrap_plots(plot_list, ncol = ncol, guides = 'collect')

    return(panel)
}

#' plot mutational history
#' @param G igraph Mutation history graph
#' @param clone_post dataframe Clone assignment posteriors
#' @param edge_label_size numeric Size of edge label
#' @param node_label_size numeric Size of node label
#' @param node_size numeric Size of nodes
#' @param arrow_size numeric Size of arrows
#' @param edge_label logical Whether to label edges
#' @param node_label logical Whether to label nodes
#' @param horizontal logical Whether to use horizontal layout 
#' @param show_clone_size logical Whether to show clone size
#' @param show_distance logical Whether to show evolutionary distance between clones
#' @param pal named vector Node colors
#' @return ggplot object
#' @export
plot_mut_history = function(
        G, clone_post = NULL, 
        edge_label_size = 4, node_label_size = 6, node_size = 10, arrow_size = 2,
        show_clone_size = TRUE, show_distance = TRUE, legend = TRUE,
        edge_label = TRUE, node_label = TRUE, horizontal = TRUE, pal = NULL
    ) {

    G = label_genotype(G)

    if (!is.null(clone_post)) {
        clone_sizes = clone_post %>%
            count(clone_opt) %>%
            {setNames(.$n, .$clone_opt)}
        V(G)$size = unname(clone_sizes[V(G)$clone])
    } else {
        V(G)$size = 1
        show_clone_size = FALSE
    }

    if (is.null(pal)) {
        getPalette = colorRampPalette(RColorBrewer::brewer.pal(n = 8, 'Set1'))
        pal = c('gray', getPalette(length(V(G))))
    }

    G_df = G %>% as_tbl_graph() %>% mutate(clone = factor(clone)) 
    
    if (!'superclone' %in% colnames(as.data.frame(activate(G_df, 'nodes')))) {
        G_df = G_df %>% mutate(superclone = clone)
    }

    # add edge length
    if (show_distance) {

        if ((!'length' %in% colnames(as.data.frame(activate(G_df, 'edges'))))) {
            G_df = G_df %>% activate(edges) %>%
                mutate(n_mut = unlist(purrr::map(str_split(to_label, ','), length))) %>%
                mutate(length = n_mut)
        }

    } else {
        G_df = G_df %>% activate(edges) %>% mutate(length = 1)
    }

    if (!edge_label) {
        G_df = G_df %>% activate(edges) %>% mutate(to_label = '')
    }

    p = G_df %>% 
        ggraph(
            layout = 'dendrogram',
            length = length
        ) + 
        geom_edge_elbow(
            aes(label = str_trunc(to_label, 20, side = 'center')),
            vjust = -1,
            hjust = 0,
            arrow = arrow(length = unit(arrow_size, "mm")),
            end_cap = circle(4, 'mm'),
            start_cap = circle(4, 'mm'),
            label_size = edge_label_size
        ) + 
        theme_void() +
        scale_x_continuous(expand = expansion(0.2)) +
        scale_color_manual(values = pal, limits = force) +
        guides(color = 'none')

    if (!legend) {
        p = p + guides(size = 'none')
    }

    size_max = G_df %>% activate(nodes) %>% pull(size) %>% max

    if (show_clone_size) {
        p = p + geom_node_point(
                aes(color = as.factor(superclone), size = size)
            ) + 
            scale_size_continuous(
                range = c(0, node_size),
                limits = c(0, size_max),
                name = 'Cells'
            )
    } else {
        p = p + geom_node_point(
                aes(color = as.factor(superclone)),
                size = node_size
            )
    }

    if (node_label) {
        p = p + geom_node_text(aes(label = clone), size = node_label_size)
    }

    if (horizontal) {
        p = p + coord_flip() + scale_y_reverse(expand = expansion(0.2))
    } else {
        p = p + scale_y_continuous(expand = expansion(0.2))
    }

    return(p)
}

#' plot single-cell CNV calls along with the clonal phylogeny
#' @param gtree tbl_graph The single-cell phylogeny
#' @param joint_post dataframe Joint single cell CNV posteriors
#' @param segs_consensus datatframe Consensus segment dataframe
#' @param annot named vector Cell annotations, keys are cell names
#' @param pal_annot named vector Colors for cell annotations
#' @param annot_title character Legend title for the annotation bar 
#' @param annot_scale ggplot scale Color scale for the annotation bar 
#' @param clone_dict named vector Clone annotations, mapping from cell name to clones
#' @param clone_bar logical Whether to display clone bar plot
#' @param clone_title character Legend title for the clone bar
#' @param clone_legend logical Whether to display the clone legend
#' @param tree_height numeric Relative height of the phylogeny plot
#' @param p_min numeric Probability threshold to display CNV calls
#' @param line_width numeric Line width for CNV heatmap
#' @param branch_width numeric Line width in the phylogeny
#' @param tip_length numeric Length of tips in the phylogeny
#' @param clone_line logical Whether to display borders for clones in the heatmap
#' @param superclone logical Wehether to display superclones in the clone bar
#' @return ggplot panel
#' @export
plot_phylo_heatmap = function(
        gtree, joint_post, segs_consensus, 
        annot = NULL, pal_annot = NULL, annot_title = 'Annotation', annot_scale = NULL,
        clone_dict = NULL, clone_bar = FALSE, pal_clone = NULL, clone_title = 'Genotype', clone_legend = TRUE, 
        p_min = 0.9, line_width = 0.1, tree_height = 1, branch_width = 0.2, tip_length = 0.2, 
        clone_line = FALSE, superclone = FALSE, segs_loh = NULL
    ) {
    
    # make sure chromosomes are in order
    joint_post = joint_post %>% mutate(CHROM = as.integer(as.character(CHROM)))
    segs_consensus = segs_consensus %>% mutate(CHROM = as.integer(as.character(CHROM)))

    # if no multi allelic CNVs
    if (!'n_states' %in% colnames(joint_post)) {
        joint_post = joint_post %>% mutate(
            n_states = ifelse(cnv_state == 'neu', 1, 0), 
            cnv_states = cnv_state
        )
    } else {
        # only keep one record per CNV and color by most likely state
        joint_post = joint_post %>% 
            group_by(cell, CHROM, seg_end, seg_start) %>%
            mutate(p_cnv = sum(p_cnv)) %>%
            ungroup() %>%
            distinct(cell, CHROM, seg_end, seg_start, .keep_all = T) %>%
            mutate(cnv_state = ifelse(n_states > 1, cnv_state_map, cnv_state))
    }

    if (!'clone' %in% colnames(data.frame(activate(gtree, 'nodes')))) {
        gtree = gtree %>%
            activate(nodes) %>%
            mutate(clone = as.integer(factor(GT)))
    }

    if (superclone) {
        gtree = gtree %>% activate(nodes) %>% mutate(clone = superclone)
    }

    gtree = mark_tumor_lineage(gtree)

    gtree = gtree %>% activate(edges) %>% mutate(length = ifelse(leaf, pmax(length, tip_length), length))

    # plot phylogeny 
    p_tree = gtree %>% 
            to_phylo() %>%
            ggtree(ladderize = T, size = branch_width) +
            # geom_rootedge(size = branch_width) +
            theme(
                plot.margin = margin(0,1,0,0, unit = 'mm'),
                axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank(),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                # axis.text.y = element_text(size = 5)
                axis.text.y = element_blank(),
                panel.background = element_rect(fill = "transparent",colour = NA),
                plot.background = element_rect(fill = "transparent", color = NA)
            ) +
            guides(color = 'none')

    # order the cells
    cell_order = p_tree$data %>% filter(isTip) %>% arrange(y) %>% pull(label)

    joint_post = joint_post %>% 
        mutate(cell = factor(cell, cell_order)) %>%
        mutate(cell_index = as.integer(droplevels(cell))) 

    # add clone lines
    if (clone_line) {

        leafs = res[[sample]]$gtree %>% 
                activate(nodes) %>% 
                filter(leaf) %>%
                as.data.frame()

        clones = unique(leafs$clone)
        clones = clones[clones != 1]

        clone_indices = sapply(
            clones,
            function(c) {                
                
                clone_cells = leafs %>% filter(clone == c) %>% pull(name)
                
                first_clone_index = which(cell_order %in% clone_cells)[1]

                return(first_clone_index)
                
            }
        )
    } else {
        clone_indices = c()
    }

    # add tumor vs normal line
    tumor_cells = gtree %>% 
        activate(nodes) %>% filter(leaf) %>%
        as.data.frame() %>% 
        filter(compartment == 'tumor') %>%
        pull(name)

    first_tumor_index = which(cell_order %in% tumor_cells)[1]

    chrom_labeller <- function(chr){
        chr[chr %in% c(19, 21, 22)] = ''
        return(chr)
    }

    # plot CNVs
    p_segs = ggplot(
            joint_post %>% mutate(
                cnv_state = ifelse(cnv_state == 'neu', NA, cnv_state)
            )
        ) +
        theme_classic() +
        geom_segment(
            aes(x = seg_start, xend = seg_end, y = cell_index, yend = cell_index, color = cnv_state, alpha = p_cnv),
            size = line_width
        ) +
        geom_segment(
            inherit.aes = F,
            aes(x = seg_start, xend = seg_end, y = 1, yend = 1),
            data = segs_consensus, size = 0, color = 'white', alpha = 0
        ) +
        geom_hline(yintercept = c(first_tumor_index, clone_indices), color = 'royalblue', size = 0.5, linetype = 'dashed') +
        # geom_hline(yintercept = c(first_tumor_index, clone_indices), color = 'gray', size = 0.5, linetype = 'solid') +
        theme(
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
            strip.background = element_blank(),
            axis.text = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_text(size = 10),
            axis.ticks = element_blank(),
            plot.margin = margin(0,0,5,0, unit = 'mm'),
            axis.line = element_blank(),
            legend.box.background = element_blank(),
            legend.background = element_blank(),
            # panel.background = element_rect(fill = "transparent",colour = NA),
            # plot.background = element_rect(fill = "transparent", color = NA)
        ) +
        scale_x_continuous(expand = expansion(0)) +
        scale_y_continuous(expand = expansion(0)) +
        facet_grid(.~CHROM, space = 'free', scale = 'free', labeller = labeller(CHROM = chrom_labeller)) +
        scale_alpha_continuous(range = c(0,1), limits = c(p_min, 1), oob = scales::squish) +
        guides(
            alpha = 'none',
            # alpha = guide_legend(),
            color = guide_legend(override.aes = c(size = 1), title = 'CNV state')
        ) +
        scale_color_manual(
            values = c('amp' = 'darkred', 'del' = 'darkblue', 'bamp' = cnv_colors[['bamp']], 'loh' = 'darkgreen', 'bdel' = 'blue'),
            labels = c('amp' = 'AMP', 'del' = 'DEL', 'bamp' = 'BAMP', 'loh' = 'CNLoH', 'bdel' = 'BDEL'),
            limits = force,
            na.translate = F
        ) +
        xlab('Genomic position')

    # excluded regions
    if (!is.null(segs_loh)) {

        segs_loh = segs_loh %>% mutate(CHROM = as.integer(as.character(CHROM)))

        p_segs = p_segs + 
            geom_rect(
                inherit.aes = F,
                data = segs_loh, 
                aes(xmin = seg_start,
                    xmax = seg_end,
                    ymin = min(joint_post$cell_index),
                    ymax = max(joint_post$cell_index)
                ),
                fill = 'darkblue'
            )

    }

    # clone annotation
    if (is.null(clone_dict)) {
        clone_dict = gtree %>%
            activate(nodes) %>%
            data.frame %>%
            mutate(
                GT = ifelse(compartment == 'normal', '', GT),
                GT = factor(GT),
                clone = ifelse(compartment == 'normal', 1, clone),
                clone = as.factor(clone)
            ) %>%
            {setNames(.$clone, .$name)}
    }

    if (is.null(pal_clone)) {

        getPalette = colorRampPalette(RColorBrewer::brewer.pal(n = 10, 'Spectral'))

        pal_clone = getPalette(length(unique(clone_dict)))

        if ('1' %in% unique(clone_dict)) {
            pal_clone = c('gray', pal_clone)
        }
        
    }

    p_clone = data.frame(
            cell = names(clone_dict),
            annot = unname(clone_dict)
        ) %>%
        mutate(cell = factor(cell, cell_order)) %>%
        filter(!is.na(cell)) %>%
        annot_bar(transpose = TRUE, legend = clone_legend, pal_annot = pal_clone, legend_title = clone_title, size = size)
    
    # external annotation
    if (!is.null(annot)) {
        
        p_annot = data.frame(
                cell = names(annot),
                annot = unname(annot)
            ) %>%
            filter(cell %in% cell_order) %>%
            mutate(cell = factor(cell, cell_order)) %>%
            annot_bar(transpose = T, pal_annot = pal_annot, legend_title = annot_title, size = size, annot_scale = annot_scale)

        if (clone_bar) {
            (p_tree | p_clone | p_annot | p_segs) + plot_layout(widths = c(tree_height, 0.25, 0.25, 15), guides = 'collect')
        } else {
            (p_tree | p_annot | p_segs) + plot_layout(widths = c(tree_height, 0.25, 15), guides = 'collect')
        }
    } else if (clone_bar) {
        (p_tree | p_clone | p_segs) + plot_layout(widths = c(tree_height, 0.25, 15), guides = 'collect')
    } else {
        (p_tree | p_segs) + plot_layout(widths = c(tree_height, 15), guides = 'collect')
    }
}


#' plot consensus CNVs
#' @param segs dataframe Consensus segments
#' @return ggplot object
#' @export
plot_consensus = function(segs) {
  
    chrom_labeller <- function(chr){
        chr[chr %in% c(19, 21, 22)] = ''
        return(chr)
    }
    
    ggplot(
        segs
    ) +
    geom_rect(
        aes(xmin = seg_start, xmax = seg_end, ymin = -0.5, ymax = 0.5, fill = cnv_state_post)
    ) +
    theme_void() +
    theme(
        panel.spacing = unit(1, 'mm'),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = 'top'
    ) +
    facet_grid(~CHROM, space = 'free_x', scale = 'free', labeller = labeller(CHROM = chrom_labeller)) +
    scale_fill_manual(
      values = cnv_colors,
      labels = cnv_labels,
      name = 'CN states'
    ) +
    ggrepel::geom_text_repel(
        aes(x = (seg_start+seg_end)/2, y = -0.5,
            label = str_remove(seg_cons, '\\d+')),
        min.segment.length = 0,
        vjust = 1, 
        hjust = 0,
        direction = 'x',
        segment.curvature = -0.2,
        segment.ncp = 3,
        segment.angle = 30,
        segment.inflect = TRUE,
        max.overlaps = 3
    ) +
    scale_y_continuous(
        expand = expansion(add = c(0.5, 0))
    ) +
    scale_x_continuous(
        expand = expansion(mult = 0.05)
    ) +
    guides(fill = 'none')
}

#' Plot single-cell smoothed expression magnitude heatmap
#' @param gexp_roll_wide matrix Cell x gene smoothed expression magnitudes
#' @param hc hclust Hierarchical clustring result
#' @param gtf dataframe Transcript GTF
#' @param lim numeric Limit for expression magnitudes
#' @param k integer Number of clusters
#' @param n_sample integer Number of cells to subsample
#' @param reverse logical Whether to reverse the cell order
#' @param plot_tree logical Whether to plot the dendrogram
#' @export
plot_exp_roll = function(gexp_roll_wide, hc, k, gtf, lim = 0.8, n_sample = 300, reverse = TRUE, plot_tree = TRUE) {

    gexp_norm_long = gexp_roll_wide %>% 
        as.data.frame() %>%
        tibble::rownames_to_column('cell') %>%
        reshape2::melt(id.var = 'cell', variable.name = 'gene', value.name = 'exp_rollmean') %>%
        left_join(gtf, by = 'gene') %>%
        mutate(gene_index = as.integer(factor(gene, unique(gene))))

    cells = unique(gexp_norm_long$cell)
    
    cell_sample = sample(cells, min(n_sample, length(cells)), replace = FALSE)
        
    p_tree = ggtree(hc, size = 0.2)

    cell_order = p_tree$data %>% filter(isTip) %>% arrange(y) %>% pull(label)

    chrom_labeller <- function(chr){
        chr[chr %in% c(18, 21)] = ''
        return(chr)
    }
    
    p_heatmap = gexp_norm_long %>%
        filter(cell %in% cell_sample) %>%
        mutate(cell = factor(cell, cell_order)) %>%
        mutate(cluster = cutree(hc, k = k)[as.character(cell)]) %>%
        arrange(cell) %>%
        mutate(cluster = factor(cluster, rev(unique(cluster)))) %>%
        mutate(cell_index = as.integer(cell)) %>%
        ggplot() +
        geom_tile(aes(x = gene_index, y = cell, fill = exp_rollmean)) +
        # geom_rect(aes(xmin = gene_start, xmax = gene_end, ymin = cell_index - 0.5, ymax = cell_index + 0.5, fill = exp_rollmean)) +
        scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = 0, limits = c(-lim,lim), oob = scales::squish) +
        theme_void() +
        scale_x_discrete(expand = expansion(0)) +
        scale_y_discrete(expand = expansion(0)) +
        theme(
            axis.text.y = element_blank(),
            legend.position = 'top',
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'wheat4', fill = NA),
            strip.text.y = element_blank(),
            # axis.title.y = element_text(angle = 90),
            # axis.title.x = element_text()
        ) +
        facet_grid(cluster~CHROM, scale = 'free', space = 'free', labeller = labeller(CHROM = chrom_labeller)) +
        guides(fill = guide_colorbar(title = 'Expression\nmagnitude')) +
        xlab('Gene index') +
        ylab('Cell')
    
    if (plot_tree) {
        (p_tree | p_heatmap) + plot_layout(widths = c(1,10))
    } else {
        p_heatmap
    }
}

#' Plot single-cell smoothed expression magnitude heatmap
#' @param gtree tbl_graph The single-cell phylogeny
#' @param label_mut logical Whether to label mutations
#' @param label_size numeric Size of mutation labels
#' @param dot_size numeric Size of mutation nodes
#' @param branch_width numeric Width of branches in tree
#' @param tip logical Whether to plot tip point
#' @param tip_length numeric Length of the tips
#' @param pal_clone named vector Clone colors
#' @export
plot_sc_tree = function(gtree, label_mut = TRUE, label_size = 3, dot_size = 2, branch_width = 0.5, tip = TRUE, tip_length = 0.5, pal_clone = NULL) {

    mut_nodes = gtree %>% activate(nodes) %>%
      filter(!is.na(site)) %>% data.frame() %>%
      select(name, site)
    
    gtree = gtree %>% activate(edges) %>% mutate(length = ifelse(leaf, pmax(length, tip_length), length))
    
    clone_dict = gtree %>%
        activate(nodes) %>%
        data.frame %>%
        mutate(
            GT = ifelse(compartment == 'normal', '', GT),
            GT = factor(GT),
            clone = as.factor(clone)
        ) %>%
      {setNames(.$clone, .$name)}
    
    OTU_dict = lapply(levels(clone_dict), function(x) names(clone_dict[clone_dict == x])) %>% setNames(levels(clone_dict))

    p_tree = gtree %>% 
        to_phylo() %>%
        groupOTU(
            OTU_dict,
            'clone'
        ) %>%
        ggtree(ladderize = T, size = branch_width) %<+%
        mut_nodes +
        layout_dendrogram() +
        geom_rootedge(size = branch_width) +
        theme(
            plot.margin = margin(0,0,0,0),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line.y = element_line(size = 0.2),
            axis.ticks.y = element_line(size = 0.2),
            axis.text.y = element_text(size = 8)
        ) +
        guides(color = 'none') +
        xlab('Number of mutations')

    if (label_mut) {
        p_tree = p_tree + 
            geom_point2(aes(subset = !is.na(site), x = branch), shape = 21, size = dot_size, fill = 'red') +
            geom_text2(
                aes(x = branch, label = str_trunc(site, 20, side = 'center')),
                size = label_size, hjust = 0, vjust = -0.5, nudge_y = 1, color = 'darkred'
            )
    }

    if (tip) {

        if (is.null(pal_clone)) {
            getPalette = colorRampPalette(pal)
            pal_clone = getPalette(nrow(mut_nodes) + 1)
        }

        p_tree = p_tree + 
            geom_tippoint(aes(color = as.factor(clone)), size=dot_size, stroke = 0.2) +
            scale_color_manual(values = pal_clone, limits = force)
    }
    
    return(p_tree)
    
}


########################### Functions for internal use ############################

plot_sc_exp = function(exp_post, segs_consensus, size = 0.05, censor = 0) {
    
    # cell_order = exp_post %>% 
    #     filter(!cnv_state %in% c('neu', 'loh')) %>%
    #     left_join(cell_annot, by = 'cell') %>%
    #     group_by(cell_group) %>%
    #     do(
    #         reshape2::dcast(., cell ~ seg, value.var = 'phi_mle') %>%
    #         tibble::column_to_rownames('cell') %>%
    #         dist() %>%
    #         hclust %>%
    #         {.$labels[.$order]} %>%
    #         as.data.frame()
    #     ) %>%
    #     set_names(c('cell_group', 'cell'))

    exp_post = exp_post %>% filter(n > 15)

    exp_post = exp_post %>% 
            inner_join(
                segs_consensus %>% select(seg = seg_cons, CHROM, seg_start, seg_end),
                by = 'seg'
            ) %>%
            mutate(phi_mle = ifelse(phi_mle > 1-censor & phi_mle < 1+censor, 1, phi_mle))
            # mutate(cell = factor(cell, cell_order$cell))

    ggplot(
        exp_post,
        aes(x = seg_start, xend = seg_end, y = cell, yend = cell, color = phi_mle)
    ) +
    theme_classic() +
    geom_segment(size = size) +
    theme(
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
        strip.background = element_blank(),
        axis.text.y = element_blank()
    ) +
    scale_x_continuous(expand = expansion(0)) +
    facet_grid(group~CHROM, space = 'free', scale = 'free') +
    scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 1, limits = c(0.5, 2), oob = scales::oob_squish)
}

plot_sc_allele = function(df_allele, bulk_subtrees, clone_post) {
    
    snp_seg = bulk_subtrees %>%
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

    snp_neu_haplo = bulk_subtrees %>% filter(sample == 1) %>%
        mutate(haplo = ifelse(pBAF > 0.5, 'major', 'minor')) %>% 
        filter(!is.na(haplo)) %>%
        {setNames(.$haplo, .$snp_id)}
    
    
    pal = RColorBrewer::brewer.pal(n = 8, 'Set1')

    p = df_allele %>% 
        left_join(
            snp_seg %>% select(snp_id, haplo, cnv_state),
            by = 'snp_id'
        ) %>% 
        mutate(haplo = ifelse(is.na(haplo), snp_neu_haplo[snp_id], haplo)) %>%
        mutate(cnv_state = ifelse(is.na(cnv_state), 'neu', cnv_state)) %>%
        mutate(pBAF = ifelse(GT == '1|0', AR, 1-AR)) %>%
        filter(!is.na(haplo)) %>%
        mutate(MAF = ifelse(haplo == 'major', pBAF, 1-pBAF)) %>%
        left_join(clone_post, by = 'cell') %>%
        arrange(clone) %>%
        arrange(CHROM, POS) %>%
        mutate(snp_index = as.integer(factor(snp_id, unique(snp_id)))) %>%
        ggplot(
            aes(x = snp_index, y = cell, color = MAF)
        ) +
        theme_classic() +
        geom_point(alpha = 0.5, pch = 16, size = 1) +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'white', fill = NA),
        ) +
        facet_grid(clone~CHROM, space = 'free', scale = 'free') +
        scale_x_discrete(expand = expansion(0)) +
        scale_color_gradient(low = pal[1], high = pal[2])
    
    return(p)
}

plot_markers = function(sample, count_mat, cell_annot, markers, clone_post, pal_annot = NULL) {

    if (is.null(pal_annot)) {
        pal_annot = getPalette(length(unique(cell_annot$annot)))
    }
    
    D = as.matrix(count_mat[,markers$gene]) %>%
        scale %>%
        reshape2::melt() %>%
        magrittr::set_colnames(c('cell', 'gene', 'exp')) %>%
        inner_join(
            cell_annot, by = 'cell'
        ) %>%
        mutate(exp = ifelse(is.na(exp), 0, exp)) %>%
        inner_join(
            clone_post, by = 'cell'
        ) %>%
        left_join(markers, by = 'gene') %>%
        arrange(p_1) %>%
        mutate(cell = factor(cell, unique(cell)))

    p_markers = ggplot(
            D,
            aes(x = cell, y = gene, fill = exp)
        ) +
        geom_tile() +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 7),
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.2, fill = NA),
            strip.background.x = element_blank(),
            strip.text.x = element_blank(),
            strip.background.y = element_rect(size = 0, fill = NA),
            strip.text.y.left = element_text(size = 6, angle = 0),
        ) +
        ylab('marker') +
        facet_grid(marker_type ~ cell_group, space = 'free_y', scale = 'free', switch="y") +
        scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', limits = c(-1.5,1.5), oob = scales::oob_squish)

    p_annot = ggplot(
            D,
            aes(x = cell, y = 'annot', fill = annot)
        ) +
        geom_tile() +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 7),
            axis.ticks.x = element_blank(),
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.2, fill = NA),
            strip.background = element_rect(size = 0, fill = NA),
            strip.text = element_text(size = 6),
            axis.title.x = element_blank(),
            strip.text.x = element_blank()
        ) +
        ylab('') +
        facet_grid(~ cell_group, space = 'free_y', scale = 'free', switch="y") +
        scale_fill_manual(values = pal_annot) +
        guides(fill = guide_legend())

    p_cnv = ggplot(
            D %>% mutate(p_cnv = 1-p_1),
            aes(x = cell, y = 'cnv', fill = p_cnv)
        ) +
        geom_tile() +
        theme_classic() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 7),
            axis.ticks.x = element_blank(),
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.2, fill = NA),
            strip.background = element_rect(size = 0, fill = NA),
            strip.text = element_text(size = 6),
            axis.title.x = element_blank()
        ) +
        ylab('') +
        facet_grid(~cell_group, space = 'free_y', scale = 'free', switch="y") +
        scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0.5, limits = c(0,1), oob = scales::oob_squish) + 
        ggtitle(sample) +
        guides(fill = 'none')

    p_cnv/p_annot/p_markers + plot_layout(heights = c(0.5,0.5,10), guides = 'collect')
    
}

do_plot = function(p, f, w, h, out_dir = '~/figures') {
    ggsave(filename = paste0(out_dir, '/', f, '.png'), plot = p, width = w, height = h, device = 'png', dpi = 300)
    options(repr.plot.width = w, repr.plot.height = h, repr.plot.res = 300)
    print(p)
}


annot_bar = function(D, transpose = FALSE, legend = TRUE, legend_title = '', size = 0.05, pal_annot = NULL, annot_scale = NULL) {

    D = D %>% mutate(cell_index = as.integer(cell))

    index_max = length(levels(D$cell))

    p = ggplot(
        D,
        aes(x = cell_index, y = '', fill = annot)
    ) +
    geom_tile(width=1, height=0.9, size = 0) +
    # geom_segment(
    #     aes(x = cell, xend = cell, y = -0.5, yend = 0.5, color = annot),
    #     size = size
    # ) +
    theme_void() +
    scale_y_discrete(expand = expansion(0)) +
    scale_x_continuous(expand = expansion(0), limits = c(1,index_max)) +
    theme(
        panel.spacing = unit(0.1, 'mm'),
        panel.border = element_rect(size = 0, color = 'black', fill = NA),
        panel.background = element_rect(fill = 'gray90'),
        strip.background = element_blank(),
        strip.text = element_blank(),
        # axis.text = element_text(size = 8),
        axis.text = element_blank(),
        plot.margin = margin(0.5,0,0.5,0, unit = 'mm')
    )

    if (!is.null(annot_scale)) {
        p = p + annot_scale
    } else {
        if (is.null(pal_annot)) {
            pal_annot = getPalette(length(unique(D$annot)))
        }
        p = p + scale_fill_manual(values = pal_annot, na.value = 'gray90')
    }

    if (transpose) {
        p = p + coord_flip() +
            theme(plot.margin = margin(0,0.5,0,0.5, unit = 'mm'))
    }

    if (legend) {
        p = p + guides(fill = guide_legend(keywidth = unit(3, 'mm'), keyheight = unit(1, 'mm'), title = legend_title))
    } else {
        p = p + guides(fill = 'none')
    }

    return(p)
}

#' @keywords internal
show_phasing = function(bulk, min_depth = 8, dot_size = 0.5, h = 50) {

    D = bulk %>% 
        filter(!is.na(pAD)) %>%
        group_by(CHROM) %>%
        mutate(pBAF = 1-pBAF, pAD = DP - pAD) %>%
        mutate(theta_ar_roll = theta_hat_roll(AD, DP-AD, h = h)) %>%
        mutate(theta_hat_roll = theta_hat_roll(pAD, DP-pAD, h = h)) %>%
        filter(DP >= min_depth) %>%
        mutate(snp_index = 1:n()) %>%
        mutate(dAR = 0.5+abs(AR-0.5)) %>%
        mutate(dAR_roll = caTools::runmean(abs(AR-0.5), align = 'center', k = 30))

    boundary = D %>% filter(boundary == 1) %>% pull(snp_index)
    
    p1 = D %>%
        mutate(state_post = 'neu') %>%
        ggplot(
            aes(x = snp_index, y = AR),
            na.rm=TRUE
        ) +
        geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'gray') +
        geom_point(
            aes(color = state_post),
            size = dot_size,
        ) +
        # geom_line(
        #     aes(x = snp_index, y = 0.5 + dAR_roll), color = 'red'
        # ) +
        # geom_line(
        #     aes(x = snp_index, y = 0.5 - dAR_roll), color = 'red'
        # ) +
        # geom_line(
        #     aes(x = snp_index, y = 0.5 + theta_ar_roll), color = 'red'
        # ) +
        theme_classic() +
        theme(
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
            strip.background = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank()
        ) +
        scale_color_manual(values = cnv_colors, limits = force) +
        ylim(0,1) +
        facet_grid(.~CHROM, space = 'free_x', scale = 'free_x') +
        geom_vline(xintercept = boundary - 1, color = 'red', size = 0.5, linetype = 'dashed') +
        guides(color = 'none')

    p2 = D %>%
        mutate(state_post = 'neu') %>%
        ggplot(
            aes(x = snp_index, y = pBAF),
            na.rm=TRUE
        ) +
        geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'gray') +
        geom_point(
            aes(color = ifelse(theta_hat_roll > 0, 'loh_1_up', 'loh_1_down')),
            size = dot_size,
        ) +
        geom_line(
            aes(x = snp_index, y = 0.5 + theta_hat_roll), color = 'red'
        ) +
        theme_classic() +
        theme(
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
            strip.background = element_blank()
        ) +
        scale_color_manual(values = cnv_colors) +
        ylim(0,1) +
        facet_grid(.~CHROM, space = 'free_x', scale = 'free_x') +
        geom_vline(xintercept = boundary - 1, color = 'red', size = 0.5, linetype = 'dashed') +
        guides(color = 'none')

    (p1 / p2) + plot_layout(guides = 'auto')
}

plot_exp = function(gexp_bulk, exp_limit = 3) {
    ggplot(
        gexp_bulk,
        aes(x = gene_index, y = logFC)
    ) +
    theme_classic() +
    geom_point(size = 1, color = 'gray', alpha = 0.8, pch = 16) +
    theme(
        panel.spacing = unit(0, 'mm'),
        panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
        strip.background = element_blank()
    ) +
    facet_grid(~CHROM, space = 'free_x', scale = 'free_x') +
    geom_hline(yintercept = 0, color = 'gray30', linetype = 'dashed') +
    geom_line(
        inherit.aes = FALSE,
        aes(x = gene_index, y = log2(phi_hat_roll), group = '1'),
        color = 'darkred',
        size = 0.35
    ) +
    ylim(-exp_limit, exp_limit)
}

plot_segs_post = function(segs_consensus) {
    segs_consensus %>% 
        filter(cnv_state != 'neu') %>%
        mutate(seg_label = paste0(seg_cons, '_', cnv_state_post)) %>%
        mutate(seg_label = factor(seg_label, unique(seg_label))) %>%
        reshape2::melt(measure.vars = c('p_loh', 'p_amp', 'p_del', 'p_bamp', 'p_bdel'), value.name = 'p') %>%
        ggplot(
            aes(x = seg_label, y = variable, fill = p, label = round(p, 2))
        ) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
        geom_tile() +
        geom_text(color = 'white')
}

# model diagnostics
plot_exp_post = function(exp_post, jitter = TRUE) {
    if (!'annot' %in% colnames(exp_post)) {
        exp_post$annot = '0'
    }
    p = exp_post %>%
        filter(n > 20) %>%
        mutate(seg_label = paste0(seg, '(', cnv_state, ')')) %>%
        mutate(seg_label = factor(seg_label, gtools::mixedsort(unique(seg_label)))) %>%
        ggplot(
            aes(x = seg_label, y = log2(phi_mle), fill = cnv_state, color = p_cnv)
        ) +
        geom_violin(size = 0) +
        geom_hline(yintercept = 0, color = 'green', linetype = 'dashed') +
        geom_hline(yintercept = log2(1.5), color = 'red', linetype = 'dashed') +
        geom_hline(yintercept = -1, color = 'blue', linetype = 'dashed') +
        facet_grid(annot~cnv_state, scale = 'free', space = 'free') +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
        scale_fill_manual(values = cnv_colors)

    if (jitter) {
        p = p + geom_jitter(size = 0.1)
    }

    return(p)
}

plot_clones = function(p_matrix, gtree, annot = TRUE, n_sample = 1e4, bar_ratio = 0.1, pal_clone = NULL, pal_annot = NULL) {

    if (is.null(pal_clone)) {
        pal_clone = c('gray', RColorBrewer::brewer.pal(n = 8, 'Set1'))
    }

    if (is.null(pal_annot)) {
        pal_annot = c('gray', RColorBrewer::brewer.pal(n = 8, 'Set1'))
    }

    p_matrix = p_matrix %>% 
        group_by(group) %>%
        filter(cell %in% sample(unique(cell), min(n_sample, length(unique(cell))))) %>%
        ungroup() %>%
        filter(cnv_state != 'neu')

    # ordering cells
    if (annot) {
        p_matrix = p_matrix %>% 
            arrange(group, annot) %>%
            mutate(cell = factor(cell, unique(cell)))
    } else {
        set.seed(0)
        p_matrix = p_matrix %>% 
            group_by(group) %>%
            sample_frac(1) %>%
            ungroup() %>%
            mutate(cell = factor(cell, unique(cell)))
    }
    
    # ordering cnvs
    cnv_order = gtree %>% 
            activate(nodes) %>%
            mutate(rank = dfs_rank(root = node_is_root())) %>%
            data.frame() %>%
            filter(!is.na(site)) %>%
            arrange(-rank) %>%
            pull(site) %>%
            map(function(x){rev(unlist(str_split(x, ',')))}) %>%
            unlist

    p_matrix = p_matrix %>%
        mutate(seg = factor(seg, cnv_order)) %>%
        arrange(seg) %>%
        mutate(seg_label = factor(seg_label, unique(seg_label)))  %>%
        mutate(group = factor(group))
    
    p = ggplot(
            p_matrix,
            aes(x = cell, y = seg_label, fill = logBF)
        ) +
        geom_tile(width=0.1, height=0.9) +
        theme_classic() +
        scale_y_discrete(expand = expansion(0)) +
        scale_x_discrete(expand = expansion(0)) +
        theme(
            panel.spacing = unit(0.1, 'mm'),
            panel.border = element_rect(size = 0.2, color = 'black', fill = NA),
            panel.background = element_rect(fill = 'white'),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            plot.margin = margin(3,0,0,0, unit = 'pt')
            # strip.text = element_text(angle = 0, size = 8, vjust = 0.5)
        ) +
        facet_grid(.~group, scale = 'free', space = 'free') +
        scale_fill_gradient2(low = pal[2], high = pal[1], midpoint = 0, limits = c(-5, 5), oob = scales::oob_squish) +
        xlab('') +
        ylab('') +
        guides(fill = guide_colorbar(barwidth = unit(3, 'mm'), barheight = unit(15, 'mm')))

    p_clones = ggplot(
            p_matrix %>% distinct(cell, group),
            aes(x = cell, y = 'clone', fill = group)
        ) +
        geom_tile(width=1, height=0.9) +
        theme_void() +
        scale_y_discrete(expand = expansion(0)) +
        scale_x_discrete(expand = expansion(0)) +
        theme(
            panel.spacing = unit(0.1, 'mm'),
            panel.border = element_rect(size = 0, color = 'black', fill = NA),
            panel.background = element_rect(fill = 'white'),
            strip.background = element_blank(),
            strip.text = element_text(angle = 0, size = 10, vjust = 0.5),
            axis.text.y = element_text(size = 8)
        ) +
        facet_grid(.~group, scale = 'free', space = 'free') +
        xlab('') +
        ylab('') + 
        scale_fill_manual(values = pal_clone) +
        guides(fill = 'none')

    if (annot) {

        p_annot = ggplot(
                p_matrix %>% distinct(cell, group, annot),
                aes(x = cell, y = '', fill = annot)
            ) +
            geom_tile(width=1, height=0.9) +
            theme_void() +
            scale_y_discrete(expand = expansion(0)) +
            scale_x_discrete(expand = expansion(0)) +
            theme(
                panel.spacing = unit(0.1, 'mm'),
                panel.border = element_rect(size = 0, color = 'black', fill = NA),
                panel.background = element_rect(fill = 'white'),
                strip.background = element_blank(),
                strip.text = element_blank(),
                axis.text.y = element_text(size = 8),
                plot.margin = margin(3,0,0,0, unit = 'pt')
            ) +
            facet_grid(.~group, scale = 'free', space = 'free') +
            xlab('') +
            ylab('') +
            scale_fill_manual(values = pal_annot) +
            guides(fill = guide_legend(keywidth = unit(2, 'mm'), keyheight = unit(2, 'mm'), title = ''))

        return((p_clones / p_annot / p) + plot_layout(height = c(bar_ratio, bar_ratio, 1), guides = 'collect'))
        
    } else {
        return((p_clones / p) + plot_layout(height = c(1,10)))
    }
    
}

plot_clone_panel = function(res, label = NULL, cell_annot = NULL, type = 'joint', ratio = 1, tvn = FALSE, tree = TRUE, p_min = 0.5, bar_ratio = 0.1, pal_clone = NULL, pal_annot = NULL) {

    if (is.null(pal_clone)) {
        n_clones = length(unique(res$clone_post$clone_opt))
        pal_clone = getPalette(max(V(res$G_m)-1, 8)) %>% c('gray', .) %>% setNames(1:n_clones)
    } 
    
    if (is.null(pal_annot) & !is.null(cell_annot)) {
        pal_annot = getPalette(length(unique(cell_annot$annot)))
    }
    

    if (type == 'joint') {
        p_matrix = res$joint_post
    } else if (type == 'allele') {
        p_matrix = res$allele_post
    } else {
        p_matrix = res$exp_post
    }

    if (!is.null(cell_annot)) {
        p_matrix = p_matrix %>% left_join(cell_annot, by = 'cell')
        annot = TRUE
    } else {
        annot = FALSE
    }

    if (!'p_opt' %in% colnames(res$clone_post)) {
        res$clone_post = res$clone_post %>% 
            rowwise() %>%
            mutate(p_opt = get(paste0('p_', clone_opt))) %>%
            ungroup()
    }

    if (tvn) {
        res$clone_post = res$clone_post %>%
            mutate(
                clone_opt = ifelse(clone_opt == 1, 'normal', 'tumor'),
                p_opt = ifelse(clone_opt == 'normal', p_1, 1-p_1)
            )
    }

    p_clones = p_matrix %>% 
        filter(seg %in% colnames(res$geno)) %>%
        inner_join(
            res$clone_post %>% filter(p_opt > p_min),
            by = 'cell'
        ) %>%
        mutate(group = clone_opt) %>%
        plot_clones(res$gtree, pal_clone = pal_clone, pal_annot = pal_annot, annot = annot, bar_ratio = bar_ratio)

    plot_title = plot_annotation(title = label, theme = theme(plot.title = element_text(hjust = 0.1)))

    if (tvn | (!tree)) {
        return(p_clones + plot_title)
    }

    p_mut = res$G_m %>% plot_mut_history(pal = pal_clone) 

    (p_mut / p_clones) + plot_layout(heights = c(ratio, 1)) + plot_title
}

cnv_heatmap = function(segs, var = 'group', label_group = TRUE) {
    
    gaps_hg38_filtered = gaps_hg38 %>% filter(end - start > 1e6)

    if (!'p_cnv' %in% colnames(segs)) {
        segs$p_cnv = 1
    }
    
    p = ggplot(
            segs
        ) +
        geom_rect(
            aes(xmin = seg_start, xmax = seg_end, ymin = -0.5, ymax = 0.5, fill = cnv_state, alpha = p_cnv)
        ) +
        geom_rect(
            data = chrom_sizes_hg38,
            aes(xmin = 0, xmax = size, ymin = -0.5, ymax = 0.5, fill = NA)
        ) +
        geom_rect(
            inherit.aes = F,
            data = gaps_hg38_filtered,
            aes(xmin = start, 
                xmax = end,
                ymin = -0.5,
                ymax = 0.5),
            fill = 'white'
        ) +
        geom_rect(
            inherit.aes = F,
            data = gaps_hg38_filtered,
            aes(xmin = start, 
                xmax = end,
                ymin = -0.5,
                ymax = 0.5),
            fill = 'gray',
            alpha = 0.5
        ) +
        theme_classic() +
        theme(
            panel.spacing = unit(0, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
            strip.background = element_blank(),
            strip.text.y = element_text(angle = 0),
            axis.text = element_blank(),
            plot.margin = margin(0, 0, 0, 0),
            axis.title.x = element_blank(),
            axis.ticks = element_blank(),
            axis.line.x = element_blank()
        ) +
        scale_fill_manual(
            values = c('neu' = 'white', cnv_colors[names(cnv_colors) != 'neu']),
            na.value = 'white',
            limits = force,
            labels = cnv_labels,
            na.translate = F,
            name = 'States'
        ) +
        scale_alpha_continuous(range = c(0,1), limits = c(0.5,1), oob = scales::squish) +
        guides(alpha = 'none') +
        scale_x_continuous(expand = expansion(add = 0)) +
        facet_grid(get(var)~CHROM, space = 'free_x', scale = 'free', drop = TRUE)

        if (!label_group) {
            p = p + theme(strip.text.y = element_blank())
        }

        return(p)

}
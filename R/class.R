
#' @title Numbat R6 class
#' @description 
#' @param out_dir the Numbat run output directory
#' @export
Numbat <- R6::R6Class("Numbat", lock_objects=FALSE,
  public = list(
    #' @field label
    label = 'sample',

    initialize = function(out_dir, verbose=TRUE) {

        self$out_dir = out_dir

        self$fetch_results(out_dir)

    },

    fetch_results = function(out_dir, i = 2) {

        self$joint_post = fread(glue('{out_dir}/joint_post_{i}.tsv'))
        self$exp_post = fread(glue('{out_dir}/exp_post_{i}.tsv'))
        self$allele_post = fread(glue('{out_dir}/allele_post_{i}.tsv'))
        self$bulk_init = fread(glue('{out_dir}/bulk_subtrees_0.tsv.gz'))
        self$bulk_subtrees = fread(glue('{out_dir}/bulk_subtrees_{i}.tsv.gz'))
        self$segs_consensus = fread(glue('{out_dir}/segs_consensus_{i-1}.tsv'))
        self$bulk_clones = fread(glue('{out_dir}/bulk_clones_{i}.tsv.gz'))
        self$geno = fread(glue('{out_dir}/geno_{i}.tsv')) %>% tibble::column_to_rownames('V1')
        self$tree_post = readRDS(glue('{out_dir}/tree_post_{i}.rds'))
        self$mut_graph = readRDS(glue('{out_dir}/mut_graph_{i}.rds'))
        self$gtree = readRDS(glue('{out_dir}/tree_final_{i}.rds'))
        self$clone_post = fread(glue('{out_dir}/clone_post_{i}.tsv'))

    },

    plot_heatmap = function(annot = NULL, multi_allelic = FALSE, clone_bar = FALSE, tip_length = 1, p_min = 0.5, branch_width = 0.2, size = 0.1) {
        plot_sc_joint(  
            self$gtree,
            self$joint_post,
            self$segs_consensus,
            tip_length = 1,
            branch_width = 0.2,
            size = 0.075,
            p_min = p_min,
            clone_bar = clone_bar,
            multi_allelic = multi_allelic,
            cell_dict = annot
        )
    },
    
    plot_mut_history = function(horizontal = TRUE, label = TRUE) {
        plot_mut_history(self$mut_graph, horizontal = horizontal, label = label)
    },

    plot_bulks = function(what = 'clones', min_depth = 8, phi_mle = TRUE, ncol = 1, legend = FALSE) {

        if (!what %in% c('clones', 'subtrees')) {
            stop('what = clones or subtrees')
        }

        if (what == 'clones') {
            bulks = self$bulk_clones
        } else {
            bulks = self$bulk_subtrees
        }

        plot_bulks(
            bulks,
            min_depth = min_depth,
            ncol = ncol,
            phi_mle = phi_mle,
            legend = legend
        )
    }
  )
)

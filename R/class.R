
#' @title numbat R6 class
#' @description used to summarize results
#' @param out_dir the numbat run output directory
#' @export
numbat <- R6::R6Class("numbat", lock_objects=FALSE,
  public = list(
    #' @field label
    label = 'sample',
    gtf = NULL,

    #' @description initialize numbat class
    #' @param out_dir output directory
    #' @param i get results from which iteration
    #' @param gtf transcript gtf dataframe
    #' @param verbose verbosity
    #' @return a new 'numbat' object
    initialize = function(out_dir, i = 2, gtf = gtf_hg38, verbose=TRUE) {

        self$out_dir = out_dir

        self$fetch_results(out_dir, i = i)

        self$gtf = gtf

    },

    #' @description fetch results from the numbat output directory
    #' @param out_dir output directory
    #' @param i get results from which iteration
    #' @return NULL
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
        self$gexp_roll_wide = fread(glue('{out_dir}/gexp_roll_wide.tsv.gz'))
        self$hc = readRDS(glue('{out_dir}/hc.rds'))

    },

    #' @description plot the single-cell CNV calls in a heatmap and the corresponding phylogeny
    #' @param annot named list of cell annotation
    #' @param geno_bar where to show genotype annotation
    #' @param pal_clone clone color palette
    #' @param p_min minimum probability threshold to filter the calls
    #' @param tip_length tree tip length
    #' @param branch_width tree branch width
    #' @param line_width heatmap line width
    #' @param multi_allelic whether to show different allelic states for same CNV
    #' @return a ggplot object
    plot_phylo_heatmap = function(annot = NULL, geno_bar = FALSE, pal_clone = NULL, multi_allelic = FALSE, p_min = 0.9, tip_length = 1, branch_width = 0.2, line_width = 0.1) {
        plot_phylo_heatmap(  
            self$gtree,
            self$joint_post,
            self$segs_consensus,
            tip_length = tip_length,
            branch_width = branch_width,
            size = line_width,
            p_min = p_min,
            clone_bar = geno_bar,
            multi_allelic = multi_allelic,
            cell_dict = annot,
            pal_clone = pal_clone
        )
    },

    #' @description plot window-smoothed expression profiles
    #' @param k number of clusters
    #' @param hc where to show clone annotation
    #' @param n_sample number of cells to subsample
    #' @return a ggplot object
    plot_exp_roll = function(k = 3, n_sample = 300) {
        
        plot_sc_roll(
            self$gexp_roll_wide,
            self$hc,
            self$gtf,
            k = k,
            n_sample = n_sample
        )
    },
    
    #' @description plot the mutation history of the tumor
    #' @param horizontal horizontal layout or vertical
    #' @param label whether to label the mutations on edges
    #' @param pal_clone clone color palette
    #' @return a ggplot object
    plot_mut_history = function(horizontal = TRUE, label = TRUE, pal_clone = NULL) {
        p = plot_mut_history(self$mut_graph, horizontal = horizontal, label = label, pal_clone)
        return(p)
    },

    #' @description plot the bulk CNV profiles
    #' @param what whether to visualize clones or subtrees
    #' @param min_depth minimum allele coverage to filter SNPs
    #' @param ncol number of columns in the plot panel
    #' @param legend whether to display CNV state legend
    #' @param phi_mle whether to plot expression fold change estimates
    #' @return a ggplot object
    plot_bulks = function(what = 'clones', min_depth = 8, phi_mle = TRUE, ncol = 1, legend = FALSE) {

        if (!what %in% c('clones', 'subtrees')) {
            stop("The parameter 'what' must match one of these values: 'clones' or 'subtrees'")
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

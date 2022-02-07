
#' @title Numbat R6 class
#' @description Used to allow users to plot results
#' @export
Numbat <- R6::R6Class("Numbat", lock_objects=FALSE,
  public = list(

    #' @field label something something
    label = 'sample',
    
    #' @field gtf something something
    gtf = NULL,

    #' @field joint_post something something
    joint_post = NULL,

    #' @field exp_post something something
    exp_post = NULL,

    #' @field allele_post something something
    allele_post = NULL,

    #' @field bulk_subtrees something something
    bulk_subtrees = NULL,

    #' @field bulk_clones something something
    bulk_clones = NULL,

    #' @field segs_consensus something something
    segs_consensus = NULL,

    #' @field tree_post something something
    tree_post = NULL,

    #' @field mut_graph something something
    mut_graph = NULL,

    #' @field gtree something something
    gtree = NULL,

    #' @field clone_post something something
    clone_post = NULL,

    #' @field gexp_roll_wide something something
    gexp_roll_wide = NULL,

    #' @field hc something something
    hc = NULL,   

    #' @description initialize Numbat class
    #' @param out_dir character string Output directory
    #' @param i integer Get results from which iteration (either 1 or 2) (default=2)
    #' @param gtf transcript gtf dataframe (default=gtf_hg38)
    #' @param verbose boolean Whether to output verbose results (default=TRUE)
    #' @return a new 'numbat' object
    initialize = function(out_dir, i = 2, gtf = gtf_hg38, verbose=TRUE) {

        self$out_dir = out_dir

        private$fetch_results(out_dir, i = i)

        self$gtf = gtf

    },

    #' @description Plot the single-cell CNV calls in a heatmap and the corresponding phylogeny
    #' @param annot named list of cell annotation (default=NULL)
    #' @param geno_bar where to show genotype annotation (default=NULL)
    #' @param pal_clone clone color palette (default=NULL)
    #' @param p_min minimum probability threshold to filter the calls (default=0.9)
    #' @param tip_length tree tip length (default=1)
    #' @param branch_width tree branch width (default=0.2)
    #' @param line_width heatmap line width (default=0.1)
    #' @param tree_height plotting height of the tree (default=1)
    #' @return a ggplot object
    plot_phylo_heatmap = function(annot = NULL, geno_bar = FALSE, pal_clone = NULL, p_min = 0.9, tip_length = 1, branch_width = 0.2, line_width = 0.1, tree_height = 1) {

        p = plot_phylo_heatmap(  
            self$gtree,
            self$joint_post,
            self$segs_consensus,
            tip_length = tip_length,
            branch_width = branch_width,
            line_width = line_width,
            p_min = p_min,
            geno_bar = geno_bar,
            annot = annot,
            pal_clone = pal_clone,
            tree_height = tree_height
        )
        return(p)
    },

    #' @description Plot window-smoothed expression profiles
    #' @param k integer Number of clusters (default=3)
    #' @param n_sample integer Number of cells to subsample (default=300)
    #' @return a ggplot object
    plot_exp_roll = function(k = 3, n_sample = 300) {
        
        p = plot_sc_roll(
            self$gexp_roll_wide,
            self$hc,
            self$gtf,
            k = k,
            n_sample = n_sample
        )
        return(p)
    },
    
    #' @description Plot the mutation history of the tumor
    #' @param horizontal horizontal layout or vertical
    #' @param label whether to label the mutations on edges
    #' @param pal_clone clone color palette
    #' @return a ggplot object
    plot_mut_history = function(horizontal = TRUE, label = TRUE, pal_clone = NULL) {
        p = plot_mut_history(self$mut_graph, horizontal = horizontal, label = label, pal_clone)
        return(p)
    },

    #' @description Plot the bulk CNV profiles
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
  ), 

    ## Private functions
    private = list(

        #  Fetch results from the numbat output directory. The function
        #  will read the files expected in the output directory in the format
        #  "{out_dir}/filename_{i}.tsv"
        #
        #  param: out_dir output directory
        #  param: i get results from which iteration
        #  return: NULL
        # 
        fetch_results = function(out_dir, i = 2) {

            self$joint_post = fread(glue('{out_dir}/joint_post_{i}.tsv'))
            self$exp_post = fread(glue('{out_dir}/exp_post_{i}.tsv'))
            self$allele_post = fread(glue('{out_dir}/allele_post_{i}.tsv'))
            self$bulk_subtrees = fread(glue('{out_dir}/bulk_subtrees_{i}.tsv.gz'))
            self$bulk_clones = fread(glue('{out_dir}/bulk_clones_{i}.tsv.gz'))
            self$segs_consensus = fread(glue('{out_dir}/segs_consensus_{i}.tsv'))
            self$tree_post = readRDS(glue('{out_dir}/tree_post_{i}.rds'))
            self$mut_graph = readRDS(glue('{out_dir}/mut_graph_{i}.rds'))
            self$gtree = readRDS(glue('{out_dir}/tree_final_{i}.rds'))
            self$clone_post = fread(glue('{out_dir}/clone_post_{i}.tsv'))
            self$gexp_roll_wide = fread(glue('{out_dir}/gexp_roll_wide.tsv.gz'))
            self$hc = readRDS(glue('{out_dir}/hc.rds'))

        }
    )
)

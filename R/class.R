
#' @title Numbat R6 class
#' @description Used to allow users to plot results
#' @return a new 'Numbat' object
#' @export
Numbat <- R6::R6Class("Numbat", lock_objects=FALSE,
  public = list(

    #' @field label character Sample name
    label = 'sample',
    
    #' @field gtf dataframe Transcript annotation
    gtf = NULL,

    #' @field joint_post dataframe Joint posterior
    joint_post = NULL,

    #' @field exp_post dataframe Expression posterior
    exp_post = NULL,

    #' @field allele_post dataframe Allele posetrior
    allele_post = NULL,

    #' @field bulk_subtrees dataframe Bulk profiles of lineage subtrees
    bulk_subtrees = NULL,

    #' @field bulk_clones dataframe Bulk profiles of clones
    bulk_clones = NULL,

    #' @field segs_consensus dataframe Consensus segments
    segs_consensus = NULL,

    #' @field tree_post list Tree posterior
    tree_post = NULL,

    #' @field mut_graph igraph Mutation history graph
    mut_graph = NULL,

    #' @field gtree tbl_graph Single-cell phylogeny
    gtree = NULL,

    #' @field clone_post dataframe Clone posteriors
    clone_post = NULL,

    #' @field gexp_roll_wide matrix Smoothed expression of single cells
    gexp_roll_wide = NULL,

    #' @field P matrix Genotype probability matrix
    P = NULL,

    #' @field treeML matrix Maximum likelihood tree as phylo object
    treeML = NULL,

    #' @field hc hclust Initial hierarchical clustering
    hc = NULL,

    #' @description initialize Numbat class
    #' @param out_dir character string Output directory
    #' @param i integer Get results from which iteration (default=2)
    #' @param gtf dataframe Transcript gtf (default=gtf_hg38)
    #' @param verbose logical Whether to output verbose results (default=TRUE)
    #' @return a new 'Numbat' object
    initialize = function(out_dir, i = 2, gtf = gtf_hg38, verbose=TRUE) {

        self$out_dir = out_dir

        private$fetch_results(out_dir, i = i)

        self$gtf = gtf

    },

    #' @description Plot the single-cell CNV calls in a heatmap and the corresponding phylogeny
    #' @param ... additional parameters passed to plot_phylo_heatmap()
    plot_phylo_heatmap = function(...) {

        p = plot_phylo_heatmap(  
            gtree = self$gtree,
            joint_post = self$joint_post,
            segs_consensus = self$segs_consensus,
            clone_post = self$clone_post,
            ...
        )
        return(p)
    },

    #' @description Plot window-smoothed expression profiles
    #' @param k integer Number of clusters
    #' @param n_sample integer Number of cells to subsample
    #' @param ... additional parameters passed to plot_exp_roll()
    plot_exp_roll = function(k = 3, n_sample = 300, ...) {

        if (is.null(self$gexp_roll_wide)) {

            self$gexp_roll_wide = read_file(inputfile=glue('{self$out_dir}/gexp_roll_wide.tsv.gz'), filetype="tsv")
            
            if (!is.null(self$gexp_roll_wide)) {
                if ('V1' %in% colnames(self$gexp_roll_wide)) {
                    self$gexp_roll_wide = self$gexp_roll_wide %>% rename(cell = V1)
                }
                self$gexp_roll_wide = self$gexp_roll_wide %>% tibble::column_to_rownames('cell')
            }
        }
        
        p = plot_exp_roll(
            gexp_roll_wide = self$gexp_roll_wide,
            hc = self$hc,
            gtf = self$gtf,
            k = k,
            n_sample = n_sample,
            ...
        )
        return(p)
    },
    
    #' @description Plot the mutation history of the tumor
    #' @param ... additional parameters passed to plot_mut_history()
    plot_mut_history = function(...) {
        p = plot_mut_history(
                G = self$mut_graph,
                clone_post = self$clone_post,
                ...
            )
        return(p)
    },

    #' @description Plot the single cell phylogeny
    #' @param ... additional parameters passed to plot_sc_tree()
    plot_sc_tree = function(...) {
        p = plot_sc_tree(
                gtree = self$gtree,
                ...
            )
        return(p)
    },

    #' @description Plot consensus segments
    #' @param ... additional parameters passed to plot_sc_tree()
    plot_consensus = function(...) {
        p = plot_consensus(
                segs = self$segs_consensus,
                ...
            )
        return(p)
    },

    #' @description Plot clone cnv profiles
    #' @param ... additional parameters passed to plot_clone_profile()
    plot_clone_profile = function(...) {
        p = plot_clone_profile(
                joint_post = self$joint_post,
                clone_post = self$clone_post,
                ...
            )
        return(p)
    },

    #' @description Re-define subclones on the phylogeny. 
    #' @param max_cost numeric Likelihood threshold to collapse internal branches
    #' @param n_cut integer Number of cuts on the phylogeny to define subclones
    cutree = function(max_cost = 0, n_cut = 0) {
        
        self$gtree = get_gtree(self$treeML, self$P, max_cost = max_cost, n_cut = n_cut)
        self$mut_graph = scistreer::get_mut_graph(self$gtree)
        self$clone_post = get_clone_post(self$gtree, self$exp_post, self$allele_post)

    }),

    ## Private functions
    private = list(

        #  Fetch results from the numbat output directory. The function
        #  will read the files expected in the output directory in the format
        #  "{out_dir}/filename_{i}.tsv"
        #
        #  param: out_dir output directory
        #  param: i get results from which iteration
        #  return: NULL

        fetch_results = function(out_dir, i = 2) {
 
            self$joint_post = read_file(inputfile=glue('{out_dir}/joint_post_{i}.tsv'), filetype="tsv")
            self$exp_post = read_file(inputfile=glue('{out_dir}/exp_post_{i}.tsv'), filetype="tsv")
            self$allele_post = read_file(inputfile=glue('{out_dir}/allele_post_{i}.tsv'), filetype="tsv")
            self$bulk_clones = read_file(inputfile=glue('{out_dir}/bulk_clones_final.tsv.gz'), filetype="tsv")
            self$segs_consensus = read_file(inputfile=glue('{out_dir}/segs_consensus_{i}.tsv'), filetype="tsv")

            self$segs_consensus = self$segs_consensus %>% relevel_chrom()
            self$joint_post = self$joint_post %>% relevel_chrom()
            self$exp_post = self$exp_post %>% relevel_chrom()
            self$allele_post = self$allele_post %>% relevel_chrom()
            self$bulk_clones = self$bulk_clones %>% relevel_chrom()
            
            self$P = read_file(inputfile=glue('{out_dir}/geno_{i}.tsv'), header = TRUE) %>% tibble::column_to_rownames('cell') %>% as.matrix
            self$treeML = read_file(inputfile=glue('{out_dir}/treeML_{i}.rds'), filetype="rds")
            self$mut_graph = read_file(inputfile=glue('{out_dir}/mut_graph_{i}.rds'), filetype="rds")
            self$gtree = read_file(inputfile=glue('{out_dir}/tree_final_{i}.rds'), filetype="rds")
            self$clone_post = read_file(inputfile=glue('{out_dir}/clone_post_{i}.tsv'), filetype="tsv")
            self$hc = read_file(inputfile=glue('{out_dir}/hc.rds'), filetype="rds")

    })
)
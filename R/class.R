
#' @title Numbat R6 class
#' @description Used to allow users to plot results
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

    #' @field hc hclust Initial hierarchical clustering
    hc = NULL,

    #' @description initialize Numbat class
    #' @param out_dir character string Output directory
    #' @param i integer Get results from which iteration (either 1 or 2) (default=2)
    #' @param gtf dataframe Transcript gtf (default=gtf_hg38)
    #' @param verbose logical Whether to output verbose results (default=TRUE)
    #' @return a new 'Numbat' object
    initialize = function(out_dir, i = 2, gtf = gtf_hg38, verbose=TRUE, old_index = FALSE) {

        self$out_dir = out_dir

        private$fetch_results(out_dir, i = i, old_index = old_index)

        self$gtf = gtf

    },

    #' @description Plot the single-cell CNV calls in a heatmap and the corresponding phylogeny
    #' @param ... additional parameters passed to plot_phylo_heatmap()
    plot_phylo_heatmap = function(...) {

        p = plot_phylo_heatmap(  
            gtree = self$gtree,
            joint_post = self$joint_post,
            segs_consensus = self$segs_consensus,
            ...
        )
        return(p)
    },

    #' @description Plot window-smoothed expression profiles
    #' @param k integer Number of clusters
    #' @param n_sample integer Number of cells to subsample
    #' @param ... additional parameters passed to plot_exp_roll()
    plot_exp_roll = function(k = 3, n_sample = 300, ...) {
        
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

        fetch_results = function(out_dir, i = 2, old_index = FALSE) {
 
            # joint_post_colnames = c("cell", "CHROM",  "seg",  "cnv_state")
            self$joint_post = read_file(inputfile=glue('{out_dir}/joint_post_{i}.tsv'), filetype="tsv")
            # exp_post_colnames = c("seg", "cnv_state", "n", "phi_mle")
            self$exp_post = read_file(inputfile=glue('{out_dir}/exp_post_{i}.tsv'), filetype="tsv")
            # allele_post_colnames = c("cell", "CHROM", "seg", "cnv_state", "major", "minor", "total", "MAF", "seg_start", "seg_end", "prior_loh", "prior_amp", "prior_del", "prior_bamp", "prior_bdel")
            self$allele_post = read_file(inputfile=glue('{out_dir}/allele_post_{i}.tsv'), filetype="tsv")
            # bulk_subtrees_colnames = c("cnv_state_post", "state_post", "p_up", "haplo_post", "haplo_naive", "theta_hat_roll", "phi_mle_roll", "lambda", "gamma")
            if (old_index) {
                self$bulk_subtrees = read_file(inputfile=glue('{out_dir}/bulk_subtrees_{i-1}.tsv.gz'), filetype="tsv")
            } else {
                self$bulk_subtrees = read_file(inputfile=glue('{out_dir}/bulk_subtrees_{i}.tsv.gz'), filetype="tsv")
            }
            # bulk_clones_colnames = c("n_genes", "n_snps", "seg_start", "seg_end", "theta_hat", "theta_mle", "theta_sigma")
            self$bulk_clones = read_file(inputfile=glue('{out_dir}/bulk_clones_{i}.tsv.gz'), filetype="tsv")
            # segs_consensus_colnames = c("sample", "CHROM", "seg", "cnv_state", "cnv_state_post", "seg_start", "seg_end", "seg_start_index", 
            #                                 "seg_end_index", "theta_mle", "theta_sigma", "phi_mle", "phi_sigma", "p_loh", "p_del", "p_amp",           
            #                                 "p_bamp", "p_bdel", "LLR", "LLR_y", "LLR_x", "n_genes", "n_snps", "component","LLR_sample", "seg_length", "seg_cons", "n_states", "cnv_states")
            if (old_index) {
                self$segs_consensus = read_file(inputfile=glue('{out_dir}/segs_consensus_{i-1}.tsv'), filetype="tsv")
            } else {
                self$segs_consensus = read_file(inputfile=glue('{out_dir}/segs_consensus_{i}.tsv'), filetype="tsv")
            }

            self$segs_consensus = self$segs_consensus %>% mutate(CHROM = factor(CHROM))
            
            # tree_post_colnames =  c("mut_nodes", "gtree", "l_matrix")
            self$tree_post = read_file(inputfile=glue('{out_dir}/tree_post_{i}.rds'), filetype="rds")
            self$mut_graph = read_file(inputfile=glue('{out_dir}/mut_graph_{i}.rds'), filetype="rds")
            self$gtree = read_file(inputfile=glue('{out_dir}/tree_final_{i}.rds'), filetype="rds")
            # clone_post_colnames = c("cell", "clone_opt", "GT_opt", "p_opt")
            self$clone_post = read_file(inputfile=glue('{out_dir}/clone_post_{i}.tsv'), filetype="tsv")
            ## gene names are the column names
            self$gexp_roll_wide = read_file(inputfile=glue('{out_dir}/gexp_roll_wide.tsv.gz'), filetype="tsv")
            
            if (!is.null(self$gexp_roll_wide)) {
                if ('V1' %in% colnames(self$gexp_roll_wide)) {
                    self$gexp_roll_wide = self$gexp_roll_wide %>% rename(cell = V1)
                }
                self$gexp_roll_wide = self$gexp_roll_wide %>% tibble::column_to_rownames('cell')
            }
            
            self$hc = read_hc_rds(inputfile=glue('{out_dir}/hc.rds'))
            self$geno = read_file(inputfile=glue('{out_dir}/geno_{i}.tsv'), filetype="tsv")
                
            if (!is.null(self$geno)) {
                if ('V1' %in% colnames(self$geno)) {
                    self$geno = self$geno %>% rename(cell = V1)
                }
                self$geno = self$geno %>% tibble::column_to_rownames('cell') %>% as.matrix
            }
    })
)

#' @keywords internal
check_fread_works = function(input) {
    tryCatch({
        return(data.table::fread(input))
    },
    error = function(e){
        message(paste0("Could not read the input file ", input, " with data.table::fread(). Please check that the file is valid."))
        return(NULL)
    })
}

#' @keywords internal
check_rds_works = function(input) {
    tryCatch({
        return(readRDS(input))
    },
    error = function(e){
        message(paste0("Could not read the input file ", input, " with readRDS(). Please check that the file is valid."))
        return(NULL)
    })
}



#' @keywords internal
return_missing_columns = function(file, expected_colnames = NULL) {
    ## if user sets expected_colnames = NULL, return NULL
    if (is.null(expected_colnames)) {
        return(NULL)
    }
    if (!is.vector(expected_colnames) || !is.character(expected_colnames)) {
        stop("The parameter 'expected_colnames' needs to be a character vector")
    }
    '%ni%' <- Negate('%in%')
    if (any(expected_colnames %ni% colnames(file))) {
        missing_columns = expected_colnames[!(expected_colnames %in% colnames(file))]
        if (length(missing_columns) == 0) {
            stop("Some mismatch exists between the expected columns and the columns in the file. This error shouldn't happen. Check and fix.")
        }
        return(missing_columns)
    } else {
        return(NULL)
    }
}


#' @keywords internal
read_file = function(inputfile, expected_colnames = NULL, filetype="tsv") {
    if (filetype == "tsv") {
        file = check_fread_works(inputfile)
    } else if (filetype == "rds") {
        file = check_rds_works(inputfile)       
    } else {
        stop("The parameter 'filetype' must be either 'tsv' or 'rds'. Please fix.")
    }
    potential_missing_columns = return_missing_columns(file, expected_colnames)
    if (!is.null(potential_missing_columns)) {
        stop(paste0("The file ", inputfile, " appears to be malformed; expected column names: ", potential_missing_columns, ". Please fix."))
    } else {
        return(file)
    }
}

#' @keywords internal
read_hc_rds = function(inputfile) {
    file = check_rds_works(inputfile)
    if (!is.list(file)) {
        stop(paste0("The file: ", inputfile, " is malformed; should be a list. Please fix."))
    }        
    hc_colnames = c("merge", "height", "order", "labels", "method", "call", "dist.method")
    '%ni%' <- Negate('%in%')
    if (any(hc_colnames %ni% names(file))) {
        missing_columns = expected_colnames[!(expected_colnames %in% names(file))]
        stop(paste0("The file ", inputfile, " appears to be malformed; expected column names: ", potential_missing_columns, ". Please fix."))
    } else {
        return(file)
    }
}


            



#' @title Numbat R6 class
#' @description Used to allow users to plot results
#' @export
Numbat <- R6::R6Class("Numbat", lock_objects=FALSE,
  public = list(

    #' @field label sample name
    label = 'sample',
    
    #' @field gtf transcript annotation
    gtf = NULL,

    #' @field joint_post joint posterior
    joint_post = NULL,

    #' @field exp_post expression posterior
    exp_post = NULL,

    #' @field allele_post allele posetrior
    allele_post = NULL,

    #' @field bulk_subtrees bulk profiles of lineage subtrees
    bulk_subtrees = NULL,

    #' @field bulk_clones bulk profiles of clones
    bulk_clones = NULL,

    #' @field segs_consensus consensus segments
    segs_consensus = NULL,

    #' @field tree_post tree posterior
    tree_post = NULL,

    #' @field mut_graph mutation history graph
    mut_graph = NULL,

    #' @field gtree single-cell phylogeny
    gtree = NULL,

    #' @field clone_post clone posteriors
    clone_post = NULL,

    #' @field gexp_roll_wide smoothed expression of single cells
    gexp_roll_wide = NULL,

    #' @field hc hclust object for initial clustering
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

        fetch_results = function(out_dir, i = 2) {
            # joint_post_colnames = c("cell", "CHROM",  "seg",  "cnv_state")
            self$joint_post = read_file(inputfile=glue('{out_dir}/joint_post_{i}.tsv'), filetype="tsv")
            # exp_post_colnames = c("seg", "cnv_state", "n", "phi_mle")
            self$exp_post = read_file(inputfile=glue('{out_dir}/exp_post_{i}.tsv'), filetype="tsv")
            # allele_post_colnames = c("cell", "CHROM", "seg", "cnv_state", "major", "minor", "total", "MAF", "seg_start", "seg_end", "prior_loh", "prior_amp", "prior_del", "prior_bamp", "prior_bdel")
            self$allele_post = read_file(inputfile=glue('{out_dir}/allele_post_{i}.tsv'), filetype="tsv")
            # bulk_subtrees_colnames = c("cnv_state_post", "state_post", "p_up", "haplo_post", "haplo_naive", "theta_hat_roll", "phi_mle_roll", "lambda", "gamma")
            self$bulk_subtrees = read_file(inputfile=glue('{out_dir}/bulk_subtrees_{i}.tsv.gz'), filetype="tsv")
            # bulk_clones_colnames = c("n_genes", "n_snps", "seg_start", "seg_end", "theta_hat", "theta_mle", "theta_sigma")
            self$bulk_clones = read_file(inputfile=glue('{out_dir}/bulk_clones_{i}.tsv.gz'), filetype="tsv")
            # segs_consensus_colnames = c("sample", "CHROM", "seg", "cnv_state", "cnv_state_post", "seg_start", "seg_end", "seg_start_index", 
            #                                 "seg_end_index", "theta_mle", "theta_sigma", "phi_mle", "phi_sigma", "p_loh", "p_del", "p_amp",           
            #                                 "p_bamp", "p_bdel", "LLR", "LLR_y", "LLR_x", "n_genes", "n_snps", "component","LLR_sample", "seg_length", "seg_cons", "n_states", "cnv_states")
            self$segs_consensus = read_file(inputfile=glue('{out_dir}/segs_consensus_{i}.tsv'), filetype="tsv")
            # tree_post_colnames =  c("mut_nodes", "gtree", "l_matrix")
            self$tree_post = read_file(inputfile=glue('{out_dir}/tree_post_{i}.rds'), filetype="rds")
            self$mut_graph = read_file(inputfile=glue('{out_dir}/mut_graph_{i}.rds'), filetype="rds")
            self$gtree = read_file(inputfile=glue('{out_dir}/tree_final_{i}.rds'), filetype="rds")
            # clone_post_colnames = c("cell", "clone_opt", "GT_opt", "p_opt")
            self$clone_post = read_file(inputfile=glue('{out_dir}/clone_post_{i}.tsv'), filetype="tsv")
            ## gene names are the column names
            self$gexp_roll_wide = read_file(inputfile=glue('{out_dir}/gexp_roll_wide.tsv.gz'), filetype="tsv")
            self$hc = read_hc_rds(inputfile=glue('{out_dir}/hc.rds'))
        }
    )
)



#' @keywords internal
check_fread_works = function(input) {
    tryCatch({
        return(data.table::fread(input))
    },
    error = function(e){
        stop(paste0("Could not read the input file ", input, " with data.table::fread(). Please check that the file is valid."))
    })
}

#' @keywords internal
check_rds_works = function(input) {
    tryCatch({
        return(readRDS(input))
    },
    error = function(e){
        stop(paste0("Could not read the input file ", input, " with readRDS(). Please check that the file is valid."))
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
read_file = function(inputfile, expected_colnames, filetype="tsv") {
    if (filetype == "tsv") {
        file = check_fread_works(inputfile)
    } else if (filetype == "rds") {
        file = check_rds_works(inputfile)
        ## all *rds files here should be lists
        if (!is.list(file)) {
            stop(paste0("The file: ", inputfile, " is malformed; should be a list. Please fix."))
        }        
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


            


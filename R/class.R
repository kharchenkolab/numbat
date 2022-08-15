## for magrittr and dplyr functions 
if (getRversion() >= "2.15.1"){
  utils::globalVariables(c(".", "AD", "AD_all", "ALT", "AR", "CHROM", "DP",  
        "DP_all", "FILTER", "FP", "GT", "ID", "LLR", "LLR_sample", 
        "LLR_x", "LLR_y", "L_x_a", "L_x_d", "L_x_n", "L_y_a", 
        "L_y_d", "L_y_n", "MAF", "OTH", "OTH_all POS",
        "QUAL", "REF", "TP", "UQ", "Y_obs", "Z", "Z_amp", "Z_bamp", 
        "Z_bdel", "Z_clone", "Z_clone_x",
        "Z_clone_y", "Z_cnv", "Z_del", "Z_loh", "Z_n", "acen_hg38", 
        "annot", "avg_entropy", "bkp",
        "branch", "cM", "cell", "cell_index", "chrom_sizes_hg38", 
        "clone", "clone_opt", "clone_size",
        "cluster", "cnv_state", "cnv_state_expand",
        "cnv_state_map", "cnv_state_post",
        "cnv_states", "compartment", "component", "cost", 
        "d_obs", "diploid", "down", "edges", "end_x",
        "end_y", "exp_rollmean", "expected_colnames", "extract", "frac_overlap_x",
        "frac_overlap_y", "from", "from_label", "from_node", "gaps_hg38", "gene gene_end",
        "gene_index", "gene_length", "gene_snps", "gene_start", 
        "genetic_map", "get_desc", "group", "gtf", "haplo", "haplo_naive", "haplo_post", 
        "haplo_theta_min", "het", "hom_alt", "i",
        "inter_snp_cm", "inter_snp_dist", "isTip", "is_desc", 
        "j", "keep", "l", "l00", "l00_x", "l00_y",
        "l10", "l10_x", "l10_y", "l11", "l11_x", "l11_y", "l20", "l20_x", 
        "l20_y", "l21", "l21_x", "l21_y", "l22",
        "l22_x", "l22_y", "l31", "l31_x", "l31_y", "l_adj", "l_clone", 
        "l_clone_x", "l_clone_y", "l_x", "l_y",
        "label", "lambda_obs", "lambda_ref", "lambdas_ref", 
        "last_mut", "leaf", "len_overlap",
        "len_x", "len_y", "lnFC", "lnFC_i", "lnFC_j", "lnFC_max_i", "lnFC_max_j",
        "logBF", "logBF_x","logBF_y", "logFC", "loh", "major", "major_count", "marker_index", 
        "min_depth", "minor", "minor_count", "mu", "mut_burden", "n_chrom_snp", "n_genes", "n_mut", 
        "n_sibling", "n_snps","n_states", "n_x", "n_y", "name", "node", "node_phylo", "nodes", "p", "pAD", 
        "pBAF", "p_0", "p_1", "p_amp", "p_bamp", "p_bdel", "p_cnv", "p_cnv_expand", "p_del", "p_loh", "p_max", 
        "p_n", "p_neu", "p_s", "p_up", "phi_mle", "phi_mle_roll", "phi_sigma", "pnorm.range", 
        "potential_missing_columns","precision", "prior_amp", "prior_bamp", "prior_bdel",
        "prior_clone prior_del",
        "prior_loh", "root", "s", "seg", "seg_cons", "seg_end", 
        "seg_end_index", "seg_label",
        "seg_length", "seg_start", "seg_start_index", 
        "segs_consensus", "seqnames",
        "set_colnames", "sig", "site", "size", "snp_id", "snp_index", 
        "snp_rate", "start_x", "start_y",
        "state", "state_post", "superclone", "theta_hat", 
        "theta_level", "theta_mle",
        "theta_sigma", "to", "to_label", "to_node", "total", "value", 
        "variable", "vcf_meta", "vertex",
        "vp", "w", "width", "write.vcf", "x", "y"))
}






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
            
            self$tree_post = read_file(inputfile=glue('{out_dir}/tree_post_{i}.rds'), filetype="rds")
            self$mut_graph = read_file(inputfile=glue('{out_dir}/mut_graph_{i}.rds'), filetype="rds")
            self$gtree = read_file(inputfile=glue('{out_dir}/tree_final_{i}.rds'), filetype="rds")
            self$clone_post = read_file(inputfile=glue('{out_dir}/clone_post_{i}.tsv'), filetype="tsv")
            self$gexp_roll_wide = read_file(inputfile=glue('{out_dir}/gexp_roll_wide.tsv.gz'), filetype="tsv")
            
            if (!is.null(self$gexp_roll_wide)) {
                if ('V1' %in% colnames(self$gexp_roll_wide)) {
                    self$gexp_roll_wide = self$gexp_roll_wide %>% rename(cell = V1)
                }
                self$gexp_roll_wide = self$gexp_roll_wide %>% tibble::column_to_rownames('cell')
            }
            
            self$hc = read_hc_rds(inputfile=glue('{out_dir}/hc.rds'))
            
    })
)

relevel_chrom = function(df) {
    if (!is.null(df)) {
        df = df %>% mutate(CHROM = factor(CHROM, 1:22))
    }
    return(df)
}

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


            


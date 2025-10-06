#####################################################################
# generate binned RNA matrix from a seurat object with RNA assay
#####################################################################

# Set seed for reproducibility
set.seed(123)

# Load required libraries
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(optparse)
  library(Seurat)
  library(dplyr)
})
options(stringsAsFactors = FALSE)

run_as_script <- sys.nframe() == 0
args <- NULL

if (run_as_script) {
  option_list <- list(
    make_option('--rnaCountsFile', type = "character", default = NULL,
                help = "rds file with RNA counts; if a seurat object is provided, it will retrieve seurat@assays[['RNA']]@counts. Also handles 10x .h5"),
    make_option('--geneBinMapCSVFile', type = "character", default = NULL,
                help = "gene bin map csv file, generated using get_gene_binned_intersections.R"),
    make_option('--outFile', type = "character", default = NULL,
                help = "output file for binned counts"),
    make_option('--barcodesKeep', type = "character", default = NULL,
                help = "barcodes to keep in supplied object; if None it will keep all. If used for generating aggregated references, then it needs to have 'cell' and 'group' as column names for barcode and cell annotation"),
    make_option('--generateAggRef', action = "store_true", default = FALSE,
                help = "Indicator for whether we are generating aggregated reference or just a count matrix")
  )

  args <- parse_args(OptionParser(option_list = option_list))
}

#####################################################################
# gene2bin converts gene-level count data to genomic bin-level data
#####################################################################
#' Select barcodes from a count matrix for downstream processing.
#'
#' @param count_mat Matrix-like object of counts with cell barcodes as column names.
#' @param barcodesKeep Either `NULL`, a character path to a tabular file, or a
#'   data.frame whose first column contains the barcodes to keep.
#'
#' @return Character vector of barcodes present in `count_mat` to retain. If
#'   `barcodesKeep` is `NULL`, all barcodes in `count_mat` are returned.
select_barcodes <- function(count_mat, barcodesKeep) {
  if (is.null(barcodesKeep)) {
    message("No barcodesKeep provided; keeping all cells in rnaCountsFile")
    return(colnames(count_mat))
  }

  if (is.character(barcodesKeep)) {
    message("barcodesKeep provided as file path; reading barcodes via data.table::fread")
    barcodes <- as.list(data.table::fread(barcodesKeep)[[1]])
  } else if (inherits(barcodesKeep, "data.frame")) {
    message("barcodesKeep provided as data.frame; using first column as barcodes")
    barcodes <- as.list(barcodesKeep[[1]])
  } else {
    stop("barcodesKeep must be NULL, a character path, or a data.frame")
  }

  count_colnames <- colnames(count_mat)
  matched_cols <- character(0)

  for (bc in barcodes) {
    if (bc %in% count_colnames) {
      matched_cols <- c(matched_cols, bc)
    } else {
      suffix_match <- grep(paste0("_", bc, "$"), count_colnames, value = TRUE)
      if (length(suffix_match) > 0) {
        matched_cols <- c(matched_cols, suffix_match)
      } else {
        warning(paste("Barcode", bc, "not found in count matrix"))
      }
    }
  }

  if (length(matched_cols) == 0) {
    stop("No matching barcodes found in the count matrix")
  }

  message(paste("Found", length(matched_cols), "matching cells out of", length(barcodes), "requested"))
  unique(matched_cols)
}

#' Convert gene-level counts to genomic bin-level counts.
#'
#' @param cnt Sparse or dense matrix of gene-by-cell counts.
#' @param gene_binmap_file Path to a CSV containing the gene-to-bin mapping with
#'   at least `gene_name` and `bingrp` columns.
#'
#' @return Sparse matrix (`dgCMatrix`) of bin-by-cell counts.
gene2bin <- function(cnt, gene_binmap_file) {
  if (missing(gene_binmap_file) || is.null(gene_binmap_file)) {
    stop("gene_binmap_file must be provided")
  }
  # Load presto package for efficient sparse matrix operations
  if (!requireNamespace("presto", quietly = TRUE)) {
    stop("Package 'presto' is needed for this function to work. Please install it.")
  }
  # load gene to bin mapping
  gene_binmap <- read.csv(gene_binmap_file)
  print(head(gene_binmap))
  print(head(rownames(cnt)))
  
  # filter the mapping to include only genes in the count matrix
  gene_binmap <- gene_binmap %>% 
    filter(gene_name %in% rownames(cnt)) %>%
    mutate(bingrp = forcats::fct_drop(bingrp))
  
  # ensure count matrix is in dgCMatrix format and keep only genes in the mapping
  if (class(cnt) != "dgCMatrix") {
    cnt <- as(cnt, "dgCMatrix")
  }
  cnt <- cnt[gene_binmap$gene_name, ]
  
  # sum counts by grouped bin 
  sumGroups <- function(X, y, MARGIN = 2) {
    require(presto)
    
    cpp_sumGroups_dgc <- function(x, p, i, ncol, groups, ngroups) {
      .Call('_presto_cpp_sumGroups_dgc', PACKAGE = 'presto', 
            x, p, i, ncol, groups, ngroups)
    }
    
    cpp_sumGroups_dgc_T <- function(x, p, i, ncol, nrow, groups, ngroups) {
      .Call('_presto_cpp_sumGroups_dgc_T', PACKAGE = 'presto', 
            x, p, i, ncol, nrow, groups, ngroups)
    }
    
    if (MARGIN == 1) {
      cpp_sumGroups_dgc_T(X@x, X@p, X@i, ncol(X), nrow(X), 
                          as.integer(y) - 1, length(unique(y)))
    } else {
      cpp_sumGroups_dgc(X@x, X@p, X@i, ncol(X), 
                         as.integer(y) - 1, length(unique(y)))
    }
  }
  
  # sum gene counts by genomic bin
  binsum <- sumGroups(cnt, gene_binmap$bingrp)
  colnames(binsum) <- colnames(cnt)
  rownames(binsum) <- levels(gene_binmap$bingrp)
  
  # convert to sparse matrix format
  bin_sparse <- as(binsum, "dgCMatrix")
  return(bin_sparse)
}

#' Load RNA counts from an RDS Seurat object or a 10x HDF5 file.
#'
#' @param filepath Path to an `.rds` or `.h5` file containing RNA counts.
#'
#' @return Count matrix or Seurat object extracted from `filepath`.
read_rna_counts <- function(filepath){
	# use tools to get the file extension
	ext = tolower(tools::file_ext(filepath))
	if (ext == 'rds'){
		message('RDS file detected')
		return(readRDS(filepath))
	} else if (ext == 'h5') {
		message("Handling HDF5 file with Seurat's Read10X_h5...\n")
		return(Read10X_h5(filepath))
	} else {
		message('Unknown file type? Exiting')
		stop()
	}
}

#####################################################################
# sample_RNAbin porocesses RNA-seq data, bin gene counts, and save files
#####################################################################
#' Generate binned RNA counts for a single sample.
#'
#' @param rnaCountsFile Path to the RNA count source (RDS Seurat object/matrix or
#'   10x HDF5 file).
#' @param outFileName Destination file for the binned counts RDS.
#' @param gene_binmap_file CSV file describing gene-to-bin mapping.
#' @param barcodesKeep Optional barcodes filter (file path or data.frame).
#' @param save Logical indicating whether to write the binned counts to disk.
#'
#' @return Invisible `dgCMatrix` of binned counts when `save = FALSE`.
sample_RNAbin <- function(rnaCountsFile, outFileName, gene_binmap_file, barcodesKeep = NULL, save = TRUE) {
  if (missing(rnaCountsFile) || is.null(rnaCountsFile)) {
    stop("rnaCountsFile must be provided")
  }
  if (missing(outFileName) || is.null(outFileName)) {
    stop("outFileName must be provided")
  }
  if (missing(gene_binmap_file) || is.null(gene_binmap_file)) {
    stop("gene_binmap_file must be provided")
  }
  # Load RNA counts file
  message("Reading RNA counts file...")
  counts_data = read_rna_counts(rnaCountsFile)
  # Check if the input is a Seurat object or a counts matrix
  if (inherits(counts_data, "Seurat")) {
    message("Detected Seurat object, extracting RNA counts...")
    count_mat <- counts_data@assays[[1]]@counts
  } else {
    message("Using provided counts matrix directly...")
    count_mat <- counts_data
  }
  
  # Derive the set of columns to keep
  selected_cols <- select_barcodes(count_mat, barcodesKeep)
  count_mat <- count_mat[, selected_cols, drop = FALSE]
  colnames(count_mat) <- selected_cols
  binned_cnts <- gene2bin(count_mat, gene_binmap_file)
  if (isTRUE(save)) {
	  saveRDS(binned_cnts, outFileName)
	  message("Binned counts saved to ", outFileName)
  } else {
    return(binned_cnts)
  }
}


#' Generate aggregated binned RNA reference counts by group.
#'
#' @param rnaCountsFile Path to RNA count source (as in `sample_RNAbin`).
#' @param outFileName Destination file for the aggregated reference counts.
#' @param gene_binmap_file CSV gene-to-bin mapping file.
#' @param barcodesKeep File path or data.frame with barcode-to-group annotations.
#' @param save Logical indicating whether to write the aggregated reference to disk.
#'
#' @return Invisible aggregated counts matrix when `save = TRUE`, otherwise the
#'   aggregated counts object.
sample_RNAbin_reference <- function(rnaCountsFile, outFileName, gene_binmap_file, barcodesKeep, save = TRUE) {
    if (is.null(barcodesKeep)) {
		stop("barcodesKeep is required when generateAggRef is TRUE; supply a file with 'cell' and 'group' columns")
	}
    if (missing(rnaCountsFile) || is.null(rnaCountsFile)) {
        stop("rnaCountsFile must be provided")
    }
    if (missing(outFileName) || is.null(outFileName)) {
        stop("outFileName must be provided")
    }
    if (missing(gene_binmap_file) || is.null(gene_binmap_file)) {
        stop("gene_binmap_file must be provided")
    }

    if (is.character(barcodesKeep)) {
		CB_annot_df <- data.table::fread(barcodesKeep)
	} else if (inherits(barcodesKeep, "data.frame")) {
		CB_annot_df <- barcodesKeep
	} else {
		stop("barcodesKeep must be a path to a file or a data.frame with 'cell' and 'group' columns")
	}

    # uses a barcodes cell x group matrix to then pseudobulk the counts 
    binned_cnts <- sample_RNAbin(rnaCountsFile, outFileName, gene_binmap_file, CB_annot_df, save = FALSE)
    #generate an aggregated matrix for numbat reference
    agg_ref_counts <- numbat::aggregate_counts(as.matrix(binned_cnts), 
					  CB_annot_df, 
					  normalized = TRUE, 
					  verbose = TRUE)
    if (isTRUE(save)) {
	  saveRDS(agg_ref_counts, outFileName)
	  message("Binned aggregated reference counts saved to ", outFileName)
	  return(invisible(agg_ref_counts))
    }
    agg_ref_counts
}

#####################################################################
# Main Execution
#####################################################################
if (run_as_script) {
	if (isTRUE(args$generateAggRef)) {
		sample_RNAbin_reference(args$rnaCountsFile,
		                      args$outFile,
		                      args$geneBinMapCSVFile,
		                      args$barcodesKeep)
	} else {
		sample_RNAbin(args$rnaCountsFile,
		             args$outFile,
		             args$geneBinMapCSVFile,
		             args$barcodesKeep)
	}
}

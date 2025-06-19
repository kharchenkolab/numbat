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

if (sys.nframe() == 0) {
  option_list = list(
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

  args = parse_args(OptionParser(option_list = option_list))
  invisible(list2env(args,environment()))
}

#####################################################################
# gene2bin converts gene-level count data to genomic bin-level data
#####################################################################
gene2bin <- function(cnt, gene_binmap = geneBinMapCSVFile) {
  # Load presto package for efficient sparse matrix operations
  if (!requireNamespace("presto", quietly = TRUE)) {
    stop("Package 'presto' is needed for this function to work. Please install it.")
  }
  # load gene to bin mapping
  gene_binmap <- read.csv(geneBinMapCSVFile)
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
sample_RNAbin <- function(rnaCountsFile, outFileName, barcodesKeep, save=TRUE) {
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
  
  # take only the barcodes in the barcodesKeep vector
  if (is.character(barcodesKeep)) {
     cat("It's a character string. Proceeding with string logic...\n")
     barcodes = as.list(data.table::fread(barcodesKeep)[[1]])
  } else if (inherits(barcodesKeep, "data.frame")) {
     cat("Supplied barcodes is a data frame type. Taking the first column as barcode col \n")
  	barcodes <- as.list(barcodesKeep[[1]])
  } else {
	  stop('Cant figure out the type of barcode names supplied...')
  }
  
  # match barcodes to column names of the count matrix
  # in case the count matrix has {sample_name}_barcode in the column name
  count_colnames <- colnames(count_mat)
  matched_cols <- c()
  
  for (bc in barcodes) {
    # Check if barcode exists exactly as is
    if (bc %in% count_colnames) {
      matched_cols <- c(matched_cols, bc)
    } else {
      # Check for suffix match (e.g., sample_name_barcode format)
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
  count_mat = count_mat[, matched_cols]
  colnames(count_mat) = matched_cols 
  binned_cnts <- gene2bin(count_mat)
  if(save==TRUE){
	  saveRDS(binned_cnts, outFileName)
	  message("Binned counts saved to ", outFileName)
  }else{return(binned_cnts)}
}


sample_RNAbin_reference <- function(seuratFileName, outFileName, barcodesKeep, save = TRUE) {
    CB_annot_df = data.table::fread(barcodesKeep)

    # uses a barcodes cell x group matrix to then pseudobulk the counts 
    binned_cnts = sample_RNAbin(rnaCountsFile, outFileName, CB_annot_df, save=FALSE)
    #generate an aggregated matrix for numbat reference
    agg_ref_counts = numbat::aggregate_counts(as.matrix(binned_cnts), 
						  CB_annot_df, 
						  normalized = TRUE, 
						  verbose = TRUE)
    if(save==TRUE){
	  saveRDS(agg_ref_counts, outFileName)
	  message("Binned aggregated reference counts saved to ", outFileName)
    }
}

#####################################################################
# Main Execution
#####################################################################
if (sys.nframe() == 0) {
	if (generateAggRef) {
		sample_RNAbin_reference(rnaCountsFile, outFile, barcodesKeep)
	} else {
		sample_RNAbin(rnaCountsFile, outFile, barcodesKeep)
	}
}

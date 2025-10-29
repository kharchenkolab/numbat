pacman::p_load(Matrix)

ensure_matrix <- function(x) {
  if (inherits(x, "table")) {
    x <- unclass(as.matrix(x))
  } else if (inherits(x, "array")) {
    x <- as.matrix(x)
  }
  x
}

loadcnt <- function(cfile){
  if(endsWith(cfile,".rds")){
    df_count <- readRDS(cfile)
    if (!inherits(df_count, "dgCMatrix")) {
    df_count <- as(ensure_matrix(df_count), "dgCMatrix")
   }
  }else{
    suppressWarnings({
      df_count =data.table::fread(cfile) |>
        tibble::column_to_rownames("V1") 
      df_count = as(as.matrix(df_count),"dgCMatrix")
    })
  }
  return(df_count)
}
agg_ref<- function(cnt,maxn,x,s){
  set.seed(s)
  sampledCB <- sample(colnames(cnt),maxn,replace = T)
  agg_cnt <- numbat::aggregate_counts(cnt[,sampledCB],data.frame(cell=sampledCB,group=x))
  return(agg_cnt)
}
agg_refs <- function(cfiles,sampleN,s){
  aggcntL <- imap(cfiles,\(f,i){
    df_count <- loadcnt(f)
    return(agg_ref(df_count,sampleN,i,s))
  })
  return(bind_cols(aggcntL) %>% as.data.frame() %>% 
           set_rownames(
             rownames(aggcntL[[1]])
           ))
}
cntsubN <- function(cnt,f,maxn,seed){
  if(ncol(cnt) < maxn){
   return(cnt[f,])
  }else{
  set.seed(seed)
  return(cnt[f,sample(colnames(cnt),size=min(maxn,ncol(cnt)))])
  }
}
binCnt <- function(bincntF,seed,maxCB=10000,add_suffix=FALSE,suffix_tag=c("_RNA","_ATAC")){
  cntmat <- map(bincntF,loadcnt)
  sharedFeature <- intersect(rownames(cntmat[[1]]),rownames(cntmat[[2]]))
  
  if(add_suffix){
    if(length(suffix_tag) != 2){
      stop("suffix_tag must be a character vector of length 2")
    }
    
    # get all unique column names across both matrices
    all_cols_1 <- colnames(cntmat[[1]])
    all_cols_2 <- colnames(cntmat[[2]])
    all_unique_cols <- unique(c(all_cols_1, all_cols_2))
    
    # subsample from all barcodes
    set.seed(seed)
    n_to_sample <- min(maxCB, length(all_unique_cols))
    sampled_cells <- sample(all_unique_cols, size = n_to_sample)
    
    result_cols <- character()
    result_mat <- NULL
    
    for(cell in sampled_cells){
      cell_in_1 <- cell %in% all_cols_1
      cell_in_2 <- cell %in% all_cols_2
      
      if(cell_in_1){
        cell_data_1 <- cntmat[[1]][sharedFeature, cell, drop=FALSE]
        colnames(cell_data_1) <- paste0(cell, suffix_tag[1])
        if(is.null(result_mat)){
          result_mat <- cell_data_1
        } else {
          result_mat <- cbind(result_mat, cell_data_1)
        }
        result_cols <- c(result_cols, colnames(cell_data_1))
      }
      
      if(cell_in_2){
        cell_data_2 <- cntmat[[2]][sharedFeature, cell, drop=FALSE]
        colnames(cell_data_2) <- paste0(cell, suffix_tag[2])
        if(is.null(result_mat)){
          result_mat <- cell_data_2
        } else {
          result_mat <- cbind(result_mat, cell_data_2)
        }
        result_cols <- c(result_cols, colnames(cell_data_2))
      }
    }
    
    colnames(result_mat) <- result_cols
    return(result_mat)
    
  } else {
    # no suffix needed, subsample independently then combine
    cnt1_sub <- cntsubN(cntmat[[1]],sharedFeature,maxCB,seed)
    cnt2_sub <- cntsubN(cntmat[[2]],sharedFeature,maxCB,seed)
    comb_cnt <- cbind(cnt1_sub,cnt2_sub)
    return(comb_cnt)
  }
}

library(purrr)
library(Matrix)

binCnt_union <- function(bincntF, seed, maxCB = 10000, 
                         add_suffix = FALSE, suffix_tag = c("_RNA","_ATAC"), 
                         runMultiomeAsSummed = FALSE) {
  # 1. Load the two count matrices
  cnt1 <- loadcnt(bincntF[1])
  cnt2 <- loadcnt(bincntF[2])

  allFeat <- union(rownames(cnt1), rownames(cnt2))
  
  # 3. For each matrix, add zero‐rows for any features it's missing
  add_missing <- function(mat, allFeat) {
    missing <- setdiff(allFeat, rownames(mat))
    if (length(missing)) {
      # create a sparse zero matrix with those features
      zero_block <- Matrix(0,
                           nrow = length(missing),
                           ncol = ncol(mat),
                           dimnames = list(missing, colnames(mat)))
      mat <- rbind(mat, zero_block)
    }
    mat[allFeat, , drop = FALSE]
  }
  cnt1_u <- add_missing(cnt1, allFeat)
  cnt2_u <- add_missing(cnt2, allFeat)

  if(runMultiomeAsSummed){
    print('Summing read counts for each modality, resulting one summed feature per cell barcode. Use the --add_suffix and --suffix_tag to run multiome as paired. Note: this is only applicable if you are running multiome data as paired and want each barcode to be counted once.  If you are looking for the true Multiome version, please follow the appropriate steps in the vignette.')
    if(add_suffix){
      print('Multiome can either be used in paired mode or as summed mode for this algorithm. Either --runMultiomeAsSummed can be true or --add_suffix. ')
      stop('Please specify either --runMultiomeAsSummed or --add_suffix, but not both.')
    }
    intersect_cells <- intersect(colnames(cnt1), colnames(cnt2))
    set.seed(seed)
    n_to_sample <- min(maxCB, length(intersect_cells))
    sampled_cells <- sample(intersect_cells, size = n_to_sample)
    cnt <- cnt1[,sampled_cells] + cnt2[,sampled_cells]
    return(cnt)
  }
  
  if(add_suffix){
    print('Adding suffix tags to the barcodes used in the count matrices. Note: this is only applicable if you are running multiome data as paired, since we want to ensure each barcode is unique for the algorithm. If you are looking for the true Multiome version, please follow the appropriate steps in the vignette.')
    if(length(suffix_tag) != 2){
      stop("suffix_tag must be a character vector of length 2")
    }
    
    all_cols_1 <- colnames(cnt1_u)
    all_cols_2 <- colnames(cnt2_u)
    all_unique_cols <- unique(c(all_cols_1, all_cols_2))
    
    # subsample from the union of ALL columns (without suffixes)
    set.seed(seed)
    n_to_sample <- min(maxCB, length(all_unique_cols))
    sampled_cells <- sample(all_unique_cols, size = n_to_sample)
    
    # for each sampled cell, check which modalities it appears in and add appropriate columns
    result_cols <- character()
    result_mat <- NULL
    
    for(cell in sampled_cells){
      cell_in_1 <- cell %in% all_cols_1
      cell_in_2 <- cell %in% all_cols_2
      
      if(cell_in_1){
        cell_data_1 <- cnt1_u[, cell, drop=FALSE]
        colnames(cell_data_1) <- paste0(cell, suffix_tag[1])
        if(is.null(result_mat)){
          result_mat <- cell_data_1
        } else {
          result_mat <- cbind(result_mat, cell_data_1)
        }
        result_cols <- c(result_cols, colnames(cell_data_1))
      }
      
      if(cell_in_2){
        cell_data_2 <- cnt2_u[, cell, drop=FALSE]
        colnames(cell_data_2) <- paste0(cell, suffix_tag[2])
        if(is.null(result_mat)){
          result_mat <- cell_data_2
        } else {
          result_mat <- cbind(result_mat, cell_data_2)
        }
        result_cols <- c(result_cols, colnames(cell_data_2))
      }
    }
    
    colnames(result_mat) <- result_cols
    return(result_mat)
    
  } else {
    # if no suffix addition is needed, thensubsample independently then combine
    cnt1_sub <- cntsubN(cnt1_u, allFeat, maxCB, seed)
    cnt2_sub <- cntsubN(cnt2_u, allFeat, maxCB, seed)
    comb_cnt <- cbind(cnt1_sub, cnt2_sub)
    return(comb_cnt)
  }
}




# input prep for allele dataframe if running multiome as paired
#' Combine and harmonize allele counts from RNA and ATAC modalities
#' @param alleleCounts_RNA Character. Path to the RNA-based allele counts file
#' @param alleleCounts_ATAC Character. Path to the ATAC-based allele counts file
#' @param output_file Character. File for combined allele outputs 
#' @param compress Logical. Whether to compress the output file (default: TRUE)
#' @return Character. Path to the output file
combine_allele_counts <- function(alleleCounts_RNA,
                                  alleleCounts_ATAC,
				  outAlleleCountsFile_RNA, 
				  outAlleleCountsFile_ATAC, 
				  output_file, 
				  addBarcodeSuff = FALSE,
                                  compress = TRUE) {
  geno_files <- c(alleleCounts_RNA, alleleCounts_ATAC)
  
  if (!all(file.exists(geno_files))) {
    missing_files <- geno_files[!file.exists(geno_files)]
    stop("The following genotype files do not exist: ", 
         paste(missing_files, collapse = ", "))
  }
  
  
#  message("Harmonizing allele counts for sample: ", sample_id)
  
  # Read first allele counts file (RNA)
  df1 <- tryCatch({
    data.table::fread(alleleCounts_RNA, header = TRUE, sep = "\t") %>% 
      filter(CHROM %in% 1:22) %>%
      mutate(CHROM = factor(CHROM, 1:22))
  }, error = function(e) {
    stop("Error reading RNA allele counts file: ", e$message)
  })
  
  # Extract unique SNPs and genotypes from RNA allele counts file
  geno1 <- df1 %>% distinct(snp_id, GT)
  
  # Read second allele counts file (ATAC)
  df2 <- tryCatch({
    data.table::fread(alleleCounts_ATAC, header = TRUE, sep = "\t") %>% 
      filter(CHROM %in% 1:22) %>%
      mutate(CHROM = factor(CHROM, 1:22))
  }, error = function(e) {
    stop("Error reading ATAC allele counts file: ", e$message)
  })
  
  # Extract unique SNPs and genotypes from ATAC allele counts file
  geno2 <- df2 %>% distinct(snp_id, GT)
  
  # Find shared SNPs between the two datasets
  sharedSNP <- intersect(geno1$snp_id, geno2$snp_id)
  message("Found ", length(sharedSNP), " shared SNPs between modalities")
  
  # Create reference dataframe from second file's genotypes
  snp_df <- geno2 %>% 
    filter(snp_id %in% sharedSNP) %>% 
    column_to_rownames("snp_id")
  
  # Update genotypes in first file to match second file for shared SNPs
  geno1 <- geno1 %>%
    mutate(GT = case_when(
      snp_id %in% sharedSNP ~ snp_df[snp_id, "GT"],
      TRUE ~ GT
    ))
  
  geno1_df <- geno1 %>% column_to_rownames("snp_id")
  df1$GT <- geno1_df[df1$snp_id, "GT"]

  if(addBarcodeSuff == TRUE){
	  df1$cell = paste0(df1$cell, '_RNA')
	  df2$cell = paste0(df2$cell, '_ATAC')
  } 
  combined_df <- rbind(df1, df2)
  
#  output_file <- file.path(output_dir, paste0(sample_id, "_allele_counts_consistent.tsv"))
  
  tryCatch({
    data.table::fwrite(
      combined_df, 
      output_file, 
      quote = FALSE, 
      row.names = FALSE, 
      #nThread = 3, 
      sep = "\t"
    )
    message("Successfully wrote combined allele counts to: ", output_file)
  }, error = function(e) {
    stop("Error writing combined file: ", e$message)
  })
  
  if (compress) {
    tryCatch({
      system2("gzip", args = output_file)
      output_file <- paste0(output_file, ".gz")
      message("Successfully compressed output file")
    }, error = function(e) {
      warning("Failed to compress output file: ", e$message)
    })
  }
  
  return(output_file)
}


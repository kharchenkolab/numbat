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
binCnt <- function(bincntF,seed,maxCB=10000){
  cntmat <- map(bincntF,loadcnt)
  sharedFeature <- intersect(rownames(cntmat[[1]]),rownames(cntmat[[2]]))
  comb_cnt <- cbind(cntsubN(cntmat[[1]],sharedFeature,maxCB,seed),
                    cntsubN(cntmat[[2]],sharedFeature,maxCB,seed))
  return(comb_cnt)
}

library(purrr)
library(Matrix)

binCnt_union <- function(bincntF, seed, maxCB = 10000) {
  # 1. Load the two count matrices
  cnt1 <- loadcnt(bincntF[1])
  cnt2 <- loadcnt(bincntF[2])

  allFeat <- union(rownames(cnt1), rownames(cnt2))
  
  # 3. For each matrix, add zero‐rows for any features it’s missing
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
    # reorder to match allFeat
    mat[allFeat, , drop = FALSE]
  }
  cnt1_u <- add_missing(cnt1, allFeat)
  cnt2_u <- add_missing(cnt2, allFeat)
  
  # 4. Subsample & combine
  comb_cnt <- cbind(
    cntsubN(cnt1_u, allFeat, maxCB, seed),
    cntsubN(cnt2_u, allFeat, maxCB, seed)
  )
  comb_cnt
}

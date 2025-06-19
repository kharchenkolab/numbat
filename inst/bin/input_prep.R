pacman::p_load(Matrix)
loadcnt <- function(cfile){
  if(endsWith(cfile,".rds")){
    df_count <- readRDS(cfile)
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
  set.seed(seed)
  return(cnt[f,sample(colnames(cnt),size=min(maxn,ncol(cnt)))])
}
binCnt <- function(bincntF,seed,maxCB=10000){
  cntmat <- map(bincntF,loadcnt)
  sharedFeature <- intersect(rownames(cntmat[[1]]),rownames(cntmat[[2]]))
  comb_cnt <- cbind(cntsubN(cntmat[[1]],sharedFeature,maxCB,seed),
                    cntsubN(cntmat[[2]],sharedFeature,maxCB,seed))
  return(comb_cnt)
}


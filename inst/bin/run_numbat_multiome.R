#==== Parse arguments ====
#packageVersion('numbat') ‘1.3.2.1’
library(optparse)
options(stringsAsFactors = F)
option_list = list(make_option("--countmat", default = NULL),
                   make_option("--alleledf", default = NULL),
                   make_option("--gtf", default = NULL),
                   make_option("--ref", default = NULL),
                   make_option("--out_dir", default = NULL),
                   make_option("--parL", default = NULL)
                   )
args = parse_args(OptionParser(option_list = option_list))

library(numbat)
library(dplyr)
print(packageVersion('numbat'))

#==== Load GRanges file and turn to GTF ====
# Load GTF and if it's IRanges, make it Numbat readable:
gr_bins <- readRDS(args$gtf)
if (inherits(gr_bins, "GRanges")) {
	gr_bins$gene = with(as.data.frame(gr_bins),paste0(seqnames,":",start,"-",end))
	gr_bins$gene_length = width(gr_bins)
	gr_bins$gene_start = start(gr_bins)
	gr_bins$gene_end = end(gr_bins)
	gr_bins$CHROM = gsub("chr","",seqnames(gr_bins))

	gr_bins_df = as.data.frame(gr_bins) %>% filter(CHROM!="X") |> 
	  dplyr::select(-strand) |>
	  dplyr::mutate_at(c("gene_start","gene_end","gene_length"),as.integer)
	rownames(gr_bins_df) = gr_bins_df$gene
	gtf <- gr_bins_df[,colnames(numbat::gtf_hg38)]
} else {
  tryCatch({
    if (!identical(colnames(gr_bins), colnames(numbat::gtf_hg38))) {
      stop("Supplied GTF column names do not match the expected Numbat GTF.")
    }}, 
    error = \(e){
      stop(sprintf("An error occurred during GTF validation: %s", e$message))
  })
}

#####################################
frac <- 0.6
#==== Load Count Matrix ====
cat('Loading data\n')
if(endsWith(args$countmat,".rds")){
  df_count <- readRDS(args$countmat) |> as.matrix()
}else{
  suppressWarnings({df_count =data.table::fread(args$countmat) |>
    tibble::column_to_rownames("V1") |>
    as.matrix()})
}

## load phasing and allele information
df_allele = read.table(args$alleledf, header=TRUE, sep="\t")
#==== Diagnosis Plots ====
alleleCov <- df_allele %>% group_by(cell) %>% 
  summarise(snp_idN = length(unique(snp_id)))
pdf("allele_CB_coverage_hist.pdf",width=5,height=5)
hist(alleleCov$snp_idN,
     main = "Histogram of SNP coverage per CB", 
     xlab = "",
     breaks = 100,
     xlim = quantile(alleleCov$snp_idN,c(0.01,0.99)),
     col = "white")
abline(v=quantile(alleleCov$snp_idN),col="darkred",lwd=2)
dev.off()
# write.table(alleleCov,"allele_CB_coverage.txt",row.names = F,quote = F)
#==== Subset for initial debug ====
sharedCB <- intersect(df_allele$cell,colnames(df_count))
testMax= 10000
if(length(sharedCB)>testMax){
  set.seed(0)
  CBs <- sample(sharedCB,testMax)
  df_allele <- df_allele %>% filter(cell %in% CBs)
  df_count <- df_count[,CBs]
}
cat('Supplied number of cells: ', ncol(df_count), '\n')
#==== Run Numbat ====
numbatPars <- formals(numbat::run_numbat)
if(!is.null(args$parL)){
    selfPar <- readRDS(args$parL)
    numbatPars[names(selfPar)] <- selfPar
}
source("/build/scripts/gtf_fixed_Numbat.R")

covs <- alleleCov %>% filter(cell %in% df_allele$cell) %>% pull(snp_idN)
numbatPars$ncores <- min(numbatPars$ncores,13)
if(mean(covs)>150){
  df_allele <-  df_allele%>%sample_frac(frac)
}

#==== Parameter change ====
cat('Running Numbat\n')
numbatPars$count_mat = df_count
numbatPars$lambdas_ref =readRDS(args$ref) %>% as.matrix()
numbatPars$df_allele =df_allele
numbatPars$gtf = gtf
numbatPars$out_dir <- args$out_dir
numbatPars$ncores_nni = numbatPars$ncores
numbatPars$p_multi = 1-numbatPars$alpha
numbatPars$max_cost = ncol(numbatPars$count_mat) * numbatPars$tau
do.call(numbat:::run_numbat,numbatPars)

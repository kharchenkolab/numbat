# gets intersections of any numbat gtf file with a supplied binned genome
suppressPackageStartupMessages({
library(GenomicRanges)
library(optparse)
options(stringsAsFactors = F)
library(dplyr)
})


option_list = list(make_option('--numbatGTFname', type = "character", default = NULL,
		           help = paste("one of numbat based GTF files such as hg38, hg19, or mm10,",
			              " used for obtaining genomic locations for gene row names",
						  " later used in the GEX count matrix")),
                   make_option('--binGR', type = "character", default = NULL, 
		           help = "binned genome genomic ranges file (200kb by default)"),
                   make_option('--outfile', type = "character", default = NULL, 
				   help = "output file")
                   )

args = parse_args(OptionParser(option_list = option_list))


numbat_gtf = args$numbatGTFname
binGR = args$binGR
outfile = args$outfile


if(numbat_gtf=='hg19'){
	gtf = numbat::gtf_hg19
} else if(numbat_gtf=='hg38'){
	gtf = numbat::gtf_hg38
} else if(numbat_gtf=='mm10'){
	gtf = numbat::gtf_mm10
}else{
  gtf = read.table(numbat_gtf, header=TRUE, stringsAsFactors = FALSE)
}

numbat_to_gr <- function(numbat_gtf){
	gr_gtf = GRanges(
			 seqnames = paste0('chr', numbat_gtf$CHROM),
			 ranges = IRanges(start = numbat_gtf$gene_start, 
					  end = numbat_gtf$gene_end),
			 gene = numbat_gtf$gene,
			 gene_length = numbat_gtf$gene_length
			)
	return(gr_gtf)
}

gtf_gr = numbat_to_gr(gtf)
bin_gr = readRDS(binGR)

merged = mergeByOverlaps(gtf_gr, bin_gr)
merged = as.data.frame(merged)
# duplicated genes are 12,133! that's a large number.
#merged[merged$gene %in% merged$gene[duplicated(merged$gene)], ]

# get which gene belongs to which binned group

ov <- findOverlaps(gtf_gr, bin_gr)
int_ranges <- pintersect(gtf_gr[queryHits(ov)], bin_gr[subjectHits(ov)])
int_widths <- width(int_ranges)
# For convenience, put it all in a data frame
ov_df <- data.frame(
		gene_index = queryHits(ov),
		bin_index  = subjectHits(ov),
		gene_name = gtf_gr$gene[queryHits(ov)],
		bin_id     = bin_gr[subjectHits(ov)],
		overlap_bp = int_widths,
		overlap_prop = int_widths/gtf_gr$gene_length[queryHits(ov)],
		stringsAsFactors = FALSE
)

# obtain regions with >50% overlap
subset_ov = subset(ov_df,  overlap_prop >= 0.5)

subset_ov$bingrp = paste0(subset_ov$bin_id.seqnames, ':', 
			  subset_ov$bin_id.start, '-', 
			  subset_ov$bin_id.end)
write.csv(subset_ov %>% select(gene_name, bingrp), outfile, row.names=FALSE)


# Write it as numbat read-able gtf as well

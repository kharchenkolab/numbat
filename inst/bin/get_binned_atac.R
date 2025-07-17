#!/usr/bin/env Rscript

# Script: get_binned_atac.R
# Description: Generate cell x counts matrix per genomic tiles from ATAC-seq fragments
# Usage: Rscript get_binned_atac.R --CB barcodes.txt --frag fragments.tsv.gz --binGR bins.rds --outFile output.tsv

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(rtracklayer)
})

ensure_matrix <- function(x) {
  if (inherits(x, "table")) {
    x <- unclass(as.matrix(x))
  } else if (inherits(x, "array")) {
    x <- as.matrix(x)
  }
  x
}

message('Get cell x counts per supplied GRanges tiles. Saves as tsv file.')
option_list = list(
  make_option('--CB', type = "character", default = NULL,
              help = "Cell barcodes file (txt). First column must be barcodes."),
  make_option('--frag', type = "character", default = NULL,
              help = "Path to fragments.tsv[.gz] file"),
  make_option('--binGR', type = "character", default = NULL,
              help = "Binned genome GRanges file (RDS format)"),
  make_option('--outFile', type = "character", default = NULL,
              help = "Output file path (.tsv)"),
  make_option('--generateAggRef', action = "store_true", default = FALSE,
              help = "Indicator for whether we are generating aggregated reference or just a count matrix")
)

args = parse_args(OptionParser(option_list = option_list))
invisible(list2env(args,environment()))

# Main function to process data and generate bin-by-cell matrix
generate_bin_cell_matrix <- function(barcode_file, fragment_file, bin_file, output_file) {
  # Read cell barcodes
  cat("Reading cell barcodes...\n")
  barcodes <- as.list(data.table::fread(barcode_file)[[1]])
  cat(paste0("Found ", length(barcodes), " cell barcodes.\n"))
  
  # Read fragment file
  cat("Reading fragment file...\n")
  fragments <- import.bed(fragment_file, extraCols = c("type" = "character", "type" = "integer"))
  colnames(mcols(fragments)) <- c("barcode", "dup_counts")
  cat(paste0("Total ", length(fragments), " fragments.\n"))
  
  # Filter fragments by cell barcodes
  cat("Filtering fragments by cell barcodes...\n")
  fragments_in_cell <- fragments[fragments$barcode %in% barcodes]
  cat(paste0("Total ", length(fragments_in_cell), " fragments in cells.\n"))
  
  # Read genomic bins
  cat("Reading genomic bins...\n")
  query <- readRDS(bin_file)
  cat(paste0("Using ", length(query), " genomic bins.\n"))
  
  # Find overlaps between bins and fragments
  cat("Finding overlaps between bins and fragments...\n")
  ov <- findOverlaps(query, fragments_in_cell)
  ov <- as.matrix(ov)
  tmp <- fragments_in_cell$barcode[ov[,2]]
  ov <- cbind(ov, match(tmp, barcodes))
  
  # Generate bin-by-cell matrix
  cat("Generating bin-by-cell matrix...\n")
  mm <- table(ov[,1], ov[,3])
  colnames(mm) <- barcodes[as.numeric(colnames(mm))]
  bins <- paste0(as.character(seqnames(query)), ":", as.character(ranges(query)))
  rownames(mm) <- bins[as.numeric(rownames(mm))]
  mm <- as(ensure_matrix(mm), "dgCMatrix")
  if (args$generateAggRef) {
	  message("Generating aggregated reference...\n")
	   agg_ref_counts = numbat::aggregate_counts(as.matrix(mm), 
						  data.table::fread(barcode_file), 
						  normalized = TRUE, 
						  verbose = TRUE)
	  message("Saving aggregated count matrix.")
	   if(endsWith(output_file, ".rds")) {
	     saveRDS(agg_ref_counts, output_file)
	   }else if( endsWith(output_file, ".tsv")) {
	     write.table(agg_ref_counts, output_file)
	   }
	  
  } else {
	  message("Proceeding with count matrix (not aggregated)\n")
	  # Write output
	  cat("Writing output to file...\n")
	  if(endsWith(output_file, ".rds")) {
	    saveRDS(mm, output_file)
	  }else if( endsWith(output_file, ".tsv")) {
	    write.table(mm, output_file, sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
	  }
	  
	  
	  cat(paste0("Success! Bin-by-cell matrix saved to: ", output_file, "\n"))
  }
  return(mm)
}

# run function
mat <- generate_bin_cell_matrix(CB,frag,binGR,outFile)

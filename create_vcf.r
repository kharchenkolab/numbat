#!/usr/bin/env Rscript
library(dplyr)
library(Matrix)
library(data.table)
library(stringr)
library(glue)
library(vcfR)
library(argparse)
devtools::load_all('~/Numbat')

parser <- ArgumentParser(description='Create VCF from scRNA dataset for phasing')
parser$add_argument('--samples', type = "character", help = "sample names")
parser$add_argument('--label', type = "character", help = "output sample label")

args <- parser$parse_args()

samples = str_split(args$samples, ',')[[1]]
label = args$label

###### Main ######

vcfs = lapply(samples, function(sample){read.vcfR(glue('/home/tenggao/pileup/{sample}/cellSNP.base.vcf'), verbose = F)})

genotype(label, samples, vcfs)

  



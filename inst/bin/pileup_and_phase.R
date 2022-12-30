#!/usr/bin/env Rscript

library(optparse)

parser = OptionParser(description='Run SNP pileup and phasing with 1000G')
parser = add_option(parser, '--label', default = 'subject', type = "character", help = "Individual label. One per run.")
parser = add_option(parser, '--samples', default = 'sample', type = "character", help = "Sample name(s); comma delimited if multiple. All samples must belong to the same individual.")
parser = add_option(parser, '--bams', default = NULL, type = "character",  help = "BAM file(s); one per sample, comma delimited if multiple.")
parser = add_option(parser, '--barcodes', default = NULL, type = "character", help = "Cell barcode file(s); one per sample, comma delimited if multiple.")
parser = add_option(parser, '--gmap', default = NULL, type = "character", help = "Path to genetic map provided by Eagle2 (e.g. Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz)")
parser = add_option(parser, '--eagle', type = "character", default = 'eagle', help = "Path to Eagle2 binary file")
parser = add_option(parser, '--snpvcf', default = NULL, type = "character", help = "SNP VCF for pileup")
parser = add_option(parser, '--paneldir', default = NULL, type = "character", help = "Directory to phasing reference panel (BCF files)")
parser = add_option(parser, '--outdir', default = './pileup_and_phase', type = "character", help = "Output directory")
parser = add_option(parser, '--ncores', default = 1, type = "integer", help = "Number of cores")
parser = add_option(parser, '--UMItag', default = "Auto", type = "character", help = "UMI tag in bam; should be Auto for 10x and XM for Slide-seq")
parser = add_option(parser, '--cellTAG', default = "CB", type = "character", help = "Cell tag in bam; should be CB for 10x and XC for Slide-seq")
parser = add_option(parser, '--smartseq', default = FALSE, action = 'store_true', help = "Running with SMART-seq mode; Supply a txt file containing directories of BAM files to --bams and a txt file containing cell names to --barcodes (each entry on its own line for both; ordering must match).")
parser = add_option(parser, '--bulk', default = FALSE, action = 'store_true', help = "Running with bulk RNA-seq mode")
args <- parse_args(parser)

suppressPackageStartupMessages({
    library(glue)
    library(stringr)
    library(data.table)
    library(dplyr)
    library(vcfR)
    library(Matrix)
    library(numbat)
})

# required args
if (any(is.null(c(args$bams, args$barcodes, args$gmap, args$snpvcf, args$paneldir)))) {
    stop('Missing one or more required arguments: bams, barcodes, gmap, snpvcf, paneldir')
}

label = args$label
samples = str_split(args$samples, ',')[[1]]
outdir = args$outdir
bams = str_split(args$bams, ',')[[1]]
if (!is.null(args$barcodes)) {
    barcodes = str_split(args$barcodes, ',')[[1]]
}
n_samples = length(samples)
ncores = args$ncores
gmap = args$gmap
eagle = args$eagle
snpvcf = args$snpvcf
paneldir = args$paneldir
UMItag = args$UMItag
cellTAG = args$cellTAG
smartseq = args$smartseq
bulk = args$bulk
genome = ifelse(str_detect(args$gmap, 'hg19'), 'hg19', 'hg38')
message(paste0('Using genome version: ', genome))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(glue('{outdir}/pileup'), showWarnings = FALSE)
dir.create(glue('{outdir}/phasing'), showWarnings = FALSE)
for (sample in samples) {
    dir.create(glue('{outdir}/pileup/{sample}'), showWarnings = FALSE)
}

## pileup

cmds = c()

if (bulk) {

    for (i in 1:n_samples) {

        bam_file = glue('{outdir}/pileup/{sample}/bam_path.tsv')
        sample_file = glue('{outdir}/pileup/{sample}/sample.tsv')

        fwrite(list(bams[i]), bam_file)
        fwrite(list(samples[i]), sample_file)
        
        cmd = glue(
            'cellsnp-lite', 
            '-S {bam_file}',
            '-i {sample_file}',
            '-O {outdir}/pileup/{samples[i]}',
            '-R {snpvcf}', 
            '-p {ncores}',
            '--minMAF 0',
            '--minCOUNT 2',
            '--UMItag None',
            '--cellTAG None',
            .sep = ' ')

        cmds = c(cmds, cmd)

    }

} else if (smartseq) {

    cmd = glue(
            'cellsnp-lite', 
            '-S {bams}',
            '-i {barcodes}',
            '-O {outdir}/pileup/{samples}',
            '-R {snpvcf}', 
            '-p {ncores}',
            '--minMAF 0',
            '--minCOUNT 2',
            '--UMItag None',
            '--cellTAG None',
            .sep = ' ')

    cmds = c(cmd)

} else {
    
    for (i in 1:n_samples) {
        
        cmd = glue(
            'cellsnp-lite', 
            '-s {bams[i]}',
            '-b {barcodes[i]}',
            '-O {outdir}/pileup/{samples[i]}',
            '-R {snpvcf}', 
            '-p {ncores}',
            '--minMAF 0',
            '--minCOUNT 2',
            '--UMItag {UMItag}',
            '--cellTAG {cellTAG}',
            .sep = ' ')

        cmds = c(cmds, cmd)

    }

}

cat('Running pileup\n')

script = glue('{outdir}/run_pileup.sh')

list(cmds) %>% fwrite(script, sep = '\n')

system(glue('chmod +x {script}'))

tryCatch({
    system(glue('sh {script}'), intern = TRUE)
},
warning = function(w){
    stop('Pileup failed')
})

## VCF creation
cat('Creating VCFs\n')

# read in the pileup VCF
vcfs = lapply(samples, function(sample) {
    vcf_file = glue('{outdir}/pileup/{sample}/cellSNP.base.vcf')
    if (file.exists(vcf_file)) {
        if (file.size(vcf_file) != 0) {
            vcf = vcfR::read.vcfR(vcf_file, verbose = F)
            return(vcf)
        } else {
            stop('Pileup VCF is empty')
        }
    } else {
        stop('Pileup VCF not found')
    }
})

# Remove chr prefix if present
vcfs = lapply(vcfs, function(vcf){
    vcf@fix[,1] <- gsub("chr", "", vcf@fix[,1])
    return(vcf)
})

numbat:::genotype(label, samples, vcfs, glue('{outdir}/phasing'), chr_prefix = TRUE)

## phasing
cat('Running phasing\n')
eagle_cmd = function(chr) {
    paste(eagle, 
        glue('--numThreads {ncores}'), 
        glue('--vcfTarget {outdir}/phasing/{label}_chr{chr}.vcf.gz'), 
        glue('--vcfRef {paneldir}/chr{chr}.genotypes.bcf'), 
        glue('--geneticMapFile={gmap}'), 
        glue('--outPrefix {outdir}/phasing/{label}_chr{chr}.phased'),
    sep = ' ')
}

cmds = lapply(1:22, function(chr){eagle_cmd(chr)})

script = glue('{outdir}/run_phasing.sh')

list(cmds) %>% fwrite(script, sep = '\n')

system(glue('chmod +x {script}'))

tryCatch({
    system2(script, stdout = glue("{outdir}/phasing.log"))
},
warning = function(w){
    stop('Phasing failed')
})

## Generate allele count dataframe
cat('Generating allele count dataframes\n')

if (genome == 'hg19') {
    gtf = gtf_hg19
} else {
    gtf = gtf_hg38
}

genetic_map = fread(gmap) %>% 
    setNames(c('CHROM', 'POS', 'rate', 'cM')) %>%
    group_by(CHROM) %>%
    mutate(
        start = POS,
        end = c(POS[2:length(POS)], POS[length(POS)])
    ) %>%
    ungroup()

for (sample in samples) {
    
    # read in phased VCF
    vcf_phased = lapply(1:22, function(chr) {
            vcf_file = glue('{outdir}/phasing/{label}_chr{chr}.phased.vcf.gz')
            if (file.exists(vcf_file)) {
                fread(vcf_file, skip = '#CHROM') %>%
                    rename(CHROM = `#CHROM`) %>%   
                    mutate(CHROM = str_remove(CHROM, 'chr'))
            } else {
                stop('Phased VCF not found')
            }
        }) %>%
        Reduce(rbind, .) %>%
        mutate(CHROM = factor(CHROM, unique(CHROM)))

    pu_dir = glue('{outdir}/pileup/{sample}')

    # pileup VCF
    vcf_pu = fread(glue('{pu_dir}/cellSNP.base.vcf'), skip = '#CHROM') %>% 
        rename(CHROM = `#CHROM`) %>%
        mutate(CHROM = str_remove(CHROM, 'chr'))

    # count matrices
    AD = readMM(glue('{pu_dir}/cellSNP.tag.AD.mtx'))
    DP = readMM(glue('{pu_dir}/cellSNP.tag.DP.mtx'))

    cell_barcodes = fread(glue('{pu_dir}/cellSNP.samples.tsv'), header = F) %>% pull(V1)

    df = numbat:::preprocess_allele(
        sample = label,
        vcf_pu = vcf_pu,
        vcf_phased = vcf_phased,
        AD = AD,
        DP = DP,
        barcodes = cell_barcodes,
        gtf = gtf,
        gmap = genetic_map
    ) %>%
    filter(GT %in% c('1|0', '0|1'))
    
    fwrite(df, glue('{outdir}/{sample}_allele_counts.tsv.gz'), sep = '\t')
    
}

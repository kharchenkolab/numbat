library(logger, quietly = T)
library(glue, quietly = T)
library(stringr, quietly = T)
library(argparse, quietly = T)
library(data.table, quietly = T)
library(dplyr, quietly = T)
library(vcfR, quietly = T)
library(Matrix, quietly = T)
devtools::load_all('~/Numbat')

parser <- ArgumentParser(description='Run SNP pileup and phasing with 1000G')
parser$add_argument('--label', type = "character", help = "individual label")
parser$add_argument('--samples', type = "character", help = "sample names, comma delimited")
parser$add_argument('--bams', type = "character", help = "bam files, one per sample, comma delimited")
parser$add_argument('--barcodes', type = "character", help = "cell barcodes, one per sample, comma delimited")
parser$add_argument('--gmap', type = "character", help = "path to genetic map provided by Eagle2")
parser$add_argument('--snpvcf', type = "character", help = "SNP VCF for pileup")
parser$add_argument('--paneldir', type = "character", help = "directory to phasing reference panel (BCF files)")
parser$add_argument('--outdir', type = "character", help = "output directory")
parser$add_argument('--ncores', type = "integer", help = "number of cores")

args <- parser$parse_args()

label = args$label
samples = str_split(args$samples, ',')[[1]]
outdir = args$outdir
bams = str_split(args$bams, ',')[[1]]
barcodes = str_split(args$barcodes, ',')[[1]]
n_samples = length(samples)
label = args$label
ncores = args$ncores
gmap = args$gmap
snpvcf = args$snpvcf
paneldir = args$paneldir

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (sample in samples) {
    dir.create(glue('{outdir}/pileup'), showWarnings = FALSE)
    dir.create(glue('{outdir}/phasing'), showWarnings = FALSE)
    dir.create(glue('{outdir}/pileup/{sample}'), showWarnings = FALSE)
}

## pileup

cmds = c()

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
        .sep = ' ')

    cmds = c(cmds, cmd)

}

cat('Running pileup')

script = glue('{outdir}/run_pileup.sh')

list(cmds) %>% fwrite(script, sep = '\n')

system(glue('chmod +x {script}'))
system2(script, stdout = glue("{outdir}/pileup.log"))

## VCF creation
cat('Creating VCFs')
vcfs = lapply(samples, function(sample){read.vcfR(glue('{outdir}/pileup/{sample}/cellSNP.base.vcf'), verbose = F)})

genotype(label, samples, vcfs, glue('{outdir}/phasing'))

## phasing
eagle_cmd = function(chr, sample) {
    paste('eagle', 
        glue('--numThreads {ncores}'), 
        glue('--vcfTarget {outdir}/phasing/{label}_chr{chr}.vcf.gz'), 
        glue('--vcfRef {paneldir}/ALL.chr{chr}.shapeit2_integrated_v1a.GRCh38.20181129.phased.bcf'), 
        glue('--geneticMapFile={gmap}'), 
        glue('--outPrefix {outdir}/phasing/{label}_chr{chr}.phased'),
    sep = ' ')
}

cmds = c()

for (sample in samples) {
    cmds = c(cmds, lapply(1:22, function(chr){eagle_cmd(chr, sample)}))
}

script = glue('{outdir}/run_phasing.sh')

list(cmds) %>% fwrite(script, sep = '\n')

system(glue('chmod +x {script}'))
system2(script, stdout = glue("{outdir}/phasing.log"))

## Generate allele count dataframe
cat('Generating allele count dataframes')

for (sample in samples) {
    
    # read in phased VCF
    vcf_phased = lapply(1:22, function(chr) {
        fread(glue('{outdir}/phasing/{label}_chr{chr}.phased.vcf.gz')) %>%
            rename(CHROM = `#CHROM`) %>%
            mutate(CHROM = str_remove(CHROM, 'chr'))   
        }) %>% Reduce(rbind, .) %>%
        mutate(CHROM = factor(CHROM, unique(CHROM)))

    pu_dir = glue('{outdir}/pileup/{sample}')

    # pileup VCF
    vcf_pu = fread(glue('{pu_dir}/cellSNP.base.vcf')) %>% rename(CHROM = `#CHROM`)

    # count matrices
    AD = readMM(glue('{pu_dir}/cellSNP.tag.AD.mtx'))
    DP = readMM(glue('{pu_dir}/cellSNP.tag.DP.mtx'))

    cell_barcodes = fread(glue('{pu_dir}/cellSNP.samples.tsv'), header = F) %>% pull(V1)
    cell_barcodes = paste0(sample, '_', cell_barcodes)

    df = preprocess_allele(
        sample = label,
        vcf_pu = vcf_pu,
        vcf_phased = vcf_phased,
        AD = AD,
        DP = DP,
        barcodes = cell_barcodes,
        gtf_transcript = gtf_transcript
    )
    
    fwrite(df, glue('{outdir}/{sample}_allele_counts.tsv.gz'), sep = '\t')
    
}
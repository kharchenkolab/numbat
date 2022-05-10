library(glue, quietly = T)
library(stringr, quietly = T)
library(argparse, quietly = T)
library(data.table, quietly = T)
library(dplyr, quietly = T)
library(vcfR, quietly = T)
library(Matrix, quietly = T)
library(numbat, quietly = T)

parser <- ArgumentParser(description='Run SNP pileup and phasing with 1000G')
parser$add_argument('--label', type = "character", required = TRUE, help = "Individual label")
parser$add_argument('--samples', type = "character", required = TRUE, help = "Sample names, comma delimited")
parser$add_argument('--bams', type = "character", required = TRUE, help = "BAM files, one per sample, comma delimited")
parser$add_argument('--barcodes', type = "character", required = TRUE, help = "Cell barcode files, one per sample, comma delimited")
parser$add_argument('--gmap', type = "character", required = TRUE, help = "Path to genetic map provided by Eagle2 (e.g. Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz)")
parser$add_argument('--eagle', type = "character", required = FALSE, default = 'eagle', help = "Path to Eagle2 binary file")
parser$add_argument('--snpvcf', type = "character", required = TRUE, help = "SNP VCF for pileup")
parser$add_argument('--paneldir', type = "character", required = TRUE, help = "Directory to phasing reference panel (BCF files)")
parser$add_argument('--outdir', type = "character", required = TRUE, help = "Output directory")
parser$add_argument('--ncores', type = "integer", required = TRUE, help = "Number of cores")
parser$add_argument('--UMItag', default = "Auto", required = FALSE, type = "character", help = "UMI tag in bam. Should be Auto for 10x and XM for Slide-seq")
parser$add_argument('--cellTAG', default = "CB", required = FALSE, type = "character", help = "Cell tag in bam. Should be CB for 10x and XC for Slide-seq")
parser$add_argument('--smartseq', action = 'store_true', help = "running with smart-seq mode")
parser$add_argument('--bulk', action = 'store_true', help = "running with bulk mode")

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

for (sample in samples) {
    dir.create(glue('{outdir}/pileup'), showWarnings = FALSE)
    dir.create(glue('{outdir}/phasing'), showWarnings = FALSE)
    dir.create(glue('{outdir}/pileup/{sample}'), showWarnings = FALSE)
}

## pileup

cmds = c()

if (bulk) {

    for (i in 1:n_samples) {

        bam_file = glue('{outdir}/pileup/{sample}/bam_path.tsv')
        sample_file = glue('{outdir}/pileup/{sample}/sample.tsv')

        fwrite(list(bams[i]), bam_file)
        fwrite(list(barcodes[i]), sample_file)
        
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
vcfs = lapply(samples, function(sample){vcfR::read.vcfR(glue('{outdir}/pileup/{sample}/cellSNP.base.vcf'), verbose = F)})

numbat:::genotype(label, samples, vcfs, glue('{outdir}/phasing'))

## phasing
eagle_cmd = function(chr, sample) {
    paste(eagle, 
        glue('--numThreads {ncores}'), 
        glue('--vcfTarget {outdir}/phasing/{label}_chr{chr}.vcf.gz'), 
        glue('--vcfRef {paneldir}/chr{chr}.genotypes.bcf'), 
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

tryCatch({
    system2(script, stdout = glue("{outdir}/phasing.log"))
},
warning = function(w){
    stop('Phasing failed')
})

## Generate allele count dataframe
cat('Generating allele count dataframes\n')

if (genome == 'hg19') {
    gtf_transcript = gtf_hg19
} else {
    gtf_transcript = gtf_hg38
}

for (sample in samples) {
    
    # read in phased VCF
    vcf_phased = lapply(1:22, function(chr) {
            vcf_file = glue('{outdir}/phasing/{label}_chr{chr}.phased.vcf.gz')
            if (file.exists(vcf_file)) {
                fread(vcf_file) %>%
                    rename(CHROM = `#CHROM`) %>%
                    mutate(CHROM = str_remove(CHROM, 'chr'))   
            }
        }) %>%
        Reduce(rbind, .) %>%
        mutate(CHROM = factor(CHROM, unique(CHROM)))

    pu_dir = glue('{outdir}/pileup/{sample}')

    # pileup VCF
    vcf_pu = fread(glue('{pu_dir}/cellSNP.base.vcf')) %>% rename(CHROM = `#CHROM`)

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
        gtf_transcript = gtf_transcript
    )
    
    fwrite(df, glue('{outdir}/{sample}_allele_counts.tsv.gz'), sep = '\t')
    
}

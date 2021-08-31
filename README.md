# Armadillo
<img src="armadillo.png" width="500">
Haplotype-aware detection of copy number variations in scRNA-Seq

1. SNP pileup
```
cmd = glue(
      'cellsnp-lite', 
      '-s {bam}.bam',
      '-b {barcodes}.tsv.gz',
      '-O /home/tenggao/pileup/{sample}',
      '-R /home/tenggao/ref/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf', 
      '-p 30',
      '--minMAF 0',
      '--minCOUNT 2',
      .sep = ' ')
```

2. Create VCF
```
Rscript /home/tenggao/Armadillo/create_vcf.r --sample sample1,sample2 --label patient 
```

3. Phasing
```
eagle_cmd = function(chr, sample) {
    paste('eagle', 
        '--numThreads 20', 
        glue('--vcfTarget /home/tenggao/phasing/{sample}_chr{chr}.vcf.gz'), 
        glue('--vcfRef /home/tenggao/ref/ALL.chr{chr}.shapeit2_integrated_v1a.GRCh38.20181129.phased.bcf'), 
        '--geneticMapFile=/home/tenggao/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz', 
        glue('--outPrefix /home/tenggao/phasing/{sample}_chr{chr}_phased'),
    sep = ' ')
}

cmds = c()

for (sample in samples) {
    cmds = c(cmds, lapply(1:22, function(chr){eagle_cmd(chr, sample)}))
}

list(cmds) %>% fwrite('~/external/WASHU/run_phasing.sh', sep = '\n')
```

4. Run Armadillo
```
source('~/Armadillo/armadillo.r')

# gtf
gtf = fread('~/ref/hg38.refGene.gtf')
cols = c('CHROM', 'source', 'biotype', 'start', 'end', 'a', 'strand', 'b', 'info')
colnames(gtf) = cols

gtf = gtf %>% filter(CHROM %in% paste0('chr', 1:22))

# transcript GTF
gtf_transcript = gtf %>% filter(biotype == 'transcript') %>% distinct(start, end, `.keep_all` = TRUE) %>%
    mutate(gene = str_extract(info, '(?<=gene_name \\").*(?=\\";)')) %>%
    mutate(CHROM = as.integer(str_remove(CHROM, 'chr'))) %>%
    select(gene, start, end, CHROM) %>%
    mutate(region = paste0('chr', CHROM, ':', start, '-', end)) %>%
    mutate(gene_length = end-start) %>%
    arrange(CHROM, gene, -gene_length) %>%
    distinct(gene, .keep_all = TRUE) %>%
    rename(gene_start = start, gene_end = end) %>%
    arrange(CHROM, gene_start)
 
# read in phased VCF
vcf_phased = lapply(1:22, function(chr) {
 fread(glue('/home/tenggao/phasing/{sample}_chr{chr}_phased.vcf.gz')) %>%
    rename(CHROM = `#CHROM`) %>%
    mutate(CHROM = str_remove(CHROM, 'chr'))   
}) %>% Reduce(rbind, .) %>%
mutate(CHROM = factor(CHROM, unique(CHROM)))

# pileup VCF
vcf_pu = fread(glue('/home/tenggao/pileup/{sample}/cellSNP.base.vcf')) %>% rename(CHROM = `#CHROM`)

# count matrices
AD = readMM(glue('/home/tenggao/pileup/{sample}/cellSNP.tag.AD.mtx'))
DP = readMM(glue('/home/tenggao/pileup/{sample}/cellSNP.tag.DP.mtx'))

count_mat = as.matrix(t(con$samples[[sample]]$misc$rawCounts))

# cell annotations
cell_barcodes = fread(glue('/home/tenggao/pileup/{sample}/cellSNP.samples.tsv'), header = F) %>% pull(V1)

# prepare allele count dataframe
df = preprocess_data(
    sample = sample,
    vcf_pu = vcf_pu,
    vcf_phased = vcf_phased,
    AD = AD,
    DP = DP,
    barcodes = cell_barcodes,
    gtf_transcript = gtf_transcript
)$df_obs

 # run
 out = armadillo_subclone(
    count_mat_obs,
    lambdas_ref,
    df,
    gtf_transcript,
    t = 1e-8,
    sample_size = 500,
    out_dir = glue('~/results/{sample}')
)
```


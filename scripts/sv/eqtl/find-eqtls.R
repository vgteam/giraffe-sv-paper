library(dplyr)
library(rtracklayer)
library(MatrixEQTL)
library(GenomicRanges)

##
## sample information
##

## PCs
pcs.df = read.table('1kgp.4pcs.tsv', as.is=TRUE)

## pedigree with population/gender information
if(!file.exists('20130606_g1k.ped')){
  download.file('ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped', '20130606_g1k.ped')
}

ped = read.table('20130606_g1k.ped', as.is=TRUE, header=TRUE, sep='\t')
info.df = ped %>% mutate(sample=Individual.ID) %>% select(sample, Gender, Population)

## save YRI samples for EUR/YRI analysis later
yri.samples = info.df %>% filter(Population=='YRI') %>% .$sample

## format into a matrix
rownames(info.df) = info.df$sample
info.df = cbind(info.df[rownames(pcs.df),], pcs.df)
info.df$sample = NULL
info.df$Gender = info.df$Gender - 1
info.df$Population = NULL
info.mat = as.matrix(t(info.df))

##
## Gene expression
##

## Gene annotation
if(!file.exists('gencode.v35.annotation.gtf.gz')){
  download.file('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz', 'gencode.v35.annotation.gtf.gz')
}
genc = import('gencode.v35.annotation.gtf.gz')

## Gene expression data
if(!file.exists('GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz')){
  download.file('https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz', 'GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz')
}
ge = read.table('GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz', as.is=TRUE, header=TRUE)

## remove .* from gene id names
genc.g = subset(genc, type=='gene')
genc.g$gene_id_short = gsub('(.+)\\..*', '\\1', genc.g$gene_id)
ge$gene_id_short = gsub('(.+)\\..*', '\\1', ge$Gene_Symbol)

## genes
mean(ge$Gene_Symbol %in% genc.g$gene_id)
mean(ge$gene_id_short %in% genc.g$gene_id_short)
com.genes = intersect(ge$gene_id_short, genc.g$gene_id_short)

## median gene expression used for the enrichment in highly expressed genes
gene.median.exp = apply(ge, 1, median)
tibble(gene=names(gene.median.exp), median.ge=gene.median.exp) %>%
  write.table(file='geuvadis.median.ge.tsv', sep='\t', row.names=FALSE, quote=FALSE)

##
## SV genotypes
##

ac.mat = read.table('2504kgp.svsite80al.ac.tsv.gz', as.is=TRUE)
ac.mat = as.matrix(ac.mat)

svs.df = read.table('svs.2504kgp.svsite80al.tsv.gz', as.is=TRUE, header=TRUE)

## subset to common samples
sum(colnames(ge) %in% colnames(ac.mat))
com.samples = intersect(colnames(ge), colnames(ac.mat))
rownames(ge) = ge$gene_id_short
ge = as.matrix(ge[, com.samples])
ac.mat = ac.mat[, com.samples]
info.mat = info.mat[, com.samples]

write.table(rownames(ge), file='eqtl-genes.txt', row.names=FALSE, col.names=FALSE, quote=FALSE)

##
## Matrix-eQTL
##

## Variant and gene positions
varspos = tibble(snpid=svs.df$svid, chr=svs.df$seqnames, pos=round((svs.df$start+svs.df$end)/2)) %>%
  filter(chr %in% paste0('chr', 1:22)) %>% 
  arrange(chr, pos) %>% as.data.frame
genepos = tibble(geneid=genc.g$gene_id_short, chr=as.character(seqnames(genc.g)),
                 s1=start(genc.g), s2=end(genc.g)) %>%
  arrange(chr, s1, s2) %>% as.data.frame

save(varspos, genepos, ac.mat, info.mat, ge, file='eqtl-test.RData')
## load('eqtl-test.RData')

## quantile normalization
quantnorm <- function(x){
  rks = rank(x, ties.method = "average")
  qnorm(rks / (length(rks)+1))
}

## proportion of samples with minor alleles
nonmajor <- function(ac){
  ac.mode = unique(ac)
  ac.mode = ac.mode[which.max(tabulate(match(ac, ac.mode)))]
  mean(ac-ac.mode>0)
}

## run eQTL discovery with MatrixEQTL
run_matrixeqtl <- function(ge.norm='none', samples=NULL, min.nonref.prop=0.05){
  ## allele counts, gene expression and covariates
  vars = SlicedData$new();
  vars$fileSliceSize = 5000
  ac = ac.mat[varspos$snpid,]
  if(!is.null(samples)){
    ac = ac[,samples]
  }
  ## remove SVs with not enough allele count variance
  nonmaj.prop = apply(ac, 1, nonmajor)
  ## at least 2 individuals with non-ref genotypes to mitigate effect by single outliers
  min.nonref.prop = ifelse(min.nonref.prop<2/ncol(ac), 2/ncol(ac), min.nonref.prop)
  ac = ac[which(nonmaj.prop>=min.nonref.prop),]
  vars$CreateFromMatrix(ac)
  vars$ResliceCombined()
  gene = SlicedData$new();
  gene$fileSliceSize = 5000
  if(!is.null(samples)){
    ge = ge[,samples]
  }
  if(ge.norm == 'quantile'){
    message('Quantile normalization...')
    ge = t(apply(ge, 1, quantnorm))
  }
  gene$CreateFromMatrix(ge)
  gene$ResliceCombined()
  if(ge.norm == 'normal'){
    message('Standard normalization...')
    gene$RowStandardizeCentered()
  }
  cvrt = SlicedData$new();
  cvrt$fileSliceSize = 5000
  if(!is.null(samples)){
    info.mat = info.mat[,samples]
  }
  cvrt$CreateFromMatrix(info.mat)
  cvrt$ResliceCombined()
  ## pvalue breaks to make both the histogram and the qqplot
  pvs.bks = sort(c(0, 10^(-200:-3), 5*10^(-200:-3), seq(.01, 1, .01)))
  ## run MatrixEQTL
  me.o = Matrix_eQTL_main(
    vars, gene, cvrt, 
    useModel=modelLINEAR, errorCovariance=numeric(), 
    output_file_name="", pvOutputThreshold=0,
    output_file_name.cis="", pvOutputThreshold.cis=1e-2,
    snpspos=varspos, genepos=genepos, cisDist=1e6,
    verbose=TRUE, pvalue.hist=pvs.bks, min.pv.by.genesnp=FALSE, noFDRsaveMemory=FALSE)
  message(sum(me.o$cis$eqtls$FDR<.01), ' cis-eQTLs at FDR<0.01')
  message(length(unique(subset(me.o$cis$eqtls, FDR<.01)$gene)), ' genes with cis-eQTLs at FDR<0.01')
  me.o$input = list(samples=ncol(ac), variants=nrow(ac))
  return(me.o)
}

me.lin.norm.all = run_matrixeqtl(ge.norm='normal', min.nonref.prop=.01)
me.lin.norm.eur = run_matrixeqtl(ge.norm='normal', samples=setdiff(colnames(ac.mat), yri.samples), min.nonref.prop=.01)
me.lin.norm.yri = run_matrixeqtl(ge.norm='normal', samples=intersect(colnames(ac.mat), yri.samples), min.nonref.prop=.01)
me.lin.quant.all = run_matrixeqtl(ge.norm='quantile', min.nonref.prop=.01)
me.lin.quant.eur = run_matrixeqtl(ge.norm='quantile', samples=setdiff(colnames(ac.mat), yri.samples), min.nonref.prop=.01)
me.lin.quant.yri = run_matrixeqtl(ge.norm='quantile', samples=intersect(colnames(ac.mat), yri.samples), min.nonref.prop=.01)

ll = list(me.lin.norm.all, me.lin.norm.eur, me.lin.norm.yri,
          me.lin.quant.all, me.lin.quant.eur, me.lin.quant.yri)
names(ll) = c('me.lin.norm.all', 'me.lin.norm.eur', 'me.lin.norm.yri',
              'me.lin.quant.all', 'me.lin.quant.eur', 'me.lin.quant.yri')

save(ll, file='eqtl-test-results-maf01.RData')

##
## Examples
##

load('eqtl-test-results-maf01.RData')

eqtl.df = rbind(
  ll$me.lin.norm.all$cis$eqtls %>% mutate(pop='all'),
  ll$me.lin.norm.eur$cis$eqtls %>% mutate(pop='eur'),
  ll$me.lin.norm.yri$cis$eqtls %>% mutate(pop='yri')) %>%
  mutate(snps=as.character(snps), gene=as.character(gene))

set.seed(123)

## ex: eqtl in eur+yri
ex.all = eqtl.df %>% filter(FDR<=.01, pop=='all') %>%
  sample_n(100)

## ex: eqtl in eur but not yri or eur+yri
ex.eur = eqtl.df %>% group_by(snps, gene) %>%
  filter(n()==1, pop=='eur', FDR<=.01)

## ex: eqtl in yri but not eur or eur+yri
ex.yri = eqtl.df %>% group_by(snps, gene) %>%
  filter(n()==1, pop=='yri', FDR<=.01)

ge.ex = ge[unique(c(ex.all$gene, ex.eur$gene, ex.yri$gene)),]
ac.ex = ac.mat[unique(c(ex.all$snps, ex.eur$snps, ex.yri$snps)),]

save(ge.ex, ac.ex, yri.samples, ex.all, ex.eur, ex.yri, file='eqtl-examples.RData')

library(dplyr)
library(rtracklayer)
library(MatrixEQTL)
library(GenomicRanges)

##
## sample information
##

## pedigree with population/gender information
if(!file.exists('20130606_g1k.ped')){
  download.file('ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped', '20130606_g1k.ped')
}

ped = read.table('20130606_g1k.ped', as.is=TRUE, header=TRUE, sep='\t')
info.df = ped %>% mutate(sample=Individual.ID) %>% select(sample, Gender, Population)

## save YRI samples for EUR/YRI analysis later
yri.samples = info.df %>% filter(Population=='YRI') %>% .$sample

## load information previously prepared for the SV-only analysis
## e.g. sample information, gene expression
load('eqtl-test.RData')

## proportion of samples with minor alleles
nonmajor <- function(ac){
  ac.mode = unique(ac)
  ac.mode = ac.mode[which.max(tabulate(match(ac, ac.mode)))]
  mean(ac-ac.mode>0)
}

## quantile normalization
quantnorm <- function(x){
  rks = rank(x, ties.method = "average")
  qnorm(rks / (length(rks)+1))
}

## run eQTL discovery with MatrixEQTL
## (takes up to a couple of hours)
run_matrixeqtl_all <- function(pop=c('all', 'eur', 'yri'), ge.norm='none', min.nonref.prop=0.05, use.cache=TRUE){
  if(pop[1] == 'eur'){
    ac.mat = ac.mat[, setdiff(colnames(ac.mat), yri.samples)]
  }
  if(pop[1] == 'yri'){
    ac.mat = ac.mat[, intersect(colnames(ac.mat), yri.samples)]
  }
  ## prepare positions and frequencies for SVs
  svs.varspos = tibble(snpid=rownames(ac.mat)) %>% 
    mutate(prop.nonref=apply(ac.mat, 1, nonmajor))
  svs.varspos = merge(varspos, svs.varspos)
  ## prepare genotype file with SNVs, indels and SVs
  ac.file = paste0('ac.012.', pop[1], '-maf', min.nonref.prop, '.tsv')
  if(!file.exists(paste0(ac.file, '.gz')) | !use.cache){
    file.remove(c(ac.file, paste0(ac.file, '.gz')))
    all.pos = lapply(1:22, function(chrn){
      message('Reading SNVs/indels in chr', chrn)
      df = read.table(paste0('1000gp-nygc-variants/CCDG_14151_B01_GRM_WGS_2020-08-05_chr', chrn, '.filtered.shapeit2-duohmm-phased.012.tsv.gz'), as.is=TRUE, header=TRUE)
      ## positions and frequencies in EUR+YRI, EUR, YRI
      all.varspos = tibble(snpid=df$id) %>%
        mutate(chr=paste0('chr', gsub('(.*):.*:.*:.*', '\\1', snpid)), 
               pos=as.integer(gsub('.*:(.*):.*:.*', '\\1', snpid)),
               prop.nonref=rowSums(df[,colnames(ac.mat)]>0)/ncol(ac.mat),
               prop.nonref2=rowSums(df[,colnames(ac.mat)]<2)/ncol(ac.mat),
               prop.nonref=ifelse(prop.nonref<prop.nonref2, prop.nonref, prop.nonref2)) %>%
        select(-prop.nonref2)
      ## format data.frame and add SV genotypes
      newdf = df[,c('id', colnames(ac.mat))]
      sv.iis = svs.varspos %>% filter(chr==paste0('chr', chrn)) %>% .$snpid
      newdf = rbind(newdf, cbind(tibble(id=sv.iis), ac.mat[sv.iis,]))
      ## position of the all variants
      min.nonref = ifelse(min.nonref.prop<2/ncol(ac.mat), 2/ncol(ac.mat), min.nonref.prop)
      all.varspos = svs.varspos %>% filter(chr==paste0('chr', chrn)) %>% rbind(all.varspos) %>%
        filter(prop.nonref >= min.nonref) %>% arrange(chr, pos)  
      ## reorder genotype matrix by position
      ord.v = 1:nrow(newdf)
      names(ord.v) = newdf$id
      newdf = newdf[ord.v[all.varspos$snpid], ]
      ## write genotypes
      write.table(newdf, file=ac.file, sep='\t', quote=FALSE, row.names=FALSE,
                  col.names=!file.exists(ac.file), append=file.exists(ac.file))
      ## return variant positions
      return(all.varspos)
    })
    all.pos = do.call(rbind, all.pos)
    save(all.pos, file=paste0('eqtl-test-allvars-', pop[1], '-maf', min.nonref.prop, '.RData'))
    system2('gzip', ac.file)
  } else {
    load(paste0('eqtl-test-allvars-', pop[1], '-maf', min.nonref.prop, '.RData'))
  }
  ## load variants
  vars = SlicedData$new();
  vars$fileDelimiter = "\t"
  vars$fileOmitCharacters = "NA"
  vars$fileSkipRows = 1
  vars$fileSkipColumns = 1
  vars$fileSliceSize = 20000
  vars$LoadFile(paste0(ac.file, '.gz'))
  all.pos = all.pos[,c('snpid', 'chr', 'pos')]
  ## gene expression
  gene = SlicedData$new();
  gene$fileSliceSize = 5000
  ge = ge[, colnames(ac.mat)]
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
  ## covariates
  cvrt = SlicedData$new();
  cvrt$fileSliceSize = 5000
  info.mat = info.mat[,colnames(ac.mat)]
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
    snpspos=all.pos, genepos=genepos, cisDist=1e6,
    verbose=TRUE, pvalue.hist=pvs.bks, min.pv.by.genesnp=FALSE, noFDRsaveMemory=FALSE)
  message(sum(me.o$cis$eqtls$FDR<.01), ' cis-eQTLs at FDR<0.01')
  message(length(unique(subset(me.o$cis$eqtls, FDR<.01)$gene)), ' genes with cis-eQTLs at FDR<0.01')
  me.o$input = list(samples=ncol(ac.mat), variants=nrow(all.pos), svs=sum(grepl('sv', all.pos$snpid)))
  return(me.o)
}

me.lin.norm.all = run_matrixeqtl_all(pop='all', ge.norm='normal', min.nonref.prop=0.01)
me.lin.norm.eur = run_matrixeqtl_all(pop='eur', ge.norm='normal', min.nonref.prop=0.01)
me.lin.norm.yri = run_matrixeqtl_all(pop='yri', ge.norm='normal', min.nonref.prop=0.01)

ll = list(me.lin.norm.all, me.lin.norm.eur, me.lin.norm.yri)
names(ll) = c('me.lin.norm.all', 'me.lin.norm.eur', 'me.lin.norm.yri')

## remove a few columns to reduce file size
for(coln in c('statistic')){
  for(nn in names(ll)){
    ll[[nn]]$cis$eqtls[,coln] = NULL
  }
}

## remove FDR>5% to reduce file size
for(nn in names(ll)){
  ll[[nn]]$cis$eqtls = subset(ll[[nn]]$cis$eqtls, FDR<=.05)
}

save(ll, file='eqtl-test-allvars-results-maf01.RData')

##
## Examples
##

## load information previously prepared for the SV-only analysis
## e.g. sample information, gene expression
load('eqtl-test.RData')

load('eqtl-test-allvars-results-maf01.RData')

ex.lead = ll$me.lin.norm.all$cis$eqtls %>% filter(FDR<=.01) %>%
  mutate(type=ifelse(grepl('sv', snps), 'SV', 'SNV-indel'),
         snps=as.character(snps), gene=as.character(gene)) %>%
  group_by(gene) %>%
  arrange(FDR) %>% do(head(., 1)) %>% filter(type == 'SV')

## ex: eqtl leads
ge.ex = ge[unique(c(ex.lead$gene)),]
ac.ex = ac.mat[unique(c(ex.lead$snps)),]

save(ge.ex, ac.ex, yri.samples, ex.lead, file='eqtl-lead-examples.RData')

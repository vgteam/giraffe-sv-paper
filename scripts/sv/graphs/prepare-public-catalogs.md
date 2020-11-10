Prepare public SV catalogs
================

``` r
library(dplyr)
library(sveval)
library(GenomicRanges)
```

## 1000 Genomes Project phase 3

``` r
if(!file.exists('ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz')){
  download.file('http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/supporting/GRCh38_positions/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz', 'ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz')
}

kgp3 = readSVvcf.multisamps('ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz')
kgp3$ac = 1

## DEL_LINE1 etc -> DEL
kgp3$type = ifelse(grepl('DEL', kgp3$type), 'DEL', kgp3$type)
## ALU etc -> INS
kgp3$type = ifelse(kgp3$type %in% c('ALU', 'LINE1', 'SVA'), 'INS', kgp3$type)
## CNV -> DEL or DUP
kgp3$type = ifelse(kgp3$type=='CNV' & kgp3$alt %in% c('<CN0>','<CN1>'), 'DEL', kgp3$type)
kgp3$type = ifelse(kgp3$type=='CNV', 'DUP', kgp3$type)
kgp3 = subset(kgp3, alt != '<CN2>')

## rename chromosome names
seqlevels(kgp3) = paste0('chr', seqlevels(kgp3))

## convert DUP to INS?
kgp3$type = ifelse(kgp3$type=='DUP', 'INS', kgp3$type)
```

## gnomAD-SV

``` r
if(!file.exists('gnomad_v2_sv.sites.pass.lifted.vcf.gz')){
  download.file('https://storage.googleapis.com/jmonlong-vg-wdl-dev-test/gnomad_v2_sv.sites.pass.lifted.vcf.gz', 'gnomad_v2_sv.sites.pass.lifted.vcf.gz')
}

gnomad = readSVvcf('gnomad_v2_sv.sites.pass.lifted.vcf.gz', other.field='AF')
gnomad$ac = 1
gnomad$af = gnomad$AF
gnomad$AF = NULL

## convert DUP to INS?
gnomad$type = ifelse(gnomad$type=='DUP', 'INS', gnomad$type)
```

## SVPOP

``` r
if(!file.exists('EEE_SV-Pop_1.ALL.genotypes.20181204.vcf.gz')){
  download.file('http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181025_EEE_SV-Pop_1/VariantCalls_EEE_SV-Pop_1/EEE_SV-Pop_1.ALL.genotypes.20181204.vcf.gz', 'EEE_SV-Pop_1.ALL.genotypes.20181204.vcf.gz')
  ## this file was gzipped twice apparently so I gunzipped it once and renamed to '.vcf.gz'
}

## svpop.af = readSVvcf('EEE_SV-Pop_1.ALL.genotypes.20181204.vcf.gz', sample.name=NULL, other.field='AF')
svpop = readSVvcf.multisamps('EEE_SV-Pop_1.ALL.genotypes.20181204.vcf.gz')
svpop$ac = 1
```

## Save R object with the formatted catalogs

``` r
save(kgp3, gnomad, svpop, file='public-sv-catalogs.RData')
```

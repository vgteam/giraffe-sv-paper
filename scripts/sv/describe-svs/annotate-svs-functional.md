Annotate SVs with functional information
================

``` r
library(dplyr)
library(sveval)
library(GenomicRanges)
library(ggplot2)
library(gridExtra)
library(knitr)
library(rtracklayer)
## list of graphs
ggp = list()
```

## Gene annotation

``` r
if(!file.exists('gencode.v35.annotation.gtf.gz')){
  download.file('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz', 'gencode.v35.annotation.gtf.gz')
}

genc = import('gencode.v35.annotation.gtf.gz')
```

## SVs genotyped in 2,504 samples from 1000 Genomes Project

``` r
## SVs grouped by site ('svsite' and 'clique' columns)
kgp = read.table('svs.2504kgp.svsite80al.tsv.gz', as.is=TRUE, header=TRUE)

## stats for each SV locus
## use the most frequent allele (and then the largest) for ac/af/size
## also saves sum/max/min across all alleles
kgp.s = kgp %>% arrange(desc(af), desc(size)) %>%
  group_by(seqnames, svsite, type, clique) %>%
  summarize(start=start[1], end=end[1], ac=ac[1], af=af[1], size=size[1], .groups='drop') %>%
  filter(size>=50) %>% makeGRangesFromDataFrame(keep.extra.columns=TRUE)
```

## Coding SVs

``` r
genc.pc.cds = subset(genc, type=='CDS' & gene_type=='protein_coding')
ol.cds = findOverlaps(kgp.s, genc.pc.cds) %>% as.data.frame %>%
  mutate(gene=genc.pc.cds$gene_name[subjectHits]) %>%
  group_by(queryHits) %>% summarize(gene=paste(unique(sort(gene)), collapse=';'))

kgp.cds = kgp.s[ol.cds$queryHits]
kgp.cds$gene = ol.cds$gene

length(kgp.cds)
```

    ## [1] 1599

## Example: population-specific and novel

``` r
## population specific information
freq.all = read.table('2504kgp.svsite80al.superpopfreq.tsv.gz', as.is=TRUE, header=TRUE)
pop.spec = read.table('pops-freq-1kgp-med1.tsv', as.is=TRUE, header=TRUE)

## novel annotation
novel = read.table('svsite.2504kgp.novel.tsv', as.is=TRUE, header=TRUE)

ex.df = pop.spec %>% filter(svsite %in% kgp.cds$svsite,
                            svsite %in% novel$svsite, 
                            af.med<.05) %>%
  arrange(desc(abs(af-af.med))) 
```

``` r
tmp = lapply(1:min(10, nrow(ex.df)), function(ii){
  subset(kgp.cds, svsite==ex.df$svsite[ii]) %>% as.data.frame %>%
    mutate(coord=paste0('[', seqnames, ':', start, '-', end,
                        '](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=',
                        seqnames, '%3A', start, '%2D', end, ')')) %>% 
    select(coord, svsite, type, size, af, gene) %>% 
    kable %>% cat(sep='\n')
  cat('\n\n')
  freq.all %>% filter(svsite==ex.df$svsite[ii]) %>% arrange(desc(af)) %>% kable %>% cat(sep='\n')
  cat('\n\n')
})
```

| coord                                                                                                              | svsite         | type | size |        af | gene   |
| :----------------------------------------------------------------------------------------------------------------- | :------------- | :--- | ---: | --------: | :----- |
| [chr6:166585488-166586678](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr6%3A166585488%2D166586678) | sv\_1468829\_0 | DEL  | 1190 | 0.1291933 | RAMACL |

| svsite         | Superpopulation |        af |
| :------------- | :-------------- | --------: |
| sv\_1468829\_0 | AFR             | 0.4659607 |
| sv\_1468829\_0 | AMR             | 0.0403458 |
| sv\_1468829\_0 | SAS             | 0.0020450 |
| sv\_1468829\_0 | EUR             | 0.0009940 |
| sv\_1468829\_0 | EAS             | 0.0000000 |

| coord                                                                                                                | svsite        | type | size |        af | gene   |
| :------------------------------------------------------------------------------------------------------------------- | :------------ | :--- | ---: | --------: | :----- |
| [chr12:132262607-132262607](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr12%3A132262607%2D132262607) | sv\_873684\_0 | INS  |   77 | 0.0924521 | GALNT9 |

| svsite        | Superpopulation |        af |
| :------------ | :-------------- | --------: |
| sv\_873684\_0 | AFR             | 0.2919818 |
| sv\_873684\_0 | AMR             | 0.0288184 |
| sv\_873684\_0 | EUR             | 0.0268390 |
| sv\_873684\_0 | EAS             | 0.0168651 |
| sv\_873684\_0 | SAS             | 0.0132924 |

| coord                                                                                                              | svsite         | type | size |        af | gene   |
| :----------------------------------------------------------------------------------------------------------------- | :------------- | :--- | ---: | --------: | :----- |
| [chr1:145305072-145305132](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr1%3A145305072%2D145305132) | sv\_2069116\_0 | DEL  |   60 | 0.0547125 | NBPF20 |

| svsite         | Superpopulation |        af |
| :------------- | :-------------- | --------: |
| sv\_2069116\_0 | AMR             | 0.1801153 |
| sv\_2069116\_0 | AFR             | 0.0378215 |
| sv\_2069116\_0 | EUR             | 0.0377734 |
| sv\_2069116\_0 | SAS             | 0.0327198 |
| sv\_2069116\_0 | EAS             | 0.0287698 |

| coord                                                                                                            | svsite        | type | size |       af | gene   |
| :--------------------------------------------------------------------------------------------------------------- | :------------ | :--- | ---: | -------: | :----- |
| [chr16:15363751-15363819](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr16%3A15363751%2D15363819) | sv\_573421\_0 | DEL  |   68 | 0.057508 | NPIPA5 |

| svsite        | Superpopulation |        af |
| :------------ | :-------------- | --------: |
| sv\_573421\_0 | EAS             | 0.1706349 |
| sv\_573421\_0 | AMR             | 0.0763689 |
| sv\_573421\_0 | SAS             | 0.0296524 |
| sv\_573421\_0 | EUR             | 0.0228628 |
| sv\_573421\_0 | AFR             | 0.0083207 |

| coord                                                                                                        | svsite        | type | size |        af | gene |
| :----------------------------------------------------------------------------------------------------------- | :------------ | :--- | ---: | --------: | :--- |
| [chr11:1017241-1017241](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr11%3A1017241%2D1017241) | sv\_895682\_0 | INS  | 1011 | 0.0427316 | MUC6 |

| svsite        | Superpopulation |        af |
| :------------ | :-------------- | --------: |
| sv\_895682\_0 | AFR             | 0.1399395 |
| sv\_895682\_0 | AMR             | 0.0230548 |
| sv\_895682\_0 | SAS             | 0.0061350 |
| sv\_895682\_0 | EAS             | 0.0039683 |
| sv\_895682\_0 | EUR             | 0.0029821 |

| coord                                                                                                              | svsite         | type | size |       af | gene |
| :----------------------------------------------------------------------------------------------------------------- | :------------- | :--- | ---: | -------: | :--- |
| [chr3:195780789-195780885](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr3%3A195780789%2D195780885) | sv\_1820058\_0 | DEL  |   96 | 0.048722 | MUC4 |

| svsite         | Superpopulation |        af |
| :------------- | :-------------- | --------: |
| sv\_1820058\_0 | AFR             | 0.1686838 |
| sv\_1820058\_0 | EUR             | 0.0526839 |
| sv\_1820058\_0 | EAS             | 0.0357143 |
| sv\_1820058\_0 | AMR             | 0.0317003 |
| sv\_1820058\_0 | SAS             | 0.0306748 |

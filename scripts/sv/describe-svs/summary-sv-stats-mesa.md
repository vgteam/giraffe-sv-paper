Summary stats for SVs in the MESA cohort
================

  - [Read population stats for each SV
    allele](#read-population-stats-for-each-sv-allele)
  - [Allele/site numbers](#allelesite-numbers)
  - [Size](#size)
  - [Overlap with simple repeats, satellites or low-complexity
    regions](#overlap-with-simple-repeats-satellites-or-low-complexity-regions)
  - [Gene annotation](#gene-annotation)
  - [Allele frequency](#allele-frequency)
  - [Alleles per SV sites](#alleles-per-sv-sites)
  - [Multi-panel figure](#multi-panel-figure)
  - [Save TSV with SV site
    information](#save-tsv-with-sv-site-information)

``` r
library(dplyr)
library(ggplot2)
library(gridExtra)
library(knitr)
library(GenomicRanges)
library(rtracklayer)
winsor <- function(x, u){
  if(any(x>u)) x[x>u] = u
  x
}
## list of graphs
ggp = list()
```

## Read population stats for each SV allele

``` r
## SVs grouped by site ('svsite' and 'clique' columns)
svs = read.table('svs.mesa2k.svsite80al.tsv.gz', as.is=TRUE, header=TRUE)

## stats for each SV locus
## use the most frequent allele (and then the largest) for ac/af/size
## also saves sum/max/min across all alleles
locs = svs %>% arrange(desc(af), desc(size)) %>%
  group_by(seqnames, svsite, type, clique) %>%
  summarize(start=start[1], end=end[1],
            ac.tot=sum(ac), ac=ac[1],
            af.tot=sum(af), af.top2=tail(head(af, 2), 1), af=af[1],
            af.top.fc=ifelse(af.top2==0, 10, af/af.top2),
            loc.n=n(),
            size.min=min(size), size.max=max(size), size=size[1],
            .groups='drop') %>%
  filter(size>=50)

set.seed(123)
sample_n(locs, 10) %>% as.data.frame
```

    ##    seqnames       svsite type clique     start       end ac.tot   ac  af.tot
    ## 1      chr6 sv_1411495_0  DEL   TRUE 170246886 170247024      7    7 0.00175
    ## 2      chr5 sv_1528143_0  DEL   TRUE 179313456 179313710     74   74 0.01850
    ## 3      chrX   sv_18563_0  DEL   TRUE    551424    551484    949  949 0.23725
    ## 4      chr3 sv_1712852_0  DEL  FALSE 195502445 195502633     36   26 0.00900
    ## 5      chr5 sv_1532605_0  DEL   TRUE 181385141 181385295   2134 2134 0.53350
    ## 6     chr16  sv_506558_0  INS   TRUE    961040    961040     23   23 0.00575
    ## 7     chr18  sv_378775_0  INS   TRUE  79765462  79765462    118  118 0.02950
    ## 8     chr16  sv_529988_0  DEL   TRUE  14689406  14695486      1    1 0.00025
    ## 9      chr7 sv_1234346_0  DEL   TRUE 107422090 107422154      1    1 0.00025
    ## 10    chr17  sv_466009_0  INS  FALSE  80744271  80744271    204   46 0.05100
    ##    af.top2      af af.top.fc loc.n size.min size.max size
    ## 1  0.00175 0.00175      1.00     1      138      138  138
    ## 2  0.01850 0.01850      1.00     1      254      254  254
    ## 3  0.23725 0.23725      1.00     1       60       60   60
    ## 4  0.00125 0.00650      5.20     4      152      190  188
    ## 5  0.53350 0.53350      1.00     1      154      154  154
    ## 6  0.00575 0.00575      1.00     1      263      263  263
    ## 7  0.02950 0.02950      1.00     1       86       86   86
    ## 8  0.00025 0.00025      1.00     1     6080     6080 6080
    ## 9  0.00025 0.00025      1.00     1       64       64   64
    ## 10 0.00625 0.01150      1.84    47       67      162  121

## Allele/site numbers

``` r
## numbers by type
rbind(locs %>% mutate(type='all') %>% group_by(type) %>% summarize(alleles=sum(loc.n), sites=n()),
      locs %>% group_by(type) %>% summarize(alleles=sum(loc.n), sites=n())) %>%
  mutate(prop.alleles=alleles/alleles[1], prop.sites=sites/sites[1]) %>% 
  kable(digits=3, format.args=list(big.mark=','))
```

| type |   alleles |   sites | prop.alleles | prop.sites |
| :--- | --------: | ------: | -----------: | ---------: |
| all  | 1,662,301 | 166,959 |        1.000 |      1.000 |
| DEL  |   180,814 |  75,520 |        0.109 |      0.452 |
| INS  | 1,481,487 |  91,439 |        0.891 |      0.548 |

``` r
## numbers of cliques
locs %>% group_by(clique) %>% summarize(sites=n()) %>% ungroup %>% mutate(prop=sites/sum(sites)) %>%
  kable(digits=3, format.args=list(big.mark=','))
```

| clique |   sites |  prop |
| :----- | ------: | ----: |
| FALSE  |  16,381 | 0.098 |
| TRUE   | 150,578 | 0.902 |

``` r
## numbers of cliques by type
locs %>% group_by(type, clique) %>% summarize(sites=n()) %>%
  group_by(type) %>% mutate(prop.type=sites/sum(sites)) %>% 
  kable(digits=3, format.args=list(big.mark=','))
```

| type | clique |  sites | prop.type |
| :--- | :----- | -----: | --------: |
| DEL  | FALSE  |  5,162 |     0.068 |
| DEL  | TRUE   | 70,358 |     0.932 |
| INS  | FALSE  | 11,219 |     0.123 |
| INS  | TRUE   | 80,220 |     0.877 |

## Size

``` r
ggp$size = locs %>% as.data.frame %>%
  ggplot(aes(x=size, fill=type)) +
  geom_histogram(position='dodge', bins=60) +
  scale_fill_brewer(palette='Set1', name='SV type') + 
  theme_bw() +
  xlab('size (bp)') +
  scale_x_log10(breaks=c(0, 50, 100, 300, 1000, 6000, 1e4, 1e5),
                labels=c(0, 50, 100, 300, '1,000', '6,000', '10,000', '100,000')) + 
  ylab('number of variants') +
  theme(legend.title=element_blank()) + 
  ## theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  theme(legend.position=c(.99, .99), legend.justification=c(1,1))
ggp$size
```

![](summary-sv-stats-mesa_files/figure-gfm/size-1.png)<!-- -->

``` r
locs %>% mutate(type='all') %>% rbind(locs) %>%
  group_by(type) %>%
  summarize(min.size=min(size), max.size=max(size),
    size.lt1kbp=mean(size<=1000), size.lt500=mean(size<=500))
```

    ## # A tibble: 3 x 5
    ##   type  min.size max.size size.lt1kbp size.lt500
    ##   <chr>    <int>    <int>       <dbl>      <dbl>
    ## 1 all         50   125187       0.944      0.897
    ## 2 DEL         50   114201       0.945      0.908
    ## 3 INS         50   125187       0.943      0.888

## Overlap with simple repeats, satellites or low-complexity regions

``` r
## simple repeats
if(!file.exists('simpleRepeat.hg38.txt.gz')){
  download.file('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz', 'simpleRepeat.hg38.txt.gz')
}
sr = read.table('simpleRepeat.hg38.txt.gz', as.is=TRUE)
sr = reduce(GRanges(sr$V2, IRanges(sr$V3, sr$V4)))
sr$repClass = 'Simple_repeat'
## repeat masker with low-complexity regions
if(!file.exists('rmsk.txt.gz')){
  download.file('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz', 'rmsk.txt.gz')
}
rm = read.table('rmsk.txt.gz', as.is=TRUE,
                colClasses=c(rep("NULL", 5), 'character', 'integer', 'integer',
                             rep('NULL', 3), 'character',  rep('NULL', 5)))
colnames(rm) = c('chr', 'start', 'end', 'repClass')
rm = subset(rm, repClass %in% c('Low_complexity', 'Simple_repeat', 'Satellite'))
rm = makeGRangesFromDataFrame(rm, keep.extra.columns=TRUE)
rm = c(rm, sr)

locs.gr = makeGRangesFromDataFrame(locs)

olRep <- function(locs.gr, rm.r){
  rm.r = reduce(rm.r)
  findOverlaps(locs.gr, rm.r) %>% as.data.frame %>%
    mutate(sv.w=width(locs.gr[queryHits]), ol.w=width(pintersect(locs.gr[queryHits], rm.r[subjectHits])),
           ol.prop=ol.w/sv.w) %>%
    group_by(queryHits) %>% summarize(ol.prop=sum(ol.prop), .groups='drop')
}

## all: simple repeats, low complexity, satellites
ol.df = olRep(locs.gr, rm)
locs$rep.sr.lc.sat = 0
locs$rep.sr.lc.sat[ol.df$queryHits] = ol.df$ol.prop
## simple repeats + low-complexity
ol.df = olRep(locs.gr, subset(rm, repClass %in% c('Simple_repeat', 'Low_complexity')))
locs$rep.sr.lc = 0
locs$rep.sr.lc[ol.df$queryHits] = ol.df$ol.prop
## simple repeats
ol.df = olRep(locs.gr, subset(rm, repClass=='Simple_repeat'))
locs$rep.sr = 0
locs$rep.sr[ol.df$queryHits] = ol.df$ol.prop
## low complexity
ol.df = olRep(locs.gr, subset(rm, repClass=='Low_complexity'))
locs$rep.lc = 0
locs$rep.lc[ol.df$queryHits] = ol.df$ol.prop
## simple repeats
ol.df = olRep(locs.gr, subset(rm, repClass=='Satellite'))
locs$rep.sat = 0
locs$rep.sat[ol.df$queryHits] = ol.df$ol.prop

locs %>% mutate(type='all') %>% rbind(locs) %>%
  group_by(type) %>% 
  summarize(rep.sr.lc.sat.50=mean(rep.sr.lc.sat>=.50), rep.sr.lc.50=mean(rep.sr.lc>=.50),
                   rep.sr.50=mean(rep.sr>=.50),
                   rep.lc.50=mean(rep.lc>=.50), rep.sat.50=mean(rep.sat>=.50)) %>%
  kable(digits=3)
```

| type | rep.sr.lc.sat.50 | rep.sr.lc.50 | rep.sr.50 | rep.lc.50 | rep.sat.50 |
| :--- | ---------------: | -----------: | --------: | --------: | ---------: |
| all  |            0.840 |        0.839 |     0.838 |     0.016 |      0.016 |
| DEL  |            0.830 |        0.829 |     0.829 |     0.015 |      0.018 |
| INS  |            0.847 |        0.847 |     0.846 |     0.016 |      0.015 |

*sr*: simple repeat; *lc*: low-complexity; *sat*: satellite DNA. *.50*
means that at least 50% of the SV region overlaps repeats.

### Non-clique SV sites are repeat-rich

We expect most non-clique sites, i.e. with very different alleles, to be
repeat variation like short-tandem repeats variation (or VNTRs). Is it?

``` r
locs %>% mutate(type='all') %>% rbind(locs) %>%
  filter(!clique) %>% 
  group_by(type) %>% 
  summarize(rep.sr.lc.sat.50=mean(rep.sr.lc.sat>=.50), rep.sr.lc.50=mean(rep.sr.lc>=.50),
                   rep.sr.50=mean(rep.sr>=.50),
                   rep.lc.50=mean(rep.lc>=.50), rep.sat.50=mean(rep.sat>=.50)) %>%
  kable(digits=3)
```

| type | rep.sr.lc.sat.50 | rep.sr.lc.50 | rep.sr.50 | rep.lc.50 | rep.sat.50 |
| :--- | ---------------: | -----------: | --------: | --------: | ---------: |
| all  |            0.976 |        0.976 |     0.976 |     0.015 |      0.007 |
| DEL  |            0.987 |        0.987 |     0.987 |     0.008 |      0.009 |
| INS  |            0.971 |        0.971 |     0.970 |     0.018 |      0.007 |

Yes, almost all are within simple repeats. What are the ones that are
not?

``` r
locs.nc = locs %>% filter(!clique, rep.sr.lc.sat<=.5)

## distance to simple repeat
locs.nc.gr = makeGRangesFromDataFrame(locs.nc)
dd = distanceToNearest(locs.nc.gr, rm) %>% as.data.frame
locs.nc$rep.dist = NA
locs.nc$rep.dist[dd$queryHits] = dd$distance

## random subset
set.seed(123)
locs.nc %>%
  filter(size.min/size.max>.8) %>% 
  mutate(coord=paste0('[', seqnames, ':', start, '-', end,
                      '](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=',
                      seqnames, '%3A', start, '%2D', end, ')')) %>% 
  select(coord, svsite, type, size, loc.n, size.min, size.max, rep.dist, rep.sr.lc.sat) %>% sample_n(10) %>%
  kable
```

| coord                                                                                                              | svsite         | type | size | loc.n | size.min | size.max | rep.dist | rep.sr.lc.sat |
| :----------------------------------------------------------------------------------------------------------------- | :------------- | :--- | ---: | ----: | -------: | -------: | -------: | ------------: |
| [chr16:14894534-14894534](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr16%3A14894534%2D14894534)   | sv\_530096\_0  | INS  |  131 |     4 |      117 |      131 |     3838 |     0.0000000 |
| [chr22:37529117-37529117](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr22%3A37529117%2D37529117)   | sv\_100952\_0  | INS  |  290 |   237 |      280 |      306 |      919 |     0.0000000 |
| [chr17:43270206-43271985](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr17%3A43270206%2D43271985)   | sv\_448934\_0  | DEL  | 1779 |     3 |     1779 |     2125 |        0 |     0.0983146 |
| [chr11:64237124-64237124](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr11%3A64237124%2D64237124)   | sv\_859480\_0  | INS  |  144 |    66 |      121 |      144 |        0 |     0.0000000 |
| [chr2:116950767-116950767](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr2%3A116950767%2D116950767) | sv\_1797411\_0 | INS  |  105 |     7 |       86 |      105 |       15 |     0.0000000 |
| [chr17:43251137-43251546](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr17%3A43251137%2D43251546)   | sv\_447044\_0  | DEL  |  409 |     4 |      407 |      430 |        0 |     0.4073171 |
| [chr17:43268211-43268917](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr17%3A43268211%2D43268917)   | sv\_448845\_0  | DEL  |  706 |     3 |      573 |      712 |     1089 |     0.0000000 |
| [chr9:70801717-70801717](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr9%3A70801717%2D70801717)     | sv\_1022302\_0 | INS  |  297 |   292 |      271 |      331 |     1315 |     0.0000000 |
| [chr17:43254597-43254971](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr17%3A43254597%2D43254971)   | sv\_447396\_0  | DEL  |  374 |     3 |      327 |      374 |      402 |     0.0000000 |
| [chr5:25540972-25540972](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr5%3A25540972%2D25540972)     | sv\_1483244\_0 | INS  |  101 |    32 |       97 |      117 |      544 |     0.0000000 |

Either very close to repeats, or in segmental duplication or
transposons, or slightly below the 80% threshold used to define matching
alleles.

## Gene annotation

``` r
if(!file.exists('gencode.v35.annotation.gtf.gz')){
  download.file('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz', 'gencode.v35.annotation.gtf.gz')
}
genc = import('gencode.v35.annotation.gtf.gz')

genc.pc = subset(genc, type %in% c('CDS', 'UTR', 'gene') & gene_type=='protein_coding')
prom = promoters(subset(genc.pc, type=='gene'))
prom$type = 'promoter'
genc.pc = c(genc.pc, prom)
ol.gene = findOverlaps(locs.gr, genc.pc) %>% as.data.frame %>%
  mutate(gene=genc.pc$gene_name[subjectHits], type=genc.pc$type[subjectHits]) %>%
  group_by(queryHits, gene) %>%
  summarize(cds=any(type=='CDS'))

rbind(ol.gene %>% filter(cds) %>% mutate(impact='cds'),
      ol.gene %>% mutate(impact='cds.prom.utr.intron')) %>%
  group_by(impact) %>%
  summarize(sv=length(unique(queryHits)), gene=length(unique(gene))) %>%
  kable
```

| impact              |    sv | gene |
| :------------------ | ----: | ---: |
| cds                 |  1561 |  408 |
| cds.prom.utr.intron | 77940 | 7740 |

## Allele frequency

``` r
ggp$af = locs %>% as.data.frame %>%
  ggplot(aes(x=af, fill=type)) +
  geom_histogram(position='dodge') +
  scale_fill_brewer(palette='Set1', name='SV type') + 
  theme_bw() +
  xlab('allele frequency') +
  ylab('number of SV loci') +
  theme(legend.position=c(.99, .99), legend.justification=c(1,1))
ggp$af
```

![](summary-sv-stats-mesa_files/figure-gfm/af-1.png)<!-- -->

``` r
## comparing allele frequency between top 2 most frequent alleles in SV loci
locs.s.3 = locs %>% filter(loc.n>1, af.top.fc>3) %>%
  mutate(loc.n=cut(loc.n, breaks=c(1:5,100,Inf), labels=c(2:5, '5-100', '>100'))) %>% 
  group_by(loc.n) %>% summarize(n=n())
ggp$af.top = locs %>% filter(loc.n>1) %>%
  mutate(loc.n=cut(loc.n, breaks=c(1:5,100,Inf), labels=c(2:5, '5-100', '>100'))) %>% 
  ggplot(aes(x=loc.n, y=winsor(af.top.fc, 10))) +
  geom_violin(scale='width', fill='orange') +
  theme_bw() +
  scale_y_continuous(breaks=1:10, labels=c(1:9, '10+')) +
  ylab('frequency fold-change between\ntop 2 most frequent alleles') +
  xlab('number of SVs in locus') +
  theme(legend.title=element_blank()) + 
  geom_hline(yintercept=3, linetype=2) +
  geom_label(aes(label=n), y=6, data=locs.s.3, size=3)
ggp$af.top
```

![](summary-sv-stats-mesa_files/figure-gfm/af-2.png)<!-- -->

``` r
locs %>% filter(loc.n>1) %>%
  summarize(af.fc.3=sum(af.top.fc>3),
            af01.af2lt01=sum(af>=.01 & af.top2<.01),
            af01.af2lt01.fc3=sum(af>=.01 & af.top2<.01 & af.top.fc>3),
            major.al=sum(af>af.top2)) %>%
  kable
```

| af.fc.3 | af01.af2lt01 | af01.af2lt01.fc3 | major.al |
| ------: | -----------: | ---------------: | -------: |
|   14720 |         7346 |             5936 |    37667 |

  - *af.fc.3*: SV sites where the most frequent allele is at least 3
    times more frequent than the seoncd most frequent allele.
  - *af01.af2lt01*: SV sites where most frequent allele with frequency
    \>1% but other alleles with \<1% frequency.
  - *af01.af2lt01.fc3*: Same + the most frequent allele is at least 3
    times more frequent than the seoncd most frequent allele.
  - *major.al*: SV sites with one allele more frequent than the other
    alleles.

## Alleles per SV sites

``` r
## number of alleles per loci
ggp$loc.al = locs %>% mutate(loc.n=cut(loc.n, breaks=c(0:5,100,Inf), labels=c(1:5, '5-100', '>100'))) %>%
  group_by(loc.n, type) %>% summarize(n=n(), .groups='drop') %>%
  arrange(loc.n) %>% group_by(type) %>% mutate(cprop=cumsum(n)/sum(n)) %>% 
  ggplot(aes(x=loc.n, y=cprop, color=type, group=type)) +
  geom_line() + 
  geom_point() + 
  theme_bw() +
  ylim(0,1) + 
  scale_color_brewer(palette='Set1', name='SV type') + 
  ylab('cumulative proportion of SV loci') +
  xlab('number of SVs in locus') +
  theme(legend.title=element_blank()) + 
  theme(legend.position=c(.99, .01), legend.justification=c(1,0))
ggp$loc.al
```

![](summary-sv-stats-mesa_files/figure-gfm/site-al-1.png)<!-- -->

``` r
## proportion of cliques, i.e. all alleles similar, in a SV locus
ggp$loc.cl = locs %>% mutate(loc.n=cut(loc.n, breaks=c(0:5,100,Inf), labels=c(1:5, '5-100', '>100'))) %>%
  group_by(loc.n, type) %>% summarize(cl.prop=mean(clique)) %>% 
  ggplot(aes(x=loc.n, y=cl.prop, color=type, group=type)) +
  geom_line() + 
  geom_point() + 
  theme_bw() +
  ylim(0,1) + 
  scale_color_brewer(palette='Set1', name='SV type') + 
  ylab('proportion of SV loci with\nalleles in clique formation') +
  xlab('number of SVs in locus') +
  theme(legend.title=element_blank()) + 
  theme(legend.position=c(.99, .99), legend.justification=c(1,1))
ggp$loc.cl
```

![](summary-sv-stats-mesa_files/figure-gfm/site-al-2.png)<!-- -->

## Multi-panel figure

``` r
## adds a legend title: a), b), etc
plot_list <- function(ggp.l, gg.names=NULL){
  if(is.null(names(ggp.l))) names(ggp.l) = paste0('g', 1:length(ggp.l))
  if(is.null(gg.names)) gg.names = names(ggp.l)
  lapply(1:length(gg.names), function(ii) ggp.l[[gg.names[ii]]] + ggtitle(paste0('(', LETTERS[ii], ')')))
}

grid.arrange(grobs=plot_list(ggp, c("af", "size", "loc.al", "loc.cl", "af.top")),
             layout_matrix=matrix(c(1,1,1,2,2,2,3,4,5), nrow=3, byrow=TRUE))
```

![](summary-sv-stats-mesa_files/figure-gfm/fig-1.png)<!-- -->

``` r
pdf('figs/fig-sv-mesa-stats.pdf', 10, 8)
grid.arrange(grobs=plot_list(ggp, c("af", "size", "loc.al", "loc.cl", "af.top")),
             layout_matrix=matrix(c(1,1,1,2,2,2,3,4,5), nrow=3, byrow=TRUE))
dev.off()
```

    ## png 
    ##   2

## Save TSV with SV site information

``` r
outf = gzfile('locs.mesa2k.svsite80al.tsv.gz', 'w')
write.table(locs, file=outf, row.names=FALSE, quote=FALSE, sep='\t')
close(outf)
```

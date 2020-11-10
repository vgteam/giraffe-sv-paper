Compare SV genotypes with public SV catalogs
================

``` r
library(dplyr)
library(sveval)
library(GenomicRanges)
library(ggplot2)
library(gridExtra)
library(knitr)
winsor <- function(x, u){
  if(any(x>u)) x[x>u] = u
  x
}
```

## Simple repeat annotation

Can be used to rescue misplaced SVs within simple repeats.

``` r
if(!file.exists('simpleRepeat.hg38.txt.gz')){
  download.file('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz', 'simpleRepeat.hg38.txt.gz')
}
sr = read.table('simpleRepeat.hg38.txt.gz', as.is=TRUE)
sr = reduce(GRanges(sr$V2, IRanges(sr$V3, sr$V4)))
```

## Public SV catalogs

Prepared using `prepare-public-catalogs.R`

``` r
load('public-sv-catalogs.RData', verbose=TRUE)
```

    ## Loading objects:
    ##   kgp3
    ##   gnomad
    ##   svpop

### SV catalog from the 1000 Genomes Project phase 3

``` r
## samples
max(kgp3$ncalls)
```

    ## [1] 2504

``` r
## SVs
length(kgp3)
```

    ## [1] 62860

``` r
## per type
table(kgp3$type)
```

    ## 
    ##   DEL   INS   INV 
    ## 44986 17088   786

### gnomAD-SV catalog

``` r
## SVs
length(gnomad)
```

    ## [1] 342702

``` r
## per type
table(gnomad$type)
```

    ## 
    ##    CPX    DEL    INS    INV 
    ##   5106 173830 163130    636

### SVPOP

``` r
## samples
max(svpop$ncalls)
```

    ## [1] 440

``` r
## SVs
length(svpop)
```

    ## [1] 91991

``` r
## per type
table(svpop$type)
```

    ## 
    ##   DEL   INS 
    ## 38542 53449

## SVs genotyped using vg

### SVs genotyped in MESA

``` r
## SVs grouped by site ('svsite' and 'clique' columns)
mesa = read.table('svs.mesa2k.svsite80al.tsv.gz', as.is=TRUE, header=TRUE)

## stats for each SV locus
## use the most frequent allele (and then the largest) for ac/af/size
## also saves sum/max/min across all alleles
mesa.s = mesa %>% arrange(desc(af), desc(size)) %>%
  group_by(seqnames, svsite, type, clique) %>%
  summarize(start=start[1], end=end[1], ac=ac[1], af=af[1], size=size[1], .groups='drop') %>%
  filter(size>=50) %>% makeGRangesFromDataFrame(keep.extra.columns=TRUE)
```

### SVs genotyped in 2,504 samples from 1000 Genomes Project

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

## Frequency distributions

Variants with at least 1% frequency.

``` r
qplot(x=kgp3$af) + theme_bw() +
  ylab('1000GP phase 3 variants') + xlab('allele frequency') + xlim(.01,1)
```

![](compare-public-catalogs_files/figure-gfm/freq_dist-1.png)<!-- -->

``` r
qplot(x=gnomad$af) + theme_bw() +
  ylab('gnomAD-SV variants') + xlab('allele frequency') + xlim(.01,1)
```

![](compare-public-catalogs_files/figure-gfm/freq_dist-2.png)<!-- -->

``` r
qplot(x=svpop$af) + theme_bw() +
  ylab('SVPOP variants') + xlab('allele frequency') + xlim(.01,1)
```

![](compare-public-catalogs_files/figure-gfm/freq_dist-3.png)<!-- -->

``` r
qplot(x=mesa.s$af) + theme_bw() +
  ylab('MESA variants') + xlab('allele frequency') + xlim(.01,1)
```

![](compare-public-catalogs_files/figure-gfm/freq_dist-4.png)<!-- -->

``` r
qplot(x=kgp.s$af) + theme_bw() +
  ylab('vg-1000GP variants') + xlab('allele frequency') + xlim(.01,1)
```

![](compare-public-catalogs_files/figure-gfm/freq_dist-5.png)<!-- -->

## Comparison with the 1000 Genomes Project phase 3 catalog

``` r
olcatstats = function(sites.gr, cat.gr, min.ol=.1, max.ins.dist=200, use.sr=FALSE){
  if(!use.sr){
    sr = NULL
  }
  sites.gr$ac = 1
  cat.gr$ac = 1
  ol.gr = suppressWarnings(svOverlap(sites.gr, cat.gr, min.ol=min.ol, max.ins.dist=max.ins.dist, simprep=sr))
  ## compute proportions
  tibble(prop.vgsite=length(unique(ol.gr$queryHits)) / length(sites.gr),
         prop.cat=length(unique(ol.gr$subjectHits)) / length(cat.gr))
}

rbind(
  olcatstats(mesa.s, kgp3) %>% mutate(set='MESA vs 1000GP'),
  olcatstats(subset(mesa.s, af>=.05), subset(kgp3, af>=.05)) %>% mutate(set='MESA freq>=5% vs 1000GP freq>=5%'),
  olcatstats(kgp.s, kgp3) %>% mutate(set='vg-1000GP vs 1000GP'),
  olcatstats(subset(kgp.s, af>=.05), subset(kgp3, af>=.05)) %>% mutate(set='vg-1000GP freq>=5% vs 1000GP freq>=5%')
) %>% select(set, prop.vgsite, prop.cat) %>% kable(digits=3)
```

| set                                     | prop.vgsite | prop.cat |
| :-------------------------------------- | ----------: | -------: |
| MESA vs 1000GP                          |       0.068 |    0.165 |
| MESA freq\>=5% vs 1000GP freq\>=5%      |       0.087 |    0.792 |
| vg-1000GP vs 1000GP                     |       0.068 |    0.166 |
| vg-1000GP freq\>=5% vs 1000GP freq\>=5% |       0.090 |    0.819 |

``` r
ol.gr = suppressWarnings(svOverlap(mesa.s, kgp3, min.ol=.1, max.ins.dist=200))
freq.mesa.kgp3.df = ol.gr %>% as.data.frame %>%
  mutate(af=mesa.s$af[queryHits], cat.af=kgp3$af[subjectHits])

ggplot(freq.mesa.kgp3.df, aes(x=af, y=cat.af)) +
  geom_point(alpha=.3) +
  xlab('allele frequency in MESA') + ylab('allele frequency in 1000GP phase 3') + 
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/freqcomp_kgp3-1.png)<!-- -->

``` r
ggplot(freq.mesa.kgp3.df, aes(x=af-cat.af)) +
  geom_histogram() +
  xlab('frequency difference (MESA - 1000GP phase 3)') +
  ylab('SV') + 
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/freqcomp_kgp3-2.png)<!-- -->

``` r
ol.gr = suppressWarnings(svOverlap(kgp.s, kgp3, min.ol=.1, max.ins.dist=200))
freq.kgp.kgp3.df = ol.gr %>% as.data.frame %>%
  mutate(af=kgp.s$af[queryHits], cat.af=kgp3$af[subjectHits])

ggplot(freq.kgp.kgp3.df, aes(x=af, y=cat.af)) +
  geom_point(alpha=.3) +
  xlab('allele frequency in vg-1000GP') + ylab('allele frequency in 1000GP phase 3') + 
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/freqcomp_kgp3-3.png)<!-- -->

``` r
ggplot(freq.kgp.kgp3.df, aes(x=af-cat.af)) +
  geom_histogram() +
  xlab('frequency difference (vg-1000GP - 1000GP phase 3)') +
  ylab('SV') + 
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/freqcomp_kgp3-4.png)<!-- -->

## Comparison with the gnomAD-SV catalog

``` r
rbind(
  olcatstats(mesa.s, gnomad) %>% mutate(set='MESA vs gnomAD-SV'),
  olcatstats(subset(mesa.s, af>=.05), subset(gnomad, af>=.05)) %>% mutate(set='MESA freq>=5% vs gnomAD-SV freq>=5%'),
  olcatstats(kgp.s, gnomad) %>% mutate(set='vg-1000GP vs gnomAD-SV'),
  olcatstats(subset(kgp.s, af>=.05), subset(gnomad, af>=.05)) %>% mutate(set='vg-1000GP freq>=5% vs gnomAD-SV freq>=5%')
) %>% select(set, prop.vgsite, prop.cat) %>% kable(digits=3)
```

| set                                        | prop.vgsite | prop.cat |
| :----------------------------------------- | ----------: | -------: |
| MESA vs gnomAD-SV                          |       0.283 |    0.092 |
| MESA freq\>=5% vs gnomAD-SV freq\>=5%      |       0.151 |    0.588 |
| vg-1000GP vs gnomAD-SV                     |       0.284 |    0.091 |
| vg-1000GP freq\>=5% vs gnomAD-SV freq\>=5% |       0.154 |    0.599 |

``` r
ol.gr = suppressWarnings(svOverlap(mesa.s, gnomad, min.ol=.1, max.ins.dist=200))
freq.mesa.gnomad.df = ol.gr %>% as.data.frame %>%
  mutate(af=mesa.s$af[queryHits], cat.af=gnomad$af[subjectHits])

ggplot(freq.mesa.gnomad.df, aes(x=af, y=cat.af)) +
  geom_point(alpha=.3) +
  xlab('allele frequency in MESA') + ylab('allele frequency in gnomAD-SV') + 
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/freqcomp_gnomad-1.png)<!-- -->

``` r
ggplot(freq.mesa.gnomad.df, aes(x=af-cat.af)) +
  geom_histogram() +
  xlab('frequency difference (MESA - gnomAD-SV)') +
  ylab('SV') + 
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/freqcomp_gnomad-2.png)<!-- -->

``` r
ol.gr = suppressWarnings(svOverlap(kgp.s, gnomad, min.ol=.1, max.ins.dist=200))
freq.kgp.gnomad.df = ol.gr %>% as.data.frame %>%
  mutate(af=kgp.s$af[queryHits], cat.af=gnomad$af[subjectHits])

ggplot(freq.kgp.gnomad.df, aes(x=af, y=cat.af)) +
  geom_point(alpha=.3) +
  xlab('allele frequency in vg-1000GP') + ylab('allele frequency in gnomAD-SV') + 
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/freqcomp_gnomad-3.png)<!-- -->

``` r
ggplot(freq.kgp.gnomad.df, aes(x=af-cat.af)) +
  geom_histogram() +
  xlab('frequency difference (vg-1000GP - gnomAD-SV)') +
  ylab('SV') + 
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/freqcomp_gnomad-4.png)<!-- -->

## Comparison with SVPOP

``` r
rbind(
  olcatstats(mesa.s, svpop) %>% mutate(set='MESA vs SVPOP'),
  olcatstats(kgp.s, svpop) %>% mutate(set='vg-1000GP vs SVPOP')
) %>% select(set, prop.vgsite, prop.cat) %>% kable(digits=3)
```

| set                | prop.vgsite | prop.cat |
| :----------------- | ----------: | -------: |
| MESA vs SVPOP      |       0.860 |    0.945 |
| vg-1000GP vs SVPOP |       0.859 |    0.933 |

``` r
ol.gr = suppressWarnings(svOverlap(mesa.s, svpop, min.ol=.1, max.ins.dist=200))
freq.mesa.svpop.df = ol.gr %>% as.data.frame %>%
  mutate(af=mesa.s$af[queryHits], cat.af=svpop$af[subjectHits])

ggplot(freq.mesa.svpop.df, aes(x=af, y=cat.af)) +
  geom_point(alpha=.3) +
  xlab('allele frequency in MESA') + ylab('allele frequency in SVPOP') + 
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/freqcomp_svpop-1.png)<!-- -->

``` r
ggplot(freq.mesa.svpop.df, aes(x=af-cat.af)) +
  geom_histogram() +
  xlab('frequency difference (MESA - SVPOP)') +
  ylab('SV') + 
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/freqcomp_svpop-2.png)<!-- -->

``` r
ol.gr = suppressWarnings(svOverlap(kgp.s, svpop, min.ol=.1, max.ins.dist=200))
freq.kgp.svpop.df = ol.gr %>% as.data.frame %>%
  mutate(af=kgp.s$af[queryHits], cat.af=svpop$af[subjectHits])

freq.kgp.svpop.df %>% filter(af>.05, cat.af>.05) %>% 
  ggplot(aes(x=af, y=cat.af)) +
  geom_point(alpha=.5) +
  xlab('allele frequency in vg-1000GP') + ylab('allele frequency in SVPOP') + 
  ## geom_bin2d() +
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/freqcomp_svpop-3.png)<!-- -->

``` r
ggplot(freq.kgp.svpop.df, aes(x=af-cat.af)) +
  geom_histogram() +
  xlab('frequency difference (vg-1000GP - SVPOP)') +
  ylab('SV') + 
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/freqcomp_svpop-4.png)<!-- -->

### SVPOP vs gnomAD-SV/1000GP

``` r
ol.gr = suppressWarnings(svOverlap(svpop, kgp3, min.ol=.1, max.ins.dist=200))
freq.svpop.kgp3.df = ol.gr %>% as.data.frame %>%
  mutate(af=svpop$af[queryHits], cat.af=kgp3$af[subjectHits])

ggplot(freq.svpop.kgp3.df, aes(x=af, y=cat.af)) +
  geom_point(alpha=.3) +
  xlab('allele frequency in SVPOP') + ylab('allele frequency in 1000GP phase 3') +
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/svpop_vs_gnomad-1.png)<!-- -->

``` r
ggplot(freq.svpop.kgp3.df, aes(x=af-cat.af)) +
  geom_histogram() +
  xlab('frequency difference (SVPOP - 1000GP phase 3)') +
  ylab('SV') + 
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/svpop_vs_gnomad-2.png)<!-- -->

``` r
ol.gr = suppressWarnings(svOverlap(svpop, gnomad, min.ol=.1, max.ins.dist=200))
freq.svpop.gnomad.df = ol.gr %>% as.data.frame %>%
  mutate(af=svpop$af[queryHits], cat.af=gnomad$af[subjectHits])

ggplot(freq.svpop.gnomad.df, aes(x=af, y=cat.af)) +
  geom_point(alpha=.3) +
  xlab('allele frequency in SVPOP') + ylab('allele frequency in gnomAD-SV') +
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/svpop_vs_gnomad-3.png)<!-- -->

``` r
ggplot(freq.svpop.gnomad.df, aes(x=af-cat.af)) +
  geom_histogram() +
  xlab('frequency difference (SVPOP - gnomAD-SV)') +
  ylab('SV') + 
  theme_bw()
```

![](compare-public-catalogs_files/figure-gfm/svpop_vs_gnomad-4.png)<!-- -->

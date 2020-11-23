Frequency in populations from the 1000 Genomes Project
================

``` r
library(dplyr)
library(ggplot2)
library(gridExtra)
library(knitr)
winsor <- function(x, u){
  if(any(x>u)) x[x>u] = u
  x
}
## list of figures
ggp = list()
```

## Read frequencies

``` r
## frequencies in different super-populations for each SV site
freq.df = read.table('2504kgp.svsite80al.superpopfreq.tsv.gz', as.is=TRUE, header=TRUE)
```

## Read null frequencies

To compute an expected distribution, sample labels were permuted.

``` r
freq.null = read.table('2504kgp.svsite80al.superpopfreq.null.tsv.gz', as.is=TRUE, header=TRUE)

freq.df = rbind(
  freq.df %>% mutate(exp='observed'),
  freq.null %>% mutate(exp='expected'))
```

## Compute the minimum, median and maximum frequency

``` r
freq.s = freq.df %>% group_by(exp, svsite) %>%
  summarize(af.min=min(af), af.max=max(af), af.med=median(af), .groups='drop')

ggp$range = ggplot(freq.s, aes(winsor(af.max-af.min, .5), fill=exp)) +
  geom_histogram(position='dodge') + theme_bw() +
  xlab('frequency range across populations/clusters') +
  scale_fill_brewer(palette='Set1') +
  scale_x_continuous(breaks=seq(0,.5,.1), labels=c(seq(0,.4,.1), '0.5+')) + 
  theme(legend.title=element_blank(),
        legend.position=c(.99,.99), legend.justification=c(1,1)) + 
  ylab('SV site')
ggp$range
```

![](pops-freq-1kgp_files/figure-gfm/freqmmm-1.png)<!-- -->

``` r
lapply(c(.1,.25,.5), function(th){
  freq.s %>% group_by(exp) %>% summarize(svsite=sum(af.max-af.min>th), .groups='drop') %>% 
    mutate(min.af.range=th)
}) %>% bind_rows %>%
  select(min.af.range, exp, svsite) %>%
  kable()
```

| min.af.range | exp      | svsite |
| -----------: | :------- | -----: |
|         0.10 | expected |     53 |
|         0.10 | observed |  37147 |
|         0.25 | expected |      0 |
|         0.25 | observed |  10185 |
|         0.50 | expected |      0 |
|         0.50 | observed |    989 |

The table shows the number of sites with a frequency range larger than
10%, 25%, and 50%.

## SV sites with large deviation from the median frequency

``` r
freq.df = freq.df %>% group_by(exp, svsite) %>% mutate(af.med=median(af))

bar.df = lapply(c(.1,.25,.5), function(th){
  freq.df %>% group_by(Superpopulation, exp) %>% summarize(sv=sum(abs(af.med-af)>th), .groups='drop') %>% 
    mutate(af.dev=paste0('delta.af>', th))
}) %>% bind_rows

ggp$olbar = bar.df %>% group_by(af.dev) %>%
  mutate(y.label=ifelse(sv>.05*max(sv) | sv==0, NA, .1*max(sv))) %>% 
  ggplot(aes(x=factor(Superpopulation), y=sv, fill=exp)) + 
  geom_bar(stat='identity', position='dodge') + theme_bw() +
  geom_label(aes(label=sv, color=exp, y=y.label),
             fill='white', position=position_dodge(.9), show.legend=FALSE) + 
  facet_grid(af.dev~., scales='free') + 
  scale_fill_brewer(palette='Set1', name='') + 
  scale_color_brewer(palette='Set1', name='') +
  theme(strip.text.y=element_text(angle=0)) + 
  xlab('population/cluster') +
  ylab('number of SV site')
ggp$olbar
```

![](pops-freq-1kgp_files/figure-gfm/popspec-1.png)<!-- -->

``` r
lapply(c(.1,.25,.5), function(th){
  freq.df %>% filter(af.med-af>th) %>%
    group_by(exp) %>%
    summarize(svsite=length(unique(svsite)), .groups='drop') %>% 
    mutate(min.af.dev=th)
}) %>% bind_rows %>%
  select(min.af.dev, exp, svsite) %>%
  kable()
```

| min.af.dev | exp      | svsite |
| ---------: | :------- | -----: |
|       0.10 | expected |      9 |
|       0.10 | observed |  13896 |
|       0.25 | observed |   1899 |
|       0.50 | observed |     79 |

The table shows the number of sites with a population-specific frequency
pattern, defined as deviating form the median frequency by at least 10%,
25%, and 50%.

## Population-specific SVs by type and size

``` r
## SVs grouped by site ('svsite' and 'clique' columns)
svs = read.table('svs.2504kgp.svsite80al.tsv.gz', as.is=TRUE, header=TRUE)
svs = svs %>% group_by(svsite, type) %>% summarize(size=mean(size), .groups='drop')

ggp$devtype = freq.df %>% group_by(exp) %>%  filter(abs(af.med-af)>.1) %>% 
  merge(svs) %>%
  ggplot(aes(x=winsor(abs(af.med-af), 1), fill=type)) +
  geom_histogram(bins=50, position='dodge') +
  scale_x_continuous(breaks=seq(.1, 1, .1), limits=c(.1, 1)) + 
  theme_bw() +
  xlab('frequency deviation from the mean (winsorized at 1)') +
  ylab('SV site in a population/cluster')
ggp$devtype
```

![](pops-freq-1kgp_files/figure-gfm/popspec_type_size-1.png)<!-- -->

``` r
ggp$devsize = freq.df %>% group_by(exp) %>%  filter(abs(af.med-af)>.1) %>% 
  merge(svs) %>%
  ggplot(aes(x=size, fill=type)) +
  geom_histogram(bins=50, position='dodge') +
  scale_x_log10() + 
  theme_bw() +
  xlab('SV size (bp)') +
  ylab('SV site in a population/cluster')
ggp$devsize
```

![](pops-freq-1kgp_files/figure-gfm/popspec_type_size-2.png)<!-- -->

## Examples of population-specific SVs

More frequent in a population

``` r
ex.df = freq.df %>%
  filter(exp=='observed', af<1, af.med<1, af>af.med) %>% 
  arrange(desc(abs(af.med-af))) %>%
  group_by(Superpopulation) %>% do(head(.,1))
ex.df
```

    ## # A tibble: 5 x 5
    ## # Groups:   Superpopulation [5]
    ##   svsite       Superpopulation    af exp      af.med
    ##   <chr>        <chr>           <dbl> <chr>     <dbl>
    ## 1 sv_48359_0   AFR             0.906 observed 0.0798
    ## 2 sv_673793_0  AMR             0.390 observed 0.0654
    ## 3 sv_595085_0  EAS             0.822 observed 0.198 
    ## 4 sv_1899095_0 EUR             0.789 observed 0.301 
    ## 5 sv_696864_0  SAS             0.862 observed 0.490

``` r
ex.df = freq.df %>% filter(exp=='observed', svsite %in% ex.df$svsite)

ggplot(ex.df, aes(x=svsite, y=af, fill=Superpopulation)) +
  geom_bar(stat='identity', position='dodge') +
  ylab('allele frequency') + 
  theme_bw()
```

![](pops-freq-1kgp_files/figure-gfm/ex_freq-1.png)<!-- -->

Less frequent in a population

``` r
ex.df = freq.df %>%
  filter(exp=='observed', af<1, af.med<1, af<af.med) %>% 
  arrange(desc(abs(af.med-af))) %>%
  group_by(Superpopulation) %>% do(head(.,1))
ex.df
```

    ## # A tibble: 5 x 5
    ## # Groups:   Superpopulation [5]
    ##   svsite      Superpopulation     af exp      af.med
    ##   <chr>       <chr>            <dbl> <chr>     <dbl>
    ## 1 sv_49763_0  AFR             0.0416 observed  0.843
    ## 2 sv_795467_0 AMR             0.239  observed  0.599
    ## 3 sv_990881_0 EAS             0.126  observed  0.725
    ## 4 sv_939820_0 EUR             0.311  observed  0.806
    ## 5 sv_634762_0 SAS             0.444  observed  0.817

``` r
ex.df = freq.df %>% filter(exp=='observed', svsite %in% ex.df$svsite)

ggplot(ex.df, aes(x=svsite, y=af, fill=Superpopulation)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_brewer(palette='Set1') + 
  ylab('allele frequency') + 
  theme_bw()
```

![](pops-freq-1kgp_files/figure-gfm/ex_rare-1.png)<!-- -->

## Multi-panel figure

``` r
## adds a legend title: a), b), etc
plot_list <- function(ggp.l, gg.names=NULL){
  if(is.null(names(ggp.l))) names(ggp.l) = paste0('g', 1:length(ggp.l))
  if(is.null(gg.names)) gg.names = names(ggp.l)
  lapply(1:length(gg.names), function(ii) ggp.l[[gg.names[ii]]] + ggtitle(paste0(letters[ii], ')')))
}

grid.arrange(grobs=plot_list(list(ggp$range, ggp$olbar + guides(fill=FALSE))),
             layout_matrix=matrix(1:2, 2))
```

![](pops-freq-1kgp_files/figure-gfm/fig-1.png)<!-- -->

``` r
pdf('figs/fig-pops-freq-kgp.pdf', 8, 8)
grid.arrange(grobs=plot_list(list(ggp$range, ggp$olbar + guides(fill=FALSE))),
             layout_matrix=matrix(1:2, 2))
dev.off()
```

    ## png 
    ##   2

## Save SVs with population-specific patterns

``` r
freq.df %>%
  filter(exp=='observed', abs(af.med-af)>.1) %>%
  ungroup %>% 
  select(-exp) %>% 
  write.table(file='pops-freq-1kgp-med1.tsv', sep='\t', row.names=FALSE, quote=FALSE)
```

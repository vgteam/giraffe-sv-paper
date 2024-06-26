Evaluation of the SV genotypes
================

``` r
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(knitr)
library(sveval)
library(RColorBrewer)
## graph list
ggp = list()
```

## Results from sveval

``` r
pr.df = rbind(
  read.table('./hgsvc-hsvlr-sveval-prcurve.tsv', header=TRUE, as.is=TRUE),
  read.table('./hgsvc-hsvlr-40x-sveval-prcurve.tsv', header=TRUE, as.is=TRUE) %>% mutate(method='hsvlr_graph_40x'),
  read.table('./hgsvc-map-sveval-prcurve.tsv', header=TRUE, as.is=TRUE),
  read.table('./giab6_hg37-sveval-prcurve.tsv', header=TRUE, as.is=TRUE),
  read.table('./giab6_hg38-sveval-prcurve.tsv', header=TRUE, as.is=TRUE)
)

pr.df %>% select(exp, method, sample) %>% unique %>% arrange(exp, method, sample) %>% kable(row.names=FALSE)
```

| exp         | method            | sample       |
| :---------- | :---------------- | :----------- |
| giab6\_hg37 | giab\_graph       | HG002        |
| giab6\_hg37 | giab\_graph\_map  | HG002        |
| giab6\_hg38 | hsvlr\_graph      | HG002        |
| hgsvc       | hgsvc\_graph      | HG00514      |
| hgsvc       | hgsvc\_graph      | HG00733      |
| hgsvc       | hgsvc\_graph      | NA19240      |
| hgsvc       | hgsvc\_graph\_map | HG00514      |
| hgsvc       | hsvlr\_graph      | HG00514      |
| hgsvc       | hsvlr\_graph      | HG00733      |
| hgsvc       | hsvlr\_graph      | NA19240      |
| hgsvc       | hsvlr\_graph\_40x | HG00514\_40x |

``` r
qual.bestf1 = pr.df %>% filter(type=='Total') %>%
  select(exp, method, sample, region, eval, min.cov, F1, qual) %>%
  group_by(exp, method, sample, region, eval, min.cov) %>%
  arrange(desc(F1)) %>% do(head(., 1)) %>% select(-F1)

pr.df = pr.df %>%
  filter(type!='Total', type!='INV', min.cov==.5,
         region%in%c('all', 'highconf')) %>% arrange(qual) %>%
  merge(qual.bestf1)
```

## vg map vs vg giraffe

``` r
## HG00514 only for HGSVC
ggp$map = pr.df %>%
  mutate(exp=paste(method, exp, sample),
         exp=factor(exp,
                    levels=c('hgsvc_graph_map hgsvc HG00514',
                             'hgsvc_graph hgsvc HG00514',
                             'giab_graph_map giab6_hg37 HG002',
                                    'giab_graph giab6_hg37 HG002'),
                    labels=c('HGSVC graph - vg map',
                             'HGSVC graph - vg giraffe',
                             'GIAB graph - vg map',
                             'GIAB graph - vg giraffe')),
         type=factor(type, levels=c('Total', 'INS', 'DEL', 'INV')),
         region=factor(region, levels=c('all','highconf'), labels=c('all', 'high-confidence')),
         eval=factor(eval, levels=c('call','geno'), labels=c('presence', 'genotype'))) %>%
  filter(!is.na(exp)) %>% 
  ggplot(aes(x=region, y=F1, fill=exp, alpha=eval, group=exp)) +
  geom_bar(stat='identity', position=position_dodge(), color='black') +
  facet_grid(.~type, scales='free', space='free') +
  scale_fill_brewer(palette='Paired', name='SV graph - mapper') +
  scale_alpha_manual(name='SV evaluation', values=c(.5,1)) +
  scale_y_continuous(breaks=seq(0,1,.1), limits=0:1) + 
  theme_bw() + 
  guides(fill=guide_legend(keyheight=1))
ggp$map
```

![](sv-evaluation_files/figure-gfm/map-1.png)<!-- -->

``` r
ggp$map.prec = pr.df %>%
  mutate(exp=paste(method, exp, sample),
         exp=factor(exp,
                    levels=c('hgsvc_graph_map hgsvc HG00514',
                             'hgsvc_graph hgsvc HG00514',
                             'giab_graph_map giab6_hg37 HG002',
                                    'giab_graph giab6_hg37 HG002'),
                    labels=c('HGSVC graph - vg map',
                             'HGSVC graph - vg giraffe',
                             'GIAB graph - vg map',
                             'GIAB graph - vg giraffe')),
         type=factor(type, levels=c('Total', 'INS', 'DEL', 'INV')),
         region=factor(region, levels=c('all','highconf'), labels=c('all', 'high-confidence')),
         eval=factor(eval, levels=c('call','geno'), labels=c('presence', 'genotype'))) %>%
  filter(!is.na(exp)) %>% 
  ggplot(aes(x=region, y=precision, fill=exp, alpha=eval, group=exp)) +
  geom_bar(stat='identity', position=position_dodge(), color='black') +
  facet_grid(.~type, scales='free', space='free') +
  scale_fill_brewer(palette='Paired', name='SV graph - mapper') +
  scale_alpha_manual(name='SV evaluation', values=c(.5,1)) +
  scale_y_continuous(breaks=seq(0,1,.1), limits=0:1) + 
  theme_bw() + 
  guides(fill=guide_legend(keyheight=1))
ggp$map.prec
```

![](sv-evaluation_files/figure-gfm/map-2.png)<!-- -->

``` r
ggp$map.rec = pr.df %>%
  mutate(exp=paste(method, exp, sample),
         exp=factor(exp,
                    levels=c('hgsvc_graph_map hgsvc HG00514',
                             'hgsvc_graph hgsvc HG00514',
                             'giab_graph_map giab6_hg37 HG002',
                                    'giab_graph giab6_hg37 HG002'),
                    labels=c('HGSVC graph - vg map',
                             'HGSVC graph - vg giraffe',
                             'GIAB graph - vg map',
                             'GIAB graph - vg giraffe')),
         type=factor(type, levels=c('Total', 'INS', 'DEL', 'INV')),
         region=factor(region, levels=c('all','highconf'), labels=c('all', 'high-confidence')),
         eval=factor(eval, levels=c('call','geno'), labels=c('presence', 'genotype'))) %>%
  filter(!is.na(exp)) %>% 
  ggplot(aes(x=region, y=recall, fill=exp, alpha=eval, group=exp)) +
  geom_bar(stat='identity', position=position_dodge(), color='black') +
  facet_grid(.~type, scales='free', space='free') +
  scale_fill_brewer(palette='Paired', name='SV graph - mapper') +
  scale_alpha_manual(name='SV evaluation', values=c(.5,1)) +
  scale_y_continuous(breaks=seq(0,1,.1), limits=0:1) + 
  theme_bw() + 
  guides(fill=guide_legend(keyheight=1))
ggp$map.rec
```

![](sv-evaluation_files/figure-gfm/map-3.png)<!-- -->

## Individual graph vs combined graph

``` r
## aggregate across samples (hgsvc has 3 samples)
ggp$graphs = pr.df %>%
  mutate(exp=paste(method, exp),
         exp=factor(exp,
                    levels=c('hgsvc_graph hgsvc',
                             'hsvlr_graph hgsvc',
                             'giab_graph giab6_hg37',
                             'hsvlr_graph giab6_hg38'),
                    labels=c('graph: HGSVC\ntruthset: HGSVC',
                             'graph: HGSVC+SVPOP+GIAB\ntruthset: HGSVC',
                             'graph: GIAB\ntruthset: GIAB',
                             'graph: HGSVC+SVPOP+GIAB\ntruthset GIAB')),
         type=factor(type, levels=c('Total', 'INS', 'DEL', 'INV')),
         region=factor(region, levels=c('all','highconf'), labels=c('all', 'high-confidence')),
         eval=factor(eval, levels=c('call','geno'), labels=c('presence', 'genotype'))) %>%
  filter(!is.na(exp)) %>% 
  group_by(exp, type, qual, method, region, eval) %>%
  select(TP, TP.baseline, FN, FP) %>% summarize_all(sum) %>%
  prf %>% 
  ggplot(aes(x=region, y=F1, fill=exp, alpha=eval, group=exp)) +
  geom_bar(stat='identity', position=position_dodge(), color='black') +
  facet_grid(.~type, scales='free', space='free') +
  scale_fill_brewer(palette='Paired', name='SV graph / truthset') +
  scale_alpha_manual(name='SV evaluation', values=c(.5,1)) +
  theme_bw() + 
  scale_y_continuous(breaks=seq(0,1,.1), limits=0:1) + 
  guides(fill=guide_legend(keyheight=2))
ggp$graphs
```

![](sv-evaluation_files/figure-gfm/graphs-1.png)<!-- -->

``` r
ggp$graphs.prec = pr.df %>%
  mutate(exp=paste(method, exp),
         exp=factor(exp,
                    levels=c('hgsvc_graph hgsvc',
                             'hsvlr_graph hgsvc',
                             'giab_graph giab6_hg37',
                             'hsvlr_graph giab6_hg38'),
                    labels=c('graph: HGSVC\ntruthset: HGSVC',
                             'graph: HGSVC+SVPOP+GIAB\ntruthset: HGSVC',
                             'graph: GIAB\ntruthset: GIAB',
                             'graph: HGSVC+SVPOP+GIAB\ntruthset GIAB')),
         type=factor(type, levels=c('Total', 'INS', 'DEL', 'INV')),
         region=factor(region, levels=c('all','highconf'), labels=c('all', 'high-confidence')),
         eval=factor(eval, levels=c('call','geno'), labels=c('presence', 'genotype'))) %>%
  filter(!is.na(exp)) %>% 
  group_by(exp, type, qual, method, region, eval) %>%
  select(TP, TP.baseline, FN, FP) %>% summarize_all(sum) %>%
  prf %>% 
  ggplot(aes(x=region, y=precision, fill=exp, alpha=eval, group=exp)) +
  geom_bar(stat='identity', position=position_dodge(), color='black') +
  facet_grid(.~type, scales='free', space='free') +
  scale_fill_brewer(palette='Paired', name='SV graph / truthset') +
  scale_alpha_manual(name='SV evaluation', values=c(.5,1)) +
  theme_bw() + 
  scale_y_continuous(breaks=seq(0,1,.1), limits=0:1) + 
  guides(fill=guide_legend(keyheight=2))
ggp$graphs.prec
```

![](sv-evaluation_files/figure-gfm/graphs-2.png)<!-- -->

``` r
ggp$graphs.rec = pr.df %>%
  mutate(exp=paste(method, exp),
         exp=factor(exp,
                    levels=c('hgsvc_graph hgsvc',
                             'hsvlr_graph hgsvc',
                             'giab_graph giab6_hg37',
                             'hsvlr_graph giab6_hg38'),
                    labels=c('graph: HGSVC\ntruthset: HGSVC',
                             'graph: HGSVC+SVPOP+GIAB\ntruthset: HGSVC',
                             'graph: GIAB\ntruthset: GIAB',
                             'graph: HGSVC+SVPOP+GIAB\ntruthset GIAB')),
         type=factor(type, levels=c('Total', 'INS', 'DEL', 'INV')),
         region=factor(region, levels=c('all','highconf'), labels=c('all', 'high-confidence')),
         eval=factor(eval, levels=c('call','geno'), labels=c('presence', 'genotype'))) %>%
  filter(!is.na(exp)) %>% 
  group_by(exp, type, qual, method, region, eval) %>%
  select(TP, TP.baseline, FN, FP) %>% summarize_all(sum) %>%
  prf %>% 
  ggplot(aes(x=region, y=recall, fill=exp, alpha=eval, group=exp)) +
  geom_bar(stat='identity', position=position_dodge(), color='black') +
  facet_grid(.~type, scales='free', space='free') +
  scale_fill_brewer(palette='Paired', name='SV graph / truthset') +
  scale_alpha_manual(name='SV evaluation', values=c(.5,1)) +
  theme_bw() + 
  scale_y_continuous(breaks=seq(0,1,.1), limits=0:1) + 
  guides(fill=guide_legend(keyheight=2))
ggp$graphs.rec
```

![](sv-evaluation_files/figure-gfm/graphs-3.png)<!-- -->

## 20x vs 40x

``` r
## here we just compare 
ggp$depth = pr.df %>%
  mutate(exp=paste(method, exp, sample),
         exp=factor(exp,
                    levels=c('hsvlr_graph hgsvc HG00514',
                             'hsvlr_graph_40x hgsvc HG00514_40x'),
                    labels=c('20x',
                             '40x')),
         type=factor(type, levels=c('Total', 'INS', 'DEL', 'INV')),
         region=factor(region, levels=c('all','highconf'), labels=c('all', 'high-confidence')),
         eval=factor(eval, levels=c('call','geno'), labels=c('presence', 'genotype'))) %>%
  filter(!is.na(exp)) %>% 
  ggplot(aes(x=region, y=F1, fill=exp, alpha=eval, group=exp)) +
  geom_bar(stat='identity', position=position_dodge(), color='black') +
  facet_grid(.~type, scales='free', space='free') +
  scale_fill_brewer(palette='Paired', name='depth') +
  scale_alpha_manual(name='SV evaluation', values=c(.5,1)) +
  scale_y_continuous(breaks=seq(0,1,.1), limits=0:1) + 
  theme_bw() + 
  guides(fill=guide_legend(keyheight=1))
ggp$depth
```

![](sv-evaluation_files/figure-gfm/depth-1.png)<!-- -->

``` r
ggp$depth.prec = pr.df %>%
  mutate(exp=paste(method, exp, sample),
         exp=factor(exp,
                    levels=c('hsvlr_graph hgsvc HG00514',
                             'hsvlr_graph_40x hgsvc HG00514_40x'),
                    labels=c('20x',
                             '40x')),
         type=factor(type, levels=c('Total', 'INS', 'DEL', 'INV')),
         region=factor(region, levels=c('all','highconf'), labels=c('all', 'high-confidence')),
         eval=factor(eval, levels=c('call','geno'), labels=c('presence', 'genotype'))) %>%
  filter(!is.na(exp)) %>% 
  ggplot(aes(x=region, y=precision, fill=exp, alpha=eval, group=exp)) +
  geom_bar(stat='identity', position=position_dodge(), color='black') +
  facet_grid(.~type, scales='free', space='free') +
  scale_fill_brewer(palette='Paired', name='depth') +
  scale_alpha_manual(name='SV evaluation', values=c(.5,1)) +
  scale_y_continuous(breaks=seq(0,1,.1), limits=0:1) + 
  theme_bw() + 
  guides(fill=guide_legend(keyheight=1))
ggp$depth.prec
```

![](sv-evaluation_files/figure-gfm/depth-2.png)<!-- -->

``` r
ggp$depth.rec = pr.df %>%
  mutate(exp=paste(method, exp, sample),
         exp=factor(exp,
                    levels=c('hsvlr_graph hgsvc HG00514',
                             'hsvlr_graph_40x hgsvc HG00514_40x'),
                    labels=c('20x',
                             '40x')),
         type=factor(type, levels=c('Total', 'INS', 'DEL', 'INV')),
         region=factor(region, levels=c('all','highconf'), labels=c('all', 'high-confidence')),
         eval=factor(eval, levels=c('call','geno'), labels=c('presence', 'genotype'))) %>%
  filter(!is.na(exp)) %>% 
  ggplot(aes(x=region, y=recall, fill=exp, alpha=eval, group=exp)) +
  geom_bar(stat='identity', position=position_dodge(), color='black') +
  facet_grid(.~type, scales='free', space='free') +
  scale_fill_brewer(palette='Paired', name='depth') +
  scale_alpha_manual(name='SV evaluation', values=c(.5,1)) +
  scale_y_continuous(breaks=seq(0,1,.1), limits=0:1) + 
  theme_bw() + 
  guides(fill=guide_legend(keyheight=1))
ggp$depth.rec
```

![](sv-evaluation_files/figure-gfm/depth-3.png)<!-- -->

## Comparison with GraphTyper

[GraphTyper](https://github.com/DecodeGenetics/graphtyper) was not
published yet when we benchmarked vg and different existing SV
genotypers in [Hickey et
al. 2020](https://doi.org/10.1186/s13059-020-1941-7). We genotyped SVs
with GraphTyper on the same benchmarking dataset: - 3 samples from HGSVC
- simulated reads for one individual from HGSVC - HG002 from the GIAB
dataset.

Of note, GraphTyper’s output includes multiple entries for the same
variant. That leads to a drop in genotyping accuracy due to
“over-genotyping”. Hence, we also include the benchmarking results
when filtering out these duplicates, which improves the “genotyping”
evaluation.

``` r
pr.vggt.df = rbind(
  read.table('hgsvc-vg_graphtyper-sveval-prcurve.tsv', as.is=TRUE, header=TRUE),
  read.table('hgsvcsim-vg_graphtyper-sveval-prcurve.tsv', as.is=TRUE, header=TRUE),
  read.table('giab-vg_graphtyper-sveval-prcurve.tsv', as.is=TRUE, header=TRUE)
)

pr.vggt.df = pr.vggt.df %>% filter(type!='Total', min.cov==.5, 
                         region%in%c('all', 'highconf')) %>% arrange(qual)

qual.bestf1.vggt = pr.vggt.df %>% group_by(exp, type, qual, method, region, eval) %>%
  select(TP, TP.baseline, FN, FP) %>% summarize_all(sum) %>% prf %>%
  group_by(exp, method, type, region, eval) %>%
  arrange(desc(F1)) %>% do(head(., 1))

gt.pal = c(brewer.pal(3,'Set1')[1], brewer.pal(3, 'Paired')[2:1])

ggp$graphtyper = qual.bestf1.vggt %>%
  mutate(exp=factor(exp,
                    levels=c('hgsvc', 'hgsvcsim', 'giab5'),
                    labels=c('HGSVC graph', 'HGSVC graph\n(simulated reads)', 'GIAB graph')),
         method=factor(method, levels=c('vg', 'graphtyper_nodup', 'graphtyper'),
                       labels=c('vg', 'GraphTyper', 'GraphTyper (all calls)')),
         type=factor(type, levels=c('INS', 'DEL')),
         region=factor(region, levels=c('all','highconf'), labels=c('all', 'high-confidence')),
         eval=factor(eval, levels=c('call','geno'), labels=c('presence', 'genotype'))) %>%
  filter(!is.na(exp)) %>% 
  ggplot(aes(x=region, y=F1, fill=method, alpha=eval, group=method)) +
  geom_bar(stat='identity', position=position_dodge(), color='black') +
  facet_grid(type~exp, scales='free', space='free') +
  scale_fill_manual(values=gt.pal, name='method') +
  scale_alpha_manual(name='SV evaluation', values=c(.5,1)) +
  scale_y_continuous(breaks=seq(0,1,.2), limits=0:1) + 
  theme_bw() + 
  guides(fill=guide_legend(keyheight=1))
ggp$graphtyper
```

![](sv-evaluation_files/figure-gfm/vg_graphtyper-1.png)<!-- -->

``` r
ggp$graphtyper.prec = qual.bestf1.vggt %>%
  mutate(exp=factor(exp,
                    levels=c('hgsvc', 'hgsvcsim', 'giab5'),
                    labels=c('HGSVC graph', 'HGSVC graph\n(simulated reads)', 'GIAB graph')),
         method=factor(method, levels=c('vg', 'graphtyper_nodup', 'graphtyper'),
                       labels=c('vg', 'GraphTyper', 'GraphTyper (all calls)')),
         type=factor(type, levels=c('INS', 'DEL')),
         region=factor(region, levels=c('all','highconf'), labels=c('all', 'high-confidence')),
         eval=factor(eval, levels=c('call','geno'), labels=c('presence', 'genotype'))) %>%
  filter(!is.na(exp)) %>% 
  ggplot(aes(x=region, y=precision, fill=method, alpha=eval, group=method)) +
  geom_bar(stat='identity', position=position_dodge(), color='black') +
  facet_grid(type~exp, scales='free', space='free') +
  scale_fill_manual(values=gt.pal, name='method') +
  scale_alpha_manual(name='SV evaluation', values=c(.5,1)) +
  scale_y_continuous(breaks=seq(0,1,.2), limits=0:1) + 
  theme_bw() + 
  guides(fill=guide_legend(keyheight=1))
ggp$graphtyper.prec
```

![](sv-evaluation_files/figure-gfm/vg_graphtyper-2.png)<!-- -->

``` r
ggp$graphtyper.rec = qual.bestf1.vggt %>%
  mutate(exp=factor(exp,
                    levels=c('hgsvc', 'hgsvcsim', 'giab5'),
                    labels=c('HGSVC graph', 'HGSVC graph\n(simulated reads)', 'GIAB graph')),
         method=factor(method, levels=c('vg', 'graphtyper_nodup', 'graphtyper'),
                       labels=c('vg', 'GraphTyper', 'GraphTyper (all calls)')),
         type=factor(type, levels=c('INS', 'DEL')),
         region=factor(region, levels=c('all','highconf'), labels=c('all', 'high-confidence')),
         eval=factor(eval, levels=c('call','geno'), labels=c('presence', 'genotype'))) %>%
  filter(!is.na(exp)) %>% 
  ggplot(aes(x=region, y=recall, fill=method, alpha=eval, group=method)) +
  geom_bar(stat='identity', position=position_dodge(), color='black') +
  facet_grid(type~exp, scales='free', space='free') +
  scale_fill_manual(values=gt.pal, name='method') +
  scale_alpha_manual(name='SV evaluation', values=c(.5,1)) +
  scale_y_continuous(breaks=seq(0,1,.2), limits=0:1) + 
  theme_bw() + 
  guides(fill=guide_legend(keyheight=1))
ggp$graphtyper.rec
```

![](sv-evaluation_files/figure-gfm/vg_graphtyper-3.png)<!-- -->

## Figures

``` r
## adds a legend title: a), b), etc
plot_list <- function(ggp.l, gg.names=NULL){
  if(is.null(names(ggp.l))) names(ggp.l) = paste0('g', 1:length(ggp.l))
  if(is.null(gg.names)) gg.names = names(ggp.l)
  lapply(1:length(gg.names), function(ii) ggp.l[[gg.names[ii]]] + ggtitle(paste0('(', LETTERS[ii], ')')))
}

## depth comparison
ggp$depth = ggp$depth + guides(fill='legend', alpha='legend')
ggp$depth.prec = ggp$depth.prec + guides(fill=FALSE, alpha=FALSE)
ggp$depth.rec = ggp$depth.rec + guides(fill=FALSE, alpha=FALSE)
grid.arrange(grobs=plot_list(ggp, c('depth', 'depth.prec', 'depth.rec')),
             layout_matrix=matrix(c(1,1,2,3), nrow=2, byrow=TRUE),
             heights=c(10,11))
```

![](sv-evaluation_files/figure-gfm/fig-1.png)<!-- -->

``` r
## graphtyper
ggp$graphtyper = ggp$graphtyper + guides(fill='legend', alpha='legend')
ggp$graphtyper.prec = ggp$graphtyper.prec + guides(fill=FALSE, alpha=FALSE)
ggp$graphtyper.rec = ggp$graphtyper.rec + guides(fill=FALSE, alpha=FALSE)
grid.arrange(grobs=plot_list(ggp, c('graphtyper', 'graphtyper.prec', 'graphtyper.rec')),
             layout_matrix=matrix(c(1,1,2,3), nrow=2, byrow=TRUE),
             heights=c(10,11))
```

![](sv-evaluation_files/figure-gfm/fig-2.png)<!-- -->

``` r
pdf('fig-sveval.pdf', 10, 8)
grid.arrange(grobs=plot_list(ggp, c('depth', 'depth.prec', 'depth.rec')),
             layout_matrix=matrix(c(1,1,2,3), nrow=2, byrow=TRUE),
             heights=c(10,11))
grid.arrange(grobs=plot_list(ggp, c('graphtyper', 'graphtyper.prec', 'graphtyper.rec')),
             layout_matrix=matrix(c(1,1,2,3), nrow=2, byrow=TRUE),
             heights=c(10,11))
dev.off()
```

    ## png 
    ##   2

``` r
lg = as_ggplot(get_legend(ggp$graphs))
lg2 = as_ggplot(get_legend(ggp$map))

pdf('fig-sveval-ind.pdf', 6,4)
lg
ggp$graphs + guides(fill=FALSE, alpha=FALSE)
ggp$graphs.prec + guides(fill=FALSE, alpha=FALSE)
ggp$graphs.rec + guides(fill=FALSE, alpha=FALSE)
lg2
ggp$map + guides(fill=FALSE, alpha=FALSE)
ggp$map.prec + guides(fill=FALSE, alpha=FALSE)
ggp$map.rec + guides(fill=FALSE, alpha=FALSE)
dev.off()
```

    ## png 
    ##   2

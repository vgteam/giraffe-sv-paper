Principal component analysis using SV genotypes in the MESA cohort
================

``` r
library(dplyr)
library(ggplot2)
library(gridExtra)
## color palette with clear separation between consecutive groups
interk <- function(x, k=4){ # Interleaves elements in x
  idx = unlist(lapply(1:k, function(kk) seq(kk, length(x), k)))
  x[idx]
}
pal <- function(n){
  pal = interk(rainbow(n, s=.8), 5)
}
## list of figures
ggp = list()
```

## Read PC results

``` r
## PCs derived from SV genotypes or TOPMed SNVs/indels
pc.df = read.table('mesa2k-pcs-svs-topmed.tsv.gz', as.is=TRUE, header=TRUE)
```

## Color samples based on TOPMed PCs

``` r
pc.tm.d = pc.df %>% select(PC1.topmed, PC2.topmed, PC3.topmed, PC4.topmed) %>% as.matrix %>% dist
hc.tm = hclust(pc.tm.d, method='ward.D')
cl.k = 10
pc.df$topmed.cl = cutree(hc.tm, cl.k)

ggp$topmed = ggplot(pc.df, aes(x=PC1.topmed, y=PC2.topmed, color=factor(topmed.cl))) +
  geom_point(alpha=.8) +
  scale_color_manual(name='cluster', values=pal(cl.k)) + 
  theme_bw()
ggp$topmed
```

![](pca-mesa_files/figure-gfm/topmed_cl-1.png)<!-- -->

``` r
ggp$topmed.3.4 = ggplot(pc.df, aes(x=PC3.topmed, y=PC4.topmed, color=factor(topmed.cl))) +
  geom_point(alpha=.8) +
  scale_color_manual(name='cluster', values=pal(cl.k)) + 
  theme_bw()
ggp$topmed.3.4
```

![](pca-mesa_files/figure-gfm/topmed_cl-2.png)<!-- -->

## PCs from SV genotypes

``` r
ggp$sv = ggplot(pc.df, aes(x=PC1.sv, y=PC2.sv, color=factor(topmed.cl))) +
  geom_point(alpha=.8) +
  scale_color_manual(name='cluster', values=pal(cl.k)) + 
  theme_bw()
ggp$sv
```

![](pca-mesa_files/figure-gfm/sv-1.png)<!-- -->

``` r
ggp$sv.3.4 = ggplot(pc.df, aes(x=PC3.sv, y=PC4.sv, color=factor(topmed.cl))) +
  geom_point(alpha=.8) +
  scale_color_manual(name='cluster', values=pal(cl.k)) + 
  theme_bw()
ggp$sv.3.4
```

![](pca-mesa_files/figure-gfm/sv-2.png)<!-- -->

## Multi-panel figure

``` r
## adds a legend title: a), b), etc
plot_list <- function(ggp.l, gg.names=NULL){
  if(is.null(names(ggp.l))) names(ggp.l) = paste0('g', 1:length(ggp.l))
  if(is.null(gg.names)) gg.names = names(ggp.l)
  lapply(1:length(gg.names), function(ii) ggp.l[[gg.names[ii]]] + ggtitle(paste0(letters[ii], ')')))
}

ggp$sv = ggp$sv + guides(color=FALSE)
ggp$topmed = ggp$topmed + theme(legend.position='top')
grid.arrange(grobs=plot_list(ggp, c('topmed', 'sv')),
             layout_matrix=matrix(c(1,2),1))
```

![](pca-mesa_files/figure-gfm/fig-1.png)<!-- -->

``` r
ggp$sv.3.4 = ggp$sv.3.4 + guides(color=FALSE)
ggp$topmed.3.4 = ggp$topmed.3.4 + theme(legend.position='top')
grid.arrange(grobs=plot_list(ggp, c('topmed.3.4', 'sv.3.4')),
             layout_matrix=matrix(c(1,2),1))
```

![](pca-mesa_files/figure-gfm/fig-2.png)<!-- -->

``` r
pdf('fig-sv-mesa-pcs.pdf', 9, 5)
grid.arrange(grobs=plot_list(ggp, c('topmed', 'sv')),
             layout_matrix=matrix(c(1,2),1))
grid.arrange(grobs=plot_list(ggp, c('topmed.3.4', 'sv.3.4')),
             layout_matrix=matrix(c(1,2),1))
dev.off()
```

    ## png 
    ##   2

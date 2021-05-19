Principal component analysis using SV genotypes in the MESA cohort
================

``` r
library(dplyr)
library(ggplot2)
library(gridExtra)
library(Rtsne)
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

## Pre-computed PCA

Because it takes about 30 mins and the genotypes cannot be shared for
MESA, the principal components were pre-computed on Terra from the
allele counts using the following commands:

``` r
## allele counts
ac = as.matrix(read.table('mesa2k.svsite80al.ac.tsv.gz', as.is=TRUE))
## subset on autosomes
svs = read.table('svs.mesa2k.svsite80al.tsv.gz', header=TRUE, as.is=TRUE)
site.chr = unique(svs[, c('svsite', 'seqnames')])                       
auto.sites = subset(site.chr, seqnames %in% paste0('chr', 1:22))$svsite
ac = ac[auto.sites, ]
## PCA
pca.o = prcomp(t(ac))

## save results
pca.o$x = pca.o$x[,1:20]
pca.o$rotation = pca.o$rotation[,1:4]
saveRDS(pca.o, file='mesa2k.svsite80al.ac.pcs.rds')

## save results as TSVs
pc.df = tibble(sample=rownames(pca.o$x)) %>% cbind(pca.o$x[,1:20])
outf = gzfile('mesa2k-pcs-svs-topmed.tsv.gz', 'w')
write.table(pc.m, file=outf, sep='\t', row.names=FALSE, quote=FALSE)
close(outf)
pc.sdev.df = tibble(pc=1:length(pca.o$sdev), sdev=pca.o$sdev)
write.table(pc.sdev.df, file='mesa2k.svsite80al.ac.pcs.sdev.tsv', sep='\t', quote=FALSE, row.names=FALSE)
```

## Read pre-computed principal components

``` r
## PCs derived from SV genotypes or TOPMed SNVs/indels
pc.df = read.table('mesa2k-pcs-svs-topmed.tsv.gz', as.is=TRUE, header=TRUE)
pc.sdev.df = read.table('mesa2k.svsite80al.ac.pcs.sdev.tsv', as.is=TRUE, header=TRUE)
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
  theme_bw() + coord_fixed()
ggp$topmed
```

![](pca-mesa_files/figure-gfm/topmed_cl-1.png)<!-- -->

``` r
ggp$topmed.3.4 = ggplot(pc.df, aes(x=PC3.topmed, y=PC4.topmed, color=factor(topmed.cl))) +
  geom_point(alpha=.8) +
  scale_color_manual(name='cluster', values=pal(cl.k)) + 
  theme_bw() + coord_fixed()
ggp$topmed.3.4
```

![](pca-mesa_files/figure-gfm/topmed_cl-2.png)<!-- -->

## PCs from SV genotypes

``` r
pc.sdev.df %>% filter(pc<=20) %>% 
  ggplot(aes(pc, sdev)) +
  geom_bar(stat='identity') + 
  theme_bw() +
  ylab('standard deviation') + xlab('principal component')
```

![](pca-mesa_files/figure-gfm/sv-1.png)<!-- -->

``` r
ggp$sv = ggplot(pc.df, aes(x=PC1, y=PC2, color=factor(topmed.cl))) +
  geom_point(alpha=.8) +
  scale_color_manual(name='cluster', values=pal(cl.k)) + 
  theme_bw() + coord_fixed()
ggp$sv
```

![](pca-mesa_files/figure-gfm/sv-2.png)<!-- -->

``` r
ggp$sv.3.4 = ggplot(pc.df, aes(x=PC3, y=PC4, color=factor(topmed.cl))) +
  geom_point(alpha=.8) +
  scale_color_manual(name='cluster', values=pal(cl.k)) + 
  theme_bw() + coord_fixed()
ggp$sv.3.4
```

![](pca-mesa_files/figure-gfm/sv-3.png)<!-- -->

## tSNE from SV genotypes

``` r
tsne.o = Rtsne(pc.df[,2:21])

tsne.df = tibble(sample=pc.df$sample, tsne1=tsne.o$Y[,1], tsne2=tsne.o$Y[,2])

ggp$svtsne = ggplot(tsne.df, aes(x=tsne1, y=tsne2)) +
  geom_point(alpha=.5) +
  theme_bw()
ggp$svtsne
```

![](pca-mesa_files/figure-gfm/svtsne-1.png)<!-- -->

## Direct PC comparison

``` r
corLabel <- function(x, y, position='top'){
  ## place label on the the right side
  l.x = range(x)
  l.x = l.x[1] + .8 * diff(l.x)
  ## place label vertically
  l.y = range(y)
  l.y = l.y[1] + ifelse(position=='top', .9, .1) * diff(l.y)
  ## result
  list(x=l.x, y=l.y, r=cor(x, y))
}

cl = corLabel(pc.df$PC1.topmed, pc.df$PC1)
ggp$pc1 = ggplot(pc.df, aes(x=PC1.topmed, y=PC1)) + geom_point(alpha=.2) + theme_bw() +
  annotate("label", x=cl$x, y=cl$y, label=paste0("rho == ", round(cl$r, 3)), parse=TRUE)
ggp$pc1
```

![](pca-mesa_files/figure-gfm/pc-1.png)<!-- -->

``` r
cl = corLabel(pc.df$PC2.topmed, pc.df$PC2)
ggp$pc2 = ggplot(pc.df, aes(x=PC2.topmed, y=PC2)) + geom_point(alpha=.2) + theme_bw() +
  annotate("label", x=cl$x, y=cl$y, label=paste0("rho == ", round(cl$r, 3)), parse=TRUE)
ggp$pc2
```

![](pca-mesa_files/figure-gfm/pc-2.png)<!-- -->

``` r
cl = corLabel(pc.df$PC3.topmed, pc.df$PC3, position='bottom')
ggp$pc3 = ggplot(pc.df, aes(x=PC3.topmed, y=PC3)) + geom_point(alpha=.2) + theme_bw() +
  annotate("label", x=cl$x, y=cl$y, label=paste0("rho == ", round(cl$r, 3)), parse=TRUE)
ggp$pc3
```

![](pca-mesa_files/figure-gfm/pc-3.png)<!-- -->

``` r
cl = corLabel(pc.df$PC4.topmed, pc.df$PC4)
ggp$pc4 = ggplot(pc.df, aes(x=PC4.topmed, y=PC4)) + geom_point(alpha=.2) + theme_bw() +
  annotate("label", x=cl$x, y=cl$y, label=paste0("rho == ", round(cl$r, 3)), parse=TRUE)
ggp$pc4
```

![](pca-mesa_files/figure-gfm/pc-4.png)<!-- -->

## Multi-panel figure

``` r
## adds a legend title: a), b), etc
plot_list <- function(ggp.l, gg.names=NULL){
  if(is.null(names(ggp.l))) names(ggp.l) = paste0('g', 1:length(ggp.l))
  if(is.null(gg.names)) gg.names = names(ggp.l)
  lapply(1:length(gg.names), function(ii) ggp.l[[gg.names[ii]]] + ggtitle(paste0('(', LETTERS[ii], ')')))
}

ggp$sv.f = ggp$sv + guides(color=FALSE)
ggp$topmed.f = ggp$topmed
grid.arrange(grobs=plot_list(ggp, c('topmed.f', 'sv.f')),
             layout_matrix=matrix(c(1,2),1))
```

![](pca-mesa_files/figure-gfm/fig-1.png)<!-- -->

``` r
ggp$sv.3.4.f = ggp$sv.3.4 + guides(color=FALSE)
ggp$topmed.3.4.f = ggp$topmed.3.4 + theme(legend.position=c(.99, .99), legend.justification=c(1,1)) +
  guides(color=guide_legend(ncol=5))
grid.arrange(grobs=plot_list(ggp, c('topmed.3.4.f', 'sv.3.4.f')),
             layout_matrix=matrix(c(1,2),1))
```

![](pca-mesa_files/figure-gfm/fig-2.png)<!-- -->

``` r
grid.arrange(grobs=plot_list(ggp, c('pc1', 'pc2', 'pc3')),
             layout_matrix=matrix(c(1:3),1))
```

![](pca-mesa_files/figure-gfm/fig-3.png)<!-- -->

``` r
pdf('figs/fig-sv-mesa-pcs.pdf', 9, 4)
grid.arrange(grobs=plot_list(ggp, c('topmed.f', 'sv.f')),
             layout_matrix=matrix(c(1,2),1))
grid.arrange(grobs=plot_list(ggp, c('topmed.3.4.f', 'sv.3.4.f')),
             layout_matrix=matrix(c(1,2),1))
grid.arrange(grobs=plot_list(ggp, c('pc1', 'pc2', 'pc3')),
             layout_matrix=matrix(c(1:3),1))
dev.off()
```

    ## png 
    ##   2

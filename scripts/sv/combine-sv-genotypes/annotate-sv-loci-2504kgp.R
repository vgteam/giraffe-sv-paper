library(sveval)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(igraph)
library(Matrix)

is_clique <- function(amat){
  diag(amat) = 1
  all(amat==1)
}

## Function to overlap and annotate SVs
## Used on slices of the full dataset to reduce memory consumption
annotateSVs <- function(gr, min.rol=.8, sr=NULL){
  ## overlaps SVs in the set
  ol.df = svOverlap(gr, gr, min.ol=min.rol, min.del.rol=min.rol,
                    range.seq.comp=TRUE, ins.seq.comp=TRUE,
                    simprep=sr) %>%
    as.data.frame %>%
    filter(queryHits!=subjectHits)
  ## make graph
  adj.mat = sparseMatrix(ol.df$queryHits, ol.df$subjectHits, x=rep(1, nrow(ol.df)), dims=rep(length(gr), 2))
  rownames(adj.mat) = colnames(adj.mat) = gr$svid
  ## extract components
  cmp.o = components(graph_from_adjacency_matrix(adj.mat, mode='undirected'))
  ## check is components are cliques
  site.df = lapply(unique(as.numeric(cmp.o$membership)), function(cmp){
    cmp.ii = which(cmp.o$membership == cmp)
    if(length(cmp.ii) < 3){
      return(tibble(svid=rownames(adj.mat)[cmp.ii], clique=TRUE, svsite=rownames(adj.mat)[cmp.ii[1]]))
    } else {
      return(tibble(svid=rownames(adj.mat)[cmp.ii], clique=is_clique(adj.mat[cmp.ii, cmp.ii]), svsite=rownames(adj.mat)[cmp.ii[1]]))
    }
  }) %>% bind_rows
  ## merge annotation
  gr %>% as.data.frame %>%
    merge(site.df) %>% 
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)  
}
  
## read SVs and make GRanges object
svs = read.table('svs.2504kgp.seq.tsv.gz', header=TRUE, as.is=TRUE)
svs = makeGRangesFromDataFrame(svs, keep.extra.columns = TRUE)
length(svs)

## simple repeat annotation
## download.file('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz', 'simpleRepeat.hg38.txt.gz')
sr = read.table('simpleRepeat.hg38.txt.gz', as.is=TRUE)
sr = reduce(GRanges(sr$V2, IRanges(sr$V3, sr$V4)))

## quick clustering to split dataset for analysis
cl.gr = reduce(svs, min.gapwidth=100)
cl.gr$n = countOverlaps(cl.gr, svs)
length(cl.gr)
max(cl.gr$n)
## qplot(x=cl.gr$n) + scale_x_log10()

## make batches of the SV clusters to have ~500 variants
## (or more if in one cluster that can't be split)
cl.l = subset(cl.gr, n>10)
batch.maxsize = 5000
cl.l$batch = cumsum(cl.l$n) %>% cut(seq(0,sum(cl.l$n)+batch.maxsize, batch.maxsize))
length(unique(cl.l$batch))

## merge each batch of large SV cluster
date()
svs.m.l = mclapply(unique(cl.l$batch), function(bb){
  ## message('Starting ', bb)
  cl.l = subset(cl.l, batch==bb)
  gr = subsetByOverlaps(svs, cl.l)
  gg = annotateSVs(gr)  
  ## message('Finished ', bb)
  gg
}, mc.cores=16)
date()

## bind the results
svs.m.l = do.call(c, svs.m.l)

## merge all isolated SVs together in one run
date()
svs.m.s = annotateSVs(subsetByOverlaps(svs, subset(cl.gr, n<=10)))
date()

## pool the results of the small and large SV clusters
svs.m = c(svs.m.l, svs.m.s)

## number of merged SVs vs number of total SVs
cat('Total variants\n')
length(svs.m)
cat('Total SV loci\n')
length(unique(svs.m$svsite))
cat('Total clique SV loci\n')
length(unique(subset(svs.m, clique)$svsite))

## 
outgz = gzfile('svs.2504kgp.svsite80al.tsv.gz', 'w')
svs.m %>% as.data.frame %>% select(-strand, -width) %>% write.table(file=outgz, sep='\t', row.names=FALSE, quote=FALSE)
close(outgz)


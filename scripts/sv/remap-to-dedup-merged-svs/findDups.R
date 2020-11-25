library(sveval)
library(VariantAnnotation)
library(GenomicRanges)
library(dplyr)
library(igraph)

## HSVLR catalog
svs = readSVvcf('../merge-svs/hsvlr.vcf.gz', vcf.object=TRUE)
names(svs) = paste0('sv', 1:length(svs))

## Make GRanges object
svs.gr = rowRanges(svs)
svs.gr$size = info(svs)$SIZE
svs.gr$type = info(svs)$SVTYPE
end(svs.gr) = info(svs)$END

## Overlap parameters
min.rol = .8      # 80% reciprocal overlap
max.ins.gap = 50  # insertions at max 50 bp from each other

## Overlap insertions together
ins.gr = svs.gr[which(svs.gr$type=='INS')]
if(length(ins.gr)>0){
  ## Cluster insertions
  ol.ins = findOverlaps(ins.gr, ins.gr, maxgap=max.ins.gap) %>%
    as.data.frame %>% filter(queryHits<subjectHits) %>% 
    mutate(qs=ins.gr$size[queryHits], ss=ins.gr$size[subjectHits],
           qid=names(ins.gr)[queryHits], sid=names(ins.gr)[subjectHits],
           rol=ifelse(qs>ss, ss/qs, qs/ss)) %>%
    select(-queryHits, -subjectHits) %>% 
    filter(rol > min.rol)
}

## Overlap deletions together
del.gr = svs.gr[which(svs.gr$type=='DEL')]
if(length(del.gr)>0){
  ## Cluster deletions
  ol.del = findOverlaps(del.gr, del.gr) %>%
    as.data.frame %>% filter(queryHits<subjectHits) %>% 
    mutate(qs=del.gr$size[queryHits], ss=del.gr$size[subjectHits],
           qss=width(pintersect(del.gr[queryHits], del.gr[subjectHits])),
           qid=names(del.gr)[queryHits], sid=names(del.gr)[subjectHits],
           rol=ifelse(qs>ss, qss/qs, qss/ss)) %>%
    select(-queryHits, -subjectHits, -qss) %>% 
    filter(rol > min.rol)
}

## Merge overlap stats for insertions and deletions
ol.df = rbind(ol.ins, ol.del)

## Find connected components
ol.g = ol.df %>% select(qid, sid) %>% as.matrix %>% graph_from_edgelist(directed=FALSE)
ol.c = components(ol.g)
ol.c.df = tibble(id=names(ol.c$membership), cmp=as.numeric(ol.c$membership))

## keep only autosomes and sex chr
svs = svs[which(as.character(seqnames(svs)) %in% paste0('chr', c(1:22, 'X','Y')))]

## Add coordinates
svs.df = tibble(seqnames=as.character(seqnames(svs)), start=start(svs), end=end(svs))
svs.df$id = names(svs)
cl.df = merge(svs.df, ol.c.df) %>% group_by(cmp, seqnames) %>%
  summarize(start=min(start)-1000, end=max(end)+1000, ids=paste(id, collapse=','))

## Write coordinates and ids for each cluster of near-dups
write.table(cl.df, file='neardups-clusters.tsv', quote=FALSE, sep='\t', row.names=FALSE)

## Write a VCF with new SV ids for each chromosome
for(chr in unique(as.character(seqnames(svs)))){
  writeVcf(svs[which(as.character(seqnames(svs))==chr)], file=paste0('hsvlr-fordedup-', chr, '.vcf'))
}

sv.info = tibble(id=names(svs), size=info(svs)$SIZE)

library(dplyr)

idx = read.table('1000G_2504_high_coverage.sequence.index', as.is=TRUE, sep='\t')
idx = idx[,c(1,11)]
colnames(idx) = c('cram', 'population')

table(idx$population)

pop.info = read.table('20131219.populations.tsv', as.is=TRUE, header=TRUE, sep='\t') %>%
  mutate(population=Population.Code, super.population=Super.Population) %>% 
  select(population, super.population) %>%
  filter(population !='')

idx.sp = idx %>% merge(pop.info) %>% group_by(super.population) %>% sample_n(1)
write.table(idx.sp, file='1kgp-hg38-1persuperpop.tsv', sep='\t', row.names=FALSE, quote=FALSE)

## second batch from different population, 2 per super population
batch1  = read.table('1kgp-hg38-1persuperpop.tsv', as.is=TRUE, header=TRUE)
idx.sp = idx %>% merge(pop.info) %>% filter(!(population %in% batch1$population)) %>% group_by(super.population) %>% sample_n(2)
write.table(idx.sp, file='1kgp-hg38-1persuperpop-batch2.tsv', sep='\t', row.names=FALSE, quote=FALSE)

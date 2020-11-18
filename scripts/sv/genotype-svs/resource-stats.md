Summary of the resources used for SV genotyping
================

Extracted from Terra using
[terra-notebook-utils](https://github.com/DataBiosphere/terra-notebook-utils).
At the time of running this, the module for cost estimation is still
experimental. Information from some of the jobs might be missing, for
example if it needed re-running or was too old. There was also an issue
on Terra while we were genotyping the 1000GP dataset that might have
affected the reporting of these metrics. Still, the vast majority of the
submissions will have complete information that we can use to estimate
the average resource requirement.

``` r
library(dplyr)
library(ggplot2)
library(knitr)
```

## Read resource metrics for both MESA and 1000GP analysis

``` r
res.df = rbind(
  read.table('./terra-mesa-subs-resources.tsv.gz', header=TRUE, as.is=TRUE) %>% mutate(dataset='MESA'),
  read.table('./terra-1kgp-subs-resources.tsv.gz', header=TRUE, as.is=TRUE) %>% mutate(dataset='1000GP')
)
res.df = res.df %>% mutate(core.hour=cpu*duration/3600)
sample_n(res.df, 5)
```

    ##         batch                             workflow shard cpu mem  duration
    ## 1 batch_5_372 e2129e2f-b933-495a-a6f4-71025ce55f77     3  32 100 13787.018
    ## 2 batch_3_500 6c77041d-ad36-4dd9-bf7d-c50154f474b7     4   8  50  4289.536
    ## 3 batch_4_289 bd7b04e0-5406-4f29-899f-d6699aa41e00     4   8  50  3948.741
    ## 4 batch_3_435 ac517022-b9e1-42b7-ae70-d5dcbdb802b3     2  32 100    39.441
    ## 5 batch_1_597 f4ec3a28-12f2-4f55-96fa-167eb493c5a6     3  16 100  5411.498
    ##          cost dataset   core.hour
    ## 1 1.215488800  1000GP 122.5512711
    ## 2 0.122551000    MESA   9.5323022
    ## 3 0.112809767  1000GP   8.7749800
    ## 4 0.005289333  1000GP   0.3505867
    ## 5 0.309205600  1000GP  24.0511022

## Aggregate per workflow

``` r
res.a = res.df %>% group_by(workflow, dataset) %>%
  summarize(shard=n(), cost=sum(cost), core.hour=sum(core.hour), .groups='drop') %>%
  mutate(preempted=shard-4)
```

## Pre-empted jobs

We used pre-emptible instances for all the steps in the pipeline. If no
jobs were pre-empted, each workflow should contain 4 jobs/shards (CRAM
conversion, mapping two chunks in parallel, SV genotyping). More shards
mean that some were pre-empted. Note: less shards most likely mean that
it was a rerun where cached intermediate results allowed to skip jobs.

``` r
res.a %>% group_by(shard) %>% summarize(workflow=n()) %>%
  ggplot(aes(x=shard, y=workflow)) + geom_bar(stat='identity') + theme_bw() +
  scale_x_continuous(breaks=1:100)
```

![](resource-stats_files/figure-gfm/nshards_sample-1.png)<!-- -->

``` r
res.a %>% summarize(mean.preempted.job=mean(preempted), prop.notpreempted=mean(preempted==0)) %>%
  kable(digits=3)
```

| mean.preempted.job | prop.notpreempted |
| -----------------: | ----------------: |
|              0.599 |             0.592 |

## Cost and number of core.hours per sample

``` r
ggplot(res.a, aes(x=core.hour, fill=factor(shard))) + geom_histogram() + theme_bw() +
  geom_vline(xintercept=80, linetype=2)
```

![](resource-stats_files/figure-gfm/corehours_underest-1.png)<!-- -->

``` r
ggplot(res.a, aes(x=core.hour, fill=dataset)) + geom_histogram() + theme_bw() +
  geom_vline(xintercept=80, linetype=2)
```

![](resource-stats_files/figure-gfm/corehours_underest-2.png)<!-- -->

The jobs below *x=100* are likely under-reported because of the issue
Terra experienced when we were genotyping the 1000GP dataset. The logs
look normal but and the duration on Terra is 0 for all or some of the
steps in the workflow, leading to these under-estimated resources.

Filtering incomplete or problematic jobs:

``` r
res.a %>% filter(core.hour>80, shard>=4) %>% 
  ggplot(aes(x=core.hour, fill=dataset)) + geom_histogram() + theme_bw()
```

![](resource-stats_files/figure-gfm/corehours-1.png)<!-- -->

``` r
res.a %>% filter(core.hour>80, shard>=4) %>% 
  ggplot(aes(x=preempted, y=core.hour, group=paste(dataset, preempted), fill=dataset)) +
  geom_boxplot() + theme_bw() +
  scale_x_continuous(breaks=0:100)
```

![](resource-stats_files/figure-gfm/corehours-2.png)<!-- -->

``` r
res.a %>% filter(core.hour>80, preempted>=0) %>%
  summarize(mean.preempted.job=mean(preempted), prop.notpreempted=mean(preempted==0)) %>%
  kable(digits=3)
```

| mean.preempted.job | prop.notpreempted |
| -----------------: | ----------------: |
|              0.599 |             0.592 |

``` r
res.a %>% mutate(dataset='all') %>% rbind(res.a) %>%
  mutate(dataset=factor(dataset, levels=c('all', 'MESA', '1000GP'))) %>% 
  filter(core.hour>80, shard>=4) %>% 
  group_by(preempted, dataset) %>%
  summarize(workflow=n(), cost=mean(cost), core.hour=mean(core.hour)) %>%
  kable(digits=3)
```

| preempted | dataset | workflow |  cost | core.hour |
| --------: | :------ | -------: | ----: | --------: |
|         0 | all     |     2910 | 2.044 |   194.388 |
|         0 | MESA    |     1221 | 2.082 |   197.981 |
|         0 | 1000GP  |     1689 | 2.017 |   191.791 |
|         1 | all     |     1358 | 2.315 |   220.046 |
|         1 | MESA    |      517 | 2.346 |   222.205 |
|         1 | 1000GP  |      841 | 2.295 |   218.718 |
|         2 | all     |      439 | 2.545 |   240.992 |
|         2 | MESA    |      172 | 2.638 |   248.246 |
|         2 | 1000GP  |      267 | 2.485 |   236.318 |
|         3 | all     |      154 | 2.790 |   263.533 |
|         3 | MESA    |       62 | 3.019 |   283.251 |
|         3 | 1000GP  |       92 | 2.636 |   250.244 |
|         4 | all     |       36 | 3.028 |   285.260 |
|         4 | MESA    |       16 | 2.986 |   277.209 |
|         4 | 1000GP  |       20 | 3.062 |   291.701 |
|         5 | all     |       12 | 3.259 |   304.873 |
|         5 | MESA    |        7 | 3.088 |   282.514 |
|         5 | 1000GP  |        5 | 3.500 |   336.177 |
|         6 | all     |        6 | 2.841 |   263.576 |
|         6 | MESA    |        3 | 3.105 |   280.710 |
|         6 | 1000GP  |        3 | 2.577 |   246.442 |
|         8 | all     |        1 | 4.389 |   410.061 |
|         8 | 1000GP  |        1 | 4.389 |   410.061 |

As expected, when a job is pre-empted, the total cost and core.hours for
this sample tend to be higher.

## Resources per task

The SV genotyping pipeline has three steps: converting the CRAM file to
FASTQ (and chunking), mapping each chunk, genotyping SVs from the
aligned reads.

``` r
## task name is not extracted yet but we know the resource requested for each of the three tasks
## -> check that there are only 3 requested resource profiles
res.df %>% select(cpu, mem) %>% unique %>% kable
```

|   | cpu | mem |
| - | --: | --: |
| 1 |  32 | 100 |
| 5 |  16 | 100 |
| 8 |   8 |  50 |

``` r
res.df %>% mutate(dataset='all') %>% rbind(res.df) %>%
  mutate(dataset=factor(dataset, levels=c('all', 'MESA', '1000GP')),
         task=ifelse(cpu==32, 'mapping', 'CRAM conversion'),
         task=ifelse(cpu==16, 'genotyping', task),
         task=factor(task, levels=c('CRAM conversion', 'mapping', 'genotyping'))) %>%
  group_by(dataset, workflow) %>%
  mutate(preempted=n()-4) %>%
  filter(preempted==0) %>%
  group_by(workflow, task, cpu, mem, dataset) %>%
  summarize(cost=sum(cost), core.hour=sum(core.hour), .groups='drop') %>%
  group_by(dataset, workflow) %>%
  mutate(prop.cost=cost/sum(cost), prop.core.hour=core.hour/sum(core.hour)) %>% 
  ungroup() %>% 
  select(task, cpu, mem, dataset, core.hour, prop.core.hour, cost, prop.cost) %>% 
  group_by(task, dataset) %>%
  summarize_all(mean) %>%
  kable(digits=3)
```

| task            | dataset | cpu | mem | core.hour | prop.core.hour |  cost | prop.cost |
| :-------------- | :------ | --: | --: | --------: | -------------: | ----: | --------: |
| CRAM conversion | all     |   8 |  50 |    12.867 |          0.076 | 0.166 |     0.092 |
| CRAM conversion | MESA    |   8 |  50 |    13.220 |          0.076 | 0.170 |     0.090 |
| CRAM conversion | 1000GP  |   8 |  50 |    12.636 |          0.076 | 0.163 |     0.093 |
| mapping         | all     |  32 | 100 |   146.261 |          0.752 | 1.451 |     0.710 |
| mapping         | MESA    |  32 | 100 |   157.660 |          0.771 | 1.564 |     0.725 |
| mapping         | 1000GP  |  32 | 100 |   138.782 |          0.739 | 1.377 |     0.700 |
| genotyping      | all     |  16 | 100 |    24.891 |          0.172 | 0.320 |     0.198 |
| genotyping      | MESA    |  16 | 100 |    27.112 |          0.154 | 0.349 |     0.185 |
| genotyping      | 1000GP  |  16 | 100 |    23.433 |          0.185 | 0.301 |     0.207 |

## Save some tables in LaTeX format

``` r
## cost/core.hour per sample
res.a %>% mutate(dataset='all') %>% rbind(res.a) %>%
  mutate(dataset=factor(dataset, levels=c('all', 'MESA', '1000GP'))) %>% 
  filter(core.hour>80, shard>=4) %>% 
  group_by(preempted, dataset) %>%
  summarize(workflow=n(), cost=mean(cost), core.hour=mean(core.hour)) %>%
  kable(digits=3, format='latex', caption='resource per sample') %>%
  cat(file='resource-stats.tex')

cat('\n\n', file='resource-stats.tex', append=TRUE)

## resource per task

res.df %>% mutate(dataset='all') %>% rbind(res.df) %>%
  mutate(dataset=factor(dataset, levels=c('all', 'MESA', '1000GP')),
         task=ifelse(cpu==32, 'mapping', 'CRAM conversion'),
         task=ifelse(cpu==16, 'genotyping', task),
         task=factor(task, levels=c('CRAM conversion', 'mapping', 'genotyping'))) %>%
  group_by(dataset, workflow) %>%
  mutate(preempted=n()-4) %>%
  filter(preempted==0) %>%
  group_by(workflow, task, cpu, mem, dataset) %>%
  summarize(cost=sum(cost), core.hour=sum(core.hour), .groups='drop') %>%
  group_by(dataset, workflow) %>%
  mutate(prop.cost=cost/sum(cost), prop.core.hour=core.hour/sum(core.hour)) %>% 
  ungroup() %>% 
  select(task, cpu, mem, dataset, core.hour, prop.core.hour, cost, prop.cost) %>% 
  group_by(task, dataset) %>%
  summarize_all(mean) %>%
  kable(digits=3, format='latex', caption='resource per task') %>%
  cat(file='resource-stats.tex', append=TRUE)
```

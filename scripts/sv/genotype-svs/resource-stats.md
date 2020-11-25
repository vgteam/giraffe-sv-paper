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

    ##         batch                             workflow shard cpu mem duration
    ## 1 batch_2_500 ccd0c5b2-c312-4925-bb76-2f462cc790f8     1  32 100 1433.384
    ## 2 batch_4_500 6448df41-e804-44f1-a2ff-9782bc4ef8fa     3  16 100 7088.285
    ## 3 batch_3_435 55c23703-4106-462c-8861-ee4f96b53179     3  16 100 6959.052
    ## 4 batch_1_597 095b0ed7-aa04-4ffd-8ec5-558180a34c8d     3  16 100 5168.621
    ## 5 batch_3_500 83755c89-3dda-44fd-b7f6-e1e704cbf5fe     1  32 100 1453.564
    ##        cost dataset core.hour
    ## 1 0.1264151    MESA  12.74119
    ## 2 0.4050182    MESA  31.50349
    ## 3 0.3976480  1000GP  30.92912
    ## 4 0.2953222  1000GP  22.97165
    ## 5 0.1281782    MESA  12.92057

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
  mutate(preempted=cut(preempted, c(-1,0,1,2,Inf), labels=c(0,1,2,'3+'))) %>% 
  group_by(preempted, dataset) %>%
  summarize(workflow=n(), cost=mean(cost), core.hour=mean(core.hour)) %>%
  kable(digits=3)
```

| preempted | dataset | workflow |  cost | core.hour |
| :-------- | :------ | -------: | ----: | --------: |
| 0         | all     |     2910 | 2.044 |   194.388 |
| 0         | MESA    |     1221 | 2.082 |   197.981 |
| 0         | 1000GP  |     1689 | 2.017 |   191.791 |
| 1         | all     |     1358 | 2.315 |   220.046 |
| 1         | MESA    |      517 | 2.346 |   222.205 |
| 1         | 1000GP  |      841 | 2.295 |   218.718 |
| 2         | all     |      439 | 2.545 |   240.992 |
| 2         | MESA    |      172 | 2.638 |   248.246 |
| 2         | 1000GP  |      267 | 2.485 |   236.318 |
| 3+        | all     |      209 | 2.867 |   270.351 |
| 3+        | MESA    |       88 | 3.022 |   282.007 |
| 3+        | 1000GP  |      121 | 2.755 |   261.874 |

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
  mutate(core.hour=paste0(round(core.hour, 3), ' (', round(prop.core.hour, 3)*100, '%)'),
         cost=paste0(round(cost, 3), ' (', round(prop.cost, 3)*100, '%)')) %>%
  select(task, cpu, mem, dataset, core.hour, cost) %>% 
    kable(digits=3)
```

| task            | cpu | mem | dataset | core.hour       | cost          |
| :-------------- | --: | --: | :------ | :-------------- | :------------ |
| CRAM conversion |   8 |  50 | all     | 12.867 (7.6%)   | 0.166 (9.2%)  |
| CRAM conversion |   8 |  50 | MESA    | 13.22 (7.6%)    | 0.17 (9%)     |
| CRAM conversion |   8 |  50 | 1000GP  | 12.636 (7.6%)   | 0.163 (9.3%)  |
| mapping         |  32 | 100 | all     | 146.261 (75.2%) | 1.451 (71%)   |
| mapping         |  32 | 100 | MESA    | 157.66 (77.1%)  | 1.564 (72.5%) |
| mapping         |  32 | 100 | 1000GP  | 138.782 (73.9%) | 1.377 (70%)   |
| genotyping      |  16 | 100 | all     | 24.891 (17.2%)  | 0.32 (19.8%)  |
| genotyping      |  16 | 100 | MESA    | 27.112 (15.4%)  | 0.349 (18.5%) |
| genotyping      |  16 | 100 | 1000GP  | 23.433 (18.5%)  | 0.301 (20.7%) |

## Save some tables in LaTeX format

``` r
## cost/core.hour per sample
res.a %>% mutate(dataset='all') %>% rbind(res.a) %>%
  mutate(dataset=factor(dataset, levels=c('all', 'MESA', '1000GP'))) %>% 
  filter(core.hour>80, shard>=4) %>%
  mutate(preempted=cut(preempted, c(-1,0,1,2,Inf), labels=c(0,1,2,'3+'))) %>% 
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
  mutate(core.hour=paste0(round(core.hour, 3), ' (', round(prop.core.hour, 3)*100, '%)'),
         cost=paste0(round(cost, 3), ' (', round(prop.cost, 3)*100, '%)')) %>%
  select(task, cpu, mem, dataset, core.hour, cost) %>% 
  kable(digits=3, format='latex', caption='resource per task') %>%
  cat(file='resource-stats.tex', append=TRUE)
```

## Find SV clusters and prepare files

```sh
Rscript findDups.R
## Output: neardups-clusters.tsv and hsvlr-fordedup.vcf
```

## Short reads

BAM/CRAM from 17 diverse individuals where downloaded from HGSVC and 1000GP high-coverage datasets:

- HGSVC: ERR894724, ERR895347
- 1000GP: HG03458 HG01435 HG00956 HG00136 HG03697 NA19376 NA18873 HG01572 HG01098 NA18942 HG01861 HG02221 NA20509 HG04020 HG03593.final.cram

The 1000GP samples were selected to span super-populations using [select1kgp.R](select1kgp.R).

## Run re-mapping experiment for each cluster

Lots of jobs so we can use batches with snakemake

```sh
mkdir -p graph
for BATCH in `seq 1 20`
do
    echo "Batch " $BATCH
    snakemake --cores 16 --batch main=$BATCH/20
done
```

## Merge calls

```sh
grep "#" called_vcf/cluster-1-*.vcf >> called_merged.vcf
grep "#" -vh called_vcf/*vcf >> called_merged.vcf
```

## Use calls to annotate read support in full VCF

```sh
Rscript mergeRemapResults.R
## Output: hsvlr-srdedup17-support.tsv
```

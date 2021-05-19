## SV graph construction

1. The three long-read SV catalogs were first merged naively. See [merge-svs](merge-svs).
1. Re-mapping short reads was used to filter out some duplicates based on read support. See [remap-to-dedup-merged-svs](remap-to-dedup-merged-svs).
1. The graph was build using a hybrid approach: one SV in each remaining cluster, then iterative alignment and augmentation of other alleles. See [make-graph](make-graph).
1. The indexes were built using Snakemake. See [make-graph-indexes](make-graph-indexes).

## SV genotyping

The SVs were genotyped on Terra using [a WDL deposited on Dockstore](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/vg_mapgaffe_call_sv_cram:sv-giraffe-paper?tab=info).
See [genotype-svs](genotype-svs) for more information.

## Combining SV genotypes

The genotype in each sample were combined and clustered to define SV alleles and sites.
See [combine-sv-genotypes](combine-sv-genotypes).

## SV analysis

Once the SVs were genotyped in the MESA cohort of 1000GP dataset, we explored their patterns. 
See [describe-svs](describe-svs/README.md).
This is where most figures and numbers described in the SV section of the manuscript come from and can be reproduced.

## eQTL analysis

The SVs were tested for association with gene expression changes using RNA-seq data from Geuvadis on a subset of the 1000GP sampels.
See [eqtl](eqtl/README.md).

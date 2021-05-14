## Genotyping evaluation

The genotyping evaluation used the HGSVC and GIAB dataset as in [https://github.com/vgteam/sv-genotyping-paper]. 
We evaluated the performance across different graphs (new combined SV graph vs HGSVC or GIAB graphs), read depth (20x vs 40x), read mapper (*map* vs *giraffe*), tools (vg vs.
The [sv-evaluation.md report](sv-evaluation.md) contains commands to make the figures used in the manuscript.
Compile with:

```
Rscript -e "rmarkdown::render('sv-evaluation.Rmd')"
```

[GraphTyper2](https://github.com/DecodeGenetics/graphtyper) was run using the [`graphtyper_kubernetes`](graphtyper_kubernetes) script.

## Genotype evaluation using trios

```
Rscript -e "rmarkdown::render('sv-trio-evaluation.Rmd')"
```


## SV genotyping on Terra

The SVs were genotyped on Terra using [a WDL deposited on Dockstore](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/vg_mapgaffe_call_sv_cram:sv-giraffe-paper?tab=info).

See the [inputs.json](inputs.json) and [outputs.json](outputs.json) templates used on Terra to genotype each sample. 
Note: `"${this.object_id}"` represents the path to the CRAM file extracted from the data table with information about aligned reads (*submitted_aligned_reads*).

It cost on average $1.11 to genotype the MESA cohort (BDC billing account), and $1.56 to genotype the 1000 Genomes Project dataset (lab billing account).
The [resource-stats.md report](resource-stats.md) computes the average resources (time, computing cores, memory) for a sample or for each step of the pipeline.
Compile this report with:

```
Rscript -e "rmarkdown::render('resource-stats.Rmd')"
```

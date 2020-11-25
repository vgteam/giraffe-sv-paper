The SVs were genotyped on Terra using [a WDL deposited on Dockstore](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/vg_mapgaffe_call_sv_cram:sv-giraffe-paper?tab=info).

See the [inputs.json](inputs.json) and [outputs.json](outputs.json) templates used on Terra to genotype each sample. 
Note: `"${this.object_id}"` represents the path to the CRAM file extracted from the data table with information about aligned reads (*submitted_aligned_reads*).

It cost on average $1.11 to genotype the MESA cohort (BDC billing account), and $1.56 to genotype the 1000 Genomes Project dataset (lab billing account).
The [resource-stats report](resource-stats.md) computes the average resources (time, computing cores, memory) for a sample or for each step of the pipeline.
Compile this report with `Rscript -e "rmarkdown::render('resource-stats.Rmd')"`.

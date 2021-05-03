We used Matrix-eQTL to test for association between our SV genotypes and gene expression across samples from the GEUVADIS dataset (4 European populations + YRI).
See [find-eqtls.R](find-eqtls.R) for the commands used to find eQTLs

The [eqtl-summary.md](eqtl-summary.md) report shows p-value distribution, QQ plots and a summary of the eQTLs found. To compile the report: 

```
Rscript -e "rmarkdown::render('eqtl-summary.Rmd')"
Rscript -e "rmarkdown::render('eqtl-summary-snv-indel-sv.Rmd')"
```

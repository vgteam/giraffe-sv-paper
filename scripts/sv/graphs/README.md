# Reports

The `.md` files are markdown reports that show the code and the resulting graphs/tables.

## SVs in the MESA cohort from TOPMed

- [Summary stats/distributions, number of SVs alleles and sites](summary-sv-stats-mesa.md)
- [Examples of SV sites and their allelic patterns](examples-svsites-mesa.md)
- [PCA using the SV allele counts](pca-mesa.md)
- [SV frequency between populations](pops-freq-mesa.md)

## SVs from the 1000 Genomes Project

- [Summary stats/distributions, number of SVs alleles and sites](summary-sv-stats-1kgp.md)
- [PCA using the SV allele counts](pca-1kgp.md)
- [SV frequency between populations](pops-freq-1kgp.md)

# To compile the reports

```sh
## MESA
Rscript -e "rmarkdown::render('summary-sv-stats-mesa.Rmd')"
Rscript -e "rmarkdown::render('examples-svsites-mesa.Rmd')"
Rscript -e "rmarkdown::render('pca-mesa.Rmd')"
Rscript -e "rmarkdown::render('pops-freq-mesa.Rmd')"

## 1000 Genomes Project
Rscript -e "rmarkdown::render('summary-sv-stats-1kgp.Rmd')"
Rscript -e "rmarkdown::render('pca-1kgp.Rmd')"
Rscript -e "rmarkdown::render('pops-freq-1kgp.Rmd')"

## Joint analysis
Rscript -e "rmarkdown::render('compare-public-catalogs.Rmd')"
```

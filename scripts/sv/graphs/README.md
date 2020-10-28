The `.md` files are markdown reports that show the code and the resulting graphs/tables.

To compile the reports:

```
Rscript -e "rmarkdown::render('summary-sv-stats-mesa.Rmd')"
Rscript -e "rmarkdown::render('examples-svsites-mesa.Rmd')"
```

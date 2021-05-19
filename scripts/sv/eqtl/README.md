## SV-eQTL discovery

We used Matrix-eQTL to test for association between our SV genotypes and gene expression across samples from the GEUVADIS dataset (4 European populations + YRI).
The commands used to find eQTLs are in the [find-eqtls.R](find-eqtls.R) script.
The expression matrix (downloaded within the script) is available at [https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz](https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz).

## Joint eQTL discovery for SNVs, indels, and SVs

The SNVs and indels were downloaded from [http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/][http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/) and formatted using the [subset-vcf.py](subset-vcf.py) script:

```sh
## list samples to keep (those with expression data)
zcat GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz | head -1 | cut -f 5- | sed 's/\t/\n/g' > sample.list
mkdir 1000gp-nygc-variants

for CHR in `seq 1 22`
do
	## download full VCF with phased SNVs and indels
    ~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -Tr -Q -l 100M -P33001 -L- fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz .
	## subset and format into 0/1/2 genotype
	zcat CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz | python subset-vcf.py | gzip > 1000gp-nygc-variants/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.012.tsv.gz
	## remove full VCF
	rm CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz
done
```

Once the SNV/indel data downloaded the eQTL discovery was performed using the commands in the [find-eqtls-all-variants.R](find-eqtls-all-variants.R) script.
It assumes the SV-only analysis was run first as it reuses some of its temporary files.

## eQTL results exploration

The [eqtl-summary.md](eqtl-summary.md) and [eqtl-summary-snv-indel-sv.md](eqtl-summary-snv-indel-sv.md) reports show p-value distribution, QQ plots and a summary of the eQTLs found in the two analysis mentioned above. 

To compile the reports: 

```
Rscript -e "rmarkdown::render('eqtl-summary.Rmd')"
Rscript -e "rmarkdown::render('eqtl-summary-snv-indel-sv.Rmd')"
```

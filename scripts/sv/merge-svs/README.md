## Download LR SV catalogs

```sh
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
wget https://github.com/vgteam/sv-genotyping-paper/raw/master/human/hgsvc/HG00514.hap0.vcf.gz
wget https://github.com/vgteam/sv-genotyping-paper/raw/master/human/hgsvc/HG00514.hap1.vcf.gz
wget https://github.com/vgteam/sv-genotyping-paper/raw/master/human/hgsvc/HG00733.hap0.vcf.gz
wget https://github.com/vgteam/sv-genotyping-paper/raw/master/human/hgsvc/HG00733.hap1.vcf.gz
wget https://github.com/vgteam/sv-genotyping-paper/raw/master/human/hgsvc/NA19240.hap0.vcf.gz
wget https://github.com/vgteam/sv-genotyping-paper/raw/master/human/hgsvc/NA19240.hap1.vcf.gz
cp /public/groups/cgl/graph-genomes/jmonlong/hgsvc/glenndata/sv-pop-explicit.vcf.gz .
```

Note: `sv-pop-explicit.vcf.gz` was produced in a previous study, see [https://github.com/vgteam/sv-genotyping-paper/tree/master/human/svpop](https://github.com/vgteam/sv-genotyping-paper/tree/master/human/svpop).

## GIAB: select PASS variants and lift-over to GRCh38

```sh
Rscript filter.R
awk '! /\#/' HG002_SVs_Tier1_v0.6.filtered.vcf | awk '{print $1"\t"($2-1)"\t"($2+length($4)-1)"\t"$3}' > giab.filtered.bed
liftOver giab.filtered.bed hg19ToHg38.over.chain.gz giab.filtered.lifted.bed giab.filtered.unlifted.bed
python3 liftVcf.py -v HG002_SVs_Tier1_v0.6.filtered.vcf -l giab.filtered.lifted.bed -o HG002_SVs_Tier1_v0.6.filtered.lifted.vcf -r hg38.fa
```

## Clean-up the VCFs

Remove genotypes and most of the INFO fields

```sh
python3 removeInfo.py -v HG002_SVs_Tier1_v0.6.filtered.lifted.vcf -o HG002_SVs_Tier1_v0.6.filtered.lifted.nosamp.vcf -p giab_
python3 removeInfo.py -v sv-pop-explicit.vcf.gz -o sv-pop-explicit.nosamp.vcf -p svpop_
for SAMP in HG00514 HG00733 NA19240
do
    for HAP in 0 1
    do
	python3 removeInfo.py -v $SAMP.hap$HAP.vcf.gz -o $SAMP.hap$HAP.nosamp.vcf -p hgsvc${SAMP}h${HAP}_
    done
done
```

## Normalize variants using bcftools

```sh
for VCF in HG002_SVs_Tier1_v0.6.filtered.lifted sv-pop-explicit
do
    bcftools norm -m -both $VCF.nosamp.vcf | bcftools norm -d none -c x --fasta-ref hg38.fa | bcftools sort | bgzip > $VCF.nosamp.norm.vcf.gz
    tabix -f ${VCF}.nosamp.norm.vcf.gz
done
for SAMP in HG00514 HG00733 NA19240
do
    for HAP in 0 1
    do
	bcftools norm -m -both $SAMP.hap$HAP.nosamp.vcf | bcftools norm -d none -c x --fasta-ref hg38.fa | bcftools sort | bgzip > $SAMP.hap$HAP.nosamp.norm.vcf.gz
	tabix -f $SAMP.hap$HAP.nosamp.norm.vcf.gz
    done
done
```

## Merge catalogs with bcftools

```sh
VCFS=""
for VCF in HG002_SVs_Tier1_v0.6.filtered.lifted sv-pop-explicit
do
    VCFS="$VCFS $VCF.nosamp.norm.vcf.gz"
done
for SAMP in HG00514 HG00733 NA19240
do
    for HAP in 0 1
    do
	VCFS="$VCFS $SAMP.hap$HAP.nosamp.norm.vcf.gz"
    done
done

bcftools merge $VCFS | bcftools norm -m -any -N | bcftools norm -d none --fasta-ref hg38.fa | bcftools sort | bgzip > hsvlr.vcf.gz
tabix -f hsvlr.vcf.gz
```

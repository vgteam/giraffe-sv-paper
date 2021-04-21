#!/usr/bin/env bash
# This preprocesses the 1000 Genomes Project GRCh38 liftover VCFs to remove
# variants in large segmental duplications; not doing this produced a lot of
# miscalls in these regions realative to a graph without variants.

set -ex
set -o pipefail

rm -f GRCh38_segdups_gt10kb.bed.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v2.0/GRCh38/SegmentalDuplications/GRCh38_segdups_gt10kb.bed.gz
zcat GRCh38_segdups_gt10kb.bed.gz | sed 's/^chr//g' > GRCh38_segdups_gt10kb.nochr.bed
for CHROM in {1..21} X Y ; do
    (
        if [[ ! -e ALL.chr${CHROM}_GRCh38.genotypes.20170504.vcf.gz ]] ; then
            aws s3 cp s3://vg-data/1kg_GRCh38/variants/ALL.chr${CHROM}_GRCh38.genotypes.20170504.vcf.gz .
        fi
        if [[ ! -e ALL.chr${CHROM}_GRCh38.genotypes.20170504.vcf.gz.tbi ]] ; then
            aws s3 cp s3://vg-data/1kg_GRCh38/variants/ALL.chr${CHROM}_GRCh38.genotypes.20170504.vcf.gz.tbi .
        fi
        if [[ ! -e ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz ]] ; then
            pbgzip -n 3 -dc ALL.chr${CHROM}_GRCh38.genotypes.20170504.vcf.gz | bcftools view -T ^GRCh38_segdups_gt10kb.nochr.bed -O v | pbgzip -n 3 -c > ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz
        fi
        if [[ ! -e ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz.tbi ]] ; then
            tabix -p vcf ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz.tbi
        fi
    ) &
done

for CHROM in {1..21} X Y ; do
    wait
done

for CHROM in {1..21} X Y ; do
    aws s3 cp ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz s3://vg-data/1kg_GRCh38/variants/subsets/no_segdups_gt10kb/ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz
    aws s3 cp ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz.tbi s3://vg-data/1kg_GRCh38/variants/subsets/no_segdups_gt10kb/ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz.tbi
done



set -e

# Set file name prefixes 
GENOME_PREFIX="GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic"
VARIANTS_PREFIX="grch38_1kg_lo_nsd_all_af001"
OUT_PREFIX="${VARIANTS_PREFIX}_${CHR}"

# Download genome
aws s3 cp s3://vg-k8s/profiling/data/${GENOME_PREFIX}.fna.gz . --no-progress

# Download variants
aws s3 cp s3://vg-data/1kg_GRCh38/variants/subsets/no_segdups_gt10kb/ALL.chr${CHR}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz variants.vcf.gz --no-progress

# Decompress genome, convert IUPAC and remove chr prefix
/usr/bin/time -v bash -c "gunzip -c ${GENOME_PREFIX}.fna.gz | sed -e 's/Y/N/g' | sed -e 's/>chrN/>chrY/g' | sed -E 's/>chr([0-9]*|X|Y) />\1 /g' > ${GENOME_PREFIX}.fna"

# Normalize and filter low frequency variants
/usr/bin/time -v bash -c "bcftools annotate -x ^INFO/AC,^INFO/AF,^INFO/AN variants.vcf.gz | bcftools norm -m -any -f ${GENOME_PREFIX}.fna | bcftools view -q 0.001 -O z > variants_filt.vcf.gz; tabix variants_filt.vcf.gz"

if [ "${CHR}" = "X" ] || [ "${CHR}" = "Y" ]; then

	# Convert genotypes to diploid
	/usr/bin/time -v bash -c "zcat variants_filt.vcf.gz | head -n 10000 | grep '^#' > variants_filt.vcf; zcat variants_filt.vcf.gz | grep -v '^#' | awk -v OFS='\t' '{for (i = 10; i <= NF; i++) {if (("'$i !~ "\\/"'") && ("'$i !~ "\\|"'")) "'$i'" = "'$i "|" $i'"} {print}}' >> variants_filt.vcf; rm variants_filt.vcf.gz; bgzip variants_filt.vcf"
fi 

# Construct variant and haplotype lists
/usr/bin/time -v bash -c "zcat variants_filt.vcf.gz | grep -v '^#' | wc -l; hisat2_extract_snps_haplotypes_VCF.py --non-rs ${GENOME_PREFIX}.fna variants_filt.vcf.gz ${OUT_PREFIX}; wc -l ${OUT_PREFIX}.snp; wc -l ${OUT_PREFIX}.haplotype"

# Upload variants and haplotype lists
aws s3 sync . s3://vg-k8s/users/jsibbesen/giraffe/paper/hisat2/variants/${VARIANTS_PREFIX}/${CHR}/ --exclude "*" --include "${OUT_PREFIX}.snp" --include "${OUT_PREFIX}.haplotype" --no-progress

set -e

# Set file name prefixes 
GENOME_PREFIX="GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic"
VARIANTS_PREFIX="HGSVC.haps"
OUT_PREFIX="grch38_hgsvc_all_whole"

# Download genome
aws s3 cp s3://vg-k8s/profiling/data/${GENOME_PREFIX}.fna.gz . --no-progress

# Download variants
aws s3 cp s3://glennhickey/outstore/HGSVC-jan5/${VARIANTS_PREFIX}.vcf.gz . --no-progress

# Decompress genome and convert IUPAC
/usr/bin/time -v bash -c "gunzip -c ${GENOME_PREFIX}.fna.gz | sed -e 's/Y/N/g' > ${GENOME_PREFIX}.fna"

# Convert to bi-allelic and normalize
/usr/bin/time -v bash -c "bcftools norm -c x -m -any -f ${GENOME_PREFIX}.fna ${VARIANTS_PREFIX}.vcf.gz > variants_norm.vcf; grep -v '^#' variants_norm.vcf | wc -l"

# Construct variant and haplotype lists
/usr/bin/time -v bash -c "hisat2_extract_snps_haplotypes_VCF.py --non-rs ${GENOME_PREFIX}.fna variants_norm.vcf ${OUT_PREFIX}; wc -l ${OUT_PREFIX}.snp; wc -l ${OUT_PREFIX}.haplotype"

# Upload variants and haplotype lists
aws s3 sync . s3://vg-k8s/users/jsibbesen/giraffe/paper/hisat2/variants/grch38_hgsvc_all/whole/ --exclude "*" --include "${OUT_PREFIX}.snp" --include "${OUT_PREFIX}.haplotype" --no-progress

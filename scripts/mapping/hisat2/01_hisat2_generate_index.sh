set -e

# Set file name prefixes 
GENOME_PREFIX="GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic"
OUT_PREFIX="${VARIANTS}_index"

# Download genome
aws s3 cp s3://vg-k8s/profiling/data/${GENOME_PREFIX}.fna.gz . --no-progress

# Decompress genome
/usr/bin/time -v bash -c "gunzip -c ${GENOME_PREFIX}.fna.gz > ${GENOME_PREFIX}.fna"

if [ "${VARIANTS}" = "grch38" ]; then

	# Construct HISAT2 index
	/usr/bin/time -v bash -c "hisat2-build -p ${CPU} ${GENOME_PREFIX}.fna ${OUT_PREFIX}"

else

	if echo ${VARIANTS} | grep -q "grch38_1kg_lo"; then

		# Remove chr prefix
		/usr/bin/time -v bash -c "sed -i -E 's/>chr([0-9]*|X|Y) />\1 /g' ${GENOME_PREFIX}.fna"
	fi

	# Download variants
	aws s3 cp s3://vg-k8s/users/jsibbesen/giraffe/paper/hisat2/variants/${VARIANTS}/ . --recursive --exclude "*" --include "*.snp" --include "*.haplotype" --no-progress

	# Combine variants and haplotype lists
	/usr/bin/time -v bash -c "cat */*.snp > variants.snp; cat */*.haplotype > variants.haplotype; wc -l variants.snp; wc -l variants.haplotype"

	# Construct HISAT2 index
	/usr/bin/time -v bash -c "hisat2-build -p ${CPU} --snp variants.snp --haplotype variants.haplotype ${GENOME_PREFIX}.fna ${OUT_PREFIX}"
fi

# Upload HISAT2 index
aws s3 sync . s3://vg-k8s/users/jsibbesen/giraffe/paper/hisat2/indexes/${VARIANTS}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress

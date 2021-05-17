set -e

# Set file name prefixes
OUT_PREFIX="hisat2_${PARA}_${REF}_sim_${NAME}"

# Download reads
aws s3 cp ${READS} reads.fq.gz --no-progress

# Download index
aws s3 cp s3://vg-k8s/users/jsibbesen/giraffe/paper/hisat2/indexes/${REF}/ . --recursive --no-progress

# De-interleave reads (https://gist.github.com/nathanhaigh/3521724)
/usr/bin/time -v bash -c 'zcat reads.fq.gz | wc -l; zcat reads.fq.gz | paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > reads_1.fq) | cut -f 5-8 | tr "\t" "\n" > reads_2.fq; wc -l reads_1.fq; wc -l reads_2.fq; zcat reads.fq.gz | head -n 20; head reads_*.fq; zcat reads.fq.gz | tail -n 20; tail reads_*.fq'

# Compress reads
/usr/bin/time -v bash -c 'gzip reads_1.fq; gzip reads_2.fq'	

# Use default HISAT2
if [ "${PARA}" = "def" ]; then

	# Map single-end reads
	/usr/bin/time -v bash -c "hisat2 -p ${CPU} -t --maxins 1065 --no-spliced-alignment -x ${REF}_index -U reads.fq.gz -S ${OUT_PREFIX}_se.sam"

	# Map paired-end reads
	/usr/bin/time -v bash -c "hisat2 -p ${CPU} -t --maxins 1065 --no-spliced-alignment -x ${REF}_index -1 reads_1.fq.gz -2 reads_2.fq.gz -S ${OUT_PREFIX}_pe.sam"

# Use sensitive HISAT2
elif [ "${PARA}" = "sens" ]; then

	# Map single-end reads
	/usr/bin/time -v bash -c "hisat2 -p ${CPU} -t --maxins 1065 --no-spliced-alignment --sensitive -x ${REF}_index -U reads.fq.gz -S ${OUT_PREFIX}_se.sam"

	# Map paired-end reads
	/usr/bin/time -v bash -c "hisat2 -p ${CPU} -t --maxins 1065 --no-spliced-alignment --sensitive -x ${REF}_index -1 reads_1.fq.gz -2 reads_2.fq.gz -S ${OUT_PREFIX}_pe.sam"

# Use very sensitive HISAT2
elif [ "${PARA}" = "vsens" ]; then

	# Map single-end reads
	/usr/bin/time -v bash -c "hisat2 -p ${CPU} -t --maxins 1065 --no-spliced-alignment --very-sensitive -x ${REF}_index -U reads.fq.gz -S ${OUT_PREFIX}_se.sam"

	# Map paired-end reads
	/usr/bin/time -v bash -c "hisat2 -p ${CPU} -t --maxins 1065 --no-spliced-alignment --very-sensitive -x ${REF}_index -1 reads_1.fq.gz -2 reads_2.fq.gz -S ${OUT_PREFIX}_pe.sam"
fi

# Compress single-end alignments
/usr/bin/time -v bash -c "samtools view -b -O BAM --threads ${CPU} ${OUT_PREFIX}_se.sam > ${OUT_PREFIX}_se.bam"

# Compress paired-end alignments
/usr/bin/time -v bash -c "samtools view -b -O BAM --threads ${CPU} ${OUT_PREFIX}_pe.sam > ${OUT_PREFIX}_pe.bam"

# Upload read alignments 
aws s3 sync . s3://vg-k8s/users/jsibbesen/giraffe/paper/hisat2/alignments/sim_review2/${PARA}/${REF}/${NAME}/ --exclude "*" --include "${OUT_PREFIX}*" --exclude "*.sam" --no-progress

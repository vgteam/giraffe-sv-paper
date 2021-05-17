set -e

# Set file name prefixes
ALIGN_PREFIX="hisat2_${PARA}_${REF}_sim_${NAME}"
OUT_PREFIX="${ALIGN_PREFIX}_gam_d100"

# Download alignments
aws s3 cp s3://vg-k8s/users/jsibbesen/giraffe/paper/hisat2/alignments/sim/${PARA}/${REF}/${NAME}/ . --recursive --exclude "*" --include "*.bam" --no-progress

# Download graph
aws s3 cp ${GRAPH} graph.xg --no-progress

# Download simulated alignments
aws s3 cp ${SIM} sim.gam --no-progress

for TYPE in se pe; do 

	# Get alignment statistics
	/usr/bin/time -v bash -c "samtools flagstat ${ALIGN_PREFIX}_${TYPE}.bam"

	# Print number of supplementary alignments
	/usr/bin/time -v bash -c "samtools view -f 2048 ${ALIGN_PREFIX}_${TYPE}.bam | wc -l"

	# Filter secondary alignments
	/usr/bin/time -v bash -c "samtools view -F 256 -b ${ALIGN_PREFIX}_${TYPE}.bam > align_primary.bam"

	# Fix chromosome names
	if echo ${ALIGN_PREFIX} | grep -q "hgsvc"; then

		if [ "${REF}" = "grch38_snp" ]; then

			wget --no-check-certificate https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_gencode2UCSC.txt
		
			/usr/bin/time -v bash -c "samtools view -H align_primary.bam | sed -E 's/SN:MT/SN:chrM/g' | sed -f <(paste -d : <(cut -f1 GRCh38_gencode2UCSC.txt) <(cut -f2 GRCh38_gencode2UCSC.txt) | sed 's/\([^:]*\):\([^:]*\)/s%SN:\1%SN:\2%/') > new_header.sam; samtools reheader new_header.sam align_primary.bam > align_primary_fixed.bam; mv align_primary_fixed.bam align_primary.bam"

		else	

			/usr/bin/time -v bash -c "samtools view -H align_primary.bam | sed -E 's/SN:chr([0-9]*|X|Y)	/SN:\1	/g' > new_header.sam; samtools reheader new_header.sam align_primary.bam > align_primary_fixed.bam; mv align_primary_fixed.bam align_primary.bam"
		fi
	fi 

	# Inject alignments
	/usr/bin/time -v bash -c "vg inject -t ${CPU} -x graph.xg align_primary.bam | vg view -a - | sed 's/\/1\",/\",/g' | sed 's/\/2\",/\",/g' | vg view -a -G -J - > align.gam"

	# Annotate alignments
	/usr/bin/time -v bash -c "vg annotate -t ${CPU} -m -x graph.xg -a align.gam > align_anno.gam"

	# Calculate alignment statistics
	/usr/bin/time -v bash -c "vg gamcompare -t ${CPU} -T -r 100 -a hisat2-${REF}-${NAME}-${PARA}-${TYPE} align_anno.gam sim.gam | cut -f1-3 | awk -v OFS='\t' '{print "'$1'", "'$2'", 1, "'$3'"}' > ${OUT_PREFIX}_${TYPE}.txt; gzip ${OUT_PREFIX}_${TYPE}.txt"
done

# Upload alignment statistics 
aws s3 sync . s3://vg-k8s/users/jsibbesen/giraffe/paper/hisat2/stats/sim/${PARA}/${REF}/${NAME}/ --exclude "*" --include "${OUT_PREFIX}*" --no-progress

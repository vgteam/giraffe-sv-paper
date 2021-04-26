set -e
printf "graph\talgorithm\treads\tpairing\tload_time\tspeed\n" > report_speed_hisat2.tsv
CPU=16


# Download reads
#aws s3 cp s3://vg-k8s/profiling/reads/real/NA19239/novaseq6000-ERR3239454-shuffled-1m.fq.gz novaseq6000.fq.gz --no-progress
#aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseq2500-ERR309934-shuffled-1m.fq.gz hiseq2500.fq.gz --no-progress
#aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseqxten-SRR6691663-shuffled-1m.fq.gz hiseqxten.fq.gz --no-progress
#for STRAIN in DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138 ; do
#    aws s3 cp s3://vg-k8s/profiling/reads/real/yeast/${STRAIN}-shuffled.fq.gz ${STRAIN}.fq.gz --no-progress
#done

aws s3 cp s3://vg-k8s/profiling/reads/real/NA19239/novaseq6000-ERR3239454-shuffled-600m.fq.gz novaseq6000.fq.gz


# De-interleave reads (https://gist.github.com/nathanhaigh/3521724)
/usr/bin/time -v bash -c 'zcat novaseq6000.fq.gz | wc -l; zcat novaseq6000.fq.gz | paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > novaseq6000_1.fq) | cut -f 5-8 | tr "\t" "\n" > novaseq6000_2.fq; wc -l novaseq6000_1.fq; wc -l novaseq6000_2.fq'
#/usr/bin/time -v bash -c 'zcat hiseqxten.fq.gz | wc -l; zcat hiseqxten.fq.gz | paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > hiseqxten_1.fq) | cut -f 5-8 | tr "\t" "\n" > hiseqxten_2.fq; wc -l hiseqxten_1.fq; wc -l hiseqxten_2.fq'
#/usr/bin/time -v bash -c 'zcat hiseq2500.fq.gz | wc -l; zcat hiseq2500.fq.gz | paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > hiseq2500_1.fq) | cut -f 5-8 | tr "\t" "\n" > hiseq2500_2.fq; wc -l hiseq2500_1.fq; wc -l hiseq2500_2.fq'
#for STRAIN in DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138 ; do
#    /usr/bin/time -v bash -c 'zcat '${STRAIN}'.fq.gz | wc -l; zcat '${STRAIN}'.fq.gz | paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > '${STRAIN}'_1.fq) | cut -f 5-8 | tr "\t" "\n" > '${STRAIN}'_2.fq; wc -l '${STRAIN}'_1.fq; wc -l '${STRAIN}'_2.fq'
#done

# Compress reads
/usr/bin/time -v bash -c 'gzip novaseq6000_1.fq; gzip novaseq6000_2.fq'	
#/usr/bin/time -v bash -c 'gzip hiseqxten_1.fq; gzip hiseqxten_2.fq'	
#/usr/bin/time -v bash -c 'gzip hiseq2500_1.fq; gzip hiseq2500_2.fq'
#for STRAIN in DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138 ; do
#    /usr/bin/time -v bash -c 'gzip '${STRAIN}'_1.fq; gzip '${STRAIN}'_2.fq'
#done

for SPECIES in human ; do
    case "${SPECIES}" in
    yeast)
        REFS=(s288c)
        # Yeast graphs aren't VCF-derived so hisat2 can't use them.
        READSETS=(DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138)
        ;;
    human)
        REFS=(grch38_snp)
        READSETS=(novaseq6000)
        ;;
    esac
    for PARA in def ; do
        for REF in ${REFS[@]} ; do
                
            # Download index
            aws s3 cp s3://vg-k8s/users/jsibbesen/giraffe/paper/hisat2/indexes/${REF}/ . --recursive --no-progress

            for READS in ${READSETS[@]} ; do

                # Set file name prefixes
                OUT_PREFIX="hisat2_${PARA}_${REF}_${READS}"
                
                
                # Use default HISAT2
                if [ "${PARA}" = "def" ]; then
                
                    # Map single-end reads
                    /usr/bin/time -v bash -c "hisat2 -p ${CPU} -t --maxins 1065 --no-spliced-alignment -x ${REF}_index -U ${READS}.fq.gz -S ${OUT_PREFIX}_se.sam 2> log_single.txt"
                
                    # Map paired-end reads
                    /usr/bin/time -v bash -c "hisat2 -p ${CPU} -t --maxins 1065 --no-spliced-alignment -x ${REF}_index -1 ${READS}_1.fq.gz -2 ${READS}_2.fq.gz -S ${OUT_PREFIX}_pe.sam 2> log_paired.txt"
                
                # Use sensitive HISAT2
                elif [ "${PARA}" = "sens" ]; then
                
                    # Map single-end reads
                    /usr/bin/time -v bash -c "hisat2 -p ${CPU} -t --maxins 1065 --no-spliced-alignment --sensitive -x ${REF}_index -U ${READS}.fq.gz -S ${OUT_PREFIX}_se.sam 2> log_single.txt"
                
                    # Map paired-end reads
                    /usr/bin/time -v bash -c "hisat2 -p ${CPU} -t --maxins 1065 --no-spliced-alignment --sensitive -x ${REF}_index -1 ${READS}_1.fq.gz -2 ${READS}_2.fq.gz -S ${OUT_PREFIX}_pe.sam 2> log_paired.txt"
                
                # Use very sensitive HISAT2
                elif [ "${PARA}" = "vsens" ]; then
                
                    # Map single-end reads
                    /usr/bin/time -v bash -c "hisat2 -p ${CPU} -t --maxins 1065 --no-spliced-alignment --very-sensitive -x ${REF}_index -U ${READS}.fq.gz -S ${OUT_PREFIX}_se.sam 2> log_single.txt"
                
                    # Map paired-end reads
                    /usr/bin/time -v bash -c "hisat2 -p ${CPU} -t --maxins 1065 --no-spliced-alignment --very-sensitive -x ${REF}_index -1 ${READS}_1.fq.gz -2 ${READS}_2.fq.gz -S ${OUT_PREFIX}_pe.sam 2> log_paired.txt"
                fi
     
                #Get the time as reported by the mapper

                MAPPED_COUNT="$(cat log_single.txt | grep "reads" | awk '{print$1}')"
                LOAD_TIME="$(cat log_single.txt | grep "Time loading" | awk -F: '{print ($2*3600) + ($3*60) + $4}' | awk '{sum+=$1} END {print sum}')"
                RUNTIME="$(cat log_single.txt | grep "Multiseed full-index search" | awk -F: '{print ($2*3600) + ($3*60) + $4}')"
                ALL_RPS="$(echo "${MAPPED_COUNT} / ${RUNTIME}" | bc -l)"
                RPS_PER_THREAD="$(echo "${ALL_RPS} / ${THREAD_COUNT}" | bc -l)"
                printf "${GRAPH}\thiseq2_${PARA}\t${READS}\tsingle\t${LOAD_TIME}\t${RPS_PER_THREAD}\n" >> report_speed_hisat2.tsv

                MAPPED_COUNT="$(cat log_paired.txt | grep "reads" | awk '{print$1}')"
                LOAD_TIME="$(cat log_paired.txt | grep "Time loading" | awk -F: '{print ($2*3600) + ($3*60) + $4}' | awk '{sum+=$1} END {print sum}')"
                RUNTIME="$(cat log_paired.txt | grep "Multiseed full-index search" | awk -F: '{print ($2*3600) + ($3*60) + $4}')"
                ALL_RPS="$(echo "${MAPPED_COUNT} / ${RUNTIME}" | bc -l)"
                RPS_PER_THREAD="$(echo "2 * ${ALL_RPS} / ${THREAD_COUNT}" | bc -l)"
                printf "${GRAPH}\thiseq2_${PARA}\t${READS}\tpaired\t${LOAD_TIME}\t${RPS_PER_THREAD}\n" >> report_speed_hisat2.tsv


                >&2 cat log_single.txt

                # Compress single-end alignments
                #/usr/bin/time -v bash -c "samtools view -b -O BAM --threads ${CPU} ${OUT_PREFIX}_se.sam > ${OUT_PREFIX}_se.bam"
                #
                ## Compress paired-end alignments
                #/usr/bin/time -v bash -c "samtools view -b -O BAM --threads ${CPU} ${OUT_PREFIX}_pe.sam > ${OUT_PREFIX}_pe.bam"

            done
        done
    done
done

#!/usr/bin/env bash

set -e

THREAD_COUNT=16

printf "graph\talgorithm\treads\tpairing\tspeed\ttotal_cpu_sec\ttotal_wall_clock_time\tmax_resident_set_size(kbytes)\n" > report_speed.tsv

#get all real read sets
aws s3 cp s3://vg-k8s/profiling/reads/real/NA19239/novaseq6000-ERR3239454-shuffled-1m.fq.gz novaseq6000.fq.gz
aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseq2500-ERR309934-shuffled-1m.fq.gz hiseq2500.fq.gz
aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseqxten-SRR6691663-shuffled-1m.fq.gz hiseqxten.fq.gz

#Get bigger read sets
#aws s3 cp s3://vg-k8s/profiling/reads/real/NA19239/novaseq6000-ERR3239454-shuffled-600m.fq.gz novaseq6000.fq.gz
#aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseq2500-ERR309934-shuffled-600m.fq.gz hiseq2500.fq.gz
#aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseqxten-SRR6691663-shuffled-600m.fq.gz hiseqxten.fq.gz

for STRAIN in DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138 ; do
    aws s3 cp s3://vg-k8s/profiling/reads/real/yeast/${STRAIN}.fq.gz ${STRAIN}.fq.gz
done

gunzip novaseq6000.fq.gz
gunzip hiseq2500.fq.gz
gunzip hiseqxten.fq.gz
for STRAIN in DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138 ; do
    gunzip ${STRAIN}.fq.gz
done

#Take out every other read so that minimap interprets it as single end
awk 'NR%8==1 || NR%8==2 || NR%8==3 || NR%8==4' novaseq6000.fq  > novaseq6000.unpaired.fq
awk 'NR%8==1 || NR%8==2 || NR%8==3 || NR%8==4' hiseq2500.fq  > hiseq2500.unpaired.fq
awk 'NR%8==1 || NR%8==2 || NR%8==3 || NR%8==4' hiseqxten.fq  > hiseqxten.unpaired.fq
for STRAIN in DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138 ; do
    awk 'NR%8==1 || NR%8==2 || NR%8==3 || NR%8==4' ${STRAIN}.fq  > ${STRAIN}.unpaired.fq
done


# Get the reference genomes
aws s3 cp s3://vg-k8s/profiling/data/hs37d5.fa hs37d5.fa
aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.gz hs38d1.fa.gz
gunzip hs38d1.fa.gz
aws s3 cp s3://vg-k8s/profiling/graphs/v2/generic/primary/S288C/primaryS288C.fa ./S288C.fa


#Build all indexes
#bwa index hs37d5.fa
#bwa index hs38d1.fa
#bwa index S288C.fa

#bowtie2-build --large-index hs37d5.fa hs37d5
#bowtie2-build --large-index hs38d1.fa hs38d1
bowtie2-build --large-index S288C.fa S288C

#minimap2 -x sr -d hs37d5.mmi hs37d5.fa
#minimap2 -x sr -d hs38d1.mmi hs38d1.fa
minimap2 -x sr -d S288C.mmi S288C.fa

for SPECIES in yeast human ; do
    case "${SPECIES}" in
    yeast)
        GRAPHS=(S288C)
        READSETS=(DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138)
        ;;
    human)
        GRAPHS=(hs38d1 hs37d5)
        READSETS=(novaseq6000)
        ;;
    esac
    
    #run all bwa mem for species
    for GRAPH in ${GRAPHS[@]} ; do
        
        # Determine where the BWA indexes are
        case "${GRAPH}" in
        S288C)
            BWA_INDEX_BASE=s3://vg-k8s/profiling/data/hs37d5.fa
            ;;
        hs38d1)
            BWA_INDEX_BASE=s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna
            ;;
        hs37d5)
            BWA_INDEX_BASE=s3://vg-k8s/profiling/graphs/v2/generic/primary/S288C/primaryS288C.fa
            ;;
        esac
        
        # Get (or put) the BWA indexes
        for EXT in amb ann bwt fai pac sa ; do
            if [ -e ${GRAPH}.fa.${EXT} ] ; then
                # Upload index made here
                aws s3 cp ${GRAPH}.fa.${EXT} ${BWA_INDEX_BASE}.${EXT} || true
            else
                # Download index
                aws s3 cp ${BWA_INDEX_BASE}.${EXT} ${GRAPH}.fa.${EXT} || true
            fi
        done
        
        for READS in ${READSETS[@]} ; do
            for PAIRING in single paired ; do
                if [[ ${PAIRING} ==  "single" ]] ; then
                    PAIRED=""
                elif [[ ${PAIRING} == "paired" ]] ; then
                    PAIRED="-p"
                fi

                /usr/bin/time -v bash -c "bwa mem -t ${THREAD_COUNT} ${PAIRED} ${GRAPH}.fa ${READS}.fq > mapped.bam 2> log.txt" 2> time-log.txt

                MAPPING_TIME="$(cat "log.txt" | grep "Processed" | sed 's/[^0-9]*\([0-9]*\) reads in .* CPU sec, \([0-9]*\.[0-9]*\) real sec/\1/g' | tr ' ' '\t' | awk '{sum1+=$1} END {print sum1}')"
                RPS_ALL_THREADS="$(cat "log.txt" | grep "Processed" | sed 's/[^0-9]*\([0-9]*\) reads in .* CPU sec, \([0-9]*\.[0-9]*\) real sec/\1 \2/g' | tr ' ' '\t' | awk '{sum1+=$1; sum2+=$2} END {print sum1/sum2}')"
                TOTAL_TIME="$(cat "log.txt" | grep "[main] Real time" | sed 's/Real time: \([0-9]*\.[0-9]*\) sec/\1/g')"
                LOAD_TIME="$(echo "${TOTAL_TIME} - ${MAPPING_TIME}" | bc -l)"
                RPS_PER_THREAD="$(echo "${RPS_ALL_THREADS} / ${THREAD_COUNT}" | bc -l)"

                #Get time and memory use as reported by /usr/bin/time

                USER_TIME="$(cat "time-log.txt" | grep "User time" | sed 's/User\ time\ (seconds):\ \([0-9]*\.[0-9]*\)/\1/g')"
                SYS_TIME="$(cat "time-log.txt" | grep "System time" | sed 's/System\ time\ (seconds):\ \([0-9]*\.[0-9]*\)/\1/g')"
                TOTAL_TIME="$(echo "${USER_TIME} + ${SYS_TIME}" | bc -l)"
                CLOCK_TIME="$(cat "time-log.txt" | grep "Elapsed (wall clock) time" | sed 's/.*\ \([0-9,:]*\)/\1/g')"
                MEMORY="$(cat "time-log.txt" | grep "Maximum resident set" | sed 's/Maximum\ resident\ set\ size\ (kbytes):\ \([0-9]*\)/\1/g')"

                
                printf "${GRAPH}\tbwa_mem\t${READS}\t${PAIRING}\t-\t${RPS_PER_THREAD}\t${TOTAL_TIME}\t${CLOCK_TIME}\t${MEMORY}\n" >> report_speed.tsv 
                printf "${GRAPH}\tbwa_mem\t${READS}\t${PAIRING}\t-\t${RPS_PER_THREAD}\t${TOTAL_TIME}\t${CLOCK_TIME}\t${MEMORY}\n" >> bwa-log.txt 
                cat log.txt >> bwa-log.txt
                cat time-log.txt >> bwa-log.txt

            done
        done
        
        # Clean up indexes for graph
        for EXT in amb ann bwt fai pac sa ; do
            rm -f ${GRAPH}.fa.${EXT}
        done
        
    done
    
    
    #Run all minimap2 for species
    for GRAPH in ${GRAPHS[@]} ; do
        
        case "${GRAPH}" in
        S288C)
            MINIMAP2_INDEX_BASE=s3://vg-k8s/profiling/graphs/v2/generic/primary/S288C/primaryS288C
            ;;
        hs38d1)
            MINIMAP2_INDEX_BASE=s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic
            ;;
        hs37d5)
            MINIMAP2_INDEX_BASE=s3://vg-k8s/profiling/data/hs37d5
            ;;
        esac
        
        # Get (or put) the minimap2 indexes
        for EXT in mmi ; do
            if [ -e ${GRAPH}.${EXT} ] ; then
                # Upload index made here
                aws s3 cp ${GRAPH}.${EXT} ${MINIMAP2_INDEX_BASE}.${EXT} || true
            else
                # Download index
                aws s3 cp ${MINIMAP2_INDEX_BASE}.${EXT} ${GRAPH}.${EXT} || true
            fi
        done
    
        for READS in ${READSETS[@]} ; do
            for PAIRING in single paired ; do
                if [[ ${PAIRING} ==  "single" ]] ; then
                    /usr/bin/time -v bash -c "minimap2 -ax sr --secondary=no -t ${THREAD_COUNT} ${GRAPH}.mmi ${READS}.unpaired.fq > mapped.bam 2> log.txt" 2> time-log.txt
                elif [[ ${PAIRING} == "paired" ]] ; then
                    /usr/bin/time -v bash -c "minimap2 -ax sr --secondary=no -t ${THREAD_COUNT} ${GRAPH}.mmi ${READS}.fq > mapped.bam 2> log.txt" 2> time-log.txt
                fi


                MAPPED_COUNT="$(cat log.txt | grep "mapped" | awk '{sum+=$3} END {print sum}')"
                INDEX_LOAD_TIME="$(cat log.txt | grep "loaded/built the index" | sed 's/.M::main::\([0-9]*\.[0-9]*\).*/\1 /g')"
                TOTAL_TIME="$(cat log.txt | grep "\[M::main\] Real time" | sed 's/.*Real time: \([0-9]*\.[0-9]*\) sec.*/\1/g')"

                RUNTIME="$(echo "${TOTAL_TIME} - ${INDEX_LOAD_TIME}" | bc -l)"
                ALL_RPS="$(echo "${MAPPED_COUNT} / ${RUNTIME}" | bc -l)"
                RPS_PER_THREAD="$(echo "${ALL_RPS} / ${THREAD_COUNT}" | bc -l)"
                REPORTED_MEMORY="$(cat log.txt | grep "Peak RSS" | sed 's/.*Peak RSS:\ \([0-9]*\.[0-9]*\)\ GB/\1/g')"
                

                #Get time and memory use as reported by /usr/bin/time

                USER_TIME="$(cat "time-log.txt" | grep "User time" | sed 's/User\ time\ (seconds):\ \([0-9]*\.[0-9]*\)/\1/g')"
                SYS_TIME="$(cat "time-log.txt" | grep "System time" | sed 's/System\ time\ (seconds):\ \([0-9]*\.[0-9]*\)/\1/g')"
                TOTAL_TIME="$(echo "${USER_TIME} + ${SYS_TIME}" | bc -l)"
                CLOCK_TIME="$(cat "time-log.txt" | grep "Elapsed (wall clock) time" | sed 's/.*\ \([0-9,:]*\)/\1/g')"
                MEMORY="$(cat "time-log.txt" | grep "Maximum resident set" | sed 's/Maximum\ resident\ set\ size\ (kbytes):\ \([0-9]*\)/\1/g')"


                printf "${GRAPH}\tminimap2\t${READS}\t${PAIRING}\t${INDEX_LOAD_TIME}\t${RPS_PER_THREAD}\t${TOTAL_TIME}\t${CLOCK_TIME}\t${MEMORY}\t${REPORTED_MEMORY}\n" >> report_speed.tsv 
                printf "${GRAPH}\tminimap2\t${READS}\t${PAIRING}\t${INDEX_LOAD_TIME}\t${RPS_PER_THREAD}\t${TOTAL_TIME}\t${CLOCK_TIME}\t${MEMORY}\t${REPORTED_MEMORY}\n" >> minimap2-log.txt 
                cat log.txt >> minimap2-log.txt
                cat time-log.txt >> minimap2-log.txt
            done
        done
        
        # Clean up indexes for graph
        for EXT in mmi ; do
            rm -f ${GRAPH}.${EXT}
        done
    done
    

    #run all bowtie2 for species
    for GRAPH in ${GRAPHS[@]} ; do
    
        # Determine where the indexes are
        case "${GRAPH}" in
        S288C)
            BOWTIE2_INDEX_BASE=s3://vg-k8s/profiling/graphs/v2/generic/primary/S288C/primaryS288C
            ;;
        hs38d1)
            BOWTIE2_INDEX_BASE=s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic
            ;;
        hs37d5)
            BOWTIE2_INDEX_BASE=s3://vg-k8s/profiling/data/hs37d5
            ;;
        esac
        
        # Get (or put) the bowtie2 indexes
        for EXT in 1.bt2l 2.bt2l 3.bt2l 4.bt2l rev.1.bt2l rev.2.bt2l ; do
            if [ -e ${GRAPH}.${EXT} ] ; then
                # Upload index made here
                aws s3 cp ${GRAPH}.${EXT} ${BOWTIE2_INDEX_BASE}.${EXT} 
            else
                # Download index
                aws s3 cp ${BOWTIE2_INDEX_BASE}.${EXT} ${GRAPH}.${EXT}
            fi
        done
        
        for READS in ${READSETS[@]} ; do
            for PAIRING in single paired ; do
                if [[ ${PAIRING} ==  "single" ]] ; then
                    PAIRED="-U"
                elif [[ ${PAIRING} == "paired" ]] ; then
                    PAIRED="--interleaved"
                fi

                /usr/bin/time -v bash -c "bowtie2 -t -p ${THREAD_COUNT} -x ${GRAPH} ${PAIRED} ${READS}.fq > mapped.bam 2> log.txt" 2> time-log.txt

                MAPPED_COUNT="$(cat log.txt | grep "reads" | awk '{print$1}')"
                LOAD_TIME="$(cat log.txt | grep "Time loading" | awk -F: '{print ($2*3600) + ($3*60) + $4}' | awk '{sum+=$1} END {print sum}')"
                RUNTIME="$(cat log.txt | grep "Multiseed full-index search" | awk -F: '{print ($2*3600) + ($3*60) + $4}')"
                ALL_RPS="$(echo "${MAPPED_COUNT} / ${RUNTIME}" | bc -l)"
                if [[ ${PAIRING} == "paired" ]] ; then
                    RPS_PER_THREAD="$(echo "2 * ${ALL_RPS} / ${THREAD_COUNT}" | bc -l)"
                elif [[ ${PAIRING} == "single" ]] ; then
                    RPS_PER_THREAD="$(echo "${ALL_RPS} / ${THREAD_COUNT}" | bc -l)"
                fi

                
                #Get time and memory use as reported by /usr/bin/time

                USER_TIME="$(cat "time-log.txt" | grep "User time" | sed 's/User\ time\ (seconds):\ \([0-9]*\.[0-9]*\)/\1/g')"
                SYS_TIME="$(cat "time-log.txt" | grep "System time" | sed 's/System\ time\ (seconds):\ \([0-9]*\.[0-9]*\)/\1/g')"
                TOTAL_TIME="$(echo "${USER_TIME} + ${SYS_TIME}" | bc -l)"
                CLOCK_TIME="$(cat "time-log.txt" | grep "Elapsed (wall clock) time" | sed 's/.*\ \([0-9,:]*\)/\1/g')"
                MEMORY="$(cat "time-log.txt" | grep "Maximum resident set" | sed 's/Maximum\ resident\ set\ size\ (kbytes):\ \([0-9]*\)/\1/g')"


                printf "${GRAPH}\tbowtie2\t${READS}\t${PAIRING}\t${LOAD_TIME}\t${RPS_PER_THREAD}\t${TOTAL_TIME}\t${CLOCK_TIME}\t${MEMORY}\n" >> report_speed.tsv 

                printf "${GRAPH}\tbowtie2\t${READS}\t${PAIRING}\t${LOAD_TIME}\t${RPS_PER_THREAD}\t${TOTAL_TIME}\t${CLOCK_TIME}\t${MEMORY}\n" >> bowtie2-log.txt 
                cat log.txt >> bowtie2-log.txt
                cat time-log.txt >> bowtie2-log.txt

            done
        done
        
        # Clean up indexes for graph
        for EXT in 1.bt2l 2.bt2l 3.bt2l 4.bt2l rev.1.bt2l rev.2.bt2l ; do
            rm -f ${GRAPH}.${EXT}
        done
    done
done

aws s3 cp report_speed.tsv s3://vg-k8s/users/xhchang/giraffe_experiments/speed/speed_report_linear.tsv
aws s3 cp bwa-log.txt s3://vg-k8s/users/xhchang/giraffe_experiments/speed/bwa_speed_log.txt
aws s3 cp minimap2-log.txt s3://vg-k8s/users/xhchang/giraffe_experiments/speed/minimap2_speed_log.txt
aws s3 cp bowtie2-log.txt s3://vg-k8s/users/xhchang/giraffe_experiments/speed/bowtie2_speed_log.txt

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


Get the reference genomes
aws s3 cp s3://vg-k8s/profiling/data/hs37d5.fa 1kg.fa
aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.gz hgsvc.fa.gz
gunzip hgsvc.fa.gz
aws s3 cp s3://vg-k8s/profiling/graphs/v2/generic/primary/S288C/primaryS288C.fa ./S288C.fa


#Build all indexes
#bwa index 1kg.fa
#bwa index hgsvc.fa
bwa index S288C.fa

#bowtie2-build --large-index 1kg.fa 1kg
#bowtie2-build --large-index hgsvc.fa hgsvc
bowtie2-build --large-index S288C.fna S288C

#minimap2 -x sr -d 1kg.mmi 1kg.fa
#minimap2 -x sr -d hgsvc.mmi hgsvc.fa
minimap2 -x sr -d S288C.mmi S288C.fa




for SPECIES in human ; do
    case "${SPECIES}" in
    yeast)
        GRAPHS=(S288C)
        READSETS=(DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138)
        ;;
    human)
        GRAPHS=(hgsvc 1kg)
        READSETS=(novaseq6000)
        ;;
    esac
    aws s3 cp s3://vg-k8s/profiling/data/hs37d5.fa.amb 1kg.fa.amb
    aws s3 cp s3://vg-k8s/profiling/data/hs37d5.fa.ann 1kg.fa.ann
    aws s3 cp s3://vg-k8s/profiling/data/hs37d5.fa.bwt 1kg.fa.bwt
    aws s3 cp s3://vg-k8s/profiling/data/hs37d5.fa.fai 1kg.fa.fai
    aws s3 cp s3://vg-k8s/profiling/data/hs37d5.fa.pac 1kg.fa.pac
    aws s3 cp s3://vg-k8s/profiling/data/hs37d5.fa.sa 1kg.fa.sa

    aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.amb hgsvc.fa.amb
    aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.ann hgsvc.fa.ann
    aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.bwt hgsvc.fa.bwt
    aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.fai hgsvc.fa.fai
    aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.pac hgsvc.fa.pac
    aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.sa hgsvc.fa.sa

    #run all bwa mem for species
    for GRAPH in ${GRAPHS[@]} ; do
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
                printf "${GRAPH}\tbwa_mem\t${READS}\t${PAIRING}\t-\t${RPS_PER_THREAD}\t${TOTAL_TIME}\t${CLOCK_TIME}\t${MEMORY}\n" >> bwt-log.tsv 
                cat log.txt >> bwa-log.txt
                cat time-log.txt >> bwa-log.txt

            done
        done
    done
    rm 1kg.fa.amb
    rm 1kg.fa.ann
    rm 1kg.fa.bwt
    rm 1kg.fa.fai
    rm 1kg.fa.pac
    rm 1kg.fa.sa
    rm hgsvc.fa.amb
    rm hgsvc.fa.ann
    rm hgsvc.fa.bwt
    rm hgsvc.fa.fai
    rm hgsvc.fa.pac
    rm hgsvc.fa.sa

    aws s3 cp s3://vg-k8s/profiling/data/hs37d5.mmi 1kg.mmi
    aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.mmi hgsvc.mmi
    #Run all minimap2 for species
    for GRAPH in ${GRAPHS[@]} ; do
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
    done
    rm 1kg.mmi
    rm hgsvc.mmi



    aws s3 cp s3://vg-k8s/profiling/data/hs37d5.1.bt2l 1kg.1.bt2l
    aws s3 cp s3://vg-k8s/profiling/data/hs37d5.2.bt2l 1kg.2.bt2l
    aws s3 cp s3://vg-k8s/profiling/data/hs37d5.3.bt2l 1kg.3.bt2l
    aws s3 cp s3://vg-k8s/profiling/data/hs37d5.4.bt2l 1kg.4.bt2l
    aws s3 cp s3://vg-k8s/profiling/data/hs37d5.rev.1.bt2l 1kg.rev.1.bt2l
    aws s3 cp s3://vg-k8s/profiling/data/hs37d5.rev.2.bt2l 1kg.rev.2.bt2l
    aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.1.bt2l hgsvc.1.bt2l
    aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.2.bt2l hgsvc.2.bt2l
    aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.3.bt2l hgsvc.3.bt2l
    aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.4.bt2l hgsvc.4.bt2l
    aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.rev.1.bt2l hgsvc.rev.1.bt2l
    aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.rev.2.bt2l hgsvc.rev.2.bt2l
    #run all bowtie2 for species
    for GRAPH in ${GRAPHS[@]} ; do
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

                printf "${GRAPH}\tbowtie2\t${READS}\t${PAIRING}\t${LOAD_TIME}\t${RPS_PER_THREAD}\t${TOTAL_TIME}\t${CLOCK_TIME}\t${MEMORY}\n" >> bowtie2-log.tsv 
                cat log.txt >> bowtie2-log.txt
                cat time-log.txt >> bowtie2-log.txt

            done
        done
    done
done


rm 1kg.1.bt2l
rm 1kg.2.bt2l
rm 1kg.3.bt2l
rm 1kg.4.bt2l
rm 1kg.rev.1.bt2l
rm 1kg.rev.2.bt2l
 rm hgsvc.1.bt2l
 rm hgsvc.2.bt2l
 rm hgsvc.3.bt2l
 rm hgsvc.4.bt2l
 rm hgsvc.rev.1.bt2l
 rm hgsvc.rev.2.bt2l

#!/usr/bin/env bash

set -ex

THREAD_COUNT=16


#aws s3 cp s3://vg-k8s/profiling/reads/real/NA19239/novaseq6000-ERR3239454-shuffled-600m.fq.gz novaseq6000.fq.gz
#aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseq2500-ERR309934-shuffled-600m.fq.gz hiseq2500.fq.gz
#aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseqxten-SRR6691663-shuffled-600m.fq.gz hiseqxten.fq.gz

aws s3 cp s3://vg-k8s/profiling/reads/real/NA19239/novaseq6000-ERR3239454-shuffled-1m.fq.gz novaseq6000.fq.gz
aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseq2500-ERR309934-shuffled-1m.fq.gz hiseq2500.fq.gz
aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseqxten-SRR6691663-shuffled-1m.fq.gz hiseqxten.fq.gz
for STRAIN in DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138 ; do
    aws s3 cp s3://vg-k8s/profiling/reads/real/yeast/${STRAIN}-shuffled.fq.gz ${STRAIN}.fq.gz
done

printf "graph\talgorithm\treads\tpairing\tload_time\tspeed\ttotal_cpu_sec\ttotal_wall_clock_time\tmax_resident_set_size(kbytes)\n" > speed_report_map.tsv
for SPECIES in human yeast ; do
    case "${SPECIES}" in
    yeast)
        GRAPHS=(yeast_subset S288C)
        READSETS=(DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138)
        ;;
    human)
        GRAPHS=(hgsvc 1kg hs38d1 hs37d5)
        READSETS=(novaseq6000 hiseqxten hiseq2500)
        ;;
    esac
    for GRAPH in ${GRAPHS[@]} ; do
        # Find indexes
        case ${GRAPH} in
        1kg)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/for-NA19239/1kg/hs37d5/1kg_hs37d5_filter
            ;;
        hgsvc)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1
            ;;
        hs37d5)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/primary/hs37d5/primaryhs37d5
            ;;
        hs38d1)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2-2/generic/primary/hs38d1/primaryhs38d1
            ;;
        yeast_subset)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/cactus/yeast_subset/yeast_subset
            ;;
        S288C)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/primary/S288C/primaryS288C
            ;;
        esac
        # Download indexes
        for EXT in xg gcsa gcsa.lcp ; do
            aws s3 cp ${GRAPH_BASE}.${EXT} ./${GRAPH}.${EXT}
        done
        for READS in ${READSETS[@]} ; do

            for PAIRING in single paired ; do 
                if [[ ${PAIRING} == "single" ]] ; then
                    PAIRED=""
                elif [[ ${PAIRING} == "paired" ]] ; then
                    PAIRED="-i"
                fi
                /usr/bin/time -v bash -c "vg map -x ${GRAPH}.xg -g ${GRAPH}.gcsa -f ${READS}.fq.gz ${PAIRED} -t 16 -p --log-time 2>log.txt >mapped.gam" 2> time-log.txt || true
                LOAD_TIME="$(cat log.txt | grep "Index load time" | sed 's/Index load time:\ \([0-9]*\.[0-9]*\)/\1/g')"
                SPEED="$(cat log.txt | grep "Mapping speed" | sed 's/Mapping\ speed:\ \([0-9]*\.[0-9]*\)\ reads per second per thread/\1/g')"

                USER_TIME="$(cat "time-log.txt" | grep "User time" | sed 's/User\ time\ (seconds):\ \([0-9]*\.[0-9]*\)/\1/g')"
                SYS_TIME="$(cat "time-log.txt" | grep "System time" | sed 's/System\ time\ (seconds):\ \([0-9]*\.[0-9]*\)/\1/g')"
                TOTAL_TIME="$(echo "${USER_TIME} + ${SYS_TIME}" | bc -l)"
                CLOCK_TIME="$(cat "time-log.txt" | grep "Elapsed (wall clock) time" | sed 's/.*\ \([0-9,:]*\)/\1/g')"
                MEMORY="$(cat "time-log.txt" | grep "Maximum resident set" | sed 's/Maximum\ resident\ set\ size\ (kbytes):\ \([0-9]*\)/\1/g')"

                echo ${GRAPH} vg_map ${READS} ${PAIRING} ${LOAD_TIME} ${SPEED} ${TOTAL_TIME} ${CLOCK_TIME} ${MEMORY}
                printf "${GRAPH}\tvg_map\t${READS}\t${PAIRING}\t${LOAD_TIME}\t${SPEED}\t${TOTAL_TIME}\t${CLOCK_TIME}\t${MEMORY}\n" >> speed_report_map.tsv 
                printf "${GRAPH}\tvg_map\t${READS}\t${PAIRING}\t${LOAD_TIME}\t${SPEED}\t${TOTAL_TIME}\t${CLOCK_TIME}\t${MEMORY}\n" >> map_speed_log.txt 
                cat time-log.txt >> map_speed_log.txt
                cat log.txt >> map_speed_log.txt
            done
        done
        
        # Clean up indexes
        for EXT in xg gcsa gcsa.lcp ; do
            rm -f ./${GRAPH}.${EXT}
        done
        
    done
done

aws s3 cp speed_report_map.tsv  s3://vg-k8s/users/xhchang/giraffe_experiments/speed/speed_report_map.tsv 
aws s3 cp map_speed_log.txt s3://vg-k8s/users/xhchang/giraffe_experiments/speed/map_speed_log.txt

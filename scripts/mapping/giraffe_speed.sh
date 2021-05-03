#!/usr/bin/env bash

set -ex

THREAD_COUNT=16

#600million real reads
aws s3 cp s3://vg-k8s/profiling/reads/real/NA19239/novaseq6000-ERR3239454-shuffled-600m.fq.gz novaseq6000.fq.gz
#aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseq2500-ERR309934-shuffled-600m.fq.gz hiseq2500.fq.gz
#aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseqxten-SRR6691663-shuffled-600m.fq.gz hiseqxten.fq.gz

#1million real reads
#aws s3 cp s3://vg-k8s/profiling/reads/real/NA19239/novaseq6000-ERR3239454-shuffled-1m.fq.gz novaseq6000.fq.gz
#aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseq2500-ERR309934-shuffled-1m.fq.gz hiseq2500.fq.gz
#aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseqxten-SRR6691663-shuffled-1m.fq.gz hiseqxten.fq.gz
for STRAIN in DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138 ; do
    aws s3 cp s3://vg-k8s/profiling/reads/real/yeast/${STRAIN}-shuffled.fq.gz ${STRAIN}.fq.gz
done


printf "graph\talgorithm\treads\tpairing\tspeed\ttotal_cpu_sec\ttotal_wall_clock_time\tmax_resident_set_size(kbytes)\n" > speed_report_giraffe.tsv
for SPECIES in human yeast ; do
    case "${SPECIES}" in
    yeast)
        GRAPHS=(yeast_subset S288C)
        READSETS=(DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138)
        ;;
    human)
        GRAPHS=(hgsvc 1000gp hs38d1)
        READSETS=(novaseq6000 hiseqxten hiseq2500)
        ;; 
    esac
    for GRAPH in ${GRAPHS[@]} ; do
        case ${GRAPH} in
        1000gp)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter
            GBWTS=("sampled.64" "full")
            ;;
        hgsvc)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1
            GBWTS=("cover.16" "full")
            ;;
        hs37d5)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/primary/hs37d5/primaryhs37d5
            GBWTS=("cover")
            ;;
        hs38d1)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2-2/generic/primary/hs38d1/primaryhs38d1
            GBWTS=("cover")
            ;;
        yeast_subset)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/cactus/yeast_subset/yeast_subset
            GBWTS=(raw)
            ;;
        S288C)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/primary/S288C/primaryS288C
            GBWTS=(raw)
            ;;
        esac
        aws s3 cp ${GRAPH_BASE}.xg ./${GRAPH}.xg
        aws s3 cp ${GRAPH_BASE}.dist ./${GRAPH}.dist
        for GBWT in ${GBWTS[@]} ; do
            if [[ "${GBWT}" == "raw" ]] ; then
                aws s3 cp ${GRAPH_BASE}.gbwt ./${GRAPH}.${GBWT}.gbwt
                aws s3 cp ${GRAPH_BASE}.gg ./${GRAPH}.${GBWT}.gg
                aws s3 cp ${GRAPH_BASE}.min ./${GRAPH}.${GBWT}.min
            else
                aws s3 cp ${GRAPH_BASE}.${GBWT}.gbwt ./${GRAPH}.${GBWT}.gbwt
                aws s3 cp ${GRAPH_BASE}.${GBWT}.gg ./${GRAPH}.${GBWT}.gg
                aws s3 cp ${GRAPH_BASE}.${GBWT}.min ./${GRAPH}.${GBWT}.min
            fi
            for READS in ${READSETS[@]} ; do

                for PARAM_PRESET in default fast ; do

                    for PAIRING in single paired ; do 
                        if [[ ${PAIRING} == "single" ]] ; then
                            PAIRED=""
                        elif [[ ${PAIRING} == "paired" ]] ; then
                            PAIRED="-i"
                        fi
                        /usr/bin/time -v  bash -c "vg giraffe -x ${GRAPH}.xg -H ${GRAPH}.${GBWT}.gbwt -g ${GRAPH}.${GBWT}.gg -m ${GRAPH}.${GBWT}.min -d ${GRAPH}.dist -f ${READS}.fq.gz -b ${PARAM_PRESET} ${PAIRED} -t 16 -p 2>log.txt >mapped.gam" 2> time-log.txt || true
                        SPEED="$(cat log.txt | grep speed | sed 's/[^0-9\.]//g')"
                        USER_TIME="$(cat "time-log.txt" | grep "User time" | sed 's/User\ time\ (seconds):\ \([0-9]*\.[0-9]*\)/\1/g')"
                        SYS_TIME="$(cat "time-log.txt" | grep "System time" | sed 's/System\ time\ (seconds):\ \([0-9]*\.[0-9]*\)/\1/g')"
                        TOTAL_TIME="$(echo "${USER_TIME} + ${SYS_TIME}" | bc -l)"
                        CLOCK_TIME="$(cat "time-log.txt" | grep "Elapsed (wall clock) time" | sed 's/.*\ \([0-9,:]*\)/\1/g')"
                        MEMORY="$(cat "time-log.txt" | grep "Maximum resident set" | sed 's/Maximum\ resident\ set\ size\ (kbytes):\ \([0-9]*\)/\1/g')"


                        echo ${GRAPH} ${GBWT} ${READS} ${PAIRING} ${SPEED} ${TOTAL_TIME} ${CLOCK_TIME} ${MEMORY}
                        printf "${GRAPH}\t${GBWT}\t${READS}\t${PAIRING}\t${SPEED}\t${TOTAL_TIME}\t${CLOCK_TIME}\t${MEMORY}\n" >> giraffe_speed_log.txt
                        cat time-log.txt >> giraffe_speed_log.txt
                        cat log.txt >> giraffe_speed_log.txt
                        printf "${GRAPH}\t${GBWT}_${PARAM_PRESET}\t${READS}\t${PAIRING}\t${SPEED}\t${TOTAL_TIME}\t${CLOCK_TIME}\t${MEMORY}\n" >> speed_report_giraffe.tsv

                    done
                done
            done
        done
    done
done

aws s3 cp speed_report_giraffe.tsv s3://vg-k8s/users/xhchang/giraffe_experiments_2/liftover/nosegdups/speed/speed_report_giraffe.tsv
aws s3 cp giraffe_speed_log.txt s3://vg-k8s/users/xhchang/giraffe_experiments_2/liftover/nosegdups/speed/giraffe_speed_log.txt

sudo shutdown -h now

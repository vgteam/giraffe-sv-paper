#!/usr/bin/env bash

set -e

THREAD_COUNT=16

printf "graph\talgorithm\treads\tpairing\tload_time\tspeed\n" > speed_report_map.tsv

aws s3 cp s3://vg-k8s/profiling/reads/real/NA19239/novaseq6000-ERR3239454-shuffled-1m.fq.gz novaseq6000.fq.gz
aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseq2500-ERR309934-shuffled-1m.fq.gz hiseq2500.fq.gz
aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseqxten-SRR6691663-shuffled-1m.fq.gz hiseqxten.fq.gz
for STRAIN in DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138 ; do
    aws s3 cp s3://vg-k8s/profiling/reads/real/yeast/${STRAIN}.fq.gz ${STRAIN}.fq.gz
done

for SPECIES in human yeast ; do
    case "${SPECIES}" in
    yeast)
        GRAPHS=(SK1)
        READSETS=(DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138)
        ;;
    human)
        GRAPHS=(hgsvc 1kg)
        READSETS=(novaseq6000 hiseqxten hiseq2500)
        ;;
    esac
    for GRAPH in ${GRAPHS[@]} ; do
        case ${GRAPH} in
        1kg)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/primary/hs37d5/primaryhs37d5
            ;;
        hgsvc)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2-2/generic/primary/hs38d1/primaryhs38d1
            ;;
        SK1)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/primary/SK1/primarySK1
            ;;
        esac
        aws s3 cp ${GRAPH_BASE}.xg ./${GRAPH}.xg
        aws s3 cp ${GRAPH_BASE}.dist ./${GRAPH}.dist
        for READS in ${READSETS[@]} ; do

            for PAIRING in single paired ; do 
                if [[ ${PAIRING} == "single" ]] ; then
                    PAIRED=""
                elif [[ ${PAIRING} == "paired" ]] ; then
                    PAIRED="-i"
                fi
                vg map -x ${GRAPH}.xg -g ${GRAPH}.gcsa -f ${READS}.fq.gz ${PAIRED} -t 16 -p 2>log.txt >mapped.gam
                LOAD_TIME="$(cat log.txt | grep "Index load time" | sed 's/Index load time:\ \([0-9]*\.[0-9]*\)/\1/g')"
                SPEED="$(cat log.txt | grep "Mapping speed" | sed 's/Mapping\ speed:\ \([0-9]*\.[0-9]*\)\ reads per second per thread/\1/g')"
                
                echo ${GRAPH} vg_map ${READS} ${PAIRING} ${LOAD_TIME} ${SPEED}
                printf "${GRAPH}\tvg_map\t${READS}\t${PAIRING}\t${LOAD_TIME}\t${SPEED}\n" >> speed_report_map.tsv 
            done
        done
    done
done

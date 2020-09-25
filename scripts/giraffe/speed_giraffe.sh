#!/usr/bin/env bash

set -e

THREAD_COUNT=16

aws s3 cp s3://vg-k8s/profiling/reads/real/NA19239/novaseq6000-ERR3239454-shuffled-1m.fq.gz novaseq6000.fq.gz
aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseq2500-ERR309934-shuffled-1m.fq.gz hiseq2500.fq.gz
aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseqxten-SRR6691663-shuffled-1m.fq.gz hiseqxten.fq.gz
for STRAIN in DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138 ; do
    aws s3 cp s3://vg-k8s/profiling/reads/real/yeast/${STRAIN}.fq.gz ${STRAIN}.fq.gz
done


printf "graph\talgorithm\treads\tpairing\tspeed\n" > speed_report_giraffe.tsv
for SPECIES in human yeast ; do
    case "${SPECIES}" in
    yeast)
        GRAPHS=(yeast_subset)
        READSETS=(DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138)
        GBWTS=(raw)
        ;;
    human)
        GRAPHS=(hgsvc 1kg)
        READSETS=(novaseq6000 hiseqxten hiseq2500)
        GBWTS=("full" "sampled.64" "cover")
        ;; 
    esac
    for GRAPH in ${GRAPHS[@]} ; do
        case ${GRAPH} in
        1kg)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/for-NA19239/1kg/hs37d5/1kg_hs37d5_filter
            ;;
        hgsvc)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1
            ;;
        yeast_subset)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/cactus/yeast_all/yeast_subset
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
                    SPEED="$(vg giraffe -x ${GRAPH}.xg -H ${GRAPH}.${GBWT}.gbwt -g ${GRAPH}.${GBWT}.gg -g ${GRAPH}.${GBWT}.min -d ${GRAPH}.dist -f ${READS}.fq.gz -b ${PARAM_PRESET} ${PAIRED} -t 16 -p 2>&1 >mapped.gam | grep speed | sed 's/[^0-9\.]//g')"
                    
                    echo ${GRAPH} ${GBWT} ${READS} ${PARAM_PRESET} ${PAIRING} ${SPEED}
                    printf "${GRAPH}\t${GBWT}\t${READS}\t${PARAM_PRESET}\t${PAIRING}\t${SPEED}\n" >> speed_report_giraffe.tsv 
                done
            done
        done
    done
done

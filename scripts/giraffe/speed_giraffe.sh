#!/usr/bin/env bash

set -e

THREAD_COUNT=16

aws s3 cp s3://vg-k8s/profiling/reads/real/NA19239/novaseq6000-ERR3239454-shuffled-1m.fq.gz novaseq6000.fq.gz
aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseq2500-ERR309934-shuffled-1m.fq.gz hiseq2500.fq.gz
aws s3 cp s3://vg-k8s/profiling/reads/real/NA19240/hiseqxten-SRR6691663-shuffled-1m.fq.gz hiseqxten.fq.gz


printf "graph\talgorithm\treads\tpairing\tspeed\n" > speed_report_giraffe.tsv
for GRAPH in hgsvc 1kg ; do
   if [[ ${GRAPH} == "1kg" ]] ; then
        GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/for-NA19239/1kg/hs37d5/1kg_hs37d5_filter
    elif [[ ${GRAPH} == "hgsvc" ]] ; then
        GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1
    fi
    aws s3 cp ${GRAPH_BASE}.xg ./${GRAPH}.xg
    aws s3 cp ${GRAPH_BASE}.dist ./${GRAPH}.dist
    for GBWT in full sampled cover ; do
        aws s3 cp ${GRAPH_BASE}.${GBWT}.gbwt ./${GRAPH}.${GBWT}.gbwt
        aws s3 cp ${GRAPH_BASE}.${GBWT}.gg ./${GRAPH}.${GBWT}.graph.gg
        for READS in novaseq6000 hiseqxten hiseq2500 ; do

            for PAIRING in single paired ; do 
                if [[ ${PAIRING} == "single" ]] ; then
                    PAIRED=""
                elif [[ ${PAIRING} == "paired" ]] ; then
                    PAIRED="-i"
                fi
                SPEED="$(vg giraffe -x ${GRAPH}.xg -H ${GRAPH}.${GBWT}.gbwt -g ${GRAPH}.${GBWT}.gg  -d ${GRAPH}.dist -f ${READS}.fq.gz ${PAIRED} -t 16 -p 2>&1 >mapped.gam | grep speed | sed 's/[^0-9\.]//g')"
                
                echo ${GRAPH} ${GBWT} ${READS} ${PAIRING} ${SPEED}
                printf "${GRAPH}\t${GBWT}\t${READS}\t${PAIRING}\t${SPEED}\n" >> speed_report_giraffe.tsv 
            done
        done
    done
done

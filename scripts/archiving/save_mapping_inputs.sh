#!/usr/bin/env bash

# save_mapping_inputs.sh: Save graphs and simulated reads used in mapping experiments

set -x

DEST_DIR=/nanopore/cgl/data/giraffe

for SPECIES in human yeast ; do
    case "${SPECIES}" in
    yeast)
        GRAPHS=(yeast_all yeast_subset S288C)
        ;;
    human)
        GRAPHS=(hgsvc 1kg 1kg-unfiltered 1kg-primary hgsvc-primary)
        ;;
    esac
    for GRAPH in ${GRAPHS[@]} ; do
        case "${SPECIES}" in
        yeast)
            GBWTS=(raw)
            ;;
        human)
            GBWTS=("raw" "full" "sampled.64" "cover" "paths")
            ;;
        esac
    
        # Most graphs simulated no reads
        READSETS=()
    
        case ${GRAPH} in
        1kg-unfiltered)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/for-NA19239/1kg/hs37d5/1kg_hs37d5
            DEST_BASE=${DEST_DIR}/mapping/graphs/for-NA19239/1kg/hs37d5/1kg_hs37d5
            # We simulated some reads
            READSETS=(novaseq6000 hiseqxten hiseq2500)
        1kg-primary)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/primary/hs37d5/primaryhs37d5
            DEST_BASE=${DEST_DIR}/mapping/graphs/generic/primary/hs37d5/primaryhs37d5
            # Grab all the other sampled GBWTs too
            GBWTS+=("sampled.1" "sampled.2" "sampled.4" "sampled.8" "sampled.16" "sampled.32" "sampled.128")
            ;;
        hgsvc-primary)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2-2/generic/primary/hs38d1/primaryhs38d1
            DEST_BASE=${DEST_DIR}/mapping/graphs/generic/primary/hs38d1/primaryhs38d1
            ;;
        1kg)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/for-NA19239/1kg/hs37d5/1kg_hs37d5_filter
            DEST_BASE=${DEST_DIR}/mapping/graphs/for-NA19239/1kg/hs37d5/1kg_hs37d5_filter
            ;;
        hgsvc)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1
            DEST_BASE=${DEST_DIR}/mapping/graphs/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1
            # Grab all the cover GBWTs too
            GBWTS+=("cover.1" "cover.2" "cover.4" "cover.8" "cover.16" "cover.32" "cover.64" "cover.128")
            # We simulated some reads
            READSETS=(novaseq6000 hiseqxten hiseq2500)
            ;;
        S288C)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/primary/S288C/primaryS288C
            DEST_BASE=${DEST_DIR}/mapping/graphs/generic/primary/S288C/primaryS288C
            ;;
        yeast_all)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/cactus/yeast_all/yeast_all
            DEST_BASE=${DEST_DIR}/mapping/graphs/generic/cactus/yeast_all/yeast_all
            # We simulated some reads
            READSETS=(DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138)
            ;;
        yeast_subset)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/cactus/yeast_subset/yeast_subset
            DEST_BASE=${DEST_DIR}/mapping/graphs/generic/cactus/yeast_subset/yeast_subset
            ;;
        esac
        
        for EXT in vg xg snarls trivial.snarls gcsa gcsa.lcp dist ; do
            # Copy the graph and all generic indexes
            aws s3 cp ${GRAPH_BASE}.${EXT} ${DEST_BASE}.${EXT}
        done
        
        
        for GBWT in ${GBWTS[@]} ; do
            # Copy each GBWT and related file that exists
            if [[ "${GBWT}" == "raw" ]] ; then
                aws s3 cp ${GRAPH_BASE}.gbwt ${DEST_BASE}.gbwt
                aws s3 cp ${GRAPH_BASE}.gg ${DEST_BASE}.gg
                aws s3 cp ${GRAPH_BASE}.min ${DEST_BASE}.min
            else
                aws s3 cp ${GRAPH_BASE}.${GBWT}.gbwt ${DEST_BASE}.${GBWT}.gbwt
                aws s3 cp ${GRAPH_BASE}.${GBWT}.gg ${DEST_BASE}.${GBWT}.gg
                aws s3 cp ${GRAPH_BASE}.${GBWT}.min ${DEST_BASE}.${GBWT}.min
            fi
        done
        
        for READS in ${READSETS[@]} ; do
            # Copy each readset simulated from this graph
            
            case ${GRAPH} in
            1kg-unfiltered)
                READ_BASE=s3://vg-k8s/profiling/reads/sim/for-NA19239/1kg/hs37d5/${READS}/out_sim_gbwt/sim
                READDEST_BASE=${DEST_DIR}/mapping/reads/sim/for-NA19239/1kg/hs37d5/${READS}/out_sim_gbwt/sim
            1kg)
                READ_BASE=""
                READDEST_BASE=""
                ;;
            hgsvc)
                READ_BASE=s3://vg-k8s/profiling/reads/sim/for-NA19240/hgsvc/grch38/${READS}/out_sim_gbwt/sim
                READDEST_BASE=${DEST_DIR}/mapping/reads/sim/for-NA19240/hgsvc/grch38/${READS}/out_sim_gbwt/sim
                ;;
            S288C)
                READ_BASE=""
                READDEST_BASE=""
                ;;
            yeast_all)
                READ_BASE=s3://vg-k8s/profiling/reads/sim/yeast/sim-${READS}
                READDEST_BASE=${DEST_DIR}/mapping/reads/sim/yeast/sim-${READS}
                ;;
            yeast_subset)
                READ_BASE=""
                READDEST_BASE=""
                ;;
            esac
            
            # Do the read copy 
            aws s3 cp ${READ_BASE}.gam ${READDEST_BASE}.gam
        done
    done
done
            

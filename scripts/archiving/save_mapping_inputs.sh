#!/usr/bin/env bash

# save_mapping_inputs.sh: Save graphs and simulated reads used in mapping experiments

set -x

DEST_DIR=/nanopore/cgl/data/giraffe

function download() {
    if [ ! -e "${2}" ] ; then
        aws s3 cp --no-progress "${1}" "${2}"
    fi
}

function download_if_exists() {
    if [ ! -e "${2}" ] ; then
        aws s3 cp --no-progress "${1}" "${2}" || true
    fi
}

function wget_download() {
    if [ ! -e "${2}" ] ; then
        wget "${1}" -O "${2}"
    fi
}

for SPECIES in human yeast ; do
    case "${SPECIES}" in
    yeast)
        GRAPHS=(yeast_all yeast_subset S288C)
        ;;
    human)
        GRAPHS=(hgsvc 1000gplons 1000gplons-unfiltered 1000gplons-sample 1000gplons-sample-withref 1000gplons-formap primary)
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
        # Most graphs have all the indexes except GCSA/LCP and nontrivial snarls (not really used)
        EXTS=(vg xg trivial.snarls gbwt dist)
        # Some graphs have accompanied reference dictionaries for surjection 
        SEQDICTS=()
        case ${GRAPH} in
        1000gplons-unfiltered)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1
            DEST_BASE=${DEST_DIR}/mapping/graphs/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1
            GBWTS=(raw-gbwt)
            # We simulated some reads
            READSETS=(novaseq6000 hiseqxten hiseq2500)
            ;;
        1000gplons-sample)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_NA19239_sample
            DEST_BASE=${DEST_DIR}/mapping/graphs/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_NA19239_sample
            GBWTS=()
            ;;
        1000gplons-sample-withref)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_NA19239_sample_withref
            DEST_BASE=${DEST_DIR}/mapping/graphs/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_NA19239_sample_withref
            EXTS=(vg xg)
            GBWTS=(raw force force.augment)
            ;;
        1000gplons-formap)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v4/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter
            DEST_BASE=${DEST_DIR}/mapping/graphs/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter_forvgmap
            EXTS=(vg xg gcsa gcsa.lcp)
            GBWTS=(raw)
            ;;
        primary)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2-2/generic/primary/hs38d1/primaryhs38d1
            DEST_BASE=${DEST_DIR}/mapping/graphs/generic/primary/hs38d1/primaryhs38d1
            # We set up for map
            EXTS+=(gcsa gcsa.lcp)
            GBWTS=(paths cover)
            SEQDICTS=("https://storage.googleapis.com/cmarkell-vg-wdl-dev/giraffe_manuscript_data/genome_references/linear_references/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.dict")
            ;;
        1000gplons)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter
            DEST_BASE=${DEST_DIR}/mapping/graphs/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter
            # Grab all the other sampled GBWTs too
            GBWTS+=("sampled.1" "sampled.2" "sampled.4" "sampled.8" "sampled.16" "sampled.32" "sampled.128")
            SEQDICTS=("https://storage.googleapis.com/cmarkell-vg-wdl-dev/giraffe_manuscript_data/genome_references/linear_references/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.no_segdup.dict")
            ;;
        hgsvc)
            GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1
            DEST_BASE=${DEST_DIR}/mapping/graphs/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1
            # Grab all the cover GBWTs too
            GBWTS+=("cover.1" "cover.2" "cover.4" "cover.8" "cover.16" "cover.32" "cover.64" "cover.128" "onlyNA19240.augment")
            # We simulated some reads
            READSETS=(novaseq6000 hiseqxten hiseq2500)
            # We set up for map
            EXTS+=(gcsa gcsa.lcp)
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
        
        for SEQ in "${SEQDICTS[@]}" ; do
            base_filename=$(basename ${SEQ})
            wget_download ${SEQ} ${DEST_DIR}/mapping/${base_filename}
        done
        
        for EXT in "${EXTS[@]}" ; do
            # Copy the graph and all generic indexes
            download ${GRAPH_BASE}.${EXT} ${DEST_BASE}.${EXT}
        done
        
        
        for GBWT in ${GBWTS[@]} ; do
            # Copy each GBWT and related file that exists
            if [[ "${GBWT}" == "raw" ]] ; then
                download ${GRAPH_BASE}.gbwt ${DEST_BASE}.gbwt
                download_if_exists ${GRAPH_BASE}.gg ${DEST_BASE}.gg
                download_if_exists ${GRAPH_BASE}.gg ${DEST_BASE}.min
            else
                # Generated sampled GBWTs come with a min file
                download ${GRAPH_BASE}.${GBWT}.gbwt ${DEST_BASE}.${GBWT}.gbwt
                download_if_exists ${GRAPH_BASE}.${GBWT}.gg ${DEST_BASE}.${GBWT}.gg
                download_if_exists ${GRAPH_BASE}.${GBWT}.min ${DEST_BASE}.${GBWT}.min
            fi
        done
        
        for READS in ${READSETS[@]} ; do
            # Copy each readset simulated from this graph
            
            case ${GRAPH} in
            1000gplons-unfiltered)
                READ_BASE=s3://vg-k8s/profiling/reads/sim/for-NA19239/1000gp/hs38d1/liftover_nosegdups/${READS}/out_sim_gbwt/sim
                READDEST_BASE=${DEST_DIR}/mapping/reads/sim/for-NA19239/1000gplons/hs38d1/${READS}/out_sim_gbwt/sim
                ;;
            hgsvc)
                READ_BASE=s3://vg-k8s/profiling/reads/sim/for-NA19240/hgsvc/grch38/${READS}/out_sim_gbwt/sim
                READDEST_BASE=${DEST_DIR}/mapping/reads/sim/for-NA19240/hgsvc/grch38/${READS}/out_sim_gbwt/sim
                ;;
            yeast_all)
                READ_BASE=s3://vg-k8s/profiling/reads/sim/yeast/sim-${READS}
                READDEST_BASE=${DEST_DIR}/mapping/reads/sim/yeast/sim-${READS}
                ;;
            *)
                READ_BASE=""
                READDEST_BASE=""
                ;;
            esac
            
            # Do the read copy 
            download ${READ_BASE}.gam ${READDEST_BASE}.gam
        done
    done
done

chmod -R g+rw "${DEST_DIR}" 2>/dev/null || true

#!/usr/bin/env bash
# plot-rocs.sh: plot ROC curves for Giraffe and compared mappers

# Where should the ROCs end up
WORKDIR="$HOME/build/vg/trash/rocs"

# Where are the data TSVs from the Kubernetes scripts?
STAT_URL="s3://vg-k8s/users/adamnovak/giraffe_experiments"

# Where are we?
SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"

mkdir -p "${WORKDIR}"

if [ ! -e "${WORKDIR}/stats" ] ; then
    mkdir -p "${WORKDIR}/stats"
    aws s3 cp ${STAT_URL}/roc_stats_giraffe.tsv "${WORKDIR}/stats/"
    aws s3 cp ${STAT_URL}/roc_stats_giraffe_primary.tsv "${WORKDIR}/stats/"
    aws s3 cp ${STAT_URL}/roc_stats_map.tsv "${WORKDIR}/stats/"
    aws s3 cp ${STAT_URL}/roc_stats_map_primary.tsv "${WORKDIR}/stats/"
    aws s3 cp ${STAT_URL}/roc_stats_bowtie2_primary.tsv "${WORKDIR}/stats/"
    aws s3 cp ${STAT_URL}/roc_stats_minimap2_primary.tsv "${WORKDIR}/stats/" 
    aws s3 cp ${STAT_URL}/roc_stats_bwa_primary.tsv "${WORKDIR}/stats/"
fi

for SPECIES in human yeast ; do
    case "${SPECIES}" in
    yeast)
        GRAPHS=(S288C yeast_all yeast_subset)
        READSETS=(DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138)
        GBWT="raw"
        ;;
    human)
        GRAPHS=(hgsvc 1kg)
        READSETS=(novaseq6000 hiseqxten hiseq2500)
        GBWT="sampled"
        ;;
    esac
    for GRAPH in ${GRAPHS[@]} ; do
        for READS in ${READSETS[@]} ; do
            if [ ! -e "${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}.tsv" ] ; then
                echo "Extracting ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}.tsv"
                cat ${WORKDIR}/stats/roc_stats_*.tsv | head -n1 > ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}.tsv
                tail -q -n +2 ${WORKDIR}/stats/roc_stats_*.tsv | grep -P "${GRAPH}(${GBWT})?(map)?(${READS})" >> ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}.tsv
                wc -l ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}.tsv
            fi
            Rscript ${SCRIPT_DIR}/plot-roc.R ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}.tsv ${WORKDIR}/roc-plot-${SPECIES}-${GRAPH}-${READS}.png
        done
    done
done

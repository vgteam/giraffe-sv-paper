#!/usr/bin/env bash
# plot-rocs.sh: plot ROC curves for Giraffe and compared mappers

set -ex

# Where should the ROCs end up
WORKDIR="$HOME/build/vg/trash/rocs"

# Where are the data TSVs from the Kubernetes scripts?
STAT_URL="s3://vg-k8s/users/adamnovak/giraffe_experiments"

# Where are we?
SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"

mkdir -p "${WORKDIR}"

mkdir -p "${WORKDIR}/stats"
for STAT_FILE in roc_stats_giraffe.tsv roc_stats_giraffe_primary.tsv roc_stats_map.tsv roc_stats_map_primary.tsv roc_stats_bowtie2_primary.tsv roc_stats_minimap2_primary.tsv roc_stats_bwa_primary.tsv ; do
    if [ ! -e "${WORKDIR}/stats/${STAT_FILE}" ] ; then
        aws s3 cp "${STAT_URL}/${STAT_FILE}" "${WORKDIR}/stats/${STAT_FILE}"
    fi
done

for SPECIES in yeast human ; do
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
    for READS in ${READSETS[@]} ; do
        # Do human-designed plots
        if [ "${SPECIES}" == "yeast" ] ; then
            for PAIRING in single paired ; do
                if [ "${PAIRING}" == "paired" ] ; then
                    PE_OPTS="-- -pe"
                else
                    PE_OPTS="-v -- -pe"
                fi
                if [ ! -e "${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv" ] ; then
                    echo "Extracting ${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv"
                    cat ${WORKDIR}/stats/roc_stats_*.tsv | head -n1 > ${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv
                    tail -q -n +2 ${WORKDIR}/stats/roc_stats_*.tsv | grep -P "(yeast_subset(${GBWT})?${READS}|S288C${READS})" | grep ${PE_OPTS} >> ${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv
                    wc -l ${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv
                fi
                if [ ! -e "${WORKDIR}/roc-plot-${SPECIES}-overall-${READS}.png" ] ; then
                    Rscript ${SCRIPT_DIR}/plot-roc.R ${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv ${WORKDIR}/roc-plot-${SPECIES}-overall-${READS}-${PAIRING}.png
                fi
            done
        fi
        
        # Do boring plots
        for GRAPH in ${GRAPHS[@]} ; do
            if [ ! -e "${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}.tsv" ] ; then
                echo "Extracting ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}.tsv"
                cat ${WORKDIR}/stats/roc_stats_*.tsv | head -n1 > ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}.tsv
                tail -q -n +2 ${WORKDIR}/stats/roc_stats_*.tsv | grep -P "${GRAPH}(${GBWT})?(${READS})" >> ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}.tsv
                wc -l ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}.tsv
            fi
            if [ ! -e "${WORKDIR}/roc-plot-${SPECIES}-${GRAPH}-${READS}.png" ] ; then
                Rscript ${SCRIPT_DIR}/plot-roc.R ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}.tsv ${WORKDIR}/roc-plot-${SPECIES}-${GRAPH}-${READS}.png
            fi
        done
    done
done

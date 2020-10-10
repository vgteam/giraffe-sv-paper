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

# Replace all names of mappers with human-readable ones
function humanize_names() {
    sed -e 's/[a-zA-Z0-9_]*bwa_mem[a-zA-Z0-9_]*/BWA/' -e 's/[a-zA-Z0-9_]*giraffe_default[a-zA-Z0-9_]*/Giraffe/' -e 's/[a-zA-Z0-9_]*giraffe_fast[a-zA-Z0-9_]*/Giraffe-Fast/' -e 's/[a-zA-Z0-9_]*map_[a-zA-Z0-9_]*/Map/'
}

for SPECIES in yeast ; do
    case "${SPECIES}" in
    yeast)
        GRAPHS=(S288C yeast_all yeast_subset)
        #READSETS=(DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138)
        READSETS=(DBVPG6044)
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
            # Yeast plots of all the subset and linear graphs together
            for PAIRING in single paired ; do
                if [ "${PAIRING}" == "paired" ] ; then
                    PE_OPTS="-- -pe"
                else
                    PE_OPTS="-v -- -pe"
                fi
                if [ ! -e "${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv" ] ; then
                    echo "Extracting ${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv"
                    cat ${WORKDIR}/stats/roc_stats_*.tsv | head -n1 > ${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv
                    # Grab all the subset and linear graph reads
                    tail -q -n +2 ${WORKDIR}/stats/roc_stats_*.tsv | grep -P "(yeast_subset(${GBWT})?${READS}|S288C(${GBWT})?${READS})" | grep ${PE_OPTS} >> ${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv
                    wc -l ${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv
                fi
                if [ ! -e "${WORKDIR}/roc-plot-${SPECIES}-overall-${READS}-${PAIRING}.png" ] ; then
                    Rscript ${SCRIPT_DIR}/plot-roc.R ${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv ${WORKDIR}/roc-plot-${SPECIES}-overall-${READS}-${PAIRING}.png
                fi
            done
        fi
        
        # Plots of giraffe normal and fast graph, map graph, and bwa linear
        for PAIRING in single paired ; do
            if [ "${PAIRING}" == "paired" ] ; then
                PE_OPTS="-- -pe"
            else
                PE_OPTS="-v -- -pe"
            fi
            if [ ! -e "${WORKDIR}/toplot-${SPECIES}-headline-${READS}-${PAIRING}.tsv" ] ; then
                echo "Extracting ${WORKDIR}/toplot-${SPECIES}-headline-${READS}-${PAIRING}.tsv"
                cat ${WORKDIR}/stats/roc_stats_*.tsv | head -n1 > ${WORKDIR}/toplot-${SPECIES}-headline-${READS}-${PAIRING}.tsv
                # Grab linear BWA
                tail -q -n +2 ${WORKDIR}/stats/roc_stats_*.tsv | grep -P "bwa" | grep "${READS}" | grep ${PE_OPTS} | humanize_names >> ${WORKDIR}/toplot-${SPECIES}-headline-${READS}-${PAIRING}.tsv
                # Grab giraffe and map non-linear
                tail -q -n +2 ${WORKDIR}/stats/roc_stats_*.tsv | grep -P "(giraffe_default|giraffe_fast|bwa_mem)" | grep "${READS}" | grep -v "_primary" | grep ${PE_OPTS} | humanize_names >> ${WORKDIR}/toplot-${SPECIES}-headline-${READS}-${PAIRING}.tsv
                
                wc -l ${WORKDIR}/toplot-${SPECIES}-headline-${READS}-${PAIRING}.tsv
            fi
            if [ ! -e "${WORKDIR}/roc-plot-${SPECIES}-headline-${READS}-${PAIRING}.png" ] ; then
                Rscript ${SCRIPT_DIR}/plot-roc.R ${WORKDIR}/toplot-${SPECIES}-headline-${READS}-${PAIRING}.tsv ${WORKDIR}/roc-plot-${SPECIES}-headline-${READS}-${PAIRING}.png
            fi
        done
        
    done
    
    continue
    
    for READS in ${READSETS[@]} ; do
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

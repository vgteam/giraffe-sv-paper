#!/usr/bin/env bash
# plot-rocs.sh: plot ROC curves for Giraffe and compared mappers

set -ex

# Where should the ROCs end up
WORKDIR="$HOME/build/vg/trash/rocs"

# Where are the data TSVs from the Kubernetes scripts?
STAT_URL="s3://vg-k8s/users/adamnovak/giraffe_experiments"
# If it's not there, look here instead
BACKUP_STAT_URL="s3://vg-k8s/users/xhchang/giraffe_experiments"

# Where are we?
SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"

mkdir -p "${WORKDIR}"

mkdir -p "${WORKDIR}/stats"
for STAT_FILE in roc_stats_giraffe.tsv roc_stats_giraffe_primary.tsv roc_stats_map.tsv roc_stats_map_primary.tsv roc_stats_bowtie2_primary.tsv roc_stats_minimap2_primary.tsv roc_stats_bwa_primary.tsv roc_stats_hisat2_1kg_novaseq6000.tsv roc_stats_hisat2_hgsvc_novaseq6000.tsv roc_stats_hisat2_1kg_hiseqxten.tsv roc_stats_hisat2_hgsvc_hiseqxten.tsv roc_stats_hisat2_1kg_hiseq2500.tsv roc_stats_hisat2_hgsvc_hiseq2500.tsv roc_stats_graphaligner.tsv ; do
    if [ ! -e "${WORKDIR}/stats/${STAT_FILE}" ] ; then
        aws s3 cp "${STAT_URL}/${STAT_FILE}" "${WORKDIR}/stats/${STAT_FILE}" || \
        aws s3 cp "${BACKUP_STAT_URL}/${STAT_FILE}" "${WORKDIR}/stats/${STAT_FILE}" || \
        (aws s3 cp "${BACKUP_STAT_URL}/${STAT_FILE}.gz" "${WORKDIR}/stats/${STAT_FILE}.gz" && gunzip "${WORKDIR}/stats/${STAT_FILE}.gz")
    fi
done

# Replace all names of mappers with human-readable ones
function humanize_names() {
    sed -e 's/[a-zA-Z0-9_.]*bwa_mem[a-zA-Z0-9_.]*/BWA/' -e 's/[a-zA-Z0-9_.]*bowtie2[a-zA-Z0-9_.]*/Bowtie2/' -e 's/[a-zA-Z0-9_.]*minimap2[a-zA-Z0-9_.]*/Minimap2/' -e 's/[a-zA-Z0-9_.]*hisat2[a-zA-Z0-9_.-]*-/Hisat2-/' -e 's/[a-zA-Z0-9_.]*giraffe_default[a-zA-Z0-9_.]*/Giraffe/' -e 's/[a-zA-Z0-9_.]*giraffe_fast[a-zA-Z0-9_.]*/GiraffeFast/' -e 's/[a-zA-Z0-9_.]*giraffe_primary[a-zA-Z0-9_.]*/GiraffePrimary/' -e 's/[a-zA-Z0-9_.]*map_primary[a-zA-Z0-9_.]*/MapPrimary/' -e 's/[a-zA-Z0-9_.]*map_[a-zA-Z0-9_.]*/Map/' -e 's/[a-zA-Z0-9_.]*graphaligner[a-zA-Z0-9_.]*/GraphAligner/'
}

for SPECIES in human yeast ; do
    case "${SPECIES}" in
    yeast)
        GRAPHS=(S288C yeast_all yeast_subset)
        HEADLINE_GRAPHS=(yeast_subset)
        LINEAR_GRAPH="S288C"
        #READSETS=(DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138)
        READSETS=(DBVPG6044)
        GBWT="raw"
        ;;
    human)
        GRAPHS=(hgsvc 1kg)
        HEADLINE_GRAPHS=(hgsvc 1kg)
        LINEAR_GRAPH="NOTAPPLICABLE"
        READSETS=(novaseq6000 hiseqxten hiseq2500)
        # Must be in regex form
        GBWT='sampled\.64'
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
                    # Grab all the subset and linear graph reads, except map linear and hisat2 which won't exist
                    tail -q -n +2 ${WORKDIR}/stats/roc_stats_bwa*.tsv ${WORKDIR}/stats/roc_stats_giraffe*.tsv ${WORKDIR}/stats/roc_stats_bowtie2*.tsv ${WORKDIR}/stats/roc_stats_map.tsv ${WORKDIR}/stats/roc_stats_minimap2*.tsv | grep ${PE_OPTS} | grep -v "map_primary" | grep -P "(yeast_subset(${GBWT})?${READS}|S288C(${GBWT})?${READS})" | humanize_names >> ${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv
                fi
                
                cat ${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv | cut -f4 | uniq -c
                
                if [ ! -e "${WORKDIR}/roc-plot-${SPECIES}-overall-${READS}-${PAIRING}.png" ] ; then
                    Rscript ${SCRIPT_DIR}/plot-roc-comparing-aligners.R ${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv ${WORKDIR}/roc-plot-${SPECIES}-overall-${READS}-${PAIRING}.png
                fi
                if [ ! -e "${WORKDIR}/qq-plot-${SPECIES}-overall-${READS}-${PAIRING}.png" ] ; then
                    Rscript ${SCRIPT_DIR}/plot-qq-comparing-aligners.R ${WORKDIR}/toplot-${SPECIES}-overall-${READS}-${PAIRING}.tsv ${WORKDIR}/qq-plot-${SPECIES}-overall-${READS}-${PAIRING}.png
                fi
            done
        fi
        
        # Plots of giraffe normal and fast graph, map graph, and common linear mappers
        for GRAPH in ${HEADLINE_GRAPHS[@]} ; do
            # Hisat2 names its graphs weird
            case "${GRAPH}" in
            hgsvc)
                HISAT_COND="grch38_hgsvc_all-def"
                ;;
            1kg)
                HISAT_COND="grch37_snp-def"
                ;;
            *)
                HISAT_COND="NOTAPPLICABLE"
                ;;
            esac
            for PAIRING in single paired ; do
                if [ "${PAIRING}" == "paired" ] ; then
                    PE_OPTS="-- -pe"
                else
                    PE_OPTS="-v -- -pe"
                fi
                if [ ! -e "${WORKDIR}/toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv" ] ; then
                    echo "Extracting ${WORKDIR}/toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv"
                    cat ${WORKDIR}/stats/roc_stats_*.tsv | head -n1 > ${WORKDIR}/toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv
                    # Grab linear BWA for this graph and these reads
                    tail -q -n +2 ${WORKDIR}/stats/roc_stats_bwa*.tsv | grep -P "(${GRAPH}(${GBWT})?${READS}|${LINEAR_GRAPH}(${GBWT})?${READS})" | grep ${PE_OPTS} | humanize_names >> ${WORKDIR}/toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv
                    
                    # Grab giraffe and map and graphaligner non-linear, and other linear mappers
                    tail -q -n +2 ${WORKDIR}/stats/roc_stats_giraffe*.tsv ${WORKDIR}/stats/roc_stats_bowtie2*.tsv ${WORKDIR}/stats/roc_stats_map.tsv ${WORKDIR}/stats/roc_stats_minimap2*.tsv ${WORKDIR}/stats/roc_stats_graphaligner.tsv | grep ${PE_OPTS} | grep -v "_primary" | grep -P "(${GRAPH}(${GBWT})?${READS}|${LINEAR_GRAPH}(${GBWT})?${READS})" | sed 's/null/0/g' | humanize_names >> ${WORKDIR}/toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv
                    
                    # hisat2 doesn't have proper graph or read set names in its entries; it is split up by file instead
                    if [ -e "${WORKDIR}/stats/roc_stats_hisat2_${GRAPH}_${READS}.tsv" ] ; then
                        # Grab it too
                        tail -q -n +2 "${WORKDIR}/stats/roc_stats_hisat2_${GRAPH}_${READS}.tsv" | grep ${PE_OPTS} | grep "${HISAT_COND}" | humanize_names | sed 's/-se//' >> ${WORKDIR}/toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv
                    fi
                    
                    cat ${WORKDIR}/toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv | cut -f4 | uniq -c
                    
                fi
                if [ ! -e "${WORKDIR}/roc-plot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.png" ] ; then
                    Rscript ${SCRIPT_DIR}/plot-roc-comparing-aligners.R ${WORKDIR}/toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv ${WORKDIR}/roc-plot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.png
                fi
                if [ ! -e "${WORKDIR}/qq-plot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.png" ] ; then
                    Rscript ${SCRIPT_DIR}/plot-qq-comparing-aligners.R ${WORKDIR}/toplot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.tsv ${WORKDIR}/qq-plot-${SPECIES}-headline_${GRAPH}-${READS}-${PAIRING}.png
                fi
            done
        done
        
    done
    
    for READS in ${READSETS[@]} ; do
        # Do boring plots
        for GRAPH in ${GRAPHS[@]} ; do
            # Hisat2 names its graphs weird
            case "${GRAPH}" in
            hgsvc)
                HISAT_COND="grch38_hgsvc_all-def"
                ;;
            1kg)
                HISAT_COND="grch37_snp-def"
                ;;
            *)
                HISAT_COND="NOTAPPLICABLE"
                ;;
            esac
            for PAIRING in single paired ; do
                if [ "${PAIRING}" == "paired" ] ; then
                    PE_OPTS="-- -pe"
                else
                    PE_OPTS="-v -- -pe"
                fi
                if [ ! -e "${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}.tsv" ] ; then
                    echo "Extracting ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}-${PAIRING}.tsv"
                    cat ${WORKDIR}/stats/roc_stats_*.tsv | head -n1 > ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}-${PAIRING}.tsv
                    # Grab everything but map linear and hisat2
                    tail -q -n +2 ${WORKDIR}/stats/roc_stats_bwa*.tsv ${WORKDIR}/stats/roc_stats_giraffe*.tsv ${WORKDIR}/stats/roc_stats_bowtie2*.tsv ${WORKDIR}/stats/roc_stats_map.tsv ${WORKDIR}/stats/roc_stats_minimap2*.tsv  ${WORKDIR}/stats/roc_stats_graphaligner.tsv | grep ${PE_OPTS} | grep -P "${GRAPH}(${GBWT})?(${READS})" | grep -v "map_primary" | sed 's/null/0/g' | humanize_names >> ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}-${PAIRING}.tsv
                    
                    # hisat2 doesn't have proper graph or read set names in its entries; it is split up by file instead
                    if [ -e "${WORKDIR}/stats/roc_stats_hisat2_${GRAPH}_${READS}.tsv" ] ; then
                        # Grab it too
                        tail -q -n +2 "${WORKDIR}/stats/roc_stats_hisat2_${GRAPH}_${READS}.tsv" | grep ${PE_OPTS} | grep "${HISAT_COND}" | humanize_names | sed 's/-se//' >> ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}-${PAIRING}.tsv
                    fi
                    
                    cat ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}-${PAIRING}.tsv | cut -f4 | uniq -c
                    
                fi
                if [ ! -e "${WORKDIR}/roc-plot-${SPECIES}-${GRAPH}-${READS}-${PAIRING}.png" ] ; then
                    Rscript ${SCRIPT_DIR}/plot-roc-comparing-aligners.R ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}-${PAIRING}.tsv ${WORKDIR}/roc-plot-${SPECIES}-${GRAPH}-${READS}-${PAIRING}.png
                fi
                if [ ! -e "${WORKDIR}/qq-plot-${SPECIES}-${GRAPH}-${READS}-${PAIRING}.png" ] ; then
                    Rscript ${SCRIPT_DIR}/plot-qq-comparing-aligners.R ${WORKDIR}/toplot-${SPECIES}-${GRAPH}-${READS}-${PAIRING}.tsv ${WORKDIR}/qq-plot-${SPECIES}-${GRAPH}-${READS}-${PAIRING}.png
                fi
            done
        done
    done
done


#Plot all mappers for human graphs
#
#MAPPERS=(bowtie2 bwa giraffe giraffe_fast giraffe_primary graphaligner hisat2 map minimap2)
#
#for GRAPH in hgsvc ; do
#    if [ "${GRAPH}" == "1kg" ] ; then
#        PRIMARY_GRAPH="hs37d5"
#        HISAT2_GRAPH="grch37_snp"
#    elif [ "${GRAPH}" == "hgsvc" ] ; then
#        PRIMARY_GRAPH="hs38d1"
#        HISAT2_GRAPH="grch38_hgsvc_all"
#    fi
#    for READS in hiseqxten hiseq2500 ; do
#        for PAIRING in single paired ; do
#            printf "correct\tmq\tscore\taligner\n" > roc_stats_${GRAPH}_${READS}_${PAIRING}.tsv
#
#            if [ "${PAIRING}" == "single" ] ; then
#                PAIRING_TAG=""
#                PAIRING_TAG_HISAT="-se"
#            else
#                PAIRING_TAG="-pe"
#                PAIRING_TAG_HISAT="-pe"
#            fi
#            for MAPPER in ${MAPPERS[@]} ; do
#                case "${MAPPER}" in
#                bowtie2)
#                    INPUT_FILE="./roc_stats_bowtie2_primary.tsv"
#                    INPUT_LABEL="${MAPPER}_${GRAPH}${READS}${PAIRING_TAG}"
#                    ;;
#                bwa)
#                    INPUT_FILE="./roc_stats_bwa_primary.tsv"
#                    INPUT_LABEL="${MAPPER}_${GRAPH}${READS}${PAIRING_TAG}"
#                    ;;
#                giraffe)
#                    INPUT_FILE="../roc_stats_giraffe.tsv"
#                    if [ "${GRAPH}" == "1kg" ] ; then
#                        INPUT_LABEL="giraffe_${GRAPH}sampled\.64${READS}${PAIRING_TAG}"
#                    elif [ "${GRAPH}" == "hgsvc" ] ; then
#                        INPUT_LABEL="giraffe_${GRAPH}full${READS}${PAIRING_TAG}"
#                    fi
#                    ;;
#                giraffe_fast)
#                    INPUT_FILE="../roc_stats_giraffe_fast.tsv"
#                    INPUT_LABEL="giraffe_${GRAPH}sampled\.64fast${READS}${PAIRING_TAG}"
#                    ;;
#                giraffe_primary)
#                    INPUT_FILE="../roc_stats_giraffe_primary.tsv"
#                    INPUT_LABEL="giraffe_${PRIMARY_GRAPH}cover${READS}${PAIRING_TAG}"
#                    ;;
#                graphaligner)
#                    INPUT_FILE="../roc_stats_graphaligner.tsv"
#                    INPUT_LABEL="${MAPPER}_${GRAPH}${READS}${PAIRING_TAG}"
#                    ;;
#                hisat2)
#                    INPUT_FILE="../hisat2_${GRAPH}_${READS}.txt"
#                    INPUT_LABEL="${MAPPER}-${HISAT2_GRAPH}-def${PAIRING_TAG_HISAT}"
#                    ;;
#                map)
#                    INPUT_FILE="../roc_stats_map_${READS}.tsv"
#                    INPUT_LABEL="${MAPPER}_${GRAPH}${READS}${PAIRING_TAG}"
#                    ;;
#                minimap2)
#                    INPUT_FILE="./roc_stats_minimap2_primary.tsv"
#                    INPUT_LABEL="${MAPPER}_${GRAPH}${READS}${PAIRING_TAG}"
#                    ;;
#                esac
#                echo  ${INPUT_LABEL}  ${INPUT_FILE}
#                if [ "${PAIRING}" == "single" ] ; then
#                    grep "${INPUT_LABEL}" ${INPUT_FILE} | grep -v "pe" |  sed "s/${INPUT_LABEL}/${MAPPER}/g" | sed 's/null/0/g' >> roc_stats_${GRAPH}_${READS}_${PAIRING}.tsv
#                else
#                    grep "${INPUT_LABEL}" ${INPUT_FILE} | sed "s/${INPUT_LABEL}/${MAPPER}/g" | sed 's/null/0/g' >> roc_stats_${GRAPH}_${READS}_${PAIRING}.tsv
#                fi
#            done
#            if [ "${PAIRING}" == "single" ] ; then
#                ./plot-roc-single.R roc_stats_${GRAPH}_${READS}_${PAIRING}.tsv plot-roc-${GRAPH}-${READS}-${PAIRING}.svg 
#            else
#                ./plot-roc-paired.R roc_stats_${GRAPH}_${READS}_${PAIRING}.tsv plot-roc-${GRAPH}-${READS}-${PAIRING}.svg
#            fi
#        done
#    done
#done
    


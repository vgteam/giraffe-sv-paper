#!/usr/bin/env bash
# giraffe-yeast-experiment.sh: run the yeast graph analysis for the Giraffe paper
set -ex
set -o pipefail

# Where should intermediate and output files go?
WORKDIR="${HOME}/build/vg/trash/yeast"
# Where should input graphs come from?
IN_ALL_STRAINS="/public/groups/cgl/users/daheller/yeast_graph/graphs/cactus_all/yeast.chop32.vg"
IN_SUBSET="/public/groups/cgl/users/daheller/yeast_graph/graphs/cactus_four/yeast.chop32.vg"
# What training FASTQ should be used for simulating reads?
IN_TRAINING_FASTQ="/public/groups/cgl/users/daheller/yeast_graph/illumina_reads/SRR4074257.fastq.gz"

mkdir -p "${WORKDIR}"

if [ ! -e ${WORKDIR}/yeast_subset.vg ] ; then
    # Import graphs to modern formats
    vg convert -p "${IN_SUBSET}" > ${WORKDIR}/yeast_subset.vg
fi

if [ ! -e ${WORKDIR}/yeast_all.vg ] ; then
    vg convert -p "${IN_ALL_STRAINS}" > ${WORKDIR}/yeast_all.vg
fi

if [ ! -e ${WORKDIR}/yeast_subset.snarls ] ; then
    vg snarls --include-trivial ${WORKDIR}/yeast_subset.vg >${WORKDIR}/yeast_subset.snarls
fi
if [ ! -e ${WORKDIR}/yeast_subset.dist ] ; then
    vg index -s ${WORKDIR}/yeast_subset.snarls -j ${WORKDIR}/yeast_subset.dist
fi
if [ ! -e  ${WORKDIR}/yeast_subset.xg ] ; then
    vg index -x ${WORKDIR}/yeast_subset.xg ${WORKDIR}/yeast_subset.vg
fi
if [ ! -e ${WORKDIR}/yeast_subset.gbwt ] ; then
    vg index -T -G ${WORKDIR}/yeast_subset.gbwt ${WORKDIR}/yeast_subset.xg
fi
if [ ! -e ${WORKDIR}/yeast_subset.gg ] ; then 
    vg gbwt -g ${WORKDIR}/yeast_subset.gg -x ${WORKDIR}/yeast_subset.xg ${WORKDIR}/yeast_subset.gbwt
fi
if [ ! -e ${WORKDIR}/yeast_subset.min ] ; then
    vg minimizer -i ${WORKDIR}/yeast_subset.min -g ${WORKDIR}/yeast_subset.gbwt -G ${WORKDIR}/yeast_subset.gbwt
fi
if [ ! -e ${WORKDIR}/yeast_subset.pruned.vg ] || [ ! -e ${WORKDIR}/yeast_subset.mapping ] ; then
    vg prune -u -p ${WORKDIR}/yeast_subset.vg -m ${WORKDIR}/yeast_subset.mapping >${WORKDIR}/yeast_subset.pruned.vg
fi
if [ ! -e ${WORKDIR}/yeast_subset.gcsa ] || [ ! -e ${WORKDIR}/yeast_subset.gcsa.lcp] ; then
    vg index -g ${WORKDIR}/yeast_subset.gcsa -f ${WORKDIR}/yeast_subset.mapping ${WORKDIR}/yeast_subset.pruned.vg
fi

if [ ! -e ${WORKDIR}/all_contigs.txt ] ; then
    vg paths -L -v ${WORKDIR}/yeast_all.vg > ${WORKDIR}/all_contigs.txt
fi
if [ ! -e ${WORKDIR}/subset_contigs.txt ] ; then
    vg paths -L -v ${WORKDIR}/yeast_subset.vg > ${WORKDIR}/subset_contigs.txt
fi

if [ ! -e ${WORKDIR}/all_strains.txt ] ; then
    cat ${WORKDIR}/all_contigs.txt | cut -f1 -d'.' | sort | uniq > ${WORKDIR}/all_strains.txt
fi
if [ ! -e ${WORKDIR}/subset_strains.txt ] ; then
    cat ${WORKDIR}/subset_contigs.txt | cut -f1 -d'.' | sort | uniq > ${WORKDIR}/subset_strains.txt
fi
if [ ! -e ${WORKDIR}/heldout_strains.txt ] ; then
    comm -1 -3 ${WORKDIR}/subset_strains.txt ${WORKDIR}/all_strains.txt > ${WORKDIR}/heldout_strains.txt
fi

for STRAIN in $(cat ${WORKDIR}/heldout_strains.txt) ; do
    if [ ! -e  "${WORKDIR}/sim-${STRAIN}.gam" ] ; then
        echo "Simulating reads for strain ${STRAIN}"
        vg sim -t 10 -a -x ${WORKDIR}/yeast_all.vg -F "${IN_TRAINING_FASTQ}" -p 570 -v 165 -i 0.00029 -n 500000 $(cat ${WORKDIR}/all_contigs.txt | grep "${STRAIN}" | sed 's/^/--path /g' | tr '\n' ' ') | vg annotate -m -a - -x ${WORKDIR}/yeast_all.vg > "${WORKDIR}/sim-${STRAIN}.gam"
    fi
    if [ ! -e "${WORKDIR}/mapped-giraffe-${STRAIN}.gam" ] ; then
        echo "Mapping strain ${STRAIN} with Giraffe"
        vg giraffe -x ${WORKDIR}/yeast_subset.xg -g ${WORKDIR}/yeast_subset.gg -H ${WORKDIR}/yeast_subset.gbwt -m ${WORKDIR}/yeast_subset.min -d ${WORKDIR}/yeast_subset.dist -i -G "${WORKDIR}/sim-${STRAIN}.gam" | vg annotate -m -a - -x ${WORKDIR}/yeast_subset.xg > "${WORKDIR}/mapped-giraffe-${STRAIN}.gam"
    fi
    if [ ! -e "${WORKDIR}/mapped-giraffe-${STRAIN}.gam" ] ; then
        echo "Mapping strain ${STRAIN} with map"
        vg map -x ${WORKDIR}/yeast_subset.xg -g ${WORKDIR}/yeast_subset.gcsa -i -G "${WORKDIR}/sim-${STRAIN}.gam" | vg annotate -m -a - -x ${WORKDIR}/yeast_subset.xg > "${WORKDIR}/mapped-map-${STRAIN}.gam"
    fi
done

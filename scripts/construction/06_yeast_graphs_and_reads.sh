#!/usr/bin/env bash
# giraffe-yeast-experiment.sh: run the yeast graph analysis for the Giraffe paper
set -ex
set -o pipefail

# Where should intermediate and output files go? No whitespace!
WORKDIR="${HOME}/build/vg/trash/yeast2"
# Where should input graphs come from?
IN_ALL_STRAINS_HAL="/public/groups/cgl/users/daheller/yeast_graph/graphs/cactus_all/cactusoutput.hal"
IN_SUBSET_HAL="/public/groups/cgl/users/daheller/yeast_graph/graphs/cactus_four/cactusoutput.hal"
# What training FASTQ should be used for simulating reads?
IN_TRAINING_FASTQ="/public/groups/cgl/users/daheller/yeast_graph/illumina_reads/SRR4074257.fastq.gz"
# Where's the input reference FASTA for the linear control?
IN_REFERENCE_FASTA="/public/groups/cgl/users/daheller/yeast_graph/assemblies/assemblies_raw/S288C.genome.fa"
# What strain is it?
REFERENCE_STRAIN="S288C"

# Where are we?
SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"

mkdir -p "${WORKDIR}"

# Where should temp files go?
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

barrier() {
    for PID in $(jobs -p) ; do
        wait
        if [ "${?}" -ne "0" ] ; then
            echo "Job failed!"
            kill $(jobs -p)
            exit 1
        fi
    done
}

index_graph() {
    # Build Giraffe and vg map indexes next to the given ${BASENAME}.vg
    BASENAME="${1}"
    shift
    
    if [ ! -e ${BASENAME}.snarls ] ; then
        vg snarls --include-trivial ${BASENAME}.vg >${BASENAME}.snarls &
    fi
    if [ ! -e  ${BASENAME}.xg ] ; then
        vg index -x ${BASENAME}.xg ${BASENAME}.vg &
    fi
    if [ ! -e ${BASENAME}.pruned.vg ] || [ ! -e ${BASENAME}.mapping ] ; then
        vg prune -u -p ${BASENAME}.vg -m ${BASENAME}.mapping >${BASENAME}.pruned.vg &
    fi
    
    barrier
    
    if [ ! -e ${BASENAME}.dist ] ; then
        vg index -s ${BASENAME}.snarls -j ${BASENAME}.dist ${BASENAME}.vg &
    fi
    if [ ! -e ${BASENAME}.gbwt ] ; then
        vg index -T -G ${BASENAME}.gbwt ${BASENAME}.xg &
    fi
    if [ ! -e ${BASENAME}.gcsa ] || [ ! -e ${BASENAME}.gcsa.lcp ] ; then
        vg index -g ${BASENAME}.gcsa -f ${BASENAME}.mapping ${BASENAME}.pruned.vg &
    fi
    
    barrier
    
    if [ ! -e ${BASENAME}.gg ] ; then 
        vg gbwt -g ${BASENAME}.gg -x ${BASENAME}.xg ${BASENAME}.gbwt
    fi
    if [ ! -e ${BASENAME}.min ] ; then
        vg minimizer -i ${BASENAME}.min -g ${BASENAME}.gbwt -G ${BASENAME}.gg
    fi
}

# Import graphs from HAL, excluding wandering components
if [ ! -e ${WORKDIR}/yeast_subset.vg ] ; then
    hal2vg "${IN_SUBSET_HAL}" --inMemory --noAncestors --chop 32 --progress > ${WORKDIR}/yeast_subset.tmp1.vg
    # This doesn't work when piping in for some reason.
    vg ids -s ${WORKDIR}/yeast_subset.tmp1.vg > ${WORKDIR}/yeast_subset.tmp2.vg
    rm ${WORKDIR}/yeast_subset.tmp1.vg
    mv ${WORKDIR}/yeast_subset.tmp2.vg ${WORKDIR}/yeast_subset.vg
fi
if [ ! -e ${WORKDIR}/yeast_all.vg ] ; then
    hal2vg "${IN_ALL_STRAINS_HAL}" --inMemory --noAncestors --chop 32 --progress > ${WORKDIR}/yeast_all.tmp1.vg
    vg ids -s ${WORKDIR}/yeast_all.tmp1.vg > ${WORKDIR}/yeast_all.tmp2.vg
    rm ${WORKDIR}/yeast_all.tmp1.vg
    mv ${WORKDIR}/yeast_all.tmp2.vg ${WORKDIR}/yeast_all.vg
fi

# Find the contig names in the full and some-held-out graphs
if [ ! -e ${WORKDIR}/all_contigs.txt ] ; then
    vg paths -L -v ${WORKDIR}/yeast_all.vg > ${WORKDIR}/all_contigs.txt
fi
if [ ! -e ${WORKDIR}/subset_contigs.txt ] ; then
    vg paths -L -v ${WORKDIR}/yeast_subset.vg > ${WORKDIR}/subset_contigs.txt
fi

# Find the strain names in both graphs
if [ ! -e ${WORKDIR}/all_strains.txt ] ; then
    cat ${WORKDIR}/all_contigs.txt | cut -f1 -d'.' | sort | uniq > ${WORKDIR}/all_strains.txt
fi
if [ ! -e ${WORKDIR}/subset_strains.txt ] ; then
    cat ${WORKDIR}/subset_contigs.txt | cut -f1 -d'.' | sort | uniq > ${WORKDIR}/subset_strains.txt
fi
if [ ! -e ${WORKDIR}/heldout_strains.txt ] ; then
    comm -1 -3 ${WORKDIR}/subset_strains.txt ${WORKDIR}/all_strains.txt > ${WORKDIR}/heldout_strains.txt
fi

if [ ! -e "${WORKDIR}/yeast_reference.fa" ] ; then
    # Get a reference strain FASTA with Cactus-style names
    cat "${IN_REFERENCE_FASTA}" | sed "s/>/>${REFERENCE_STRAIN}./g" > ${WORKDIR}/yeast_reference.fa
fi
if [ ! -e "${WORKDIR}/yeast_reference.vg" ] ; then
    # Make linear graph for the reference strain
    vg construct -p -r ${WORKDIR}/yeast_reference.fa  | vg convert -p - >${WORKDIR}/yeast_reference.vg
fi

# Index the linerar reference graph we're mapping to
index_graph ${WORKDIR}/yeast_reference

# Index the pangenome graph we're mapping to
index_graph ${WORKDIR}/yeast_subset

# Index the source graph to map back to for positive control
index_graph ${WORKDIR}/yeast_all

# Index for BWA
if [ ! -e ${WORKDIR}/yeast_reference.fa.amb  ] ; then
    bwa index ${WORKDIR}/yeast_reference.fa
fi

for FILENAME in $(cd ${WORKDIR} && ls yeast_reference.*) ; do
    if [[ "${FILENAME#*.}" == ".mapping" || "${FILENAME#*.}" == ".pruned.vg" ]] ; then
        continue
    fi
    aws s3 cp "${WORKDIR}/${FILENAME}" "s3://vg-k8s/profiling/graphs/v2/generic/primary/${REFERENCE_STRAIN}/primary${REFERENCE_STRAIN}.${FILENAME#*.}"
done
for FILENAME in $(cd ${WORKDIR} && ls yeast_all.*) ; do
    if [[ "${FILENAME#*.}" == ".mapping" || "${FILENAME#*.}" == ".pruned.vg" ]] ; then
        continue
    fi
    aws s3 cp "${WORKDIR}/${FILENAME}" "s3://vg-k8s/profiling/graphs/v2/generic/cactus/yeast_all/${FILENAME}"
done
for FILENAME in $(cd ${WORKDIR} && ls yeast_subset.*) ; do
    if [[ "${FILENAME#*.}" == ".mapping" || "${FILENAME#*.}" == ".pruned.vg" ]] ; then
        continue
    fi
    aws s3 cp "${WORKDIR}/${FILENAME}" "s3://vg-k8s/profiling/graphs/v2/generic/cactus/yeast_subset/${FILENAME}"
done

for STRAIN in $(cat ${WORKDIR}/heldout_strains.txt) ; do
    # Prepare reads
    if [ ! -e  "${WORKDIR}/sim-${STRAIN}.gam" ] ; then
        echo "Simulating reads for strain ${STRAIN}"
        vg sim -t 10 -a -x ${WORKDIR}/yeast_all.vg -F "${IN_TRAINING_FASTQ}" -p 570 -v 165 -i 0.00029 -n 500000 $(cat ${WORKDIR}/all_contigs.txt | grep "${STRAIN}" | sed 's/^/--path /g' | tr '\n' ' ') | vg annotate -m -a - -x ${WORKDIR}/yeast_all.vg > "${WORKDIR}/sim-${STRAIN}.gam" &
    fi
done

barrier

for STRAIN in $(cat ${WORKDIR}/heldout_strains.txt) ; do
    if [ ! -e  "${WORKDIR}/sim-${STRAIN}.fq.gz" ] ; then
        # Keep vg pairing style
        (vg view -i -X "${WORKDIR}/sim-${STRAIN}.gam" > "${WORKDIR}/sim-${STRAIN}.fq" && gzip "${WORKDIR}/sim-${STRAIN}.fq") &
    fi
done

barrier

for FILENAME in $(cd ${WORKDIR} && ls sim-*) ; do
    aws s3 cp "${WORKDIR}/${FILENAME}" "s3://vg-k8s/profiling/reads/sim/yeast/${FILENAME}"
done

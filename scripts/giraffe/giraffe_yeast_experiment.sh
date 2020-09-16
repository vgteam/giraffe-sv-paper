#!/usr/bin/env bash
# giraffe-yeast-experiment.sh: run the yeast graph analysis for the Giraffe paper
set -ex
set -o pipefail

# Where should intermediate and output files go? No whitespace!
WORKDIR="${HOME}/build/vg/trash/yeast"
# Where should input graphs come from?
IN_ALL_STRAINS="/public/groups/cgl/users/daheller/yeast_graph/graphs/cactus_all/yeast.chop32.vg"
IN_SUBSET="/public/groups/cgl/users/daheller/yeast_graph/graphs/cactus_four/yeast.chop32.vg"
# What training FASTQ should be used for simulating reads?
IN_TRAINING_FASTQ="/public/groups/cgl/users/daheller/yeast_graph/illumina_reads/SRR4074257.fastq.gz"
# Where's the input reference FASTA for the linear control?
IN_REFERENCE_FASTA="/public/groups/cgl/users/daheller/yeast_graph/assemblies/assemblies_raw/SK1.genome.fa"
# What strain is it?
REFERENCE_STRAIN="SK1"

# Where are we?
SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"

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

mkdir -p "${WORKDIR}"

# Import graphs to modern formats
if [ ! -e ${WORKDIR}/yeast_subset.vg ] ; then
    vg convert -p "${IN_SUBSET}" > ${WORKDIR}/yeast_subset.vg
fi
if [ ! -e ${WORKDIR}/yeast_all.vg ] ; then
    vg convert -p "${IN_ALL_STRAINS}" > ${WORKDIR}/yeast_all.vg
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

# Index for BWA
if [ ! -e ${WORKDIR}/yeast_reference.fa.amb  ] ; then
    bwa index ${WORKDIR}/yeast_reference.fa
fi

for STRAIN in $(cat ${WORKDIR}/heldout_strains.txt) ; do
    # Prepare reads
    if [ ! -e  "${WORKDIR}/sim-${STRAIN}.gam" ] ; then
        echo "Simulating reads for strain ${STRAIN}"
        vg sim -t 10 -a -x ${WORKDIR}/yeast_all.vg -F "${IN_TRAINING_FASTQ}" -p 570 -v 165 -i 0.00029 -n 500000 $(cat ${WORKDIR}/all_contigs.txt | grep "${STRAIN}" | sed 's/^/--path /g' | tr '\n' ' ') | vg annotate -m -a - -x ${WORKDIR}/yeast_all.vg > "${WORKDIR}/sim-${STRAIN}.gam"
    fi
    if [ ! -e  "${WORKDIR}/sim-${STRAIN}.fq" ] ; then
        # Make sure to drop pair number from read name because BWA needs exact match (fragment names)
        vg view -i -X "${WORKDIR}/sim-${STRAIN}.gam" | sed 's/_1$//g' | sed 's/_2$//g' > "${WORKDIR}/sim-${STRAIN}.fq"
    fi
    
    # Do mapping
    if [ ! -e "${WORKDIR}/mapped-subset-giraffe-${STRAIN}.gam" ] ; then
        echo "Mapping strain ${STRAIN} to graph with Giraffe"
        vg giraffe -x ${WORKDIR}/yeast_subset.xg -g ${WORKDIR}/yeast_subset.gg -H ${WORKDIR}/yeast_subset.gbwt -m ${WORKDIR}/yeast_subset.min -d ${WORKDIR}/yeast_subset.dist -i -G "${WORKDIR}/sim-${STRAIN}.gam" | vg annotate -m -a - -x ${WORKDIR}/yeast_subset.xg > "${WORKDIR}/mapped-subset-giraffe-${STRAIN}.gam" &
    fi
    if [ ! -e "${WORKDIR}/mapped-reference-giraffe-${STRAIN}.gam" ] ; then
        echo "Mapping strain ${STRAIN} to linear reference with Giraffe"
        vg giraffe -x ${WORKDIR}/yeast_reference.xg -g ${WORKDIR}/yeast_reference.gg -H ${WORKDIR}/yeast_reference.gbwt -m ${WORKDIR}/yeast_reference.min -d ${WORKDIR}/yeast_reference.dist -i -G "${WORKDIR}/sim-${STRAIN}.gam" | vg annotate -m -a - -x ${WORKDIR}/yeast_reference.xg > "${WORKDIR}/mapped-reference-giraffe-${STRAIN}.gam" &
    fi
    if [ ! -e "${WORKDIR}/mapped-subset-map-${STRAIN}.gam" ] ; then
        echo "Mapping strain ${STRAIN} to graph with map"
        vg map -x ${WORKDIR}/yeast_subset.xg -g ${WORKDIR}/yeast_subset.gcsa -i -G "${WORKDIR}/sim-${STRAIN}.gam" | vg annotate -m -a - -x ${WORKDIR}/yeast_subset.xg > "${WORKDIR}/mapped-subset-map-${STRAIN}.gam" &
    fi
    if [ ! -e "${WORKDIR}/mapped-reference-map-${STRAIN}.gam" ] ; then
        echo "Mapping strain ${STRAIN} to linear reference with map"
        vg map -x ${WORKDIR}/yeast_reference.xg -g ${WORKDIR}/yeast_reference.gcsa -i -G "${WORKDIR}/sim-${STRAIN}.gam" | vg annotate -m -a - -x ${WORKDIR}/yeast_reference.xg > "${WORKDIR}/mapped-reference-map-${STRAIN}.gam" &
    fi
    if [ ! -e "${WORKDIR}/mapped-reference-bwa-${STRAIN}.bam" ] ; then
        # Drop supplementary alignments and convert to BAM
        bwa mem ${WORKDIR}/yeast_reference.fa -p "${WORKDIR}/sim-${STRAIN}.fq" | samtools view -b /dev/stdin > "${WORKDIR}/mapped-reference-bwa-${STRAIN}.bam" &
    fi
    
    barrier
    
    if [ ! -e "${WORKDIR}/mapped-reference-bwa-${STRAIN}.gam" ] ; then
        # Inject and annotate and also hack read names while we go through JSON.
        # TODO: make inject do this??? Or do it at the SAM stage?
        samtools view -F 2048 -b "${WORKDIR}/mapped-reference-bwa-${STRAIN}.bam" | vg inject -x ${WORKDIR}/yeast_reference.xg - | vg view -aj - | sed 's/\/1/_1/g' | sed 's/\/2/_2/g' | vg view -JGa - | vg annotate -m -x ${WORKDIR}/yeast_reference.xg -a - > "${WORKDIR}/mapped-reference-bwa-${STRAIN}.gam" &
    fi
    if [ ! -e "${WORKDIR}/mapped-reference-bwa-${STRAIN}.sup.gam" ] ; then
        # Also bring along supplementary alignments
        samtools view -f 2048 -b "${WORKDIR}/mapped-reference-bwa-${STRAIN}.bam" | vg inject -x ${WORKDIR}/yeast_reference.xg - | vg view -aj - | sed 's/\/1/_1/g' | sed 's/\/2/_2/g' | vg view -JGa - | vg annotate -m -x ${WORKDIR}/yeast_reference.xg -a - > "${WORKDIR}/mapped-reference-bwa-${STRAIN}.sup.gam" &
    fi
    
    barrier
    
    # Do all the comparisons
    if [ ! -e "${WORKDIR}/mapped-subset-giraffe-${STRAIN}.compared.gam" ] ; then
        vg gamcompare -s -r 100 "${WORKDIR}/mapped-subset-giraffe-${STRAIN}.gam" "${WORKDIR}/sim-${STRAIN}.gam" > "${WORKDIR}/mapped-subset-giraffe-${STRAIN}.compared.gam" &
    fi
    if [ ! -e "${WORKDIR}/mapped-reference-giraffe-${STRAIN}.compared.gam" ] ; then
        vg gamcompare -s -r 100 "${WORKDIR}/mapped-reference-giraffe-${STRAIN}.gam" "${WORKDIR}/sim-${STRAIN}.gam" > "${WORKDIR}/mapped-reference-giraffe-${STRAIN}.compared.gam" &
    fi
    if [ ! -e "${WORKDIR}/mapped-subset-map-${STRAIN}.compared.gam" ] ; then
        vg gamcompare -s -r 100 "${WORKDIR}/mapped-subset-map-${STRAIN}.gam" "${WORKDIR}/sim-${STRAIN}.gam" > "${WORKDIR}/mapped-subset-map-${STRAIN}.compared.gam" &
    fi
    if [ ! -e "${WORKDIR}/mapped-reference-map-${STRAIN}.compared.gam" ] ; then
        vg gamcompare -s -r 100 "${WORKDIR}/mapped-reference-map-${STRAIN}.gam" "${WORKDIR}/sim-${STRAIN}.gam" > "${WORKDIR}/mapped-reference-map-${STRAIN}.compared.gam" &
    fi
    if [ ! -e "${WORKDIR}/mapped-reference-bwa-${STRAIN}.compared.gam" ] ; then
        # We need to do two gamcompares and combine for linear mappers.
        python3 "${SCRIPT_DIR}/../linear_mappers/combine_reads.py" <(vg gamcompare -s -r 100 ${WORKDIR}/mapped-reference-bwa-${STRAIN}.gam "${WORKDIR}/sim-${STRAIN}.gam" | vg view -aj -) <(vg gamcompare -s -r 100 ${WORKDIR}/mapped-reference-bwa-${STRAIN}.sup.gam "${WORKDIR}/sim-${STRAIN}.gam" | vg view -aj -) /dev/stdout | vg view -JGa - > "${WORKDIR}/mapped-reference-bwa-${STRAIN}.compared.gam" &
    fi
    
    barrier
    
done

# TODO: make ROCs and reports

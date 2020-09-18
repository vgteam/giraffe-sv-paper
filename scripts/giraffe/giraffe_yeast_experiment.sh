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
IN_REFERENCE_FASTA="/public/groups/cgl/users/daheller/yeast_graph/assemblies/assemblies_raw/SK1.genome.fa"
# What strain is it?
REFERENCE_STRAIN="SK1"

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
    
    # Map reads
    for GRAPH in reference subset all ; do
        for MAPPER in map giraffe bwa ; do
            if [ ! -e "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.gam" ] ; then
                if [ "${MAPPER}" == giraffe ] ; then
                     vg giraffe -x ${WORKDIR}/yeast_${GRAPH}.xg -g ${WORKDIR}/yeast_${GRAPH}.gg -H ${WORKDIR}/yeast_${GRAPH}.gbwt -m ${WORKDIR}/yeast_${GRAPH}.min -d ${WORKDIR}/yeast_${GRAPH}.dist -i -G "${WORKDIR}/sim-${STRAIN}.gam" | vg annotate -m -a - -x ${WORKDIR}/yeast_${GRAPH}.xg > "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.gam" &
                elif [ "${MAPPER}" == map ] ; then
                    vg map -x ${WORKDIR}/yeast_${GRAPH}.xg -g ${WORKDIR}/yeast_${GRAPH}.gcsa -i -G "${WORKDIR}/sim-${STRAIN}.gam" | vg annotate -m -a - -x ${WORKDIR}/yeast_${GRAPH}.xg > "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.gam" &
                elif [ "${MAPPER}" == bwa ] && [ "${GRAPH}" == reference ] ; then
                    # Run the BWA mapping and injecting all as one backgorund job.
                    # Inject and annotate and also hack read names while we go through JSON.
                    # TODO: make inject do this??? Or do it at the SAM stage?
                    ((bwa mem ${WORKDIR}/yeast_${GRAPH}.fa -p "${WORKDIR}/sim-${STRAIN}.fq" | samtools view -b /dev/stdin > "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.bam") && \
                    (samtools view -f 2048 -b "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.bam" | vg inject -x ${WORKDIR}/yeast_${GRAPH}.xg - | vg view -aj - | sed 's/\/1/_1/g' | sed 's/\/2/_2/g' | vg view -JGa - | vg annotate -m -x ${WORKDIR}/yeast_${GRAPH}.xg -a - > "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.sup.gam") && \
                    (samtools view -F 2048 -b "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.bam" | vg inject -x ${WORKDIR}/yeast_${GRAPH}.xg - | vg view -aj - | sed 's/\/1/_1/g' | sed 's/\/2/_2/g' | vg view -JGa - | vg annotate -m -x ${WORKDIR}/yeast_${GRAPH}.xg -a - > "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.gam")) &
                fi
            fi
        done
    done
    
    barrier
    
    # Compare reads
    for GRAPH in reference subset all ; do
        for MAPPER in map giraffe bwa ; do
            if [ ! -e "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.compared.gam" ] ; then
                if [ "${MAPPER}" == bwa ] ; then
                    if [ "${GRAPH}" == reference ] ; then
                        # We need to do two gamcompares and combine for linear mappers.
                        python3 "${SCRIPT_DIR}/../linear_mappers/combine_reads.py" <(vg gamcompare -s -r 100 ${WORKDIR}/mapped-reference-bwa-${STRAIN}.gam "${WORKDIR}/sim-${STRAIN}.gam" | vg view -aj -) <(vg gamcompare -s -r 100 ${WORKDIR}/mapped-reference-bwa-${STRAIN}.sup.gam "${WORKDIR}/sim-${STRAIN}.gam" | vg view -aj -) /dev/stdout | vg view -JGa - > "${WORKDIR}/mapped-reference-bwa-${STRAIN}.compared.gam" &
                    fi
                else
                    # Do it the normal way
                    vg gamcompare -s -r 100 "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.gam" "${WORKDIR}/sim-${STRAIN}.gam" > "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.compared.gam" &
                fi
            fi
        done
    done
    
    barrier
    
    # Dump to JSON
    for GRAPH in reference subset all ; do
        for MAPPER in map giraffe bwa ; do
            if [ "${MAPPER}" == bwa ] && [ "${GRAPH}" != reference ] ; then
                # BWA only works on reference
                continue
            fi
            if [ ! -e "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.compared.json" ] ; then
                vg view -aj "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.compared.gam" > "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.compared.json" &
            fi
        done
    done
    
    barrier
    
    for GRAPH in reference subset all ; do
        for MAPPER in map giraffe bwa ; do
            if [ "${MAPPER}" == bwa ] && [ "${GRAPH}" != reference ] ; then
                # BWA only works on reference
                continue
            fi
            if [ ! -e "${WORKDIR}/report-${GRAPH}-${MAPPER}-${STRAIN}.tsv" ] ; then
                MAPPED_COUNT="$(grep path "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.compared.json" | wc -l)"
                CORRECT_COUNT="$(grep correctly_mapped "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.compared.json" | wc -l)"
                MAPQ60_TOTAL="$(grep mapping_quality\":\ 60 "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.compared.json" | wc -l)"
                MAPQ60_WRONG="$(grep -v correctly_mapped "${WORKDIR}/mapped-${GRAPH}-${MAPPER}-${STRAIN}.compared.json" | grep mapping_quality\":\ 60 | wc -l)"
                printf "${GRAPH}\t${MAPPER}\t${STRAIN}\t${MAPPED_COUNT}\t${CORRECT_COUNT}\t${MAPQ60_TOTAL}\t${MAPQ60_WRONG}\n" > "${WORKDIR}/report-${GRAPH}-${MAPPER}-${STRAIN}.tsv"
            fi
        done
    done
done

printf "#GRAPH\tMAPPER\tSTRAIN\tMAPPED_COUNT\tCORRECT_COUNT\tMAPQ60_TOTAL\tMAPQ60_WRONG\n" > ${WORKDIR}/report.tsv
cat ${WORKDIR}/report-*.tsv >> ${WORKDIR}/report.tsv

# TODO: make ROCs

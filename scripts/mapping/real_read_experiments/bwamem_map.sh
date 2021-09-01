#!/usr/bin/env bash
# bwamem_map.sh: map HG002, HG003, 150bp and 250bp paired reads with BWAMEM to the HS38d1 reference

set -ex
set -o pipefail

function download() {
    if [ ! -e "${2}" ] ; then
        aws s3 cp --no-progress "${1}" "${2}"
    fi
}

function wget_download() {
    if [ ! -e "${2}" ] ; then
        wget "${1}" -O "${2}"
    fi
}

function copy() {
    if [ ! -e "${2}" ] ; then
        cp "${1}" "${2}"
    fi
}

SAMPLE_NAMES=("HG002" "HG003")
READSETS=("150bp" "250bp")
WORK_DIR=${HOME}/run_bwamem_mapping
REF_FASTA="GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.no_segdup.fna"

# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORK_DIR
download s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.renamed.fna.gz "${WORK_DIR}/${REF_FASTA}.gz"
gzip -d "${WORK_DIR}/${REF_FASTA}.gz"

for SAMPLE_NAME in "${SAMPLE_NAMES[@]}" ; do
  for READSET in "${READSETS[@]}" ; do
    if [[ ${SAMPLE_NAME} == *"HG002"* ]]; then
      if [[ ${READSET} == *"150bp"* ]]; then
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R1.fastq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R1.fastq.gz"
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R2.fastq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R2.fastq.gz"
      else
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002_read_pair_1.fq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R1.fastq.gz"
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002_read_pair_2.fq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R2.fastq.gz"
      fi
    elif [[ ${SAMPLE_NAME} == *"HG003"* ]]; then
      if [[ ${READSET} == *"150bp"* ]]; then  
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG003.novaseq.pcr-free.35x.R1.fastq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R1.fastq.gz"
        wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG003.novaseq.pcr-free.35x.R2.fastq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R2.fastq.gz"
      else
        continue
      fi
    fi
    
    # run bwa mem mapper
    cd $WORK_DIR
    docker run \
    -e REF_FASTA=${REF_FASTA} \
    -e READ1="${SAMPLE_NAME}.${READSET}.R1.fastq.gz" \
    -e READ2="${SAMPLE_NAME}.${READSET}.R2.fastq.gz" \
    -v ${PWD}:${HOME} -w ${HOME} biocontainers/bwa:v0.7.17_cv1 \
    bwa mem \
    -t 16 \
    ${REF_FASTA} \
    ${READ1} \
    ${READ2} \
    > bwamem_${SAMPLE_NAME}_${READSET}.sam
  done
done



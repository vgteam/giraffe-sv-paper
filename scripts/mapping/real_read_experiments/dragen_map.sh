#!/usr/bin/env bash
# dragen_map.sh: map HG002, HG003, 150bp and 250bp paired reads with Illumina's DRAGEN module to the HS38d1 reference

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
WORK_DIR=${HOME}/run_dragen_mapping
REF_FASTA="GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.no_segdup.fna"

# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORK_DIR
download s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.renamed.fna.gz "${WORK_DIR}/${REF_FASTA}.gz"
gzip -d "${WORK_DIR}/${REF_FASTA}.gz"

# Generate the DRAGEN reference index
mkdir -p ${WORK_DIR}/dragen_index ${WORK_DIR}/tmp && cd ${WORK_DIR}
dragen --build-hash-table true \
  --output-directory ${WORK_DIR}/dragen_index \
  --ht-reference ${WORK_DIR}/${REF_FASTA}

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
    mkdir -p ${WORK_DIR}/dragen_output_${SAMPLE_NAME}_${READSET}
    dragen -f \
      -r ${WORK_DIR}/dragen_index \
      -1 "${SAMPLE_NAME}.${READSET}.R1.fastq.gz" \
      -2 "${SAMPLE_NAME}.${READSET}.R2.fastq.gz" \
      --RGID 1 \
      --RGSM ${SAMPLE_NAME} \
      --verbose --bin_memory=50000000000 --enable-map-align true --enable-variant-caller true \
      --pair-by-name=true \
      --enable-map-align-output=true \
      --intermediate-results-dir ${WORK_DIR}/tmp \
      --output-directory ${WORK_DIR}/dragen_output_${SAMPLE_NAME}_${READSET} \
      --output-file-prefix dragen_output_${SAMPLE_NAME}_${READSET} 2> dragen_output_${SAMPLE_NAME}_${READSET}.stderr
  done
done





#!/usr/bin/env bash
# giraffe_map.sh: map HG002, HG003, 150bp and 250bp paired reads with VG GIRAFFE, VG GIRAFFE FAST, and VG GIRAFFE PRIMARY to the HS38d1-based graph references

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
MAPPER_TYPES=("DEFAULT" "FAST" "PRIMARY")
WORK_DIR=${HOME}/run_giraffe_mapping
REF_FASTA="GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.no_segdup.fna"

# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORK_DIR
download s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.xg "${WORK_DIR}/1000GPlons_hs38d1_filter.xg"
download s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.sampled.64.gbwt "${WORK_DIR}/1000GPlons_hs38d1_filter.sampled.64.gbwt"
download s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.sampled.64.min "${WORK_DIR}/1000GPlons_hs38d1_filter.sampled.64.min"
download s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.sampled.64.gg "${WORK_DIR}/1000GPlons_hs38d1_filter.sampled.64.gg"
download s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.dist "${WORK_DIR}/1000GPlons_hs38d1_filter.dist"
download s3://vg-k8s/profiling/graphs/v2-2/generic/primary/hs38d1/primaryhs38d1.xg "${WORK_DIR}/primaryhs38d1.xg"
download s3://vg-k8s/profiling/graphs/v2-2/generic/primary/hs38d1/primaryhs38d1.dist "${WORK_DIR}/primaryhs38d1.dist"
download s3://vg-k8s/profiling/graphs/v2-2/generic/primary/hs38d1/primaryhs38d1.cover.gbwt "${WORK_DIR}/primaryhs38d1.cover.gbwt"
download s3://vg-k8s/profiling/graphs/v2-2/generic/primary/hs38d1/primaryhs38d1.cover.gg "${WORK_DIR}/primaryhs38d1.cover.gg"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/giraffe_manuscript_data/genome_references/linear_references/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.no_segdup.dict "${WORK_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.no_segdup.dict"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/giraffe_manuscript_data/genome_references/linear_references/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.dict "${WORK_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.dict"

for SAMPLE_NAME in "${SAMPLE_NAMES[@]}" ; do
  for READSET in "${READSETS[@]}" ; do
    for MAPPER_TYPE in "${MAPPER_TYPES[@]}" ; do
      SEQ_DICT="GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.no_segdup.dict"
      XG_IDX="-x 1000GPlons_hs38d1_filter.xg"
      GBWT_IDX="-H 1000GPlons_hs38d1_filter.sampled.64.gbwt"
      MIN_IDX="-m 1000GPlons_hs38d1_filter.sampled.64.min"
      GG_IDX="-g 1000GPlons_hs38d1_filter.sampled.64.gg"
      DIST_IDX="-d 1000GPlons_hs38d1_filter.dist"
      FAST_PARAM=""
      if [[ ${SAMPLE_NAME} == *"HG002"* ]]; then
        if [[ ${READSET} == *"150bp"* ]]; then
          wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R1.fastq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R1.fastq.gz"
          wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R2.fastq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R2.fastq.gz"
        else
          wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002_read_pair_1.fq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R1.fastq.gz"
          wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002_read_pair_2.fq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R2.fastq.gz"
        fi
      elif [[ ${SAMPLE_NAME} == *"HG003"* ]]; then
        if [[ ! ${MAPPER_TYPE} == *"DEFAULT"* ]]; then
          continue
        else
          if [[ ${READSET} == *"150bp"* ]]; then
            wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG003.novaseq.pcr-free.35x.R1.fastq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R1.fastq.gz"
            wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG003.novaseq.pcr-free.35x.R2.fastq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R2.fastq.gz"
          else
            continue
          fi
        fi
      fi
      if [[ ${MAPPER_TYPE} == *"PRIMARY"* ]]; then
        SEQ_DICT="GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.dict"
        XG_IDX="-x primaryhs38d1.xg"
        GBWT_IDX="-H primaryhs38d1.cover.gbwt"
        MIN_IDX=""
        GG_IDX="-g primaryhs38d1.cover.gg"
        DIST_IDX="-d primaryhs38d1.dist"
      elif [[ ${MAPPER_TYPE} == *"FAST"* ]]; then
        FAST_PARAM="-b fast"
      fi
    
    # run giraffe mapper
    cd $WORK_DIR
    docker run \
    -e SEQ_DICT=${SEQ_DICT} \
    -e XG_IDX=${XG_IDX} \
    -e GBWT_IDX=${GBWT_IDX} \
    -e MIN_IDX=${MIN_IDX} \
    -e GG_IDX=${GG_IDX} \
    -e DIST_IDX=${DIST_IDX} \
    -e READ1="${SAMPLE_NAME}.${READSET}.R1.fastq.gz" \
    -e READ2="${SAMPLE_NAME}.${READSET}.R2.fastq.gz" \
    -v ${PWD}:${HOME} -w ${HOME} quay.io/vgteam/vg:v1.31.0 \
    vg giraffe \
    -t 16 \
    ${XG_IDX} \
    ${GBWT_IDX} \
    ${MIN_IDX} \
    ${GG_IDX} \
    ${DIST_IDX} \
    -f ${READ1} \
    -f ${READ2} \
    ${FAST_PARAM} \
    --output-format BAM \
    --ref-paths ${SEQ_DICT} > giraffe_${MAPPER_TYPE}_${SAMPLE_NAME}_${READSET}.bam
    done
  done
done


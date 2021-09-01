#!/usr/bin/env bash
# If rerunning, the Toil 6 and vg 1.33 releases are probably the right choices.
# Run in screen on a Kubernetes cluster: TOIL_APPLIANCE_SELF=quay.io/ucsc_cgl/toil:5.4.0a1-1530a9190357fc058333f3e929049ef9593a6784-py3.7 toil launch-cluster --provisioner aws -T kubernetes -z us-west-2a adamnovak-toil-vg --leaderNodeType t3a.medium --nodeTypes=t3a.medium,r5ad.24xlarge,r5d.24xlarge/r5ad.24xlarge:2.50,i3.8xlarge:1.50 --workers 1-4,0-1,0-8,0-6 --keyPairName anovak@soe.ucsc.edu
# When rerunning, use --defaultPreemptable and/or throw max 5 of the on-demand workers at it

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

SAMPLE_NAME="HG002"
READSETS=("150bp" "250bp")
WORK_DIR=${HOME}/run_vg_map_mapping

# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORK_DIR
download s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.xg "${WORK_DIR}/1000GPlons_hs38d1_filter.xg"
download s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.gcsa "${WORK_DIR}/1000GPlons_hs38d1_filter.gcsa"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/giraffe_manuscript_data/genome_references/graph_references/path.list "${WORK_DIR}/path.list"

for READSET in "${READSETS[@]}" ; do
  if [[ ${READSET} == *"150bp"* ]]; then
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R1.fastq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R1.fastq.gz"
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R2.fastq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R2.fastq.gz"
  else
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002_read_pair_1.fq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R1.fastq.gz"
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002_read_pair_2.fq.gz "${WORK_DIR}/${SAMPLE_NAME}.${READSET}.R2.fastq.gz"
  fi
cd $WORK_DIR
virtualenv --system-site-packages --python python3 venv
. venv/bin/activate
pip3 install --upgrade git+https://github.com/vgteam/toil-vg.git@b319a1b22df6dac585b7f95bc1a603577452d443#egg=toil-vg
READ1="${SAMPLE_NAME}.${READSET}.R1.fastq.gz" \
READ2="${SAMPLE_NAME}.${READSET}.R2.fastq.gz" \
JOBSTORE="${WORK_DIR}/vg_map_${SAMPLE_NAME}_${READSET}_jobstore"
OUTSTORE="${WORK_DIR}/vg_map_${SAMPLE_NAME}_${READSET}_outstore"
LOGFILE="${WORK_DIR}/vg_map_${SAMPLE_NAME}_${READSET}.log"
TMP_DIR="${WORK_DIR}/vg_map_${SAMPLE_NAME}_${READSET}.tmp"
rm -fr ${LOGFILE} ${TMP_DIR} ${OUTSTORE}
mkdir -p ${TMP_DIR} ${OUTSTORE}
toil clean ${JOBSTORE}
toil-vg map \
${JOBSTORE} \
${SAMPLE_NAME} \
${OUTSTORE} \
--batchSystem singleMachine \
--container Docker \
--logFile ${LOGFILE} \
--workDir ${TMP_DIR} \
--cleanWorkDir onSuccess \
--whole_genome_config \
--vg_docker 'quay.io/vgteam/vg:v1.31.0' \
--bam_output \
--id_ranges path.list \
--fastq ${READ1} ${READ2} \
--xg_index 1000GPlons_hs38d1_filter.xg \
--gcsa_index 1000GPlons_hs38d1_filter.gcsa \
--statePollingWait 120 \
--rescueJobsFrequency 120
done







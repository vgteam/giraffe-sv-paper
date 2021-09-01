#!/usr/bin/env bash
# real_read_deepvariant_call.sh: indel realign input BAMs and run the DeepVariant Genotyper

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

WORK_DIR=${HOME}/run_deepvariant_genotyping
SAMPLE_NAME="HG002"
BAM_FILE="${1}"
INPUT_BAM=$(basename ${BAM_FILE})
REF_FASTA="GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.no_segdup.fna.gz"
SEQ_DICT="GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.no_segdup.dict"


# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORK_DIR
copy ${BAM_FILE} "${WORK_DIR}/${INPUT_BAM}"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/giraffe_manuscript_data/genome_references/linear_references/${SEQ_DICT} "${WORK_DIR}/${SEQ_DICT}"
download s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.renamed.fna.gz "${WORK_DIR}/${REF_FASTA}"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/benchmark_data/hg002_cohort/${SAMPLE_NAME}/${SAMPLE_NAME}_GRCh38_1_22_v4.2.1_benchmark.xian_adam.vcf.gz "${WORK_DIR}/${SAMPLE_NAME}_GRCh38_1_22_v4.2.1_benchmark.xian_adam.vcf.gz"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/benchmark_data/hg002_cohort/${SAMPLE_NAME}/${SAMPLE_NAME}_GRCh38_1_22_v4.2.1_benchmark.xian_adam.vcf.gz.tbi "${WORK_DIR}/${SAMPLE_NAME}_GRCh38_1_22_v4.2.1_benchmark.xian_adam.vcf.gz.tbi"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/benchmark_data/hg002_cohort/${SAMPLE_NAME}/${SAMPLE_NAME}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.xian_adam.bed "${WORK_DIR}/${SAMPLE_NAME}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.xian_adam.bed"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/giraffe_manuscript_data/deepvariant_models/model.ckpt-25600.meta "${WORK_DIR}/model.ckpt-25600.meta"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/giraffe_manuscript_data/deepvariant_models/model.ckpt-25600.data-00000-of-00001 "${WORK_DIR}/model.ckpt-25600.data-00000-of-00001"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/giraffe_manuscript_data/deepvariant_models/model.ckpt-25600.index "${WORK_DIR}/model.ckpt-25600.index"

# Sort and reorder input BAM
cd $WORK_DIR
docker run \
-e INPUT_BAM=${INPUT_BAM} \
-v ${PWD}:${HOME} -w ${HOME} quay.io/ucsc_cgl/samtools:latest \
samtools sort -@ 32 ${INPUT_BAM} -O BAM > positionsorted.${INPUT_BAM} && rm ${INPUT_BAM}

docker run \
-e INPUT_BAM=${INPUT_BAM} \
-e SEQ_DICT=${SEQ_DICT} \
-v ${PWD}:${HOME} -w ${HOME} broadinstitute/picard:2.21.9 \
  java -Xmx20g -XX:ParallelGCThreads=16 -jar /usr/picard/picard.jar \
  ReorderSam \
  VALIDATION_STRINGENCY=SILENT \
  INPUT=positionsorted.${INPUT_BAM} \
  OUTPUT=reordered.positionsorted.${INPUT_BAM} \
  SEQUENCE_DICTIONARY=${SEQ_DICT} && rm positionsorted.${INPUT_BAM}

# Indel realign input BAM
cd $WORK_DIR
docker run \
-e INPUT_BAM=${INPUT_BAM} \
-e SAMPLE_NAME=${SAMPLE_NAME} \
-v ${PWD}:${HOME} -w ${HOME} quay.io/ucsc_cgl/samtools:latest \
samtools addreplacerg \
-@ 32 -O BAM -r ID:1 -r LB:lib1 -r SM:${SAMPLE_NAME} -r PL:illumina -r PU:unit1 \
reordered.positionsorted.${INPUT_BAM} > gatk_ready.reordered.positionsorted.${INPUT_BAM}

docker run \
-e INPUT_BAM=${INPUT_BAM} \
-v ${PWD}:${HOME} -w ${HOME} quay.io/ucsc_cgl/samtools:latest \
samtools index -@ 32 gatk_ready.reordered.positionsorted.${INPUT_BAM}

docker run \
-e SAMPLE_NAME=${SAMPLE_NAME} \
-e REF_FASTA=${REF_FASTA} \
-e INPUT_BAM=${INPUT_BAM} \
-v ${PWD}:${HOME} -w ${HOME} broadinstitute/gatk3:3.8-1 \
  java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
  -R ${REF_FASTA} \
  -I gatk_ready.reordered.positionsorted.${INPUT_BAM} -o ${SAMPLE_NAME}.intervals

awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' ${SAMPLE_NAME}.intervals > ${SAMPLE_NAME}.intervals.bed
docker run \
-e SAMPLE_NAME=${SAMPLE_NAME} \
-e INPUT_BAM=${INPUT_BAM} \
-e REF_FASTA=${REF_FASTA} \
-v ${PWD}:${HOME} -w ${HOME} quay.io/biocontainers/abra2:2.24--h7d875b9_0 \
  --targets ${SAMPLE_NAME}.intervals.bed \
  --in gatk_ready.reordered.positionsorted.${INPUT_BAM} \
  --out indel_realigned.${INPUT_BAM} \
  --ref ${REF_FASTA} \
  --threads 16

# Run DeepVariant genotyper on input BAM
cd ${WORK_DIR}
docker run \
-e REF_FASTA=${REF_FASTA} \
-e INPUT_BAM=${INPUT_BAM} \
-e SAMPLE_NAME=${SAMPLE_NAME} \
-v ${PWD}:${HOME} -w ${HOME} google/deepvariant:1.1.0 \
/opt/deepvariant/bin/run_deepvariant \
  --make_examples_extra_args 'min_mapping_quality=1' \
  --customized_model model.ckpt-25600 \
  --model_type=WGS \
  --ref=${REF_FASTA} \
  --reads=indel_realigned.${INPUT_BAM} \
  --output_vcf=${SAMPLE_NAME}.deepvariant.vcf.gz \
  --output_gvcf=${SAMPLE_NAME}.deepvariant.g.vcf.gz \
  --intermediate_results_dir=tmp_deepvariant.${SAMPLE_NAME} \
  --num_shards=16

# Evaluate Called File
cd $WORK_DIR
rm -r $WORK_DIR/happy_vcfeval_output
TRUTH_VCF="${SAMPLE_NAME}_GRCh38_1_22_v4.2.1_benchmark.xian_adam.vcf.gz" \
TRUTH_BED="${SAMPLE_NAME}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.xian_adam.bed" \
CALLED_VCF="${SAMPLE_NAME}.deepvariant.vcf.gz" \
docker run \
-e TRUTH_VCF=${TRUTH_VCF} \
-e TRUTH_BED=${TRUTH_BED} \
-e CALLED_VCF=${CALLED_VCF} \
-e REF_FASTA=${REF_FASTA} \
-v ${PWD}:${HOME} -w ${HOME} jmcdani20/hap.py:v0.3.12 \
  /opt/hap.py/bin/hap.py \
  ${VCF_TRUTH} \
  ${CALLED_VCF} \
  -f ${TRUTH_BED} \
  --reference ${REF_FASTA} \
  --threads 16 \
  --engine=vcfeval \
  -o happy_vcfeval_output





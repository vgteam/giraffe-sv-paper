#!/bin/bash
module load python/3.7 singularity
#source /data/markellocj/test_toil_vg_run/toil_vg_tools/toilvg_venv/bin/activate
source /data/Udpbinfo/usr/markellocj/test_vg_pedigree_giraffe_grch38/tiny_unittest_workdir/toilvg_test_venv/bin/activate
INPUT_DIR="/data/markellocj/raw_read_data/HG002_cohort/HG002"
OUTSTORE="/data/Udpbinfo/usr/markellocj/redo_giraffe_paper_revisions/xian_adam_grch38_liftover_experiments/vg_map_snp1kg_all_nosegdup_2x250_outstore"
JOBSTORE="/data/Udpbinfo/usr/markellocj/redo_giraffe_paper_revisions/xian_adam_grch38_liftover_experiments/vg_map_snp1kg_all_nosegdup_2x250_jobstore"
LOGFILE="/data/Udpbinfo/usr/markellocj/redo_giraffe_paper_revisions/xian_adam_grch38_liftover_experiments/vg_map_snp1kg_all_nosegdup_2x250_workflow.log"
TMP_DIR="/data/Udpbinfo/usr/markellocj/redo_giraffe_paper_revisions/xian_adam_grch38_liftover_experiments/tmp_vg_map_nosegdup_2x250"
export TOIL_SLURM_ARGS='-t 15:00:00'
export SINGULARITY_CACHEDIR=/data/markellocj/singularity_cache
cd /data/Udpbinfo/usr/markellocj/redo_giraffe_paper_revisions/xian_adam_grch38_liftover_experiments
rm -fr ${LOGFILE} ${TMP_DIR} ${OUTSTORE}
mkdir -p ${TMP_DIR} ${OUTSTORE}
toil clean ${JOBSTORE}
toil-vg map \
${JOBSTORE} \
HG002 \
${OUTSTORE} \
--batchSystem slurm \
--container Singularity \
--logFile ${LOGFILE} \
--workDir ${TMP_DIR} \
--cleanWorkDir onSuccess \
--whole_genome_config \
--vg_docker 'quay.io/vgteam/vg:v1.31.0' \
--bam_output \
--id_ranges /data/markellocj/graph_reference/snp1kg_ref_wgs/xians_adams_liftover_grch38/no_segdup/path.list \
--fastq ${INPUT_DIR}/HG002_read_pair_1.fq.gz ${INPUT_DIR}/HG002_read_pair_2.fq.gz \
--xg_index /data/markellocj/graph_reference/snp1kg_ref_wgs/xians_adams_liftover_grch38/no_segdup/1000GPlons_hs38d1_filter.xg \
--gcsa_index /data/markellocj/graph_reference/snp1kg_ref_wgs/xians_adams_liftover_grch38/no_segdup/1000GPlons_hs38d1_filter.gcsa \
--statePollingWait 120 \
--rescueJobsFrequency 120


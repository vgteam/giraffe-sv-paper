#!/bin/bash
module load singularity samtools picard
cd /data/Udpbinfo/usr/markellocj/redo_giraffe_paper_revisions/xian_adam_grch38_liftover_experiments
singularity exec -H ${PWD}:${HOME} \
-B /data/markellocj/graph_reference/primary_ref_wgs/xian_adam_grch38:${HOME}/xian_adam_grch38 \
-B /data/markellocj/raw_read_data/HG002_cohort_precision_fda_reads/HG002:${HOME}/input_files \
-B /data/markellocj/fasta_references/grch38_reference:${HOME}/grch38_reference \
--pwd ${HOME} \
docker://quay.io/vgteam/vg:v1.31.0 \
vg giraffe \
-x xian_adam_grch38/primaryhs38d1.xg \
-H xian_adam_grch38/primaryhs38d1.cover.gbwt \
-g xian_adam_grch38/primaryhs38d1.cover.gg \
-d xian_adam_grch38/primaryhs38d1.dist \
-f input_files/HG002.novaseq.pcr-free.35x.R1.fastq.gz \
-f input_files/HG002.novaseq.pcr-free.35x.R2.fastq.gz \
--output-format BAM \
--ref-paths grch38_reference/xian_adams_grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.dict \
-p \
-t 140 > giraffe_all_wgs.RAW.XIAN_ADAM_GRCH38.PRIMARY.NO_MIN.HG002.bam

samtools sort --threads 140 -O BAM giraffe_all_wgs.RAW.XIAN_ADAM_GRCH38.PRIMARY.NO_MIN.HG002.bam > giraffe_all_wgs.RAW.XIAN_ADAM_GRCH38.PRIMARY.NO_MIN.HG002.positionsorted.bam && rm giraffe_all_wgs.RAW.XIAN_ADAM_GRCH38.PRIMARY.NO_MIN.HG002.bam
samtools index -@ 140 giraffe_all_wgs.RAW.XIAN_ADAM_GRCH38.PRIMARY.NO_MIN.HG002.positionsorted.bam

#time java -Xmx100g -XX:ParallelGCThreads=130 -jar $PICARDJARPATH/picard.jar \
#ReorderSam \
#VALIDATION_STRINGENCY=SILENT \
#INPUT=giraffe_all_wgs.grch38_pgrc_decoys_mar13.HG003.positionsorted.bam \
#OUTPUT=giraffe_all_wgs.grch38_pgrc_decoys_mar13.HG003.positionsorted.reordered.bam \
#SEQUENCE_DICTIONARY=/data/markellocj/fasta_references/grch38_reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.dict


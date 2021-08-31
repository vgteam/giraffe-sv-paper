#!/bin/bash
module load singularity samtools picard
cd /data/Udpbinfo/usr/markellocj/redo_giraffe_paper_revisions/xian_adam_grch38_liftover_experiments
singularity exec -H ${PWD}:${HOME} \
-B /data/markellocj/graph_reference/snp1kg_ref_wgs/xians_adams_liftover_grch38/no_segdup:${HOME}/no_segdup \
-B /data/markellocj/raw_read_data/HG002_cohort/HG002:${HOME}/input_files \
-B /data/markellocj/fasta_references/grch38_reference/xian_adams_grch38:${HOME}/xian_adams_grch38 \
--pwd ${HOME} \
docker://quay.io/vgteam/vg:v1.31.0 \
vg giraffe \
-x no_segdup/1000GPlons_hs38d1_filter.xg \
-H no_segdup/1000GPlons_hs38d1_filter.sampled.64.gbwt \
-m no_segdup/1000GPlons_hs38d1_filter.sampled.64.min \
-g no_segdup/1000GPlons_hs38d1_filter.sampled.64.gg \
-d no_segdup/1000GPlons_hs38d1_filter.dist \
-f input_files/HG002_read_pair_1.fq.gz \
-f input_files/HG002_read_pair_2.fq.gz \
--output-format BAM \
--ref-paths xian_adams_grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.no_segdup.dict \
-b fast \
-p \
-t 100 > giraffe_all_wgs.RAW.1000GP.FAST.NO_SEGDUP.HG002.2x250.bam

samtools sort --threads 100 -O BAM giraffe_all_wgs.RAW.1000GP.FAST.NO_SEGDUP.HG002.2x250.bam > giraffe_all_wgs.RAW.1000GP.FAST.NO_SEGDUP.HG002.2x250.positionsorted.bam && rm giraffe_all_wgs.RAW.1000GP.FAST.NO_SEGDUP.HG002.2x250.bam
samtools index -@ 100 giraffe_all_wgs.RAW.1000GP.FAST.NO_SEGDUP.HG002.2x250.positionsorted.bam

#time java -Xmx100g -XX:ParallelGCThreads=130 -jar $PICARDJARPATH/picard.jar \
#ReorderSam \
#VALIDATION_STRINGENCY=SILENT \
#INPUT=giraffe_all_wgs.grch38_pgrc_decoys_mar13.HG003.positionsorted.bam \
#OUTPUT=giraffe_all_wgs.grch38_pgrc_decoys_mar13.HG003.positionsorted.reordered.bam \
#SEQUENCE_DICTIONARY=/data/markellocj/fasta_references/grch38_reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.dict


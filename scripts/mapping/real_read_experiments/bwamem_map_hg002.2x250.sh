#!/bin/bash
cd /data/Udpbinfo/usr/markellocj/redo_giraffe_paper_revisions/xian_adam_grch38_liftover_experiments
module load bwa samtools picard
time bwa mem \
-t 32 \
/data/markellocj/fasta_references/grch38_reference/xian_adams_grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna \
/data/markellocj/raw_read_data/HG002_cohort/HG002/HG002_read_pair_1.fq.gz \
/data/markellocj/raw_read_data/HG002_cohort/HG002/HG002_read_pair_2.fq.gz \
> bwamem_hg002_wgs.2x250.sam

samtools sort -@ 32 -O BAM bwamem_hg002_wgs.2x250.sam > bwamem_hg002_wgs.2x250.positionsorted.bam && rm bwamem_hg002_wgs.2x250.sam
samtools index -@ 32 bwamem_hg002_wgs.2x250.positionsorted.bam

time java -Xmx100g -XX:ParallelGCThreads=32 -jar $PICARDJARPATH/picard.jar \
ReorderSam \
VALIDATION_STRINGENCY=SILENT \
INPUT=bwamem_hg002_wgs.2x250.positionsorted.bam \
OUTPUT=bwamem_hg002_wgs.2x250.positionsorted.reordered.bam \
SEQUENCE_DICTIONARY=/data/markellocj/fasta_references/grch38_reference/xian_adams_grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.dict

samtools index -@ 32 bwamem_hg002_wgs.2x250.positionsorted.reordered.bam && rm bwamem_hg002_wgs.2x250.positionsorted.bam bwamem_hg002_wgs.2x250.positionsorted.bam.bai


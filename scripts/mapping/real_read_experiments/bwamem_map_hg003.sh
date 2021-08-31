#!/bin/bash
cd /data/Udpbinfo/usr/markellocj/redo_giraffe_paper_revisions/xian_adam_grch38_liftover_experiments
module load bwa samtools picard
time bwa mem \
-t 32 \
/data/markellocj/fasta_references/grch38_reference/xian_adams_grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna \
/data/markellocj/raw_read_data/HG002_cohort_precision_fda_reads/HG003/HG003.novaseq.pcr-free.35x.R1.fastq.gz \
/data/markellocj/raw_read_data/HG002_cohort_precision_fda_reads/HG003/HG003.novaseq.pcr-free.35x.R2.fastq.gz \
> bwamem_hg003_wgs.sam

samtools sort -@ 32 -O BAM bwamem_hg003_wgs.sam > bwamem_hg003_wgs.positionsorted.bam && rm bwamem_hg003_wgs.sam
samtools index -@ 32 bwamem_hg003_wgs.positionsorted.bam

time java -Xmx100g -XX:ParallelGCThreads=32 -jar $PICARDJARPATH/picard.jar \
ReorderSam \
VALIDATION_STRINGENCY=SILENT \
INPUT=bwamem_hg003_wgs.positionsorted.bam \
OUTPUT=bwamem_hg003_wgs.positionsorted.reordered.bam \
SEQUENCE_DICTIONARY=/data/markellocj/fasta_references/grch38_reference/xian_adams_grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.dict

samtools index -@ 32 bwamem_hg003_wgs.positionsorted.reordered.bam && rm bwamem_hg003_wgs.positionsorted.bam bwamem_hg003_wgs.positionsorted.bam.bai


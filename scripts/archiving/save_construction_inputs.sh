#!/usr/bin/env bash

# save_construction_inputs.sh: Save files needed to reconstruct graphs, and other input files that may be hard to obtain.

set -x

DEST_DIR=/nanopore/cgl/data/giraffe

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

CONSTRUCTION_DIR="${DEST_DIR}/construction"
mkdir -p "${CONSTRUCTION_DIR}"

download s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.gz "${CONSTRUCTION_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.gz"
download s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.renamed.fna.gz "${CONSTRUCTION_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.renamed.fna.gz"
wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v2.0/GRCh38/SegmentalDuplications/GRCh38_segdups_gt10kb.bed.gz "${CONSTRUCTION_DIR}/GRCh38_segdups_gt10kb.bed.gz"
wget_download ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz "${CONSTRUCTION_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
wget_download ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz "${CONSTRUCTION_DIR}/GCA_000786075.2_hs38d1_genomic.fna.gz"
for CHROM in {1..22} X Y ; do
    download s3://vg-data/1kg_GRCh38/variants/ALL.chr${CHROM}_GRCh38.genotypes.20170504.vcf.gz "${CONSTRUCTION_DIR}/ALL.chr${CHROM}_GRCh38.genotypes.20170504.vcf.gz"
    download s3://vg-data/1kg_GRCh38/variants/ALL.chr${CHROM}_GRCh38.genotypes.20170504.vcf.gz.tbi "${CONSTRUCTION_DIR}/ALL.chr${CHROM}_GRCh38.genotypes.20170504.vcf.gz.tbi"
    download s3://vg-data/1kg_GRCh38/variants/subsets/no_segdups_gt10kb/ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz "${CONSTRUCTION_DIR}/ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz"
    download s3://vg-data/1kg_GRCh38/variants/subsets/no_segdups_gt10kb/ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz.tbi "${CONSTRUCTION_DIR}/ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz.tbi"
done
download s3://glennhickey/outstore/HGSVC-jan5/HGSVC.haps.vcf.gz "${CONSTRUCTION_DIR}/HGSVC.haps.vcf.gz"

# Also grab the real reads
READS_DIR="${DEST_DIR}/mapping/reads/real"
mkdir -p "${READS_DIR}/NA19239"
download s3://vg-k8s/profiling/reads/real/NA19239/novaseq6000-ERR3239454-shuffled-1m.fq.gz "${READS_DIR}/NA19239/novaseq6000-ERR3239454-shuffled-1m.fq.gz"
mkdir -p "${READS_DIR}/NA19240"
download s3://vg-k8s/profiling/reads/real/NA19240/hiseq2500-ERR309934-shuffled-1m.fq.gz "${READS_DIR}/NA19240/hiseq2500-ERR309934-shuffled-1m.fq.gz"
download s3://vg-k8s/profiling/reads/real/NA19240/hiseqxten-SRR6691663-shuffled-1m.fq.gz "${READS_DIR}/NA19240/hiseqxten-SRR6691663-shuffled-1m.fq.gz"

# And the HAL files for the yeast experiments
mkdir -p "${CONSTRUCTION_DIR}/yeast/cactus_all"
mkdir -p "${CONSTRUCTION_DIR}/yeast/cactus_four"
copy /public/groups/cgl/users/daheller/yeast_graph/graphs/cactus_all/cactusoutput.hal "${CONSTRUCTION_DIR}/yeast/cactus_all/cactusoutput.hal"
copy /public/groups/cgl/users/daheller/yeast_graph/graphs/cactus_four/cactusoutput.hal "${CONSTRUCTION_DIR}/yeast/cactus_four/cactusoutput.hal"
mkdir -p "${READS_DIR}/yeast"
copy /public/groups/cgl/users/daheller/yeast_graph/illumina_reads/SRR4074257.fastq.gz "${READS_DIR}/yeast/SRR4074257.fastq.gz"
copy /public/groups/cgl/users/daheller/yeast_graph/assemblies/assemblies_raw/S288C.genome.fa "${CONSTRUCTION_DIR}/yeast/S288C.genome.fa"

# Put the README in place
SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"
cp "${SCRIPT_DIR}/archive-readme.md" "${DEST_DIR}/README.md"

#!/usr/bin/env bash
# Was run inside a Docker container: quay.io/ucsc_cgl/toil:5.3.0a1-1f0930b7b9ecc31ca556d15ab07ff836ba85eb23-py3.7
# When rerunning, use quay.io/vgteam/vg:ci-2859-9275d40d72da4232f07a82173da54a996ff31527 or another build of vg 9275d40d72da4232f07a82173da54a996ff31527 as --vg_docker
# If rerunning, the (currently upcoming) Toil 5.3 and vg 1.32 releases are probably the right choices.
# Should work on an AWS r5d.24xlarge, but was moved to an r5a.8xlarge with a 5000GB root volume when only the GCSA indexes were left to save money.

mkdir /tmp/work
cd /tmp/work
echo "Making graphs"
virtualenv --system-site-packages --python python3 venv
. venv/bin/activate
pip3 install --upgrade git+https://github.com/DataBiosphere/toil.git@5e6a1d6952367838a7220d7e1ef5c157f28c5813#egg=toil
reset
pip3 install --upgrade git+https://github.com/vgteam/toil-vg.git@9411d73c09b17440e5a90ed43a79a3ec4db0b354#egg=toil-vg
reset
# VCF order must be 1-22, X, Y, chrM to match FASTA
# Leaving M as "chrM" to match the SV graphs

VCFS=()
for CHROM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; do
    VCFS+=("s3://vg-k8s/users/adamnovak/projects/1000gp-giraffe/1000g/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHROM}.filtered.shapeit2-duohmm-phased.vcf.gz")
done
VCFS+=("s3://vg-k8s/users/adamnovak/projects/1000gp-giraffe/1000g/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.noME.vcf.gz" \
"s3://vg-k8s/users/adamnovak/projects/1000gp-giraffe/1000g/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrY.recalibrated_variants.vcf.gz" \
"s3://vg-k8s/users/adamnovak/projects/1000gp-giraffe/1000g/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_others.recalibrated_variants.vcf.gz")

export TOIL_KUBERNETES_HOST_PATH=/var/lib/toil
export SINGULARITY_CACHEDIR=/var/lib/toil/singularity-cache
toil-vg construct \
    aws:us-west-2:adamnovak-make-1000gp-graphs-2 \
    aws:us-west-2:vg-k8s/profiling/graphs/v3-2/for-NA19239/1000gp/hs38d1 \
    --fasta s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.gz \
    --coalesce_regions s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.minor_contigs.tsv \
    --vcf "${VCFS[@]}" \
    --vcf_phasing "${VCFS[@]}" \
    --fasta_regions \
    --remove_chr_prefix \
    --alt_paths \
    --out_name 1000GP_hs38d1 \
    --all_index \
    --force_phasing True \
    --gbwt_prune \
    --filter_samples NA19239 NA19240 \
    --pangenome \
    --merge_graphs \
    --keep_vcfs \
    --whole_genome_config \
    --realTimeLogging \
    --logInfo \
    --batchSystem single_machine \
    --disableCaching false \
    --vg_docker quay.io/adamnovak/vg:graphbuild \
    --defaultDisk 200G --defaultMemory 200G --restart 2>&1 | tee -a log.txt
                              

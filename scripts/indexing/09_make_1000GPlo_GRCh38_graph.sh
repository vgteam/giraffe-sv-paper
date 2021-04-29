#!/usr/bin/env bash
# Was run inside a Docker container: quay.io/ucsc_cgl/toil:5.3.0-py3.7
# If rerunning, the Toil 5.4 and vg v1.32.0 releases are probably the right choices.
# Should work on an AWS r5d.24xlarge with lots of EBS storage assigned, but
# would be faster to run with --batchSystem kubernetes against a Kubernetes
# cluster.

rm -Rf /tmp/work  
mkdir /tmp/work
cd /tmp/work
virtualenv --system-site-packages --python python3 venv
. venv/bin/activate
pip3 install --upgrade git+https://github.com/vgteam/toil-vg.git@d00b307f7739e3904e54d4fa440bb053d7194235#egg=toil-vg
toil-vg generate-config --whole_genome >config.cfg
sed -i'' config.cfg -e "s/gcsa-index-mem: '110G'/gcsa-index-mem: '700G'/g" -e "s/gcsa-index-disk: '2200G'/gcsa-index-disk: '4096G'/g"
toil clean aws:us-west-2:adamnovak-make-1000gp-graphs-3
# VCF order must be 1-22, X, Y to match FASTA
VCFS=()
for CHROM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; do
  VCFS+=("s3://vg-data/1kg_GRCh38/variants/ALL.chr${CHROM}_GRCh38.genotypes.20170504.vcf.gz")
done
VCFS+=("s3://vg-data/1kg_GRCh38/variants/ALL.chrX_GRCh38.genotypes.20170504.vcf.gz" \
     "s3://vg-data/1kg_GRCh38/variants/ALL.chrY_GRCh38.genotypes.20170504.vcf.gz")
export TOIL_KUBERNETES_HOST_PATH=/var/lib/toil
export SINGULARITY_CACHEDIR=/var/lib/toil/singularity-cache
toil-vg construct \
  aws:us-west-2:adamnovak-make-1000gp-graphs-3 \
  aws:us-west-2:vg-k8s/profiling/graphs/v3/for-NA19239/1000gplo/hs38d1 \
  --fasta s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.gz \
  --coalesce_regions s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.minor_contigs.tsv \
  --vcf "${VCFS[@]}" \
  --vcf_phasing "${VCFS[@]}" \
  --fasta_regions \
  --remove_chr_prefix \
  --alt_paths \
  --out_name 1000GPlo_hs38d1 \
  --all_index \
  --force_phasing True \
  --gbwt_prune \
  --filter_samples NA19239 NA19240 \
  --pangenome \
  --merge_graphs \
  --keep_vcfs \
  --config config.cfg \
  --realTimeLogging \
  --logInfo \
  --batchSystem single_machine \
  --disableCaching false \
  --vg_docker quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975 \
  --container Singularity \
  --defaultDisk 200G --defaultMemory 200G 2>&1 | tee -a log.txt
          


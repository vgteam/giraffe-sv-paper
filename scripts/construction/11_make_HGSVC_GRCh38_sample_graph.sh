#!/usr/bin/env bash
# Run in screen on a Kubernetes cluster: TOIL_APPLIANCE_SELF=quay.io/ucsc_cgl/toil:5.4.0a1-1530a9190357fc058333f3e929049ef9593a6784-py3.7 toil launch-cluster --provisioner aws -T kubernetes -z us-west-2d adamnovak-toil-vg --leaderNodeType t3a.medium --nodeTypes=t3a.medium,r5ad.24xlarge,r5d.24xlarge/r5ad.24xlarge:2.00,i3.8xlarge:0.80 --workers 1-4,0-1,0-4,0-3 --keyPairName anovak@soe.ucsc.edu
# If rerunning, the Toil 5.4 and vg v1.32.0 releases are probably the right choices.

rm -Rf /tmp/work  
mkdir /tmp/work
cd /tmp/work
virtualenv --system-site-packages --python python3 venv
. venv/bin/activate
pip3 install --upgrade git+https://github.com/vgteam/toil-vg.git@b319a1b22df6dac585b7f95bc1a603577452d443#egg=toil-vg
toil-vg generate-config --whole_genome >config.cfg
sed -s -i'' config.cfg -e "s/construct-mem: '64G'/construct-mem: '128G'/g" -e "s/construct-disk: '64G'/construct-disk: '128G'/g" -e "s/gbwt-index-mem: '35G'/gbwt-index-mem: '70G'/g" -e "s/gbwt-index-disk: '100G'/gbwt-index-disk: '200G'/g"
VCFS=(s3://glennhickey/outstore/HGSVC-jan5/HGSVC.haps.vcf.gz)
export TOIL_KUBERNETES_HOST_PATH=/var/lib/toil
export SINGULARITY_CACHEDIR=/var/lib/toil/singularity-cache
toil-vg construct \
  aws:us-west-2:adamnovak-make-hgsvc-sample-graph \
  aws:us-west-2:vg-k8s/profiling/graphs/v3/for-NA19240/hgsvc/hs38d1 \
  --fasta s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.gz \
  --coalesce_regions s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.minor_contigs.tsv \
  --vcf "${VCFS[@]}" \
  --vcf_phasing "${VCFS[@]}" \
  --fasta_regions \
  --remove_chr_prefix \
  --alt_paths \
  --out_name HGSVC_hs38d1 \
  --all_index \
  --force_phasing True \
  --gbwt_prune \
  --sample_graph NA19240 \
  --merge_graphs \
  --keep_vcfs \
  --config config.cfg \
  --realTimeLogging \
  --logInfo \
  --batchSystem kubernetes \
  --disableCaching false \
  --vg_docker quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975 \
  --container Singularity \
  --defaultDisk 400G --defaultMemory 400G 2>&1 | tee -a log.txt
              
              


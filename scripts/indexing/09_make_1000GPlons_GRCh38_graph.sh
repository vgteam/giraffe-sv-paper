#!/usr/bin/env bash
# If rerunning, the Toil 6 and vg v1.33.0 releases are probably the right choices.
# Run in screen on a Kubernetes cluster: TOIL_APPLIANCE_SELF=quay.io/ucsc_cgl/toil:5.4.0a1-1530a9190357fc058333f3e929049ef9593a6784-py3.7 toil launch-cluster --provisioner aws -T kubernetes -z us-west-2a adamnovak-toil-vg --leaderNodeType t3a.medium --nodeTypes=t3a.medium,r5ad.24xlarge,r5d.24xlarge/r5ad.24xlarge:2.50,i3.8xlarge:1.50 --workers 1-4,0-1,0-8,0-6 --keyPairName anovak@soe.ucsc.edu
# Doesn't build the GCSA indexes for vg map.
# Does build the sample graph positive control.
# When rerunning, use --defaultPreemptable and/or throw max 5 of the on-demand workers at it

rm -Rf /tmp/work  
mkdir /tmp/work
cd /tmp/work
virtualenv --system-site-packages --python python3 venv
. venv/bin/activate
pip3 install --upgrade git+https://github.com/vgteam/toil-vg.git@b319a1b22df6dac585b7f95bc1a603577452d443#egg=toil-vg
toil-vg generate-config --whole_genome >config.cfg
sed -i'' config.cfg -e "s/gcsa-index-mem: '110G'/gcsa-index-mem: '700G'/g" -e "s/gcsa-index-disk: '2200G'/gcsa-index-disk: '4096G'/g" -e "s/construct-mem: '64G'/construct-mem: '128G'/g" -e "s/construct-disk: '64G'/construct-disk: '128G'/g" -e "s/gbwt-index-mem: '35G'/gbwt-index-mem: '70G'/g" -e "s/gbwt-index-disk: '100G'/gbwt-index-disk: '200G'/g"
toil clean aws:us-west-2:adamnovak-make-1000gplons-graphs
# VCF order must be 1-22, X, Y to match FASTA
VCFS=()
for CHROM in {1..22} X Y ; do
  VCFS+=("s3://vg-data/1kg_GRCh38/variants/subsets/no_segdups_gt10kb/ALL.chr${CHROM}_GRCh38.genotypes.20170504.no_segdups_gt10kb.vcf.gz")
done
export TOIL_KUBERNETES_HOST_PATH=/var/lib/toil
export SINGULARITY_CACHEDIR=/var/lib/toil/singularity-cache
toil-vg construct \
  aws:us-west-2:adamnovak-make-1000gplons-graphs \
  aws:us-west-2:vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1 \
  --fasta s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.renamed.fna.gz \
  --coalesce_regions s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.minor_contigs.tsv \
  --vcf "${VCFS[@]}" \
  --vcf_phasing "${VCFS[@]}" \
  --fasta_regions \
  --alt_paths \
  --out_name 1000GPlons_hs38d1 \
  --xg_index --gbwt_index --trivial_snarls_index --distance_index \
  --force_phasing True \
  --filter_samples NA19239 NA19240 \
  --pangenome \
  --sample_graph NA19239 \
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
          


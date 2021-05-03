#!/usr/bin/env bash
# Run in a screen on a Toil cluster: TOIL_APPLIANCE_SELF=quay.io/ucsc_cgl/toil:4.2.0a1-edb0e583e91452e80ef6add0f9f0e8eae9dbc2d4-py3.7 toil launch-cluster -z us-west-2a adamnovakgraphbuild --leaderNodeType t2.medium --keyPairName anovak@kolossus
# Was initially started with toil-vg
# git+https://github.com/vgteam/toil-vg.git@fe90d202f6f9a9e9564c38256475b9dd4e9e368c#egg=toil-vg
# and Toil
# quay.io/ucsc_cgl/toil:4.2.0a1-99c950d3bc5aafb460e2f751decb034848c792f2-py3.7
# and then restarted (--restart) with toil-vg
# git+https://github.com/vgteam/toil-vg.git@cfe15995a8da67257af1b9aca25fa1a686a839c4#egg=toil-vg
# and
# git+https://github.com/vgteam/toil-vg.git@fbbf88ccbabf2c82a658d9d4028ae0ce9139a2c7#egg=toil-vg
# and Toil
# quay.io/ucsc_cgl/toil:4.2.0a1-8739e4a730faa078cba4e515e263ad1e60b00a8e-py3.7
# and
# quay.io/ucsc_cgl/toil:4.2.0a1-edb0e583e91452e80ef6add0f9f0e8eae9dbc2d4-py3.7
# Should not be re-run as is because it will fail at the GCSA indexing step; we used the files that were generated and built the GCSA manually.

mkdir /tmp/work
cd /tmp/work
virtualenv --system-site-packages --python python3 venv
. venv/bin/activate

pip3 install --upgrade git+https://github.com/vgteam/toil-vg.git@cfe15995a8da67257af1b9aca25fa1a686a839c4#egg=toil-vg

toil-vg generate-config --whole_genome | sed "s/xg-index-disk.*/xg-index-disk: '250G'/g" | sed "s/xg-index-mem.*/xg-index-mem: '340G'/g" > config.cfg

toil-vg construct \
aws:us-west-2:adamnovak-make-hgsvc-graphs-2 \
aws:us-west-2:vg-k8s/profiling/graphs/v2-2/for-NA19240/hgsvc/hs38d1 \
--fasta s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.gz \
--vcf s3://glennhickey/outstore/HGSVC-jan5/HGSVC.haps.vcf.gz \
--vcf_phasing s3://glennhickey/outstore/HGSVC-jan5/HGSVC.haps.vcf.gz \
--fasta_regions \
--remove_chr_prefix \
--out_name HGSVC_hs38d1 \
--flat_alts \
--all_index \
--force_phasing True \
--gbwt_prune \
--haplo_sample NA19240 \
--pangenome \
--normalize \
--merge_graphs \
--keep_vcfs \
--config config.cfg \
--realTimeLogging \
--logInfo \
'--container' Docker \
--disableCaching false \
--vg_docker quay.io/vgteam/vg:ci-1954-8ff022c3a36c5fbc2a63faf477c5bf9ac37e29d7 \
--batchSystem mesos \
--provisioner=aws \
--nodeTypes=r5.24xlarge \
--maxNodes=20 \
--minNodes=0 \
--nodeStorage 2300 \
--retryCount 2 \
--enableUnlimitedPreemptableRetries \
 2>&1 | tee -a log.txt

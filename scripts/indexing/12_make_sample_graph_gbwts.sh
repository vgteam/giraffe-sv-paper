#!/usr/bin/env bash
# Run with vg v1.32.0 on out 64 core/1TB memory machine, but this level of hardware is not expected to be necessary.

# 1000GPlo
aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplo/hs38d1/1000GPlo_hs38d1_NA19239_sample_withref.vg .

CHROM_OPTS=()
for CHROM in {1..22} X Y ; do
    aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplo/hs38d1/1000GPlo_hs38d1-vcfs/ALL.chr${CHROM}_GRCh38.genotypes.20170504_NA19239.vcf.gz .
    aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplo/hs38d1/1000GPlo_hs38d1-vcfs/ALL.chr${CHROM}_GRCh38.genotypes.20170504_NA19239.vcf.gz.tbi .
    CHROM_OPTS+=(ALL.chr${CHROM}_GRCh38.genotypes.20170504_NA19239.vcf.gz)
done

vg index -x 1000GPlo_hs38d1_NA19239_sample_withref.xg 1000GPlo_hs38d1_NA19239_sample_withref.vg

vg gbwt -p --force-phasing --discard-overlaps -x 1000GPlo_hs38d1_NA19239_sample_withref.vg -o 1000GPlo_hs38d1_NA19239_sample_withref.force.gbwt -v "${CHROM_OPTS[@]}"

vg gbwt -p -g 1000GPlo_hs38d1_NA19239_sample_withref.force.sampled64.gg -o 1000GPlo_hs38d1_NA19239_sample_withref.force.sampled64.gbwt -x 1000GPlo_hs38d1_NA19239_sample_withref.xg -l 1000GPlo_hs38d1_NA19239_sample_withref.force.gbwt -n 64

aws s3 cp 1000GPlo_hs38d1_NA19239_sample_withref.force.gbwt s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplo/hs38d1/1000GPlo_hs38d1_NA19239_sample_withref.force.gbwt
aws s3 cp 1000GPlo_hs38d1_NA19239_sample_withref.force.sampled64.gbwt s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplo/hs38d1/1000GPlo_hs38d1_NA19239_sample_withref.force.sampled64.gbwt
aws s3 cp 1000GPlo_hs38d1_NA19239_sample_withref.force.sampled64.gg s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplo/hs38d1/1000GPlo_hs38d1_NA19239_sample_withref.force.sampled64.gg
aws s3 cp 1000GPlo_hs38d1_NA19239_sample_withref.xg s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplo/hs38d1/1000GPlo_hs38d1_NA19239_sample_withref.xg

# HGSVC
aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1_NA19240_sample_withref.vg .

aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1-vcfs/HGSVC.haps_NA19240.vcf.gz .
aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1-vcfs/HGSVC.haps_NA19240.vcf.gz.tbi .
                    
vg index --force-phasing --discard-overlaps -x HGSVC_hs38d1_NA19240_sample_withref.xg -G HGSVC_hs38d1_NA19240_sample_withref.force.gbwt -v HGSVC.haps_NA19240.vcf.gz HGSVC_hs38d1_NA19240_sample_withref.vg

vg gbwt -p -g HGSVC_hs38d1_NA19240_sample_withref.force.sampled64.gg -o HGSVC_hs38d1_NA19240_sample_withref.force.sampled64.gbwt -x HGSVC_hs38d1_NA19240_sample_withref.xg -l HGSVC_hs38d1_NA19240_sample_withref.force.gbwt -n 64

aws s3 cp HGSVC_hs38d1_NA19240_sample_withref.xg s3://vg-k8s/profiling/graphs/v3/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1_NA19240_sample_withref.xg
aws s3 cp HGSVC_hs38d1_NA19240_sample_withref.force.gbwt s3://vg-k8s/profiling/graphs/v3/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1_NA19240_sample_withref.force.gbwt
aws s3 cp HGSVC_hs38d1_NA19240_sample_withref.force.sampled64.gbwt s3://vg-k8s/profiling/graphs/v3/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1_NA19240_sample_withref.force.sampled64.gbwt
aws s3 cp HGSVC_hs38d1_NA19240_sample_withref.force.sampled64.gg s3://vg-k8s/profiling/graphs/v3/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1_NA19240_sample_withref.force.sampled64.gg
              


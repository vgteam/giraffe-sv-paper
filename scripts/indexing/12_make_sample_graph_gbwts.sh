#!/usr/bin/env bash
# Run with vg v1.32.0 on our 64 core/1TB memory machine, but this level of hardware is not expected to be necessary.

set -ex

# First make an NA19240-only GBWT for the full HGSVC graph
if [[ ! -e HGSVC_hs38d1.gbwt ]] ; then
    aws s3 cp s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.gbwt .
fi
if [[ ! -e HGSVC_hs38d1.xg ]] ; then
    aws s3 cp s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.xg .
fi
if [[ ! -e HGSVC_hs38d1.onlyNA19240.gbwt ]] ; then
    # Note that vg gbwt accepts multiple -R options in one command but
    # mishandles tham and will remove the wrong samples. See
    # https://github.com/vgteam/vg/issues/3284. To work around this we remove
    # the samples in sequence.
    vg gbwt HGSVC_hs38d1.gbwt -o HGSVC_hs38d1.onlyNA19240.step1.gbwt -R HG00514
    vg gbwt HGSVC_hs38d1.onlyNA19240.step1.gbwt -o HGSVC_hs38d1.onlyNA19240.gbwt -R HG00733
    rm HGSVC_hs38d1.onlyNA19240.step1.gbwt
fi
if [[ ! -e HGSVC_hs38d1.onlyNA19240.augment.gbwt ]] ; then
    vg gbwt HGSVC_hs38d1.onlyNA19240.gbwt -o HGSVC_hs38d1.onlyNA19240.augment.gbwt -x HGSVC_hs38d1.xg -a
fi
aws s3 cp HGSVC_hs38d1.onlyNA19240.augment.gbwt s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.onlyNA19240.augment.gbwt


# Now do a sample graph GBWT for 1000GPlons

# Download everything
if [[ ! -e 1000GPlons_hs38d1_NA19239_sample_withref.vg ]] ; then
    aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_NA19239_sample_withref.vg .
fi

if [[ ! -e 1000GPlons_hs38d1_NA19239_sample_withref.xg ]] ; then
    aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_NA19239_sample_withref.xg . || true
fi

if [[ ! -e 1000GPlons_hs38d1_NA19239_sample_withref.force.gbwt ]] ; then
    aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_NA19239_sample_withref.force.gbwt . || true
fi

CHROM_OPTS=()
for CHROM in {1..22} X Y ; do
    if [[ ! -e ALL.chr${CHROM}_GRCh38.genotypes.20170504_NA19239.vcf.gz ]] ; then
        aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1-vcfs/ALL.chr${CHROM}_GRCh38.genotypes.20170504_NA19239.vcf.gz .
    fi
    if [[ ! -e ALL.chr${CHROM}_GRCh38.genotypes.20170504_NA19239.vcf.gz.tbi ]] ; then
        aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1-vcfs/ALL.chr${CHROM}_GRCh38.genotypes.20170504_NA19239.vcf.gz.tbi .
    fi
    CHROM_OPTS+=(ALL.chr${CHROM}_GRCh38.genotypes.20170504_NA19239.vcf.gz)
done

# Index 1000GPlons sample graph

if [[ ! -e 1000GPlons_hs38d1_NA19239_sample_withref.xg ]] ; then
    vg index -x 1000GPlons_hs38d1_NA19239_sample_withref.xg 1000GPlons_hs38d1_NA19239_sample_withref.vg
fi

if [[ ! -e 1000GPlons_hs38d1_NA19239_sample_withref.force.gbwt ]] ; then
    vg gbwt -p --force-phasing --discard-overlaps -x 1000GPlons_hs38d1_NA19239_sample_withref.vg -o 1000GPlons_hs38d1_NA19239_sample_withref.force.gbwt -v "${CHROM_OPTS[@]}"
fi

# If we use a local haplotype cover we will introduce switch errors; augment instead.
if [[ ! -e 1000GPlons_hs38d1_NA19239_sample_withref.force.augment.gg || ! -e 1000GPlons_hs38d1_NA19239_sample_withref.force.augment.gbwt ]] ; then
    vg gbwt -p -g 1000GPlons_hs38d1_NA19239_sample_withref.force.augment.gg -o 1000GPlons_hs38d1_NA19239_sample_withref.force.augment.gbwt -x 1000GPlons_hs38d1_NA19239_sample_withref.xg -a 1000GPlons_hs38d1_NA19239_sample_withref.force.gbwt
fi


# (Re)Upload everything
aws s3 cp 1000GPlons_hs38d1_NA19239_sample_withref.xg s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_NA19239_sample_withref.xg
aws s3 cp 1000GPlons_hs38d1_NA19239_sample_withref.force.gbwt s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_NA19239_sample_withref.force.gbwt
aws s3 cp 1000GPlons_hs38d1_NA19239_sample_withref.force.augment.gbwt s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_NA19239_sample_withref.force.augment.gbwt
aws s3 cp 1000GPlons_hs38d1_NA19239_sample_withref.force.augment.gg s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_NA19239_sample_withref.force.augment.gg

            
              


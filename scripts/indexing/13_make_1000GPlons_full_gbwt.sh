#!/usr/bin/env bash
# Run with vg v1.32.0 on out 64 core/1TB memory machine, but this level of hardware is not expected to be necessary.

set -ex

# Download everything
if [[ ! -e 1000GPlons_hs38d1_filter.xg ]] ; then
    aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.xg .
fi

if [[ ! -e 1000GPlons_hs38d1_filter.gbwt ]] ; then
    aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.gbwt .
fi

# If we use a local haplotype cover we will introduce switch errors; augment instead.
if [[ ! -e 1000GPlons_hs38d1_filter.full.gg || ! -e 1000GPlons_hs38d1_filter.full.gbwt ]] ; then
    vg gbwt -p -g 1000GPlons_hs38d1_filter.full.gg -o 1000GPlons_hs38d1_filter.full.gbwt -x 1000GPlons_hs38d1_filter.xg -a 1000GPlons_hs38d1_filter.gbwt
fi

# Upload 
aws s3 cp 1000GPlons_hs38d1_filter.full.gbwt s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.full.gbwt
aws s3 cp 1000GPlons_hs38d1_filter.full.gg s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.full.gg

              


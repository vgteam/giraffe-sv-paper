#!/usr/bin/env bash
# Run with vg v1.32.0 on our 64 core/1TB memory machine with a few TB of free disk space

set -ex

export TMPDIR="$(pwd)/tmp"
mkdir -p "${TMPDIR}"

# Download everything
if [[ ! -e 1000GPlons_hs38d1_filter.xg ]] ; then
    aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.xg .
fi
if [[ ! -e 1000GPlons_hs38d1_filter.vg ]] ; then
    aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.vg .
fi
if [[ ! -e 1000GPlons_hs38d1_filter.gbwt ]] ; then
    aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.gbwt .
fi
if [[ ! -e 1000GPlons_hs38d1_filter.maxid.txt ]] ; then
    aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.maxid.txt . || true
fi
if [[ ! -e 1000GPlons_hs38d1_filter.gcsa ]] ; then
    aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.gcsa . || true
fi
if [[ ! -e 1000GPlons_hs38d1_filter.gcsa.lcp ]] ; then
    aws s3 cp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.gcsa.lcp . || true
fi

# Make indexes
if [[ ! -e 1000GPlons_hs38d1_filter.maxid.txt ]] ; then
    vg stats --node-id-range 1000GPlons_hs38d1_filter.vg | cut -f2 -d':' >1000GPlons_hs38d1_filter.maxid.txt
fi

if [[ ! -e 1000GPlons_hs38d1_filter.gcsa || ! -e 1000GPlons_hs38d1_filter.gcsa.lcp ]] ; then
    vg autoindex -p 1000GPlons_hs38d1_filter -R "GCSA" -R "LCP" -P "VG w/ Variant Paths:1000GPlons_hs38d1_filter.vg" -P "XG:1000GPlons_hs38d1_filter.xg" -P "MaxNodeID:1000GPlons_hs38d1_filter.maxid.txt" -P "GBWT:1000GPlons_hs38d1_filter.gbwt"
fi

# Upload. But index building will take a while, so don't hang around waiting for MFA if it is expired 
aws s3 cp --quiet 1000GPlons_hs38d1_filter.maxid.txt s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.maxid.txt &
aws s3 cp --quiet 1000GPlons_hs38d1_filter.gcsa s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.gcsa &
aws s3 cp --quiet 1000GPlons_hs38d1_filter.gcsa.lcp s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplons/hs38d1/1000GPlons_hs38d1_filter.gcsa.lcp &

wait
wait
wait


              


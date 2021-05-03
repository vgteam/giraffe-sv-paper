#!/usr/bin/env bash
# Run with vg v1.32.0 on our 64 core/1TB memory machine, but this level of hardware is not expected to be necessary.

set -ex

# Download everything
if [[ ! -e HGSVC_hs38d1.vg ]] ; then
    aws s3 cp s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.vg .
fi

# Index HGSVC full graph including alt allele path names
if [[ ! -e HGSVC_hs38d1.annotation.xg ]] ; then
    vg index -L -x HGSVC_hs38d1.annotation.xg HGSVC_hs38d1.vg
fi

# Upload everything
aws s3 cp HGSVC_hs38d1.annotation.xg s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.annotation.xg

              


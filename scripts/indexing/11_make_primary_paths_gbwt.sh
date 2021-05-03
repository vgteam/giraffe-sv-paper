#!/usr/bin/env bash
# Run on our 64 core/1TB memory machine, but this level of hardware is not expected to be necessary.

set -ex

# Need a "trash" directory to rerun
mkdir -p trash

aws s3 cp s3://vg-k8s/profiling/graphs/v2-2/generic/primary/hs38d1/primaryhs38d1.xg trash/

vg index -T -G trash/primaryhs38d1.paths.gbwt trash/primaryhs38d1.xg
vg gbwt -g trash/primaryhs38d1.paths.gg -x trash/primaryhs38d1.xg trash/primaryhs38d1.paths.gbwt
vg minimizer -t 16 -p -i trash/primaryhs38d1.paths.min -g trash/primaryhs38d1.paths.gbwt -G trash/primaryhs38d1.paths.gg

aws s3 cp trash/primaryhs38d1.paths.gbwt s3://vg-k8s/profiling/graphs/v2-2/generic/primary/hs38d1/primaryhs38d1.paths.gbwt
aws s3 cp trash/primaryhs38d1.paths.gg s3://vg-k8s/profiling/graphs/v2-2/generic/primary/hs38d1/primaryhs38d1.paths.gg
aws s3 cp trash/primaryhs38d1.paths.min s3://vg-k8s/profiling/graphs/v2-2/generic/primary/hs38d1/primaryhs38d1.paths.min

              


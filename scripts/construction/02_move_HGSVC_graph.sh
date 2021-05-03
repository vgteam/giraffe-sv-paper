for EXT in dist min xg vg snarls trivial.snarls gbwt ; do
    aws s3 cp  s3://vg-k8s/profiling/graphs/v2-2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.${EXT} s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.${EXT}
done

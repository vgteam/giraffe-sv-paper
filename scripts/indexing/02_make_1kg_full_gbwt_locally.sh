ITER=2
for GRAPH_BASE in s3://vg-k8s/profiling/graphs/v2/for-NA19239/1kg/hs37d5/1kg_hs37d5_filter ; do
    for GBWT_TYPE in full ; do
        aws s3 cp --no-progress ${GRAPH_BASE}.xg input.xg
        if [[ "${GBWT_TYPE}" == "cover" ]] ; then
            vg gbwt -p -g output.gg -o output.gbwt -x input.xg -P
        elif [[ "${GBWT_TYPE}" == "sampled" ]] ; then
            aws s3 cp --no-progress ${GRAPH_BASE}.gbwt input.gbwt
            vg gbwt -p -g output.gg -o output.gbwt -x input.xg -l input.gbwt
        elif [[ "${GBWT_TYPE}" == "full" ]] ; then
            aws s3 cp --no-progress ${GRAPH_BASE}.gbwt input.gbwt
            vg gbwt -p -g output.gg -o output.gbwt -x input.xg -a input.gbwt
        fi
        aws s3 cp --no-progress output.gbwt ${GRAPH_BASE}.${GBWT_TYPE}.gbwt
        aws s3 cp --no-progress output.gg ${GRAPH_BASE}.${GBWT_TYPE}.gg
    done
done

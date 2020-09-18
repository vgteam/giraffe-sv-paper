printf "graph\tgbwt\treads\tpairing\tcorrect\tmapq60\twrong_mapq60\tidentity\tscore\n" > report.tsv
printf "correct\tmq\tscore\taligner\n" > roc_stats.tsv

#Get graphs
aws s3 cp s3://vg-k8s/profiling/graphs/v2/for-NA19239/1kg/hs37d5/1kg_hs37d5_filter.gfa ./1kg.gfa
aws s3 cp s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.gfa ./hgsvc.gfa
aws s3 cp s3://vg-k8s/profiling/graphs/v2/for-NA19239/1kg/hs37d5/1kg_hs37d5_filter.xg ./1kg.xg
aws s3 cp s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.xg ./hgsvc.xg


for GRAPH in 1kg hgsvc ; do
    for READS in novaseq6000 hiseqxten hiseq2500 ; do

        if [[ ${GRAPH} == "1kg" ]] ; then
            aws s3 cp s3://vg-k8s/profiling/reads/sim/for-NA19239/1kg/hs37d5/${READS}/out_sim/sim.gam ./sim.gam
            aws s3 cp s3://vg-k8s/profiling/reads/sim/for-NA19239/1kg/hs37d5/${READS}/out_sim/sim.fq.gz ./sim.fq.gz
        elif [[ ${GRAPH} == "hgsvc" ]] ; then
            aws s3 cp s3://vg-k8s/profiling/reads/sim/for-NA19240/hgsvc/grch/${READS}/out_sim_gbwt/sim.gam ./sim.gam
            aws s3 cp s3://vg-k8s/profiling/reads/sim/for-NA19240/hgsvc/grch/${READS}/out_sim_gbwt/sim.fq.gz ./sim.fq.gz
        fi

        GraphAligner -g ${GRAPH}.gfa -f sim.fq.gz -a mapped.gam -x vg --seeds-mxm-cache-prefix ${GRAPH}_seeds

        
        ~/vg annotate -m -x ${GRAPH}.xg -a mapped.gam | ~/vg gamcompare -s -r 100 - sim.gam 2> count | ~/vg view -aj - > compared.json

        CORRECT_COUNT="$(grep correctly_mapped compared.json | wc -l)"
        IDENTITY="$(jq '.identity' compared.json | awk '{sum+=$1} END {print sum/NR}')"
        SCORE="$(sed -n '2p' count | sed 's/[^0-9\.]//g')"

        echo ${GRAPH} ${READS} ${SPEED} ${CORRECT_COUNT} ${SCORE}
        printf "${GRAPH}\tgraphaligner\t${READS}\t-\t-\t${CORRECT_COUNT}\t-\t-\t${IDENTITY}\t${SCORE}\n" >> report.tsv

        jq -r '(if .correctly_mapped then 1 else 0 end|tostring) + "," + (.mapping_quality|tostring) + "," + (.score|tostring)' compared.json | sed 's/,/\t/g' | sed "s/$/\tgraphaligner_${GRAPH}${READS}/" >> roc_stats.tsv
        
        
    done
done
aws s3 cp report.tsv s3://vg-k8s/users/xhchang/giraffe_experiments/report_graphaligner.tsv
aws s3 cp roc_stats.tsv s3://vg-k8s/users/xhchang/giraffe_experiments/roc_stats_graphaligner.tsv

#!/usr/bin/env bash

#accuracy_giraffe.sh: Run giraffe for all conditions and report the accuracy stats and roc

set -e

printf "graph\tgbwt\treads\tpairing\tspeed\tcorrect\tmapq60\twrong_mapq60\tidentity\tscore\n" > report.tsv
printf "correct\tmq\tscore\taligner\n" > roc_stats.tsv
for SPECIES in human yeast ; do
    case "${SPECIES}" in
    yeast)
        GRAPHS=(SK1 yeast_all yeast_subset)
        READSETS=(DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138)
        GBWTS=(raw)
        ;;
    human)
        GRAPHS=(hgsvc 1kg)
        READSETS=(novaseq6000 hiseqxten hiseq2500)
        GBWTS=(full sampled cover)
        ;;
    esac
    for GRAPH in ${GRAPHS[@]} ; do
        for READS in ${READSETS[@]} ; do
            case ${GRAPH} in
            1kg)
                GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/for-NA19239/1kg/hs37d5/1kg_hs37d5_filter
                aws s3 cp s3://vg-k8s/profiling/reads/sim/for-NA19239/1kg/hs37d5/${READS}/out_sim/sim.filt1M.gam ./sim.gam
                ;;
            hgsvc)
                GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1
                aws s3 cp s3://vg-k8s/profiling/reads/sim/for-NA19240/hgsvc/grch38/${READS}/out_sim_gbwt/sim.gam ./sim.gam
                ;;
            SK1)
                GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/primary/SK1/primarySK1
                aws s3 cp s3://vg-k8s/profiling/reads/sim/yeast/sim-${READS}.gam ./sim.gam
                ;;
            yeast_all)
                GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/cactus/yeast_all/yeast_all
                aws s3 cp s3://vg-k8s/profiling/reads/sim/yeast/sim-${READS}.gam ./sim.gam
                ;;
            yeast_subset)
                GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/generic/cactus/yeast_all/yeast_subset
                aws s3 cp s3://vg-k8s/profiling/reads/sim/yeast/sim-${READS}.gam ./sim.gam
                ;;
            esac
            aws s3 cp ${GRAPH_BASE}.xg ./graph.xg
            aws s3 cp ${GRAPH_BASE}.dist ./graph.dist
            for GBWT in ${GBWTS[@]} ; do
                rm graph.gbwt
                rm graph.gg
                rm graph.min
                if [[ "${GBWT}" == "raw" ]] ; then
                    aws s3 cp ${GRAPH_BASE}.gbwt ./graph.gbwt
                    aws s3 cp ${GRAPH_BASE}.gg ./graph.gg
                    aws s3 cp ${GRAPH_BASE}.min ./graph.min
                else
                    aws s3 cp ${GRAPH_BASE}.${GBWT}.gbwt ./graph.gbwt
                    aws s3 cp ${GRAPH_BASE}.${GBWT}.gg ./graph.gg
                    aws s3 cp ${GRAPH_BASE}.${GBWT}.min ./graph.min
                fi
                for PAIRING in single paired ; do 
                    if [[ ${PAIRING} == "paired" ]] ; then
                        PAIRED="-i"
                    elif [[ ${PAIRING} == "single" ]] ; then 
                        PAIRED=""
                    fi
                    vg giraffe -x graph.xg -H graph.gbwt -g graph.gg  -d graph.dist -G sim.gam ${PAIRED} -t 24 -p > mapped.gam 
                    vg gamcompare -r 100 -s <(vg annotate -m -x graph.xg -a mapped.gam) sim.gam 2>count | vg view -aj - > compared.json
                    CORRECT_COUNT="$(sed -n '1p' count | sed 's/[^0-9]//g')"
                    SCORE="$(sed -n '2p' count | sed 's/[^0-9\.]//g')"
                    MAPQ="$(grep mapping_quality\":\ 60 compared.json | wc -l)"
                    MAPQ60="$(grep -v correctly_mapped compared.json | grep mapping_quality\":\ 60 | wc -l)"
                    IDENTITY="$(jq '.identity' compared.json | awk '{sum+=$1} END {print sum/NR}')"
                    echo ${GRAPH} ${GBWT} ${READS} ${PAIRING} ${SPEED} ${CORRECT_COUNT} ${MAPQ} ${MAPQ60} ${IDENTITY} ${SCORE}
                    printf "${GRAPH}\t${GBWT}\t${READS}\t${PAIRING}\t${SPEED}\t${CORRECT_COUNT}\t${MAPQ}\t${MAPQ60}\t${IDENTITY}\t${SCORE}\n" >> report.tsv

                    jq -r '(if .correctly_mapped then 1 else 0 end|tostring) + "," + (.mapping_quality|tostring) + "," + (.score|tostring)' compared.json | sed 's/,/\t/g' | sed "s/$/\tgiraffe_${GRAPH}${GBWT}${READS}${PAIRING}/" >> roc_stats.tsv
                done
            done
        done
    done
done
sed -i 's/single//g ; s/paired/-pe/g ; s/null/0/g' roc_stats.tsv

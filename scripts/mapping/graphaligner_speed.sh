printf "graph\tgbwt\treads\ttotal_time\twall_clock_time\tmemory\n" > speed_report_graphaligner.tsv

#Get graphs
aws s3 cp s3://vg-k8s/profiling/graphs/v2/for-NA19239/1kg/hs37d5/1kg_hs37d5_filter.gfa ./1kg.gfa
aws s3 cp s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.gfa ./hgsvc.gfa

THREAD_COUNT=16

aws s3 cp s3://vg-k8s/profiling/reads/real/NA19239/novaseq6000-ERR3239454-shuffled-600m.fq.gz novaseq6000.fq.gz
for SPECIES in human ; do
    case "${SPECIES}" in
    yeast)
        GRAPHS=(yeast_subset)
        READSETS=(DBVPG6044 DBVPG6765 N44 UWOPS034614 UWOPS919171 Y12 YPS138)
        ;;
    human)
        GRAPHS=(hgsvc 1kg)
        READSETS=(novaseq6000)
        ;;
    esac
    for GRAPH in ${GRAPHS[@]} ; do
        for READS in ${READSETS[@]} ; do

            /usr/bin/time -v timeout -k1 2h bash -c "GraphAligner -g ${GRAPH}.gfa -f ${READS}.fq.gz -a mapped.gam -x vg -t ${THREAD_COUNT} --seeds-mxm-cache-prefix ${GRAPH}_seeds" 2> time-log.txt || true

            USER_TIME="$(cat "time-log.txt" | grep "User time" | sed 's/User\ time\ (seconds):\ \([0-9]*\.[0-9]*\)/\1/g')"
            SYS_TIME="$(cat "time-log.txt" | grep "System time" | sed 's/System\ time\ (seconds):\ \([0-9]*\.[0-9]*\)/\1/g')"
            TOTAL_TIME="$(echo "${USER_TIME} + ${SYS_TIME}" | bc -l)"
            CLOCK_TIME="$(cat "time-log.txt" | grep "Elapsed (wall clock) time" | sed 's/.*\ \([0-9,:]*\)/\1/g')"
            MEMORY="$(cat "time-log.txt" | grep "Maximum resident set" | sed 's/Maximum\ resident\ set\ size\ (kbytes):\ \([0-9]*\)/\1/g')"


            echo ${GRAPH} ${READS} ${TOTAL_TIME} ${CLOCK_TIME} ${MEMORY}
            printf "${GRAPH}\tGraphAligner\t${READS}\t${TOTAL_TIME}\t${CLOCK_TIME}\t${MEMORY}\n" >> speed_report_graphaligner.tsv

        done 
    done
done

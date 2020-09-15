#!/usr/bin/env bash

set -e

printf "graph\talgorithm\treads\tpairing\tcorrect\tmapq60\twrong_mapq60\tidentity\tscore\n" > report_minimap2.tsv
printf "graph\talgorithm\treads\tpairing\tcorrect\tmapq60\twrong_mapq60\tidentity\tscore\n" > report_bwa_mem.tsv
printf "graph\talgorithm\treads\tpairing\tcorrect\tmapq60\twrong_mapq60\tidentity\tscore\n" > report_bowtie2.tsv
printf "correct\tmq\tscore\taligner\n" > roc_stats_minimap2.tsv
printf "correct\tmq\tscore\taligner\n" > roc_stats_bwa_mem.tsv
printf "correct\tmq\tscore\taligner\n" > roc_stats_bowtie2.tsv

THREADS=16

#Get reference genomes
aws s3 cp  s3://vg-k8s/profiling/data/hs37d5.fa ./1kg.fa
aws s3 cp s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.gz ./hgsvc.fa.gz
gunzip hgsvc.fa.gz

#Fix chromosome names for hgsvc so that they match the graph
sed -i -r 's/chr([0-9]*|X|Y) (\s)/\1\2/g' hgsvc.fa

#Get xgsj
aws s3 cp s3://vg-k8s/profiling/graphs/v2/for-NA19239/1kg/hs37d5/1kg_hs37d5_filter.xg ./1kg.xg
aws s3 cp s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.xg ./hgsvc.xg


#Index genomes
bwa index 1kg.fa
bwa index hgsvc.fa

bowtie2-build --large-index 1kg.fa 1kg_bowtie2
bowtie2-build --large-index hgsvc.fa hgsvc_bowtie2

minimap2 -x sr -d 1kg.mmi 1kg.fa
minimap2 -x sr -d hgsvc.mmi hgsvc.fa


for GRAPH in hgsvc 1kg ; do
    for READS in novaseq6000 hiseqxten hiseq2500; do
        if [[ ${GRAPH} == "1kg" ]] ; then
            aws s3 cp s3://vg-k8s/profiling/reads/sim/for-NA19239/1kg/hs37d5/${READS}/out_sim/sim.gam ./sim.gam
            aws s3 cp s3://vg-k8s/profiling/reads/sim/for-NA19239/1kg/hs37d5/${READS}/out_sim/sim.fq.gz ./sim.fq.gz
        elif [[ ${GRAPH} == "hgsvc" ]] ; then
            aws s3 cp s3://vg-k8s/profiling/reads/sim/for-NA19240/hgsvc/grch38/${READS}/out_sim_gbwt/sim.gam ./sim.gam
            aws s3 cp s3://vg-k8s/profiling/reads/sim/for-NA19240/hgsvc/grch38/${READS}/out_sim_gbwt/sim.fq.gz ./sim.fq.gz
        fi
        gunzip sim.fq.gz
        sed 's/_1$//g' sim.fq | sed 's/_2$//g' > sim.paired.fq

        for ALGORITHM in minimap2 bwa_mem bowtie2 ; do
            for PAIRING in single paired ; do 
                if [[ ${ALGORITHM} == "minimap2" ]] ; then
                    if [[ ${PAIRING} == "paired" ]] ; then
                        minimap2 -ax sr --secondary=no -t ${THREADS} ${graph}.fa sim.paired.fq > mapped.bam
                    elif [[ ${PAIRING} == "single" ]] ; then 
                        minimap2 -ax sr --secondary=no -t ${THREADS} ${graph}.fa sim.fq > mapped.bam
                    fi
                elif [[ ${ALGORITHM} == "bwa_mem" ]] ; then
                    if [[ ${PAIRING} == "paired" ]] ; then
                        bwa mem -t ${THREADS} -p ${GRAPH}.fa sim.paired.fq > mapped.bam
                    elif [[ ${PAIRING} == "single" ]] ; then 
                        bwa mem -t ${THREADS} -p ${GRAPH}.fa sim.fq > mapped.bam
                    fi

                elif [[ ${ALGORITHM} == "bowtie2" ]] ; then
                    if [[ ${PAIRING} == "paired" ]] ; then
                        bowtie2 -t -p ${THREADS} -X 1065 ${GRAPH}_bowtie2 --interleaved sim.paired.fq > mapped.bam
                    elif [[ ${PAIRING} == "single" ]] ; then 
                        bowtie2 -t -p ${THREADS} ${GRAPH}_bowtie2 -U sim.paired.fq > mapped.bam
                    fi
                fi
                samtools view -F 2048 -b mapped.bam > mapped.primary.bam
                samtools view -f 2048 -b mapped.bam > mapped.secondary.bam

                vg inject -x ${graph}.xg mapped.primary.bam > mapped.primary.gam
                vg inject -x ${graph}.xg mapped.secondary.bam > mapped.secondary.gam
    
                if [[ ${PAIRING} == "paired" ]] ; then
                    vg view -aj mapped.primary.gam | sed 's/\/1/_1/g' | sed 's/\/2/_2/g' | vg view -aGJ - | vg annotate -m -x ${graph}.xg -a - | vg gamcompare -r 100 -s - sim.gam 2> count | vg view -aj - > compared.primary.json
                    vg view -aj mapped.secondary.gam | sed 's/\/1/_1/g' | sed 's/\/2/_2/g' | vg view -aGJ - | vg annotate -m -x ${graph}.xg -a - | vg gamcompare -r 100 - sim.gam| vg view -aj - > compared.secondary.json
                elif [[ ${PAIRING} == "single" ]] ; then 
                     vg annotate -m -x ${graph}.xg -a mapped.primary.gam | vg gamcompare -s -r 100 - sim.gam 2> count | vg view -aj - > compared.primary.json
                     vg annotate -m -x ${graph}.xg -a mapped.secondary.gam | vg gamcompare -r 100 - sim.gam | vg view -aj - > compared.secondary.json
                fi
                python ../combine_reads.py compared.primary.json compared.secondary.json compared.json
                sed -i '/^$/d' compared.json
    
                CORRECT_COUNT="$(grep correctly_mapped compared.json | wc -l)"
                SCORE="$(sed -n '2p' count | sed 's/[^0-9\.]//g')"
                MAPQ="$(grep mapping_quality\":\ 60 compared.json | wc -l)"
                MAPQ60="$(grep -v correctly_mapped compared.json | grep mapping_quality\":\ 60 | wc -l)"
                IDENTITY="$(jq '.identity' compared.json | awk '{sum+=$1} END {print sum/NR}')"
                echo ${GRAPH} ${READS} ${PAIRING} ${SPEED} ${CORRECT_COUNT} ${MAPQ} ${MAPQ60} ${SCORE}
                printf "${GRAPH}\t${ALGORITHM}\t${READS}\t${PAIRING}\t-\t${CORRECT_COUNT}\t${MAPQ}\t${MAPQ60}\t${IDENTITY}\t${SCORE}\n" >> report_${ALGORITHM}.tsv
    
                jq -r '(if .correctly_mapped then 1 else 0 end|tostring) + "," + (.mapping_quality|tostring) + "," + (.score|tostring)' compared.json | sed 's/,/\t/g' | sed "s/$/\t${ALGORITHM}_${GRAPH}${READS}${PAIRING}/" >> roc_stats_${ALGORITHM}.tsv

            done
        done
    done
done
sed -i 's/single//g ; s/paired/-pe/g ; s/null/0/g' roc_stats_minimap2.tsv
sed -i 's/single//g ; s/paired/-pe/g ; s/null/0/g' roc_stats_bwa_mem.tsv
sed -i 's/single//g ; s/paired/-pe/g ; s/null/0/g' roc_stats_bowtie2.tsv

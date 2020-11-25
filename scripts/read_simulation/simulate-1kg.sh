echo download real reads for training
if [ $SEQ = 'novaseq6000' ]
then
    aws s3 cp --quiet s3://vg-k8s/profiling/reads/real/NA19239/novaseq6000-ERR3239454-shuffled-1m.fq.gz .
    FASTQ=novaseq6000-ERR3239454-shuffled-1m.fq.gz
fi
if [ $SEQ = 'hiseq2500' ]
then
    aws s3 cp --quiet s3://vg-k8s/profiling/reads/real/NA19240/hiseq2500-ERR309934-shuffled-1m.fq.gz .
    FASTQ=hiseq2500-ERR309934-shuffled-1m.fq.gz
fi
if [ $SEQ = 'hiseqxten' ]
then
    aws s3 cp --quiet s3://vg-k8s/profiling/reads/real/NA19240/hiseqxten-SRR6691663-shuffled-1m.fq.gz .
    FASTQ=hiseqxten-SRR6691663-shuffled-1m.fq.gz
fi

echo download graph indexes
aws s3 cp --quiet s3://vg-k8s/profiling/graphs/v2/for-NA19239/1kg/hs37d5/1kg_hs37d5.full.gbwt .
aws s3 cp --quiet s3://vg-k8s/profiling/graphs/v2/for-NA19239/1kg/hs37d5/1kg_hs37d5.xg .
aws s3 cp --quiet s3://vg-k8s/profiling/graphs/v2/for-NA19239/1kg/hs37d5/1kg_hs37d5.vg .

echo simulate reads
vg sim -r -n $NREADS -a -s 12345 -p 570 -v 165 -i 0.00029 -x 1kg_hs37d5.xg -g 1kg_hs37d5.full.gbwt --sample-name NA19239 --ploidy-regex "hs37d5:0,NC_007605:0,X:1,Y:1,chrY_.*:1,chrEBV:0,.*:2" -F $FASTQ > sim.raw.gam

echo annotate with paths
vg annotate -p -x 1kg_hs37d5.vg -a sim.raw.gam > sim.gam

echo convert to fastq
vg view -X -a sim.gam | gzip > sim.fq.gz

echo format true position information
vg view -a sim.gam | jq -c -r '[.name] + if (.annotation.features | length) > 0 then [.annotation.features | join(",")] else ["."] end + if .refpos != null then [.refpos[] | .name, if .offset != null then .offset else 0 end] else [] end + [.score] + if .mapping_quality == null then [0] else [.mapping_quality] end | @tsv' > true.pos

echo save files
aws s3 cp sim.gam  s3://vg-k8s/profiling/reads/sim/for-NA19239/1kg/hs37d5/${SEQ}/out_sim_gbwt/
aws s3 cp sim.fq.gz  s3://vg-k8s/profiling/reads/sim/for-NA19239/1kg/hs37d5/${SEQ}/out_sim_gbwt/
aws s3 cp true.pos  s3://vg-k8s/profiling/reads/sim/for-NA19239/1kg/hs37d5/${SEQ}/out_sim_gbwt/

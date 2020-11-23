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

echo download indexes
aws s3 cp --quiet s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.gbwt .
aws s3 cp --quiet s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1.xg .

echo simulate and annotate
vg sim -r -n $NREADS -a -s 12345 -p 570 -v 165 -i 0.00029 -x HGSVC_hs38d1.xg -g HGSVC_hs38d1.gbwt --sample-name NA19240 --ploidy-regex "JTFH.*:0,KN70.*:0,Y:0,chrY_.*:0,chrEBV:0,.*:2" -F $FASTQ | vg annotate -p -x HGSVC_hs38d1.xg -a - > sim.gam

echo convert to fastq
vg view -X -a sim.gam | gzip > sim.fq.gz

echo format true position information
vg view -a sim.gam | jq -c -r '[.name] + if (.annotation.features | length) > 0 then [.annotation.features | join(",")] else ["."] end + if .refpos != null then [.refpos[] | .name, if .offset != null then .offset else 0 end] else [] end + [.score] + if .mapping_quality == null then [0] else [.mapping_quality] end | @tsv' > true.pos

echo save files
aws s3 cp sim.gam  s3://vg-k8s/profiling/reads/sim/for-NA19240/hgsvc/grch38/${SEQ}/out_sim_gbwt/
aws s3 cp sim.fq.gz  s3://vg-k8s/profiling/reads/sim/for-NA19240/hgsvc/grch38/${SEQ}/out_sim_gbwt/
aws s3 cp true.pos  s3://vg-k8s/profiling/reads/sim/for-NA19240/hgsvc/grch38/${SEQ}/out_sim_gbwt/


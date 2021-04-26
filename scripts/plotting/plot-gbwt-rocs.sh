#!/usr/bin/env bash
#Get the files to make gbwt roc plots

#aws s3 cp --no-progress --recursive --exclude "*" --include "roc_stats_giraffe_gbwts_*novaseq6000.tsv"  s3://vg-k8s/users/xhchang/giraffe_experiments/ .

printf "correct\tmq\tscore\taligner\n" > roc-stats-1000gp-novaseq6000-sampled-single.tsv
printf "correct\tmq\tscore\taligner\n" > roc-stats-1000gp-novaseq6000-sampled-paired.tsv
printf "correct\tmq\tscore\taligner\n" > roc-stats-hgsvc-novaseq6000-cover-single.tsv
printf "correct\tmq\tscore\taligner\n" > roc-stats-hgsvc-novaseq6000-cover-paired.tsv

#grep hgsvcfullnovaseq6000 ../roc_stats_giraffe.tsv | grep -v pe | sed 's/hgsvcfullnovaseq6000/full/g' >> roc-stats-hgsvc-novaseq6000-cover-single.tsv
#grep hgsvcfullnovaseq6000-pe ../roc_stats_giraffe.tsv |  sed 's/hgsvcfullnovaseq6000-pe/full/g' >> roc-stats-hgsvc-novaseq6000-cover-paired.tsv
#grep hg38d1covernovaseq6000 ../roc_stats_giraffe_primary.tsv | grep -v pe | sed 's/hs38d1covernovaseq6000/primary/g' >> roc-stats-hgsvc-novaseq6000-cover-single.tsv
#grep hg38d1covernovaseq6000-pe ../roc_stats_giraffe_primary.tsv |  sed 's/hs38d1covernovaseq6000-pe/primary/g' >> roc-stats-hgsvc-novaseq6000-cover-paired.tsv
grep cover roc_stats_giraffe_gbwts_hgsvc_novaseq6000.tsv | grep -v pe |  sed 's/giraffe_hgsvccover\([0-9]*\)novaseq6000/cover_00\1/g ; s/0016/016/g ; s/0032/032/g ; s/0064/064/g ; s/00128/128/g' >> roc-stats-hgsvc-novaseq6000-cover-single.tsv
grep cover roc_stats_giraffe_gbwts_hgsvc_novaseq6000.tsv | grep pe |  sed 's/giraffe_hgsvccover\([0-9]*\)novaseq6000-pe/cover_00\1/g ; s/0016/016/g ; s/0032/032/g ; s/0064/064/g ; s/00128/128/g' >> roc-stats-hgsvc-novaseq6000-cover-paired.tsv


#grep 1kgfullnovaseq6000 ../roc_stats_giraffe.tsv | grep -v pe | sed 's/1kgfullnovaseq6000/full/g' >> roc-stats-1kg-novaseq6000-sampled-single.tsv
#grep 1kgfullnovaseq6000-pe ../roc_stats_giraffe.tsv |  sed 's/1kgfullnovaseq6000-pe/full/g' >> roc-stats-1kg-novaseq6000-sampled-paired.tsv
grep 1000gpcovernovaseq6000 ../roc_stats_giraffe_primary.tsv | grep -v pe | sed 's/1000gpcovernovaseq6000/primary/g' >> roc-stats-1000gp-novaseq6000-sampled-single.tsv

grep 1000gpcovernovaseq6000-pe ../roc_stats_giraffe_primary.tsv |  sed 's/1000gpcovernovaseq6000-pe/primary/g' >> roc-stats-1000gp-novaseq6000-sampled-paired.tsv

grep sampled roc_stats_giraffe_gbwts_1000gp_novaseq6000.tsv | grep -v pe |  sed 's/giraffe_1000gpsampled\([0-9]*\)novaseq6000/sampled_00\1/g ; s/0016/016/g ; s/0032/032/g ; s/0064/064/g ; s/00128/128/g' >> roc-stats-1000gp-novaseq6000-sampled-single.tsv

grep sampled roc_stats_giraffe_gbwts_1000gp_novaseq6000.tsv | grep pe |  sed 's/giraffe_1000gpsampled\([0-9]*\)novaseq6000-pe/sampled_00\1/g ; s/0016/016/g ; s/0032/032/g ; s/0064/064/g ; s/00128/128/g' >> roc-stats-1000gp-novaseq6000-sampled-paired.tsv

./plot-roc-gbwts.R roc-stats-1000gp-novaseq6000-sampled-single.tsv plot-roc-1000gp-novaseq6000-sampled-single.svg 
./plot-roc-gbwts.R roc-stats-1000gp-novaseq6000-sampled-paired.tsv plot-roc-1000gp-novaseq6000-sampled-paired.svg
./plot-roc-gbwts.R roc-stats-hgsvc-novaseq6000-cover-single.tsv plot-roc-hgsvc-novaseq6000-cover-single.svg
./plot-roc-gbwts.R roc-stats-hgsvc-novaseq6000-cover-paired.tsv plot-roc-hgsvc-novaseq6000-cover-paired.svg

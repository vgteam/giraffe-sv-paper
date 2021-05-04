#!/usr/bin/env bash
# Make the reference FASTA for GRCh38 with decoys

rm -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz GCA_000786075.2_hs38d1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz            
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz
zcat GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz GCA_000786075.2_hs38d1_genomic.fna.gz > GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna

cat GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna | sed -r -e 's/>chr([0-9XYMT]+) />\1 /g' >GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.renamed.fna

gzip GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna
aws s3 cp GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.gz s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.gz

pbgzip GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.renamed.fna
aws s3 cp GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.renamed.fna.gz s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.renamed.fna.gz
            


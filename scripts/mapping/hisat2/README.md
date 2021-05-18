# HISAT2 

These scripts are used for creating indexes, mapping reads and evaluating alignments for HISAT2

To run them, you will need a machine with the `aws` command configured.

You will also need `kubectl` set up to access a Kubernetes cluster which has a service account named `vg-svc` and a secret named `shared-s3-credentials` with the contents of a `~/.aws` directory that grants single-factor access to AWS. You may need to adjust the scripts to reference a different service account or secret name if you cannot provide resources with these names. The Kubernetes cluster will need to be able to fulfil requests for up to 250 Gi of memory, and 300 Gi of ephemeral storage.

These scripts expect to write to the `s3://vg-k8s` bucket. You will need to have write access to this bucket, or you will need to replace these writes with writes to a bucket where you intend to store your generated artifacts.

All the scripts include different environmental variables (e.g. `REF` or `READS`) that need to be set before running. The number of threads are set using `CPU`.

## Running the scripts

1. Parse variants and haplotypes for each chromosome (`CHR`) using either `00_hisat2_prepare_variants_1kg.sh` or `00_hisat2_prepare_variants_hgsvc.sh` (dependent on the variant set). 
   * Docker image: `jsibbesen/hisat2-s3script:2.2.0-s2`
2. Create an index of the reference using `01_hisat2_generate_index.sh`. 
   * Docker image: `jsibbesen/hisat2-s3script:2.2.0-s2`
3. Map reads using either `02_hisat2_map_reads_1kg.sh` or `02_hisat2_map_reads_hgsvc.sh` (dependent on the variant set).
   * Docker image: `quay.io/jsibbesen/hisat2-s3script:hisat2-2.2.1-s2`
4. Calculate mapping accuracy using either `03_hisat2_calc_map_stats_1kg.sh` or `03_hisat2_calc_map_stats_hgsvc.sh` (dependent on the variant set).
   * Docker image: `quay.io/jsibbesen/vg-s3script:vg-1.31.0-s1`

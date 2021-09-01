# Genotyping Scripts

`real_read_deepvariant_call.sh` runs ABRA2-based indel reanlignment and genotypes with DeepVariant v1.1.0.

To run it, you will need Docker installed and runnable from command line.
The script will write to the directory `${HOME}/run_deepvariant_genotyping` and expects as input a BAM file.

`real_read_dragen_call.sh` runs ABRA2-based indel reanlignment and genotypes with Illumina's DRAGEN module.

To run it, you will need to run it on a server that houses an Illumina DRAGEN module and have Docker installed and runnable from command line.
The script will write to the directory `${HOME}/run_dragen_genotyping` and expects as input a BAM file.

`graphtyper_kubernetes` runs genotyping with GraphTyper.

To run it, you will need `kubectl` set up to access a Kubernetes cluster which has a secret named `shared-s3-credentials` with the contents of a `~/.aws` directory that grants single-factor access to AWS. You may need to adjust the scripts to reference a different secret name if you cannot provide resources with these names. The Kubernetes cluster will need to be able to fulfil requests for up to 26 cores, 200 Gi of memory, and 200 Gi of ephemeral storage.

This script expects to write to the `s3://vg-k8s` bucket. You will need to have write access to this bucket, or you will need to replace these writes with writes to a bucket where you intend to store your generated artifacts.

## Running the scripts

Run the deepvariant script:

```
./real_read_deepvariant_call.sh ${INPUT_BAM_FILE}
```

Run the DRAGEN script:

```
./real_read_dragen_call.sh ${INPUT_BAM_FILE}
```

Run the graphtyper script, and wait for the Kubernetes jobs it generates to complete:

```
./graphtyper_kubernetes
```


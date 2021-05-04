# Genotyping Scripts

There is only one script here; it rund genotyping with GraphTyper.

To run it, you will need `kubectl` set up to access a Kubernetes cluster which has a secret named `shared-s3-credentials` with the contents of a `~/.aws` directory that grants single-factor access to AWS. You may need to adjust the scripts to reference a different secret name if you cannot provide resources with these names. The Kubernetes cluster will need to be able to fulfil requests for up to 26 cores, 200 Gi of memory, and 200 Gi of ephemeral storage.

This script expects to write to the `s3://vg-k8s` bucket. You will need to have write access to this bucket, or you will need to replace these writes with writes to a bucket where you intend to store your generated artifacts.

## Running the scripts

Run the script, and wait for the Kubernetes jobs it generates to complete:

```
./graphtyper_kubernetes
```


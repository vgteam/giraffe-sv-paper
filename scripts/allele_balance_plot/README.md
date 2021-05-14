# Allele Balance Plot Scripts

These run Giraffe, VG-MAP, and BWA-MEM and create an allele balance plot to assess mapping reference allele bias. 

To run them, you will need a machine with the `aws` command configured.

You will also need `kubectl` set up to access a Kubernetes cluster which has a service account named `vg-svc` and a secret named `shared-s3-credentials` with the contents of a `~/.aws` directory that grants single-factor access to AWS. You may need to adjust the scripts to reference a different service account or secret name if you cannot provide resources with these names. The Kubernetes cluster will need to be able to fulfil requests for up to 24 cores, 250 Gi of memory, and 300 Gi of ephemeral storage.

These scripts expect to write to the `s3://vg-k8s` bucket. You will need to have write access to this bucket, or you will need to replace these writes with writes to a bucket where you intend to store your generated artifacts.


## Running the scripts

The Kubernetes script can be run:

```
./kubernetes_pileup_call
```

This script will produce a vcf file `all_calls.vcf.gz` that is copied to an s3 bucket.

The plotting script can then be run on the unzipped vcf file:

```python plot_allele_balance.py all_calls.vcf plot.svg
```

To produce the allele balance plot in `plot.svg`

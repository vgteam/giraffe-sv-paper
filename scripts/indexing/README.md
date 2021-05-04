# Graph Indexing Scripts

These index graphs for evaluating Giraffe's read mapping ability, and for running some competing algorithms.

To run them, you will need a machine with the `aws` command configured, and a compatible version of [Toil](https://github.com/DataBiosphere/toil) installed.

You will also need `kubectl` set up to access a Kubernetes cluster which has a service account named `vg-svc` and a secret named `shared-s3-credentials` with the contents of a `~/.aws` directory that grants single-factor access to AWS. You may need to adjust the scripts to reference a different service account or secret name if you cannot provide resources with these names. The Kubernetes cluster will need to be able to fulfil requests for up to 32 cores, 350 Gi of memory, and 2200 Gi of ephemeral storage.

These scripts expect to write to the `s3://vg-k8s` bucket. You will need to have write access to this bucket, or you will need to replace these writes with writes to a bucket where you intend to store your generated artifacts.

## Running the scripts

First, run the first script, and wait for the Kubernetes jobs it starts to finish:

```
./01_make_HGSVC_processed_gbwts.sh
```

Then, run the second script, and wait for the Kubernetes jobs it starts to finish:

```
./02_make_HGSVC_processed_gbwt_indexes.sh
```

Proceed similarly with the next few scripts:

```
./03_make_HGSVC_alt_xg.sh
```

```
./04_prune_HGSVC.sh
```

```
./05_make_HGSVC_gcsa.sh
```

```
./06_make_1000GPlons_GRCh38_sampled_gbwts.sh
```

```
./07_make_1000GPlons_GRCh38_cover_gbwts.sh
```

Then, for the next script, make sure you have `vg` v1.32.0 (or a compatible version you wish to use instead) installed and on your `$PATH` as `vg`. Then run the next three scripts:

```
./08_make_positive_control_gbwts.sh
./09_make_1000GPlons_full_gbwt.sh
./10_make_HGSVC_annotation_xg.sh
```

The next script was not originally run with `vg` v1.32.0, but it is expected to work with it, so run it too:

```
./11_make_primary_paths_gbwt.sh
```

Finally, the last script uses Kubernetes again, so run it and wait for the jobs it creates to finish:

```
./12_make_primary_cover_gbwt.sh
```





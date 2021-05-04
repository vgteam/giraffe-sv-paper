# Graph Construction Scripts

These construct graphs for evaluating Giraffe read mapping.

To run them, you will need a machine with the `aws` command configured, and a compatible version of [Toil](https://github.com/DataBiosphere/toil) installed.

You will also need `kubectl` set up to access a Kubernetes cluster which has a service account named `vg-svc` and a secret named `shared-s3-credentials` with the contents of a `~/.aws` directory that grants single-factor access to AWS. You may need to adjust the scripts to reference a different service account or secret name if you cannot provide resources with these names. The Kubernetes cluster will need to be able to fulfil requests for up to 32 cores, 150 Gi of memory, and 2200 Gi of ephemeral storage.

These scripts expect to read files from the `s3://vg-data` and `s3://glennhickey` buckets. You will need to have access to this bucket, or else replace reads from it with reads from a different bucket where you have placed the archived input data from our data archive.

These scripts expect to read files from the mount point `/public/groups/cgl/users/daheller/yeast_graph`. You will need to provide this mount point, with the relevant files from the data archive, or else edit the scripts to reference a different local path for these files. 

These scripts expect to write to the `s3://vg-k8s` bucket. You will need to have write access to this bucket, or you will need to replace these writes with writes to a bucket where you intend to store your generated artifacts.

## Running the scripts

First, create a Toil Mesos cluster:

```
TOIL_APPLIANCE_SELF=quay.io/ucsc_cgl/toil:4.2.0a1-edb0e583e91452e80ef6add0f9f0e8eae9dbc2d4-py3.7 toil launch-cluster -z us-west-2a adamnovakgraphbuild --leaderNodeType t2.medium --keyPairName anovak@kolossus
```

If you do not have an SSH key set up on AWS named `anovak@kolossus`, for which you have the private key loaded into your SSH agent or decrypted on disk, you will need to replace `anovak@kolossus` with the name of such a key in your AWS account.

Then, connect to the cluster:

```
toil ssh-cluster -z us-west-2a adamnovakgraphbuild
```

Then, start a screen session:

```
screen
```

After transferring the first script over to the machine, run it:

```
./01_make_HGSVC_graph.sh
```

Note that the script is expected to try and fail to generate GCSA indexes; these will be generated later. We have better methods for building graphs with GCSA indexes that do not fail (see for example `05_make_1000GPlons_GRCh38_graph_with_gcsa.sh`). This first script reproduces the method that was actually used to generate the HGSVC graphs in the paper, not the method we recommend for use in new projects.

After running the script, disconnect from the cluster and destroy it:

```
toil destroy-cluster -z us-west-2a adamnovakgraphbuild
```

You also will want to destroy the job store, holding the resumeable workflow:

```
toil clean aws:us-west-2:adamnovak-make-hgsvc-graphs-2
```

Then, run the second script to move the generated graphs into position for the rest of the analysis:

```
./02_move_HGSVC_graph.sh
```

Then, run the third script, which will launch a Kubernetes job, and wait for the job's successful completion:

```
./03_make_primary_graph.sh
```

Then, create a Toil Kubernetes cluster. (We do **not** recommend using this version of Toil for new work; the run will be much more efficient if you use the cluster configuration as for `05_make_1000GPlons_GRCh38_graph_with_gcsa.sh` below. However, this is the version that was actually used to generate the graphs used in the paper).

```
TOIL_APPLIANCE_SELF=quay.io/ucsc_cgl/toil:5.4.0a1-1530a9190357fc058333f3e929049ef9593a6784-py3.7 toil launch-cluster --provisioner aws -T kubernetes -z us-west-2a adamnovak-toil-vg --leaderNodeType t3a.medium --nodeTypes=t3a.medium,r5ad.24xlarge,r5d.24xlarge/r5ad.24xlarge:2.50,i3.8xlarge:1.50 --workers 1-4,0-1,0-8,0-6 --keyPairName anovak@soe.ucsc.edu
```

If you do not have an SSH key set up on AWS named `anovak@soe.ucsc.edu`, for which you have the private key loaded into your SSH agent or decrypted on disk, you will need to replace `anovak@soe.ucsc.edu` with the name of such a key in your AWS account.

Then, connect to the cluster:

```
toil ssh-cluster -z us-west-2a adamnovak-toil-vg
```

Then, start a screen session:

```
screen
```

Then, after transferring it to the cluster, run the next script:

```
./04_make_1000GPlons_GRCh38_graph.sh
```

When it completes, leave the cluster.

Then, tear down the cluster:
```
toil destroy-cluster -z us-west-2a adamnovak-toil-vg
```

Re-make it with a slightly different (improved) version of Toil (you will also want to have Toil 46f0b22c4edf4dfa127ebcb2b0025dd286b5b874 or a compatible release installed for this step):
```
TOIL_APPLIANCE_SELF=quay.io/ucsc_cgl/toil:5.4.0a1-46f0b22c4edf4dfa127ebcb2b0025dd286b5b874-py3.7 toil launch-cluster --provisioner aws -T kubernetes -z us-west-2a adamnovak-toil-vg --leaderNodeType t3a.medium --nodeTypes=t3a.medium,r5ad.24xlarge,r5d.24xlarge/r5ad.24xlarge:2.50,i3.8xlarge:1.50 --workers 1-4,0-1,0-8,0-10 --keyPairName anovak@soe.ucsc.edu
```

Then, connect to the cluster:

```
toil ssh-cluster -z us-west-2a adamnovak-toil-vg
```

Then, start a screen session:

```
screen
```

Then, after transferring it to the cluster, run the next script:

```
./05_make_1000GPlons_GRCh38_graph_with_gcsa.sh
```

After that finishes, disconnect from the cluster, and then tear it down:
```
toil destroy-cluster -z us-west-2a adamnovak-toil-vg
```

Finally, run the last script locally:
```
./06_yeast_graphs_and_reads.sh
```







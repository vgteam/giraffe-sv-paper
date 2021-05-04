# Graph Indexing Scripts

These run Giraffe and competing mappers, so that speed and accuracy can be assessed.

To run them, you will need a machine with the `aws` command configured.

You will also need `kubectl` set up to access a Kubernetes cluster which has a service account named `vg-svc` and a secret named `shared-s3-credentials` with the contents of a `~/.aws` directory that grants single-factor access to AWS. You may need to adjust the scripts to reference a different service account or secret name if you cannot provide resources with these names. The Kubernetes cluster will need to be able to fulfil requests for up to 24 cores, 250 Gi of memory, and 300 Gi of ephemeral storage.

These scripts expect to write to the `s3://vg-k8s` bucket. You will need to have write access to this bucket, or you will need to replace these writes with writes to a bucket where you intend to store your generated artifacts.

The yeast graph experiment script expect to read files from the mount point `/public/groups/cgl/users/daheller/yeast_graph`. You will need to provide this mount point, with the relevant files from the data archive, or else edit the script to reference a different local path for these files. 

## Running the scripts

All the Kubernetes scripts can be run at the same time:

```
./*_kubernetes
```

Results will be available when the jobs they create finish.

The speed scripts expect to run on standard-sized AWS instances, as described in the paper, so results are comparable between them. To run them, transfer each to an AWS instance of that type with `aws` and the versions of `vg` and the other mappers to be tested installed, and with access to the appropriate output data bucket, and run them from a directory with sufficient fast storage.

These scripts test `map` and `giraffe` in `vg`:
```
./map_speed.sh
```
```
./giraffe_speed.sh
```

This script tests the speed of [GraphAligner](https://github.com/maickrau/GraphAligner):
```
./graphaligner_speed.sh
```

This script tests the speed of [hisat2](https://daehwankimlab.github.io/hisat2/):
```
./hisat2_speed.sh
```

This script tests the linear genome mappers [bowtie2](https://github.com/BenLangmead/bowtie2#readme), [minimap2](https://github.com/lh3/minimap2#readme), and [bwa mem](https://github.com/lh3/bwa#readme), all of which must be installed:
```
./linear_mappers_speed.sh
```

The yeast mapping experiment script doesn't depend on any of the others. If run after the yesst graph build script, it will re-use the constructed graphs stored in local files, but if not it will construct them itself. It is run like this:

```
./giraffe_yeast_experiment.sh
```

None of the `.py` scripts need to be run; they are used by other scripts or not at all.

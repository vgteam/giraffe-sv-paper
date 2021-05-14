# Graph Construction Scripts

These save code and data, and upload them to Zenodo.

To run them, you will need a machine with the `aws` command configured, as well as access to `docker`.

You will also need to be able to read from these AWS S3 buckets where the data to be archived was produced:

* `s3://vg-k8s`
* `s3://vg-data`
* `s3://glennhickey`

You will need the `.tex` files for the paper checked out from the Overleaf git repo into `$HOME/build/giraffe-paper/`.

You will also need a local mountpoint in which to deposit the archive, at `/nanopore/cgl/data/giraffe`.

To actually make the archive available, you will need `ZENODO_DEPOSITION` and `ZENODO_TOKEN` environment variables populated with a Zenodo deposition ID to upload the software archive and product files to, and a Zenodo access token to authenticate the upload. Note that Zenodo does not reliably support in-place update; exitsing files will need to be deleted on the Zenodo side. You will also need `ipfs` version 0.8.0 to hash and share the data.

These scripts aren't really meant to be run by anyone but the authors, but they
serve as an explaination of how the archive was assembled.

Information about the resulting archive is available in the [overall archive readme](archive-readme.md) and the [software archive readme](software-readme.md).

## Running the scripts

Run each script in turn locally:

```
./save_construction_inputs.sh
./save_mapping_inputs.sh
./save_products.sh
./archive_code_and_containers.sh
```

After running all the scripts, you will want to add the data to IPFS to get a hash to pin and publish, and for integrity chacking (as no integrity checking is done by the scripts). Note that it must be added through a symlink in your IPFS root, like this:

```
ipfs add -r --nocopy /public/home/anovak/projects/gg/data/cgl-data/giraffe
```

You will also went to publish the Zenodo deposition you uploaded the software archive to, so it becomes visible under its assigned DOI.

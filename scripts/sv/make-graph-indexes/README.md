Indexes made using a version of the [vg snakemake pipeline](https://github.com/vgteam/vg_snakemake).

```sh
CPU=8
MEMORY=200

snakemake --configfile snakemake_config.yaml --config mapper=giraffe --cores $CPU --resources mem_mb=${MEMORY}000 -p main
```

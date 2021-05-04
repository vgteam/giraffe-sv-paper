# Giraffe Mapper Evaluation and Application Scripts

This repository contains scripts used to reproduce our work with the new Giraffe short read mapper in [vg](https://github.com/vgteam/vg).

The scripts expect to be run in roughly this order:

* Giraffe mapping evaluation workflow
    * [Preparation scripts to preprocess input files](scripts/prep)
    * [Graph construction scripts, to make test graphs for Giraffe to map to](scripts/construction)
    * [Indexing scripts to prepare the constructed graphs for mapping](scripts/indexing)
    * [Read simulation scripts to produce simulated reads with known graph positions, for evaluating mapping correctness](scripts/read_simulation)
    * [Mapping scripts, for assessing the speed and accuracy of Giraffe against competing mappers](scripts/mapping)
    * [Genotyping scripts, for assessing competing genotyping methods](scripts/genotyping)
    * [Plotting scripts, for plotting the results of the mapping scripts](scripts/plotting)
    * [Dedicated allele balance plotting scripts, for producing plots of how length-changing variants affect read coverage at variable sites](scripts/allele_balance_plot)
* [Structural variant calling workflow](scripts/sv)
* [Code and data archiving](scripts/archiving)

Note that the top level workflows are not automated. Within each section, you will have to manually prepare the environment for and run each script. Some scripts expect to run locally with `vg` or `snakemake` installed and sufficient memory and scratch space, some scripts expect to run with access to a Kubernetes cluster, and some scripts expect to be launched on a Toil-managed autoscaling Mesos or Kubernetes cluster. We provide hints as to how to set up such environments, but a full tutorial is not given here. Additionally, scripts that launch asynchronous Kubernetes jobs do not include code to wait for the jobs to complete; that monitoring must be provided by you.

We provide scripts as close to what we actually ran as possible; these scripts will not be fully portable to your environment without modification. If you do not have access to UCSC's AWS storage buckets (such as `s3://vg-k8s` or `s3://vg-data`), or if you would like to avoid overwriting the original analysis artifacts, some scripts will have to be adapted to point at where you intend to keep your artifacts for your repetition of the analysis. Additionally, scripts designed to kick off Kubernetes jobs may need to be adapted to reference your Kubernetes environment's AWS credential secrets or namespace names.

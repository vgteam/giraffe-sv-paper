# Prep Scripts

These scripts prepare input data for graph construction. The prepared data is already available in our paper's data archive, but you might want to regenerate it.

To run them, you will need a machine with the `aws` command configured, as well as standard bioinformatics tools like `pbgzip` and `tabix`.

These scripts expect to read files from the `s3://vg-data` bucket. You will need to have access to this bucket, or else replace reads from it with reads from a different bucket where you have placed the archived input data from our data archive.

These scripts expect to write to the `s3://vg-data` and `s3://vg-k8s` buckets. You will need to have write access to these buckets, or you will need to replace these writes with writes to buckets where you intend to store your generated artifacts.

## Running the scripts

First, make the combined reference files including GRCh38 and decoy sequences:

```
./01_make_reference.sh
```

Then, subset the 1000 Genomes Project liftover VCFs to remove varinats in segmental duplications:

```
./02_preprocess_1000GPlons_VCFs.sh
```

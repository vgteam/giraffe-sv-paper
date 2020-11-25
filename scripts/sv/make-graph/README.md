## Prepare VCFs for hybrid construction

Remove near-duplicates using re-mapping output and separate variants for augmentation.

```sh
Rscript prepareVcf.R
```

## Make .vg chunks and augmentation using mpmap if necessary

```sh
python3 prepareAugmentFiles.py -r hg38.fa -v hsvlr_srdedup17_foraug.vcf -s SVseqs-foraugment.fa -f 5000
```

## Run iterative augmentation on each chunk and concatenate back xto chromosome-level .vg files

```sh
snakemake --config chunk_bed=hsvlr_srdedup17_augment_chunks.bed --cores $CPU -k
```

## Construct graphs for unlocalized, unplaced and decoy contigs

```sh
for CHR in `seq 1 22` X Y M
do
    echo chr$CHR >> chrs.txt
done
cut -f1 hg38.fa.fai | grep -v "_alt" | sort > noalts.txt
sort chrs.txt | comm -23 noalts.txt - > undecoy.txt
for CHR in `cat undecoy.txt`
do
    samtools faidx hg38.fa $CHR >> hg38-undecoy.fa
done
vg construct -r hg38-undecoy.fa > hg38-hsvlr_srdedup17_aug-undecoy.vg
```

## hs38d1 decoy sequences

```sh
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz
gunzip GCA_000786075.2_hs38d1_genomic.fna.gz
vg construct -r GCA_000786075.2_hs38d1_genomic.fna > hg38-hsvlr_srdedup17_aug-hs38d1.vg
```

## Construct graphs for chrM

```sh
vg construct -r hg38.fa -R chrM -C > hg38-hsvlr_srdedup17_aug-chrM.vg
aws s3 cp hg38-hsvlr_srdedup17_aug-chrM.vg $SROOT/
```

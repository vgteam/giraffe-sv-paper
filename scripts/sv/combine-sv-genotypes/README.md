## Combined VCF across samples

A combined VCF was created using a [WDL workflow deposited on Dockstore](https://dockstore.org/workflows/github.com/jmonlong/wdl-workflows/bcftools_merge:giraffe-paper?tab=info).
Briefly, it uses bcftools and custom scripts to split multi-allelic variants, extract canonical SVs from vg calls, potentially merge heterozygous variants into homozygous variants, left-align and sort each VCF.
The VCF are then merged using `bcftools merge`.

The `align_variants.py` and `merge_exact_hets.py` are used by the WDL and are also deposited here for record keeping.
The docker container used by the WDL, and containing these scripts, was built from [this repository](https://github.com/jmonlong/docker-merge-sv-vg). 

Briefly, each sample is processed with the following commands:

```sh
bcftools norm -m -both svs.vcf -O z > svs.split.vcf.gz
python3 align_variants.py -i svs.split.vcf.gz -f reference.fa -o svs.split.aligned.vcf
bcftools view --exclude 'GT="0/0" || GT="0" || GT~"\\."' svs.split.aligned.vcf | bcftools norm -f reference.fa | bcftools sort | python3 merge_exact_hets.py | bgzip > svs.split.aligned.merged.sorted.vcf.gz
```

Of note, versions of some modules/python have changed since the publication. 
To use these scripts, you might need to install specific versions of the modules, for example: `pip3 install biopython==1.77 pyvcf3`.
Otherwise, you could use the docker image used by the WDL (`jmonlong/merge-sv-vg@sha256:5536c9077f18457a9895125d38708b8286f6b8c62b9a5a4c1d6fbd140dafd803`).
Second note, [vcfwave](https://github.com/vcflib/vcflib/blob/master/doc/vcfwave.md) aims at doing the same as the `align_variants.py` script and could be a good replacement in case of issues.

Once each single-sample VCF has been decomposed, they are merged with 

```sh
bcftools merge -0 -m none -l vcf_list.txt -O z > merged.vcf.bgz
```

## Define SV sites

Information for each SV allele was read from the combined VCF using:

```r
svs = readSVvcf.multisamps('merged.vcf.bgz', keep.ids=TRUE, min.sv.size=30,
                           keep.ref.seq=TRUE, keep.ins.seq=TRUE)
outf = gzfile('svs.2504kgp.seq.tsv.gz', 'w')
svs %>% as.data.frame %>% select(-strand, -width) %>% 
    write.table(file=outf, sep='\t', quote=FALSE, row.names=FALSE)
close(outf)
```

This information was then used to cluster the SV alleles into SV loci (sites) using:

- `annotate-sv-loci-mesa2k.R` for the 2,000 samples from MESA, creating `svs.2504kgp.svsite80al.tsv.gz`
- `annotate-sv-loci-2504kgp.R` for the 1000 Genomes Project, creating `svs.mesa2k.svsite80al.tsv.gz`

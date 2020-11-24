## Combined VCF across samples

A combined VCF was created using a [WDL workflow deposited on Dockstore](https://dockstore.org/workflows/github.com/jmonlong/wdl-workflows/bcftools_merge:giraffe-paper?tab=info).
Briefly, it uses bcftools and custom scripts to split multi-allelic variants, extract canonical SVs from vg calls, potentially merge heterozygous variants into homozygous variants, left-align and sort each VCF.
The VCF are then merged using `bcftools merge`.

The `align_variants.py` and `merge_exact_hets.py` are used by the WDL and are also deposited here for record keeping.
The docker container used by the WDL, and containing these scripts, was built from [this repository](https://github.com/jmonlong/docker-merge-sv-vg). 

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

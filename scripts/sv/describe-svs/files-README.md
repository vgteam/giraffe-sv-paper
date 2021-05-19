## Description of the files used in this analysis and/or distributed

Description of the files used in the scripts in this folder. 
It also describes the columns of some files distributed or archived

### `vggiraffe-sv-2504kgp-svsites.tsv.gz`

Information about each SV site, including frequencies of the most expressed and second most expressed alleles, and frequency in the 5 super-populations.
The columns in the file correspond to:

- `seqnames`/`start`/`end` coordinates of the SV
- `svsite` ID of the SV site
- `type` SV type (DEL or INS)
- `size` SV size in bp
- `clique` are all the allele similar to each other, forming a clique? I.e. are the allelic only due to small variation between alleles?
- `af.tot` total allele frequency, considering all alleles
- `af` allele frequency of the major allele
- `af.top2` allele frequency of the second most frequent allele
- `af.top.fc` fold change between major and second most frequent allele
- `af.EUR` allele frequency in the EUR super population
- `af.EAS` allele frequency in the EAS super population
- `af.AMR` allele frequency in the AMR super population
- `af.SAS` allele frequency in the SAS super population
- `af.AFR` allele frequency in the AFR super population
- `ref` sequence of the *REF*erence allele
- `alt` sequence of the *ALT*ernate allele

### `vggiraffe-sv-2504kgp-all-bysvsites.tsv.gz` or `svs.2504kgp.svsite80al.tsv.gz`

Information about each SV allele, grouped by SV site.
The columns in the file correspond to:

- `seqnames`/`start`/`end` coordinates of the SV
- `svid` ID of the SV
- `type` SV type (DEL or INS)
- `size` SV size in bp
- `af` the allele frequency
- `ac` the total allele count
- `nrefs` the number of samples with reference genotype
- `ncalls` the number of samples with non-missing genotype
- `ref` sequence of the *REF*erence allele
- `alt` sequence of the *ALT*ernate allele
- `clique` when grouped by SV site, are all the allele similar to each other, forming a clique? I.e. are the allelic only due to small variation between alleles?
- `svsite` ID of the SV site

### `vggiraffe-sv-2504kgp-svsite-ac.tsv.gz` or `2504kgp.svsite80al.ac.tsv.gz`

The allele counts for each of the 2,504 unrelated samples across SV sites.
This file is a matrix where each row correspond to a SV site and each column to a sample (the first column contain the ID of the SV site).

### `vggiraffe-sv-2504kgp-svsite-gq.tsv.gz` or `2504kgp.svsite80al.gq.tsv.gz`

The genotype quality for each of the 2,504 unrelated samples across SV sites
This file is a matrix where each row correspond to a SV site and each column to a sample (the first column contain the ID of the SV site).

### `vggiraffe-sv-2504kgp-pcgenes.tsv.gz`

SVs around protein-coding genes (promoter, UTR, intron or exon) in the 1000 Genomes Project dataset.
The columns in the file correspond to:

- `seqnames`/`start`/`end` coordinates of the SV
- `svsite` ID of the SV site
- `type` SV type (DEL or INS)
- `size` SV size in bp
- `gene.name` the gene name
- `impact` location relative to the gene: CDS (coding), intronic, UTR, promoter.

### `vggiraffe-sv-superpop-af-diff-med10.csv.gz`

- `svsite` ID of the SV site
- `Superpopulation` name of the super population
- `af` allele frequency in this super population
- `af.med` median allele frequency across all super populations
- `type` SV type
- `size` SV size in bp

### `vggiraffe-sv-mesa-svsites.tsv.gz` or `locs.mesa2k.svsite80al.tsv.gz`

Information about each SV site, including frequencies of the most expressed and second most expressed alleles.
The columns in the file correspond to:

- `seqnames`/`start`/`end` coordinates of the SV
- `svsite` ID of the SV site
- `type` SV type (DEL or INS)
- `clique` are all the allele similar to each other, forming a clique? I.e. are the allelic only due to small variation between alleles?
- `af.tot` total allele frequency, considering all alleles
- `ac.tot` total number of alleles, considering all alleles
- `af` allele frequency of the major allele
- `af.top2` allele frequency of the second most frequent allele
- `af.top.fc` fold change between major and second most frequent allele
- `loc.n` the number of allele in the SV site
- `size.min` the size of the shortest allele
- `size.max` the size of the largest allele
- `size` SV size in bp of the major allele
- `rep.sr.lc.sat` overlap a simple repeat, low-complexity region or satellite?
- `rep.sr.lc` overlap a simple repeat, low-complexity region?
- `rep.sr` overlap a simple repeat?
- `rep.lc` overlap a low-complexity region?
- `rep.sat` overlap a satellote?

### `vggiraffe-sv-mesa-all-bysvsites.tsv.gz` or `svs.mesa2k.svsite80al.tsv.gz`

Information about each SV allele, grouped by SV site.
The columns in the file correspond to:

- `seqnames`/`start`/`end` coordinates of the SV
- `svid` ID of the SV
- `type` SV type (DEL or INS)
- `size` SV size in bp
- `af` the allele frequency
- `ac` the total allele count
- `nrefs` the number of samples with reference genotype
- `ncalls` the number of samples with non-missing genotype
- `ref` sequence of the *REF*erence allele
- `alt` sequence of the *ALT*ernate allele
- `clique` when grouped by SV site, are all the allele similar to each other, forming a clique? I.e. are the allelic only due to small variation between alleles?
- `svsite` ID of the SV site


# Giraffe Paper Products

This archive contains reusable work products for the manuscript *Pangenomics enables genotyping known structural variants in 5,202 diverse genomes*. It is organized as follows:

* `HGSVC_*`: Minimal graph indexes needed for running Giraffe against the HGSVC structural variant graph, with 6 embedded haplotypes and a path cover of graph regions without haplotypes. Giraffe will be able to relatively quickly regenerate the `.min` and `.gg` index files at runtime.
    * `HGSVC_hs38d1.full.gbwt`: Haplotype index
    * `HGSVC_hs38d1.dist`: Distance index
    * `HGSVC_hs38d1.xg`: Graph in xg format

## Structural variant products

* VCF at the SV site level. Alleles were combined if matching (>=80% reciprocal overlap or sequence similarity). The allele was counted across all alleles at each site for each sample.
	* `vggiraffe-sv-mesa-svsites.vcf.gz` and `vggiraffe-sv-mesa-svsites.vcf.gz.tbi`: VCF and index for the 2,000 MESA samples. No genotypes, just SV site definition and allele frequency estimate.
	* `vggiraffe-sv-2504kgp-svsites.vcf.gz` and `vggiraffe-sv-2504kgp-svsites.vcf.gz.tbi`: VCF and index for the 2,504 unrelated individuals in the 1000 Genomes Project. No genotypes, just SV site definition and allele frequency estimate (including for each of the super populations EUR/AFR/EAS/SAS/AMR).
	* `vggiraffe-sv-2504kgp-svsites.gt.vcf.gz` and `vggiraffe-sv-2504kgp-svsites.gt.vcf.gz.tbi`: VCF and index for the 2,504 unrelated individuals in the 1000 Genomes Project. Includes allele counts, genotypes and genotype qualities, in addition to INFO such as allele frequency in all samples or for each of the super populations EUR/AFR/EAS/SAS/AMR.
* Raw VCFs: VCF containing all the information from `vg call` (inc. GL) but at the allele level, i.e. >1M alleles. 
	* `vggiraffe-sv-2504kgp-raw.vcf.gz` and `vggiraffe-sv-2504kgp-raw.vcf.gz.tbi` VCF and index for the 2,504 unrelated individuals of the 1000 Genomes Project.
	* `vggiraffe-sv-relkgp-raw.vcf.gz` and `vggiraffe-sv-relkgp-raw.vcf.gz.tbi` VCF and index for the related individuals of the 1000 Genomes Project.
* eQTLs from the Geuvadis data
	* `vggiraffe-sv-eqtl-geuvadis.FDR01.csv`: SV-eQTLs when analyzing SVs in all the samples or separately EUR and YRI samples.
	* `vggiraffe-geuvadis-sveqtl-gene-families.csv`: Gene families significantly enriched in SV-eQTLs.
	* `vggiraffe-geuvadis-eqtl-snv-indel-svs.csv.gz`: eQTLs when jointly analyzing SNV/indel/SV in all the samples.
	* `vggiraffe-geuvadis-eqtl-svonly.csv`: SV-eQTLs in genes with no SNV/indel eQTLs.
* `vggiraffe-sv-superpop-af-diff-med10.csv.gz`: SVs with population signatures in the 1000 Genomes Project dataset (>10% difference with the median frequency across all super population)
* `vggiraffe-sv-2504kgp-pcgenes.tsv.gz`: SVs around protein-coding genes (promoter, UTR, intron or exon) in the 1000 Genomes Project dataset

More details are provided about the information in each of these files in [`https://github.com/vgteam/giraffe-sv-paper/tree/master/scripts/sv/eqtl/files-README.md`](https://github.com/vgteam/giraffe-sv-paper/tree/master/scripts/sv/eqtl/files-README.md) and [`https://github.com/vgteam/giraffe-sv-paper/tree/master/scripts/sv/describe-svs/files-README.md`](https://github.com/vgteam/giraffe-sv-paper/tree/master/scripts/sv/describe-svs/files-README.md).



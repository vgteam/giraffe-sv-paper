## SV-eQTLS

```sh
FILES="vggiraffe-sv-eqtl-geuvadis.FDR01.csv
vggiraffe-geuvadis-sveqtl-gene-families.csv
vggiraffe-geuvadis-eqtl-svonly.csv
vggiraffe-geuvadis-eqtl-snv-indel-svs.csv.gz"
for ff in $FILES
do
	aws s3 cp eqtl/$ff s3://vg-k8s/users/jmonlong/manu-giraffe-sv/products/
done
```

## SVs and frequencies

```sh
#### 1000GP unrelated individuals
## SVs
aws s3 cp describe-svs/vggiraffe-sv-2504kgp-svsites.tsv.gz s3://vg-k8s/users/jmonlong/manu-giraffe-sv/products/vggiraffe-sv-2504kgp-svsites.tsv.gz
aws s3 cp describe-svs/svs.2504kgp.svsite80al.tsv.gz s3://vg-k8s/users/jmonlong/manu-giraffe-sv/products/vggiraffe-sv-2504kgp-all-bysvsites.tsv.gz
## allele counts and genotype qualities in matrix form
aws s3 cp describe-svs/2504kgp.svsite80al.ac.tsv.gz s3://vg-k8s/users/jmonlong/manu-giraffe-sv/products/vggiraffe-sv-2504kgp-svsite-ac.tsv.gz
aws s3 cp describe-svs/2504kgp.svsite80al.gq.tsv.gz s3://vg-k8s/users/jmonlong/manu-giraffe-sv/products/vggiraffe-sv-2504kgp-svsite-gq.tsv.gz
## Super population frequencies
gzip -f describe-svs/vggiraffe-sv-superpop-af-diff-med10.csv
aws s3 cp describe-svs/vggiraffe-sv-superpop-af-diff-med10.csv.gz s3://vg-k8s/users/jmonlong/manu-giraffe-sv/products/vggiraffe-sv-superpop-af-diff-med10.csv.gz
## subset in around protein coding genes
aws s3 cp describe-svs/vggiraffe-sv-2504kgp-pcgenes.tsv.gz s3://vg-k8s/users/jmonlong/manu-giraffe-sv/products/vggiraffe-sv-2504kgp-pcgenes.tsv.gz

#### 2K MESA individuals
aws s3 cp describe-svs/locs.mesa2k.svsite80al.tsv.gz s3://vg-k8s/users/jmonlong/manu-giraffe-sv/products/vggiraffe-sv-mesa-svsites.tsv.gz
aws s3 cp describe-svs/svs.mesa2k.svsite80al.tsv.gz s3://vg-k8s/users/jmonlong/manu-giraffe-sv/products/vggiraffe-sv-mesa-all-bysvsites.tsv.gz
```

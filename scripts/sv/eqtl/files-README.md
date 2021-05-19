## Description of the files used in this analysis and/or distributed

### `vggiraffe-sv-eqtl-geuvadis.FDR01.csv`

SV-eQTLs at FDR 1%. 
The columns in the file correspond to:

- `seqnames`/`start`/`end` coordinates of the SV
- `svid` ID of the SV
- `type` SV type (DEL or INS)
- `size` SV size in bp
- `pop` population(s) analyzed: `EUR + YRI`, `EUR`, or `YRI`
- `gene` Ensembl ID of the gene
- `gene_name` gene name
- `gene_type` gene type, e.g. `protein_coding` or `pseudogene`
- `beta` the effect coefficient
- `pvalue` the raw p-value
- `FDR` the adjusted p-value, or false discovery rate

### `vggiraffe-geuvadis-sveqtl-gene-families.csv`

Gene families enriched in SV-eQTLs
The columns in the file correspond to:

- `gene_family` name of the gene family
- `egenes` number of eGenes, i.e. genes with SV-eQTLs
- `esvs` number of SV-eQTLs
- `prop` proportion of the tested genes in the family that have SV-eQTLs
- `pv.hyper` p-value of the enrichment test (hypergeometric test)
- `qv.hyper` p-value of the enrichment test (hypergeometric test) adjusted for multiple testing

### `vggiraffe-geuvadis-eqtl-snv-indel-svs.csv.gz`

eQTLs at FDR 1% for the joint SNV + indels + SV analysis.
The columns in the file correspond to:

- `seqnames`/`start`/`end` coordinates of the variant
- `id` ID of the variant
- `type` variant type (*DEL*, *INS*, or *SNV-indel*)
- `size` variant size in bp (if SV)
- `gene` Ensembl ID of the gene
- `gene_name` gene name
- `gene_type` gene type, e.g. `protein_coding` or `pseudogene`
- `beta` the effect coefficient
- `pvalue` the raw p-value
- `FDR` the adjusted p-value, or false discovery rate

### `vggiraffe-geuvadis-eqtl-svonly.csv`

eQTLs in genes with only SV-eQTLs. From the joint SNV/indel/SV eQTL analysis, these genes had no SNV or indel eQTLs.
The columns in the file correspond to:

- `gene_type` the gene type. E.g. protein_coding, pseudogene, ...
- `gene_name` the gene name
- `gene` the gene ENSEMBL ID
- `pop` the population tested: EUR+YRI, EUR only, or YRI only
- `beta` the beta coefficient of the effect
- `FDR` significance (adjusted pvalue~false discovery rate)
- `svid` the ID of the SV
- `type` the type of SV
- `size` the size of the SV
- `seqnames` the chromosome where the SV is located
- `start` the starting position of the SV
- `end` the ending position of the SV

library(VariantAnnotation)

#### GIAB
vcf = readVcf('HG002_SVs_Tier1_v0.6.vcf.gz')
## PASS filter
table(rowRanges(vcf)$FILTER)
vcf = vcf[which(rowRanges(vcf)$FILTER == 'PASS')]
seqlevels(vcf) = paste0('chr', seqlevels(vcf))
writeVcf(vcf, 'HG002_SVs_Tier1_v0.6.filtered.vcf')

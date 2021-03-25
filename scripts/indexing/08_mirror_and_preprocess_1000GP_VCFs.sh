#!/usr/bin/env bash
# This mirrors the subset of the NYGC high coverage 1000 Genomes and related samples to S3.
# It also preprocesses the chrX phased VCF to fix the missing "ME" INFO header.
# To change where it goes, replace 's3://vg-k8s/users/adamnovak/projects/1000gp-giraffe/1000g/ftp/' with another path.
# You will need pbgzip installed: https://github.com/nh13/pbgzip

VCFS=()
for CHROM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; do
  VCFS+=("ftp://ftp.ebi.ac.uk/1000g/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHROM}.filtered.shapeit2-duohmm-phased.vcf.gz")
done
VCFS+=("ftp://ftp.ebi.ac.uk/1000g/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.vcf.gz" \
"ftp://ftp.ebi.ac.uk/1000g/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrY.recalibrated_variants.vcf.gz" \
"ftp://ftp.ebi.ac.uk/1000g/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_others.recalibrated_variants.vcf.gz")
ITER=0
for VCF in "${VCFS[@]}" ; do
FILENAME="$(echo "${VCF}" | rev | cut -f1 -d'/' | rev)"
DEST="$(echo "${VCF}" | sed 's|^ftp://ftp.ebi.ac.uk/1000g/ftp/|s3://vg-k8s/users/adamnovak/projects/1000gp-giraffe/1000g/ftp/|g')"
kubectl delete job adamnovak-get-vcf-${ITER} ; kubectl apply -f - <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: adamnovak-get-vcf-${ITER}
spec:
  ttlSecondsAfterFinished: 259200
  template:
    spec:
      containers:
      - name: main
        imagePullPolicy: Always
        image: quay.io/vgteam/vg:v1.31.0
        command:
        - /bin/bash
        - -c
        - |
          set -e -o pipefail
          curl -sSL "${VCF}" | aws s3 cp - "${DEST}"
          curl -sSL "${VCF}.tbi" | aws s3 cp - "${DEST}.tbi"
        volumeMounts:
        - mountPath: /tmp
          name: scratch-volume
        - mountPath: /root/.aws
          name: s3-credentials
        resources:
          limits:
            cpu: 1
            memory: "4Gi"
            ephemeral-storage: "10Gi"
        env:
        - name: DEBIAN_FRONTEND
          value: noninteractive
      restartPolicy: Never
      volumes:
      - name: scratch-volume
        emptyDir: {}
      - name: s3-credentials
        secret:
          secretName: shared-s3-credentials
      serviceAccountName: vg-svc
  backoffLimit: 0
EOF
((ITER=ITER+1))
done

# Now apply the ME INFO header
wget ftp://ftp.ebi.ac.uk/1000g/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.vcf.gz
pbgzip -n 16 -c -d CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.vcf.gz | sed 's|^#CHROM|##INFO=<ID=ME,Number=A,Type=Float,Description="Unknown">\n#CHROM|' | pbgzip -n 16 -c /dev/stdin | tee CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.noME.vcf.gz | aws s3 cp - s3://vg-k8s/users/adamnovak/projects/1000gp-giraffe/1000g/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.noME.vcf.gz 
tabix -p vcf CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.noME.vcf.gz
aws s3 cp CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.noME.vcf.gz s3://vg-k8s/users/adamnovak/projects/1000gp-giraffe/1000g/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.noME.vcf.gz.tbi


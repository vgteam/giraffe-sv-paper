#!/usr/bin/env bash

    
kubectl delete job adamnovak-make-1000gp-graphs-minaf ; kubectl apply -f - <<'EOF'
apiVersion: batch/v1
kind: Job
metadata:
  name: adamnovak-make-1000gp-graphs-minaf
spec:
  ttlSecondsAfterFinished: 259200
  template:
    spec:
      containers:
      - name: main
        imagePullPolicy: Always
        image: quay.io/ucsc_cgl/toil:5.3.0-py3.7
        command:
        - /bin/bash
        - -c
        - |
          set -e
          mkdir /tmp/work
          cd /tmp/work
          echo "Making graphs"
          virtualenv --system-site-packages --python python3 venv
          . venv/bin/activate
          pip3 install git+https://github.com/vgteam/toil-vg.git@75b7006380d5762bd520f3ed76055454f8c478d5#egg=toil-vg
          toil-vg generate-config --whole_genome >config.cfg
          sed -i'' config.cfg -e "s/construct-mem: '64G'/construct-mem: '128G'/g" -e "s/construct-disk: '64G'/construct-disk: '128G'/g"
          # VCF order must be 1-22, X, Y, chrM to match FASTA
          # Leaving M as "chrM" to match the SV graphs
          VCFS=()
          for CHROM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; do
              VCFS+=("s3://vg-k8s/profiling/graphs/v3-2/for-NA19239/1000gp/hs38d1/1000GP_hs38d1-vcfs/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHROM}.filtered.shapeit2-duohmm-phased_filter.vcf.gz")
          done
          VCFS+=("s3://vg-k8s/profiling/graphs/v3-2/for-NA19239/1000gp/hs38d1/1000GP_hs38d1-vcfs/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.noME_filter.vcf.gz" \
            "s3://vg-k8s/profiling/graphs/v3-2/for-NA19239/1000gp/hs38d1/1000GP_hs38d1-vcfs/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrY.recalibrated_variants_filter.vcf.gz" \
            "s3://vg-k8s/profiling/graphs/v3-2/for-NA19239/1000gp/hs38d1/1000GP_hs38d1-vcfs/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_others.recalibrated_variants_filter.vcf.gz")
          toil clean aws:us-west-2:adamnovak-make-1000gp-graphs-minaf
          toil-vg construct \
              aws:us-west-2:adamnovak-make-1000gp-graphs-minaf \
              aws:us-west-2:vg-k8s/profiling/graphs/v3-5/for-NA19239/1000gp/hs38d1 \
              --fasta s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.gz \
              --coalesce_regions s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.minor_contigs.tsv \
              --vcf "${VCFS[@]}" \
              --vcf_phasing "${VCFS[@]}" \
              --fasta_regions \
              --remove_chr_prefix \
              --alt_paths \
              --out_name 1000GP_hs38d1 \
              --all_index \
              --force_phasing True \
              --gbwt_prune \
              --min_af 0.01 \
              --merge_graphs \
              --keep_vcfs \
              --config config.cfg \
              --realTimeLogging \
              --logInfo \
              --batchSystem kubernetes \
              --container Singularity \
              --disableCaching false \
              --defaultDisk 200G --defaultMemory 200G \
              --vg_docker quay.io/vgteam/vg:ci-2890-655a9622c3d60e87f14b88d943fbd8554214a975
        lifecycle:
          preStop:
            exec:
              command:
              - /bin/bash
              - -c
              - |
                pgrep toil-vg | xargs kill -2
                pgrep -P 1 | xargs kill -2
                sleep 120
        volumeMounts:
        - mountPath: /tmp
          name: scratch-volume
        - mountPath: /root/.aws
          name: s3-credentials
        resources:
          limits:
            cpu: 2
            memory: "4Gi"
            ephemeral-storage: "50Gi"
        env:
        - name: TOIL_KUBERNETES_OWNER
          value: adamnovak
        - name: TOIL_AWS_SECRET_NAME
          value: shared-s3-credentials
        - name: TOIL_KUBERNETES_HOST_PATH
          value: /data/scratch
        - name: TOIL_WORKDIR
          value: /var/lib/toil
        - name: SINGULARITY_CACHEDIR
          value: /var/lib/toil/singularity-cache
        - name: DEBIAN_FRONTEND
          value: noninteractive
      restartPolicy: Never
      terminationGracePeriodSeconds: 120
      volumes:
      - name: scratch-volume
        emptyDir: {}
      - name: s3-credentials
        secret:
          secretName: shared-s3-credentials
      serviceAccountName: vg-svc
  backoffLimit: 0
EOF

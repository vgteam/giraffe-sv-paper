#!/usr/bin/env bash

# Was originally started with smaller resource limits and restarted a couple of times.
kubectl delete job adamnovak-make-linear-38-single ; kubectl apply -f - <<'EOF'
apiVersion: batch/v1
kind: Job
metadata:
  name: adamnovak-make-linear-38-single
spec:
  ttlSecondsAfterFinished: 259200
  template:
    spec:
      containers:
      - name: main
        imagePullPolicy: Always
        image: quay.io/ucsc_cgl/toil:4.3.0a1-ab42ab9e7ee078c3d27395e2a8b0041e4d7b8c72-py3.7
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
          pip3 install --upgrade git+https://github.com/vgteam/toil-vg.git@074b6504ea7899b3fdea06b38416c17dd78d353b#egg=toil-vg
          toil-vg generate-config --whole_genome | sed "s/gcsa-index-disk.*/gcsa-index-disk: '250G'/g" | sed "s/xg-index-mem.*/xg-index-mem: '100G'/g" | sed "s/snarl-index-mem.*/snarl-index-mem: '140G'/g" | sed "s/distance-index-mem.*/distance-index-mem: '140G'/g" > config.cfg
          time toil-vg construct \
          aws:us-west-2:adamnovak-make-linear-38-single \
          aws:us-west-2:vg-k8s/profiling/graphs/v2-2/generic/primary/hs38d1 \
          --fasta s3://vg-k8s/profiling/data/GCA_000001405.15_GRCh38_no_alt_analysis_set_plus_GCA_000786075.2_hs38d1_genomic.fna.gz \
          --fasta_regions \
          --remove_chr_prefix \
          --out_name primary_hs38d1 \
          --gcsa_index --xg_index --trivial_snarls_index --distance_index --id_ranges_index \
          --primary \
          --merge_graphs \
          --config config.cfg \
          --realTimeLogging \
          --logInfo \
          --container Singularity \
          --disableCaching false \
          --vg_docker quay.io/vgteam/vg:ci-2284-dc119fa046aa7131a1a8e026be36da2d79bc2f22 \
          --batchSystem singleMachine
        lifecycle:
          preStop:
            exec:
              command:
              - /bin/bash
              - -c
              - |
                pgrep -P 1 | xargs kill -2
        volumeMounts:
        - mountPath: /tmp
          name: scratch-volume
        - mountPath: /root/.aws
          name: s3-credentials
        resources:
          limits:
            cpu: 32
            memory: "150Gi"
            ephemeral-storage: "2200Gi"
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
      volumes:
      - name: scratch-volume
        emptyDir: {}
      - name: s3-credentials
        secret:
          secretName: shared-s3-credentials
      serviceAccountName: vg-svc
  backoffLimit: 0
EOF


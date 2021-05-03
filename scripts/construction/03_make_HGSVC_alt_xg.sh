# The HGSVC graph wasn't initially made with alternate allele paths in the XG,
# so replace the XG with one with alternate allele paths, so that simulated
# reads can be annotated with positions on them.
GRAPH_BASE=s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1
kubectl delete job adamnovak-make-xg
cat <<EOF | tee /dev/stderr | kubectl apply -f -
apiVersion: batch/v1
kind: Job
metadata:
  name: adamnovak-make-xg
spec:
  ttlSecondsAfterFinished: 259200
  template:
    spec:
      containers:
      - name: main
        imagePullPolicy: Always
        image: "quay.io/vgteam/vg:ci-2035-42bb4f3123c79006f0d4ffe8e6287627c1dc50ae"
        command:
        - /bin/bash
        - -c
        - |
          set -ex
          mkdir /tmp/work
          cd /tmp/work
          aws s3 cp --no-progress ${GRAPH_BASE}.vg input.vg
          vg index -x output.xg -L input.vg
          aws s3 cp --no-progress output.xg ${GRAPH_BASE}.xg
        volumeMounts:
        - mountPath: /tmp
          name: scratch-volume
        - mountPath: /root/.aws
          name: s3-credentials
        resources:
          limits:
            cpu: 2
            memory: "180Gi"
            ephemeral-storage: "100Gi"
        env:
        - name: DEBIAN_FRONTEND
          value: noninteractive
        - name: VG_FULL_TRACEBACK
          value: "1"
      restartPolicy: Never
      priorityClassName: medium-priority
      volumes:
      - name: scratch-volume
        emptyDir: {}
      - name: s3-credentials
        secret:
          secretName: shared-s3-credentials
  backoffLimit: 0
EOF

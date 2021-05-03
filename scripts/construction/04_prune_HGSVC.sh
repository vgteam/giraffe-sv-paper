ITER=1
for GRAPH_BASE in s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1 ; do
kubectl delete job adamnovak-make-gcsa-${ITER}
cat <<EOF | tee /dev/stderr | kubectl apply -f -
apiVersion: batch/v1
kind: Job
metadata:
  name: adamnovak-make-gcsa-${ITER}
spec:
  ttlSecondsAfterFinished: 259200
  template:
    spec:
      containers:
      - name: main
        imagePullPolicy: Always
        image: "quay.io/vgteam/vg:ci-1954-8ff022c3a36c5fbc2a63faf477c5bf9ac37e29d7"
        command:
        - /bin/bash
        - -c
        - |
          set -ex
          mkdir /tmp/work
          cd /tmp/work
          aws s3 cp ${GRAPH_BASE}.vg input.vg
          aws s3 cp ${GRAPH_BASE}.gbwt input.gbwt
          vg prune input.vg --threads 32 --mapping output.pruned.mapping --unfold-paths --gbwt-name input.gbwt --progress > output.pruned.vg
          aws s3 cp output.pruned.vg ${GRAPH_BASE}.pruned.vg
          aws s3 cp output.pruned.mapping ${GRAPH_BASE}.pruned.mapping
          # These steps were in the original script but failed 
          # vg index -g output.gcsa --threads 32 --mapping output.pruned.mapping output.pruned.vg
          # aws s3 cp output.gcsa ${GRAPH_BASE}.gcsa
        volumeMounts:
        - mountPath: /tmp
          name: scratch-volume
        - mountPath: /root/.aws
          name: s3-credentials
        resources:
          limits:
            cpu: 32
            memory: "350Gi"
            ephemeral-storage: "2200G"
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
((ITER++))
done

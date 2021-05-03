ITER=1
for GRAPH_BASE in s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1 ; do
for GBWT_TYPE in full cover sampled ; do
kubectl delete job adamnovak-make-min-${ITER}
cat <<EOF | tee /dev/stderr | kubectl apply -f -
apiVersion: batch/v1
kind: Job
metadata:
  name: adamnovak-make-min-${ITER}
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
          aws s3 cp --no-progress ${GRAPH_BASE}.${GBWT_TYPE}.gbwt input.gbwt
          aws s3 cp --no-progress ${GRAPH_BASE}.${GBWT_TYPE}.gg input.gg
          vg minimizer -t 16 -p -i output.min -g input.gbwt -G input.gg
          aws s3 cp --no-progress output.min ${GRAPH_BASE}.${GBWT_TYPE}.min
        volumeMounts:
        - mountPath: /tmp
          name: scratch-volume
        - mountPath: /root/.aws
          name: s3-credentials
        resources:
          limits:
            cpu: 16
            memory: "120Gi"
            ephemeral-storage: "200Gi"
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
done

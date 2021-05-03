ITER=1
for GRAPH_BASE in s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1 ; do
for GBWT_TYPE in cover sampled ; do
kubectl delete job adamnovak-make-gbwt-${ITER}
cat <<EOF | tee /dev/stderr | kubectl apply -f -
apiVersion: batch/v1
kind: Job
metadata:
  name: adamnovak-make-gbwt-${ITER}
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
          aws s3 cp ${GRAPH_BASE}.xg input.xg
          if [[ "${GBWT_TYPE}" == "cover" ]] ; then
            vg gbwt -p -g output.gg -o output.gbwt -x input.xg -P
          else
            aws s3 cp ${GRAPH_BASE}.gbwt input.gbwt
            vg gbwt -p -g output.gg -o output.gbwt -x input.xg -l input.gbwt
          fi
          aws s3 cp output.gbwt ${GRAPH_BASE}.${GBWT_TYPE}.gbwt
          aws s3 cp output.gg ${GRAPH_BASE}.${GBWT_TYPE}.gg
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
((ITER++))
done
done

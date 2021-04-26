ITER=1
for GRAPH_BASE in s3://vg-k8s/profiling/graphs/v3/for-NA19239/1000gplo/hs38d1/1000GPlo_hs38d1_filter ; do
for GBWT_TYPE in full ; do
kubectl delete job xhchang-make-gbwt-${ITER}
cat <<EOF | tee /dev/stderr | kubectl apply -f -
apiVersion: batch/v1
kind: Job
metadata:
  name: xhchang-make-gbwt-${ITER}
spec:
  ttlSecondsAfterFinished: 259200
  template:
    spec:
      containers:
      - name: main
        imagePullPolicy: Always
        image: xhchang/vg:giraffe-paper
        command:
        - /bin/bash
        - -c
        - |
          set -ex
          mkdir /tmp/work
          cd /tmp/work
          aws s3 cp --no-progress ${GRAPH_BASE}.xg input.xg
          aws s3 cp --no-progress ${GRAPH_BASE}.dist input.dist
          if [[ "${GBWT_TYPE}" == "cover" ]] ; then
            vg gbwt -p -g output.gg -o output.gbwt -x input.xg -P
          elif [[ "${GBWT_TYPE}" == "sampled" ]] ; then
            aws s3 cp --no-progress ${GRAPH_BASE}.gbwt input.gbwt
            vg gbwt -p -g output.gg -o output.gbwt -x input.xg -l input.gbwt
          elif [[ "${GBWT_TYPE}" == "full" ]] ; then
            aws s3 cp --no-progress ${GRAPH_BASE}.gbwt input.gbwt
            vg gbwt -p -g output.gg -o output.gbwt -x input.xg -a input.gbwt
          fi
          aws s3 cp --no-progress output.gbwt ${GRAPH_BASE}.${GBWT_TYPE}.gbwt
          aws s3 cp --no-progress output.gg ${GRAPH_BASE}.${GBWT_TYPE}.gg

          vg minimizer -t 1 -p -i output.min -d input.dist -g output.gbwt -G output.gg
          aws s3 cp --no-progress output.min ${GRAPH_BASE}.${GBWT_TYPE}.min
        volumeMounts:
        - mountPath: /tmp
          name: scratch-volume
        - mountPath: /root/.aws
          name: s3-credentials
        resources:
          limits:
            cpu: 2
            memory: "200Gi"
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

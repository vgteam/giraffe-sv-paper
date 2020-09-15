ITER=1
GBWT_TYPE=cover
for GRAPH_BASE in s3://vg-k8s/profiling/graphs/v2/for-NA19240/hgsvc/hs38d1/HGSVC_hs38d1 s3://vg-k8s/profiling/graphs/v2/for-NA19239/1kg/hs37d5/1kg_hs37d5_filter ; do
for SAMPLED_PATHS in 1 2 4 8 16 32 64 128 ; do
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
        image: "quay.io/vgteam/vg:v1.26.1"
        command:
        - /bin/bash
        - -c
        - |
          set -ex
          mkdir /tmp/work
          cd /tmp/work
          aws s3 cp ${GRAPH_BASE}.xg input.xg
          if [[ "${GBWT_TYPE}" == "cover" ]] ; then
            vg gbwt -p -g output.gg -o output.gbwt -x input.xg -P -n ${SAMPLED_PATHS}
          else
            aws s3 cp ${GRAPH_BASE}.gbwt input.gbwt
            vg gbwt -p -g output.gg -o output.gbwt -x input.xg -l input.gbwt -n ${SAMPLED_PATHS}
          fi
          aws s3 cp output.gbwt ${GRAPH_BASE}.${GBWT_TYPE}.${SAMPLED_PATHS}.gbwt
          aws s3 cp output.gg ${GRAPH_BASE}.${GBWT_TYPE}.${SAMPLED_PATHS}.gg
          vg minimizer -t 16 -p -i output.min -g output.gbwt -G output.gg
          aws s3 cp --no-progress output.min ${GRAPH_BASE}.${GBWT_TYPE}.${SAMPLED_PATHS}.min
        volumeMounts:
        - mountPath: /tmp
          name: scratch-volume
        - mountPath: /root/.aws
          name: s3-credentials
        resources:
          limits:
            cpu: 16
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

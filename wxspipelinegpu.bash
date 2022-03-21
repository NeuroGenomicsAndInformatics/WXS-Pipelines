#!/bin/bash
export PATH="/opt/miniconda/bin:$PATH"
export TMP_DIR="/scratch1/fs1/cruchagac/parabricks-tmp"
[ ! -d $TMP_DIR ] && mkdir $TMP_DIR
if [[ ! -d /scratch1/fs1/cruchagac/${USER}/c1out/logs ]]; then mkdir /scratch1/fs1/cruchagac/${USER}/c1out/logs; fi
JOB_GROUP="/${USER}/compute-cruchagac"
bgadd -L 10 ${JOB_GROUP}
LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/cruchagac:/scratch1/fs1/cruchagac \
/scratch1/fs1/ris/application/parabricks-license:/opt/parabricks" \
LSF_DOCKER_NETWORK=host \
LSF_DOCKER_RUN_LOGLEVEL=DEBUG \
LSF_DOCKER_ENTRYPOINT=/bin/sh \
TMP_DIR="/scratch1/fs1/cruchagac/parabricks-tmp" \
bsub -g ${JOB_GROUP} \
-J ngi-${USER}-stage1gpu-%J \
-n 32 \
-o /scratch1/fs1/cruchagac/${USER}/c1out/logs/GPU.%J.out \
-R 'gpuhost rusage[mem=196GB] span[hosts=1]' \
-M 196GB \
-G compute-cruchagac \
-q general \
-gpu "num=4:gmodel=TeslaV100_SXM2_32GB:j_exclusive=yes" \
-a 'docker(gcr.io/ris-registry-shared/parabricks)' pbrun "$@"

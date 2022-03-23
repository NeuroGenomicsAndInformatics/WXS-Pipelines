#!/bin/bash
export PATH="/opt/miniconda/bin:$PATH"
export TMP_DIR="/scratch1/fs1/cruchagac/${USER}/parabricks-tmp"
[ ! -d $TMP_DIR ] && mkdir $TMP_DIR
if [[ ! -d /scratch1/fs1/cruchagac/${USER}/c1out/logs ]]; then mkdir /scratch1/fs1/cruchagac/${USER}/c1out/logs; fi
JOB_GROUP="/${USER}/compute-cruchagac"
bgadd -L 10 ${JOB_GROUP}
LSF_DOCKER_VOLUMES="/storage1/fs1/cruchagac/Active:/storage1/fs1/cruchagac/Active \
/scratch1/fs1/cruchagac:/scratch1/fs1/cruchagac \
/scratch1/fs1/ris/application/parabricks:/opt/parabricks \
$HOME:$HOME" \
LSF_DOCKER_NETWORK=host \
LSF_DOCKER_RUN_LOGLEVEL=DEBUG \
LSF_DOCKER_ENTRYPOINT=/bin/sh \
bsub -g ${JOB_GROUP} \
-J ngi-${USER}-stage1gpu-%J \
-N \
-o /scratch1/fs1/cruchagac/${USER}/c1out/logs/GPU.%J.out \
-R 'select[gpuhost && mem>64GB && ncpus>15] rusage[mem=64GB] span[hosts=1]' \
-M 64GB \
-G compute-cruchagac \
-q general \
-gpu "num=1:gmodel=TeslaV100_SXM2_32GB:j_exclusive=yes" \
-a 'docker(gcr.io/ris-registry-shared/parabricks)' pbrun "$@" --tmp-dir $TMP_DIR

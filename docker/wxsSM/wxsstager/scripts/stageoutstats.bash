#!/bin/bash
rsync -r ${OUTDIR}/ ${FINAL_OUTDIR}
bash /scripts/statsupdate.bash

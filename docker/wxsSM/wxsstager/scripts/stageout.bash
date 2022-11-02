#!/bin/bash
rsync -r ${OUTDIR}/ ${FINAL_OUTDIR}
rm -R ${OUTDIR}

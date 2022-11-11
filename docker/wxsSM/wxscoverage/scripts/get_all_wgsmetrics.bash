#!/bin/bash
if [ -z $FULLSMID ]; then
  bash /scripts/gatkwgsmetrics.bash $1
  bash /scripts/gatkwgsmetrics_ex.bash $1
  bash /scripts/gatkrawwgsmetrics.bash $1
else
  bash /scripts/gatkwgsmetrics_pipe.bash
  bash /scripts/gatkwgsmetrics_ex_pipe.bash
  bash /scripts/gatkrawwgsmetrics_pipe.bash
fi

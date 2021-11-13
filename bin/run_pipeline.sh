#!/bin/bash
# Process the CFHT data into a previously-created butler repo.

# exit script if any command fails
set -e
# print all the commands as they are run
set -o xtrace

REPO=repo
NUMPROC=8

COLLECTIONS="MegaPrime/raw/all,MegaPrime/calib/unbounded,refcats"

pipetask --long-log run -j $NUMPROC -b $REPO/butler.yaml -p config/'DRP.yaml#singleFrame' -i $COLLECTIONS --register-dataset-types -o singleFrame &> singleframe.log


# For testing things: just two detectors, debug logging, no multiprocessing
#pipetask run -b $REPO/butler.yaml -p config/'DRP.yaml#singleFrame' -i $COLLECTIONS --register-dataset-types -o singleFrame -d "detector=7" --debug

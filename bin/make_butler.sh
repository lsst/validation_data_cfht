#!/bin/bash
# Setup a butler repo for validation_data_cfht to run a pipeline on.

# exit script if any command fails
set -e
# print all the commands as they are run
set -o xtrace

REPO=repo

butler create $REPO
butler register-instrument $REPO lsst.obs.cfht.MegaPrime

butler register-dataset-type $REPO gaia_dr2_20200414 SimpleCatalog htm7
butler ingest-files -t symlink $REPO gaia_dr2_20200414 refcats gaia_dr2_20200414.ecsv
butler register-dataset-type $REPO ps1_pv3_3pi_20170110 SimpleCatalog htm7
butler ingest-files -t symlink $REPO ps1_pv3_3pi_20170110 refcats ps1_pv3_3pi_20170110.ecsv
butler register-dataset-type $REPO sdss_dr9_fink_v5b SimpleCatalog htm7
butler ingest-files -t symlink $REPO sdss_dr9_fink_v5b refcats sdss-dr9-fink-v5b.ecsv

butler register-skymap $REPO -C config/makeSkyMap.py -c name="discrete"
butler write-curated-calibrations $REPO lsst.obs.cfht.MegaPrime

butler ingest-raws $REPO raw/*.fz --transfer link
butler define-visits $REPO lsst.obs.cfht.MegaPrime

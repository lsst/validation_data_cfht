# You have to have source an LSST stack install already.
# This varies by system and installation and so is not included in this script.

# w_2018_37 was the latest weekly on lsst-dev the week before DM-14868 got merged.
# I also had to set up the tickets/DM-14868 version of obs_subaru
# But that won't be necessary once this ticket is merged, so it's
# commented out in these reprocessing instructions.
setup obs_cfht -t w_2018_37
setup -k -r ~wmwv/local/lsst/obs_cfht 
# To be run from within validation_data_cfht repo
# were we were at DM-14868
setup -k -r .

mkdir data
echo lsst.obs.cfht.MegacamMapper > data/_mapper
ingestImages.py data ${VALIDATION_DATA_CFHT_DIR}/raw/*.fz   --mode copy
# Link in the reference catalogs
ln -s ${VALIDATION_DATA_CFHT_DIR}/ref_cats data/ref_cats

export OMP_NUM_THREADS=1  # Suppress OMP parallelism.  We parallelize by CCD.
processCcd.py data \
     --output data \
     @${VALIDATION_DATA_CFHT_DIR}/Cfht.list \
     -j 4 \
   >& processCcd.log

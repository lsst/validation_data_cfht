validation_data_cfht
====================

Data from CFHT to validate the performance of the LSST DM stack.

Data are used by `validate_drp`.

Credit to Dominique Boutigny for providing the data and starting tests.

Create an input repo
--------------------

If you would like to create an input repo to then run `processCcd.py` on, run:

```
setup obs_cfht
setup validation_data_cfht

mkdir -p Cfht/input
echo lsst.obs.cfht.MegacamMapper > Cfht/input/_mapper
ingestImages.py Cfht/input ${VALIDATION_DATA_CFHT_DIR}/raw/*.fz --mode link
```


Git LFS
-------

To clone and use this repository, you'll need Git Large File Storage (LFS).

Our [Developer Guide](http://developer.lsst.io/en/latest/tools/git_lfs.html) explains how to setup Git LFS for LSST development.

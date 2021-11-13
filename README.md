validation_data_cfht
====================

Raw data from CFHT to validate the performance of the LSST DM stack and provide reference comparisons for development.
The data consists of two full CFHT images centered at ra/dec ``(214.884832, 52.6622199)``, plus
refcats covering the imaging region.

Data are targeted at use by LSST DM Stack developers, continuous integration tests and `faro`

The "raw" data are not actually truly "raw", but rather post-ISR processing of the CFHT data by the Terapix(?) archive.
http:terapix. [**]

These data are provided without a preset Butler repository, but with the necessary scripts and config files to setup and process the data easily.

[**] Fix this URL

Credit to Dominique Boutigny for selecting and providing these data from the public CFHT survey.

Example usage
-------------

The ``bin/`` directory contains all necessary scripts for configuring and processing this data.
To create a gen3 butler `repo/` directory and ingest the raws and refcats to run processing on, run:

```
bash bin/make_butler.sh
```

To process the data with the configurations in `config/`, run:

```
bash bin/run_pipeline.sh
```

Files
-----
path                  | description
:---------------------|:-----------------------------
`bin`                 | Scripts to process the data.
`config`              | Configurations for setting up the butler and running pipetasks.
`raw`                 | Photometrically and astrometrically calibrated data
                      |   as processed by Terapix
`data`                | Butler repo of ingested raw data and processCcd results
`ref_cats`            | HTM indexed catalog files from SDSS, Pan-Starrs, and Gaia
                      |   for astrometric and photometric calibration.


Git LFS
-------

To clone and use this repository, you'll need Git Large File Storage (LFS).

Our [Developer Guide](http://developer.lsst.io/en/latest/tools/git_lfs.html) explains how to setup Git LFS for LSST development.

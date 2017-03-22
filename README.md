validation_data_cfht
====================

Raw data from CFHT and a Butler repository of the processed outputs to validate the performance of the LSST DM stack and provide reference comparisons for development.

Data are targeted at use by LSST DM Stack developers, continuous integration tests, and `validate_drp`.  General users who wish to simply see examples of the steps and products of the LSST DM Stack processing may likely also find this useful.

The use model for reprocessing these data for comparisons and continuous integration, is to use these data in as an input repository and then write to a separate output repository (by, e.g.,  specifying `--output` on the command line call to a Task).  One might then wish to compare the outputs of the two repositories at the image, catalog, or summary statistic or visualization plot level.

The "raw" data are not actually truly "raw", but rather post-ISR processing of the CFHT data by the Terapix(?) archive.
http:terapix. [**]

These data are provided as a Butler repository in its POSIX filesystem mode*.  This provides straightforward access through both the Butler mechanism for clean programatic use, but also relatively simple access directly from the file system for quick inspection and orientation to the data.

[*] Mode isn't quite the right word here.
[**] Fix this URL

Credit to Dominique Boutigny for selecting and providing these data from the public CFHT survey.

Example usage
-------------

```
setup validation_data_cfht

NEW_OUTPUT_REPO=CFHT_data

export ASTROMETRY_NET_DATA_DIR=${VALIDATION_DATA_CFHT_DIR}/astrometry_net_data

processCcd.py ${VALIDATION_DATA_CFHT_DIR}/data \
    --output ${NEW_OUTPUT_REPO} \
    @${VALIDATION_DATA_CFHT_DIR}/Cfht.list \
    -j 4 \
    --clobber-config 
```

Notes:
 * There will be a `${NEW_OUTPUT_REPO}/_parent` link back to the input repository `${VALIDATION_DATA_CFHT_DIR}/data`.
 * The list of images (`dataIds`) to process is in `@${VALIDATION_DATA_CFHT_DIR}/Cfht.list`
 * `-j 4` specifies using 4 cores.  You may wish to change to an appropriate number on your system, but the intent is that `-j 4` should be a reasonable default in 2016.
 * We specifically use `--clobber-config` here because we're running off an already existing repository.

Analyzing the repository
------------------------
One might then choose to use `validate_drp` to analyze the peformance of the results against SRD metrics.

```
setup validate_drp

validateDrp.py CFHT_data
```

Recreating the repository
-------------------------
This repository was created using `examples/runCfhtTest.sh` from the `validate_drp` package.

To fully recreate this Butler `repo` from the `raw` data, set the `mapper` and add the `ingesetImages.py` step:

```
setup validation_data_cfht

mkdir data
echo lsst.obs.cfht.MegacamMapper > data/_mapper
ingestImages.py data ${VALIDATION_DATA_CFHT_DIR}/raw/*.fz   --mode copy

export ASTROMETRY_NET_DATA_DIR=${VALIDATION_DATA_CFHT_DIR}/astrometry_net_data
processCcd.py data \
    @${VALIDATION_DATA_CFHT_DIR}/Cfht.list \
    -j 4 \
    --logdest processCcd.log
```

Notes
 1. We use `--copy` to create a full copy of the raw images in the repo.
 2. No separate `--output` repo was specified.  In this case we intentionally wish to create
the output products in the same repo.

The packages and versions used were recorded into `eups_setup_used.txt`:

```
eups list --setup | awk '{printf "%-30s %s\n", $1, $2}' > eups_setup_used.txt
```

and the repo was made read-only to prevent accidental writes to this repo.

```
chmod -R ugo-w data
```

Files
-----
path                  | description
:---------------------|:-----------------------------
`raw`                 | Photometrically and astrometrically calibrated data
                      |   as processed by Terapix
`data`                | Butler repo of ingested raw data and processCcd results
`astrometry_net_data` | SDSS DR9 catalog files in astrometry.net format
                      |   as photometrically recalibrated by Doug Finkbeiner's group.
`processCcd.log`      | Output from the processCcd run on the `data` repo.
`eups_setup_used.txt` | EUPS setup configuration for ingestImages and processCcd run
`Cfht.list`           | List of dataIds in this repo.  For use in running Tasks.


Git LFS
-------

To clone and use this repository, you'll need Git Large File Storage (LFS).

Our [Developer Guide](http://developer.lsst.io/en/latest/tools/git_lfs.html) explains how to setup Git LFS for LSST development.

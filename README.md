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

See `examples/runCfhtTest.sh` from the `validate_drp` package for information about how to process the raw data contained in this repository.

Notes:
 * The list of images (`dataIds`) to process is in `@${VALIDATION_DATA_CFHT_DIR}/Cfht.list`

Analyzing the repository
------------------------
One might then choose to use `validate_drp` to analyze the peformance of the results against SRD metrics.

```
setup validate_drp

validateDrp.py CFHT_data
```

Recreating the repository
-------------------------
This repository was created using `examples/runCfhtTest.sh -C` from the `validate_drp` package.
See that file for how to process the raw data contained in this repository.

Notes

 1. We use `-C` to create a full copy of the raw images in the repo.
 2. The packages and versions used were recorded into `eups_setup_used.txt`:

    ```
    eups list --setup | awk '{printf "%-30s %s\n", $1, $2}' > eups_setup_used.txt
    ```

    and the repo was made read-only to prevent accidental writes to this repo.
    
    ```
    chmod -R ugo-w data
    ```

 3. The reference catalogs get linked in to the new repo by `runCfhtTest.sh`, but it is an absolute link. We replace with a relative link that by running: `rm data/input/ref_cats && ln -s ../../ref_cats/ data/input/`.

Files
-----
path                  | description
:---------------------|:-----------------------------
`raw`                 | Photometrically and astrometrically calibrated data
                      |   as processed by Terapix
`data`                | Butler repo of ingested raw data and processCcd results
`astrometry_net_data` | SDSS DR9 catalog files in astrometry.net format
                      |   as photometrically recalibrated by Doug Finkbeiner's group.
`ref_cats`            | HTM indexed catalog files from both SDSS and Pan-Starrs for
                      |   astrometric and photometric calibration.
`processCcd.log`      | Output from the processCcd run on the `data` repo.
`eups_setup_used.txt` | EUPS setup configuration for ingestImages and processCcd run
`Cfht.list`           | List of dataIds in this repo.  For use in running Tasks.


Git LFS
-------

To clone and use this repository, you'll need Git Large File Storage (LFS).

Our [Developer Guide](http://developer.lsst.io/en/latest/tools/git_lfs.html) explains how to setup Git LFS for LSST development.

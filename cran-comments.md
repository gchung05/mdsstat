## Test environments
* local Windows 10 install, R 3.6.2
* travis-ci, R 3.6.2
* win-builder devel
* win-builder release
* rhub

## R CMD check results
There were no ERRORs, WARNINGs, and no NOTEs.

## R CRAN web check results
There were no ERRORs, WARNINGs, and 1 NOTE.
NOTE: "Namespace in Imports field not imported from: ‘mds’
  All declared Imports should be used."
(Confirmed with Uwe 2019-01-23) Examples use data (mds_ts) that depend on mds package.

## r-hub check results
There were 0 ERRORs, 0 WARNINGs, and 0 NOTES.

## revdepcheck results
There were 0 PROBLEMS and 0 FAILURES.

## Downstream dependencies
There are no downstream dependencies.

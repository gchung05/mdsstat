`mdsstat` 0.2.2
---------------------------------------
UNDER DEVELOPMENT

**Implemented Updates**

- Added SPRT algorithm
- Added GPS algorithm

**Bigfixes**

- Fixed `analysis_of` attribute for all algorithms

**Potential Updates**

- Add GPS algorithm
- Add BCPNN algorithm
- Add Poisson MaxSPRT algorithm

`mdsstat` 0.2.1
---------------------------------------

**Implemented Updates**

- Added continuity adjustment option to DPA algorithms `prr()` and `ror()`
- Added CUSUM algorithm
- Added ROR algorithm
- Added support for more complete hierarchical descriptions (see mds package)
- Reorganized R scripts for easier navigation

**Bugfixes**

- `shewhart()` mu and sigma no longer calculated using the current index month

# `mdsstat` 0.2.0
---------------------------------------

**Implemented Updates**

- `run_algos()`: add parameter to skip/warn/stop when trying to run DPA on non-DPA time series
- Allow `run_algos()` to run on a list of `mds_ts` or other properly formatted objects

**Bugfixes**

- `shewhart()` now outputs the correct signal logic
- `prr()` now outputs the correct signal logic and p-value

# `mdsstat` 0.1.0
---------------------------------------

- Initial Release. Yay!

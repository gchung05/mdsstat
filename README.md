
<!-- README.md is generated from README.Rmd. Please edit that file -->
The `mdsstat` package:

-   Standardizes the output of various statistical trending algorithms
-   Allows running of multiple algorithms on the same data
-   Allows running of both **disproportionality** and **quality control** algorithms
-   Creates lightweight analysis definitions and output files for auditability, documentation, and reproducibility

**Why?**

There are many ways to trend medical device event data. Some are drawn from the [quality control](https://en.wikipedia.org/wiki/Quality_control) discipline, others from disproportionality analysis used in pharmacoepidemiology, and yet others from the general field of statistics.

There is a need to rigorously compare and contrast these various methods to more fully understand their respective performance and applicability in surveillance of medical devices.

**How?**

The `mdsstat` package aims to provide a collection of statistical trending algorithms used in medical device surveillance. Furthermore, each algorithm is written with a standardized, reusable framework philosophy. The same input data can be fed through multiple algorithms. All algorithms return **results that can be sorted, stacked, and compared.**

This package is written in tandem with the `mds` package. These are complementary in the sense that:

-   `mds` standardizes medical device event data.
-   `mdsstat` standardizes the statistical trending of medical device event data.

While `mdsstat` algorithms can run on generic R data frames, additional efficiency and traceability benefits are derived by running on data frames generated by `mds::time_series()` from the `mds` package.

Refer to the package vignette for available algorithms and guided examples.


# pema: Penalized Meta-Analysis <!--a href='https://osf.io/zcvbs/'><img src='https://github.com/cjvanlissa/pema/raw/master/docs/pema_icon.png' align="right" height="139" /></a-->

<!-- [![CRAN status](https://www.r-pkg.org/badges/version/pema)](https://cran.r-project.org/package=pema) -->
<!-- [![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/pema?color=blue)](https://r-pkg.org/pkg/pema) -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/cjvanlissa/pema/workflows/R-CMD-check/badge.svg)](https://github.com/cjvanlissa/pema/actions)
<!-- [![R-CMD-check](https://github.com/cjvanlissa/pema/workflows/R-CMD-check/badge.svg)](https://github.com/cjvanlissa/pema/actions) -->
<!-- [![codecov](https://codecov.io/gh/cjvanlissa/pema/branch/master/graph/badge.svg?token=7S9XKDRT4M)](https://codecov.io/gh/cjvanlissa/pema) -->
[![Contributor
Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](https://www.contributor-covenant.org/version/2/0/code_of_conduct.html)
<!-- [![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/3969/badge)](https://bestpractices.coreinfrastructure.org/projects/3969) -->
<!--[![DOI](http://joss.theoj.org/papers/10.21105/joss.00978/status.svg)](https://doi.org/10.21105/joss.00978)-->

Conduct *pe*nalized *m*eta-*a*nalysis (*“pema”*) In meta-analysis, there
are often between-study differences. These can be coded as moderator
variables, and controlled for using meta-regression. However, if the
number of moderators is large relative to the number of studies, such an
analysis may be overfitted. Penalized meta-regression is useful in these
cases, because it shrinks the regression slopes of irrelevant moderators
towards zero.

<!-- ## Where do I start? -->
<!-- For most users, the recommended starting point is to [read the paper](https://osf.io/zcvbs/), currently in press in [Data Science](https://content.iospress.com/journals/data-science/Pre-press/Pre-press), -->
<!-- which introduces the pema workflow, explains the underlying tools, and illustrates how the `pema` package can be used to create a new project that follows the workflow. -->

## Installing the package

Use [R-universe](https://cjvanlissa.r-universe.dev) to install the
development version of `pema` by running the following code:

``` r
options(repos = c(
    cjvanlissa = 'https://cjvanlissa.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))

install.packages('pema')
```

## Citing pema

The `pema` validation paper is currently in preprint stage. The code and
results of the validation study are publicly available. You can cite
`pema` using the following citation (please use the same citation for
either the package, or the paper); consider updating that citation once
the preprint is published:

> Van Lissa, C. J., & van Erp, S. (2021, December 9). Select relevant
> moderators using Bayesian regularized meta-regression. Retrieved from
> psyarxiv.com/6phs5

## About this repository

This repository contains the source code for the R-package called
`pema`.

## Contributing and Contact Information

We are always eager to receive user feedback and contributions to help
us improve both the workflow and the software. Major contributions
warrant coauthorship to the package. Please contact the lead author at
<c.j.vanlissa@uu.nl>, or:

-   [File a GitHub issue](https://github.com/cjvanlissa/pema) for
    feedback, bug reports or feature requests
-   [Make a pull request](https://github.com/cjvanlissa/pema/pulls) to
    contribute your code or prose

By participating in this project, you agree to abide by the [Contributor
Code of Conduct v2.0](https://www.contributor-covenant.org/).
Contributions to the package must adhere to the [tidyverse style
guide](https://style.tidyverse.org/). When contributing code, please add
tests for that contribution to the `tests/testthat` folder, and ensure
that these tests pass in the [GitHub Actions
panel](https://github.com/cjvanlissa/pema/actions/workflows/R-CMD-check).


<!-- README.md is generated from README.Rmd. Please edit that file -->

# SCENE

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

`SCENE` is a analytic tool for dissecting cellular cohabitation using single-cell and bulk RNA-seq data.

## System Requirements

### Hardware Requirements
The `SCENE` package requires computers with enough RAM to support the operations defined by a user.
Applying SCENE on 5000 samples (comprising single cells and bulk samples) requires about 12 GB of RAM and about 40 minutes to complete all calculation tasks on a workstation with 96 CPU threads. For larger cohorts, we recommend running the following functions on a high performance computer cluster and visualizing the results locally: `preprocessing2`, `determineK_runNMF`, `runcNMF`, and `identifyConsensusProgram3`.

### Software Requirements
The `SCENE` development version is developed and tested on *Linux* operating systems. The developmental version of the package has been tested on the following systems:

Linux: Ubuntu 20.04

Before setting up the `SCENE` package, users should have `R` version 4.0.5 or higher installed.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lijxug/SCENE")
```

Installation typically takes several minutes, depending on the installation of dependencies.

## Documentation
Full documentation is being drafting and will be online soon.
A demo case is currently available at [Code Ocean (Capsule DOI: 10.24433/CO.4624314.v1)](https://codeocean.com/capsule/7693770).
All the heavy calculations in this capsule have been performed previously.
Running this capsule online takes about 20 minutes to establish the environment and 5 minutes to print the figures.

## Release Notes
- version 1.0: First released version

## Lisence
The code is released under the GNU GPL-3 License.


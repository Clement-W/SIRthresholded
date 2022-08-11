<br/>
<h1 align="left">Variable selection with Sliced Inverse Regression (SIR) thresholded</h1>

<!-- badges: start -->

[![R-CMD-check](https://github.com/Clement-W/SIRthresholded/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Clement-W/SIRthresholded/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

## Description

This package offers an implementation of the SIR (Sliced Inverse Regression) method, along with a thresholded version of SIR that allows variable selection. Fore more information, you can check the vignette of the package on [that link](https://clement-w.github.io/SIRthresholded/articles/SIRthresholded.html)

The paper will soon be available.

## Install

### From CRAN

(not available yet)

### From GitHub

To install the current development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("clement-w/sirthresholded", build_vignettes = TRUE)
```

It is also possible to clone the repository, and install it manually (note: pandoc is required to build the vignette):

```sh
git clone git@github.com:Clement-W/SIRthresholded.git
cd SIRthresholded
R
```

```r
# install.packages("devtools")
devtools::install(build_vignettes = TRUE)
library(SIRthresholded)
```

### Access the vignette

Once that package is installed, you can access the vignette with

```r
vignette("SIRthresholded")
```

##

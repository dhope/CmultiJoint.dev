
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CmultiJoint.dev

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-0.0.2-blue.svg)](https://github.com/dhope/CmultiJoint.dev)
[![Project Status: Concept - Minimal or no implementation has been done
yet, or the repository is only intended to be a limited example, demo,
or
proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
<!-- badges: end -->

The goal of CmultiJoint.dev is to develop functions around bird density
modelling. See <https://github.com/davidiles/Bird_Detectability> for
details and usage cases.

## Installation

You can install the development version of CmultiJoint.dev from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("dhope/CmultiJoint.dev")
```

## Usage

Basic usage is shown below. Data is created using the
`simulate_point_counts()` function. See
[data-raw](https://github.com/dhope/CmultiJoint.dev/blob/main/data-raw/data_SimulatedPointCounts.R)
for how it was created

``` r
library(CmultiJoint.dev)

nsurvey <-  nrow(CmultiJoint.dev::SimulatedPointCounts$rarray)
fit <- cmulti_fit_joint(Yarray = CmultiJoint.dev::SimulatedPointCounts$Yarray,
                        rarray =CmultiJoint.dev::SimulatedPointCounts$rarray,
                        tarray = CmultiJoint.dev::SimulatedPointCounts$tarray,
                           X1 =NULL, # Design matrix for tau
                           X2 =NULL, # Design matrix for phi
)
log_offsetscpp <- 
  calculate_offsets(fit,
                    Yarray = CmultiJoint.dev::SimulatedPointCounts$Yarray,
                    rarray =CmultiJoint.dev::SimulatedPointCounts$rarray,
                    tarray = CmultiJoint.dev::SimulatedPointCounts$tarray,
                    X1 = NULL,
                    X2 = NULL)
```

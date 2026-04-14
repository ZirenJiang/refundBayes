---
---
---

# Package: <a href="https://zirenjiang.github.io/refundBayes/">refundBayes</a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/ZirenJiang/refundBayes/actions/workflows/rhub.yaml/badge.svg)](https://github.com/ZirenJiang/refundBayes/actions) [![CRAN status](https://www.r-pkg.org/badges/version/refundBayes)](https://cran.r-project.org/package=refundBayes) [![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/refundBayes)](https://cran.r-project.org/package=refundBayes)

<!-- badges: end -->

**refundBayes** provides a convenient interface for Bayesian functional regression using [Stan](https://mc-stan.org/). The package supports models with scalar, functional, or survival outcomes and allows both scalar and functional predictors. Its formula syntax mirrors that of [`mgcv::gam`](https://cran.r-project.org/package=mgcv), making it accessible to users familiar with the frequentist framework while providing full Bayesian posterior inference.

**Documentation website:** <https://zirenjiang.github.io/refundBayes/>

## Supported Models

| Model                                  | Function       | Response                 | Predictors          |
|:----------------------|:----------------|:----------------|:----------------|
| Scalar-on-Function Regression (SoFR)   | `sofr_bayes()` | Scalar                   | Functional / Scalar |
| Function-on-Scalar Regression (FoSR)   | `fosr_bayes()` | Functional               | Scalar              |
| Function-on-Function Regression (FoFR) | `fofr_bayes()` | Functional               | Functional / Scalar |
| Functional Cox Regression (FCox)       | `fcox_bayes()` | Survival (time-to-event) | Functional / Scalar |

## Installation

Install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("ZirenJiang/refundBayes")
```

## Quick Start

``` r
library(refundBayes)
data(example_data_sofr)
# Bayesian Scalar-on-Function Regression
fit <- sofr_bayes(
  y ~ X1 + s(tmat, by = lmat * wmat, bs = "cr", k = 10),
  data    = example_data_sofr,
  family  = gaussian(),
  niter   = 2000,
  nwarmup = 1000,
  nchain  = 1
)

summary(fit)
plot(fit)
```

## Tutorials

Detailed vignettes with full model descriptions, prior specifications, and worked examples are available on the [package website](https://zirenjiang.github.io/refundBayes/articles/):

-   [Bayesian Scalar-on-Function Regression (SoFR)](https://zirenjiang.github.io/refundBayes/articles/sofr_bayes_vignette.html)
-   [Bayesian Function-on-Scalar Regression (FoSR)](https://zirenjiang.github.io/refundBayes/articles/fosr_bayes_vignette.html)
-   [Bayesian Function-on-Function Regression (FoFR)](https://zirenjiang.github.io/refundBayes/articles/fofr_bayes_vignette.html)
-   [Bayesian Functional Cox Regression (FCox)](https://zirenjiang.github.io/refundBayes/articles/fcox_bayes_vignette.html)

## Citation

If you use **refundBayes** in your work, please cite:

> Jiang, Z., Crainiceanu, C., and Cui, E. (2025). Tutorial on Bayesian Functional Regression Using Stan. *Statistics in Medicine*, 44(20--22), e70265. <https://doi.org/10.1002/sim.70265>

``` bibtex
@article{jiang2025tutorial,
  title   = {Tutorial on {B}ayesian Functional Regression Using {S}tan},
  author  = {Jiang, Ziren and Crainiceanu, Ciprian and Cui, Erjia},
  journal = {Statistics in Medicine},
  volume  = {44},
  number  = {20--22},
  pages   = {e70265},
  year    = {2025},
  doi     = {10.1002/sim.70265}
}
```

## Related Packages

-   [refund](https://cran.r-project.org/package=refund) -- Frequentist regression with functional data
-   [mgcv](https://cran.r-project.org/package=mgcv) -- Generalized additive models
-   [rstan](https://mc-stan.org/rstan/) -- R interface to Stan
-   [brms](https://paul-buerkner.github.io/brms/) -- Bayesian regression models using Stan

## Contact

Questions, bug reports, and feature requests are welcome. For bugs and feature requests, the preferred channel is the [GitHub issue tracker](https://github.com/ZirenJiang/refundBayes/issues), so that the discussion is visible to other users. For general questions, collaboration inquiries, or anything not suitable for a public issue, please contact the maintainer:

-   **Ziren Jiang** (maintainer) --- [jian0746\@umn.edu](mailto:jian0746@umn.edu)

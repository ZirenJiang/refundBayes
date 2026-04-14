# refundBayes 0.5.1

## New features

* Added `fofr_bayes()` for Bayesian Function-on-Function Regression (FoFR),
  supporting functional responses with functional and scalar predictors. The
  bivariate coefficient surface β(s, t) is represented via a tensor-product
  basis with dual-direction smoothness (random-effect reparameterisation in
  the predictor direction and a penalty-matrix prior in the response
  direction).
* Added a dedicated vignette *Bayesian Function-on-Function Regression
  (FoFR)* with the full model specification, prior table, and a worked
  example.
* Expanded the pkgdown site to include the new FoFR reference entry and
  vignette.

## Documentation

* Rewrote `README.md`: added a supported-models table, links to per-function
  vignettes, a citation to Jiang et al. (2025, *Statistics in Medicine*),
  and CRAN status / downloads badges.
* Minor fixes to `fcox_bayes()` examples so that pkgdown can parse and
  render the reference page.

## Internal

* Added a standalone Stan program (`Simulation/StanFoFR_Gaussian.stan`) and
  a formal simulation script (`Simulation/FoFR_Simulation.R`) for
  reproducible FoFR benchmarking without recompiling Stan code via
  `refundBayes` at every run.

# refundBayes 0.5.0

Initial public release of **refundBayes**, a package providing a convenient
interface for Bayesian functional regression using Stan. The package is
designed to mirror the `mgcv::gam` formula syntax familiar to users of
[refund](https://cran.r-project.org/package=refund), while delivering full
Bayesian posterior inference via [rstan](https://mc-stan.org/rstan/).

## Supported models

* `sofr_bayes()` — Bayesian Scalar-on-Function Regression, supporting
  Gaussian, binomial, and Poisson families, with one or more functional
  predictors alongside scalar covariates.
* `fosr_bayes()` — Bayesian Function-on-Scalar Regression with FPCA-based
  residual structure for modelling subject-level functional deviations.
* `fcox_bayes()` — Bayesian Functional Cox Regression for time-to-event
  outcomes with functional and scalar predictors, including posterior
  inference on the log-hazard ratio surface.

## Modelling framework

* Functional coefficients represented via spline bases constructed through
  `mgcv::smoothCon()`, with spectral reparameterisation
  (`mgcv::smooth2random()`) into fixed and random effect components.
* Non-centered parameterisation of penalised spline coefficients for
  efficient HMC sampling.
* Inverse-gamma priors on smoothing variance components and weakly
  informative priors on fixed-effect coefficients.
* Support for multiple spline bases (cubic regression, thin-plate, P-spline,
  etc.) via the `bs` argument in the formula interface.
* Generic `summary()` and `plot()` methods for posterior summaries and
  visualisation of functional coefficients with credible bands.

## Documentation and infrastructure

* Per-function vignettes illustrating model specification, prior choices,
  fitting, and post-processing for SoFR, FoSR, and FCox.
* Example datasets (`example_data_sofr`, `example_data_FoSR`,
  `example_data_Cox`) shipped with the package.
* pkgdown documentation site published at
  <https://zirenjiang.github.io/refundBayes/>.
* Companion tutorial paper: Jiang, Crainiceanu, and Cui (2025),
  *Tutorial on Bayesian Functional Regression Using Stan*,
  Statistics in Medicine, 44(20–22), e70265,
  <doi:10.1002/sim.70265>.

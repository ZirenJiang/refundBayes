# Plot the estimated functional coefficients with the corresponding credible interval(s).

Produces coefficient plots tailored to the model family:

- [`sofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/sofr_bayes.md),
  [`fcox_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fcox_bayes.md):
  one curve plot per functional predictor coefficient \\\beta(s)\\, with
  pointwise and/or CMA credible bands.

- [`fosr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fosr_bayes.md):
  one curve plot per scalar predictor coefficient function
  \\\alpha_p(t)\\.

- [`fofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fofr_bayes.md):
  curve plots for scalar predictor coefficient functions \\\alpha_p(t)\\
  (if any), followed by heatmap plots for each bivariate coefficient
  surface \\\beta_q(s, t)\\ (posterior mean) from the functional
  predictors.

- [`fpca_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fpca_bayes.md):
  posterior mean function \\\mu(t)\\ with a pointwise credible band; a
  combined plot of the (fixed) FPC eigenfunctions \\\phi_j(t)\\; a
  point-and-error-bar plot of the posterior of the eigenvalue SDs
  \\\lambda_j\\; and a histogram of the residual-SD posterior
  \\\sigma\_\epsilon\\.

## Usage

``` r
# S3 method for class 'refundBayes'
plot(x = NULL, ..., prob = 0.95, include = "both")
```

## Arguments

- x:

  A fitted object returned by
  [`sofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/sofr_bayes.md),
  [`fosr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fosr_bayes.md),
  [`fofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fofr_bayes.md),
  [`fcox_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fcox_bayes.md),
  or
  [`fpca_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fpca_bayes.md).

- ...:

  Other parameters

- prob:

  Coverage probability for the credible interval(s). Defaults to 0.95.

- include:

  Type of interval to include. `"pointwise"` produces pointwise credible
  intervals; `"CMA"` produces the CMA credible band; `"both"` produces
  both. Defaults to `"both"`. Only used for
  [`sofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/sofr_bayes.md)
  /
  [`fcox_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fcox_bayes.md)
  curve plots.

## Value

A named list of `ggplot` objects. For FoFR, scalar-predictor curves are
named `scalar_<p>` and bivariate-predictor heatmaps are named
`bivar_<q>`. For FPCA, the plots are named `mu`, `efunctions`,
`evalues`, and `sigma`.

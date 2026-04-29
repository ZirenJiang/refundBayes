# Bayesian Functional Cox Regression

## Introduction

This vignette provides a detailed guide to the
[`fcox_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fcox_bayes.md)
function in the `refundBayes` package, which fits Bayesian Functional
Cox Regression (FCR) models for time-to-event outcomes using Stan. The
function extends the scalar-on-function regression framework to survival
analysis with right censoring, and is designed with a syntax similar to
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html).

The methodology follows the tutorial by Jiang, Crainiceanu, and Cui
(2025), *Tutorial on Bayesian Functional Regression Using Stan*,
published in *Statistics in Medicine*.

## Install the `refundBayes` Package

The `refundBayes` package can be installed from CRAN:

``` r
install.packages("refundBayes")
```

For the latest version of the `refundBayes` package, users can install
from GitHub:

``` r
library(remotes)
remotes::install_github("https://github.com/ZirenJiang/refundBayes")
```

## Statistical Model

### The Functional Cox Regression Model

The Functional Cox Regression (FCR) model extends the classical Cox
proportional hazards model to settings where one or more predictors are
functional (i.e., curves or trajectories observed over a continuum).

For subject $i = 1,\ldots,n$, let $T_{i}$ be the event time and $C_{i}$
the censoring time. The observed data consist of
$\left\lbrack Y_{i},\delta_{i},\mathbf{Z}_{i},\{ W_{i}\left( t_{m} \right),t_{m} \in \mathcal{T}\} \right\rbrack$,
where:

- $Y_{i} = \min\left( T_{i},C_{i} \right)$ is the observed follow-up
  time,
- $\delta_{i}$ is the censoring indicator with $\delta_{i} = 0$ if the
  event is observed ($Y_{i} = T_{i}$) and $\delta_{i} = 1$ if the event
  is censored ($Y_{i} < T_{i}$),
- $\mathbf{Z}_{i}$ is a $p \times 1$ vector of scalar predictors,
- $\{ W_{i}\left( t_{m} \right),t_{m} \in \mathcal{T}\}$ for
  $m = 1,\ldots,M$ is a functional predictor observed over a domain
  $\mathcal{T}$.

The domain $\mathcal{T}$ is **not** restricted to $\lbrack 0,1\rbrack$;
it is determined by the actual time points in the data (e.g.,
$\mathcal{T} = \lbrack 1,1440\rbrack$ for minute-level 24-hour activity
data). See the section on `tmat`, `lmat`, and `wmat` below for details.

The model assumes a proportional hazards structure:

$$h_{i}(t) = h_{0}(t)\exp\left( \eta_{i} \right)$$

where $h_{0}(t)$ is the baseline hazard function, and $\eta_{i}$ is the
linear predictor defined as:

$$\eta_{i} = \eta_{0} + \int_{\mathcal{T}}W_{i}(s)\beta(s)\, ds + \mathbf{Z}_{i}^{t}{\mathbf{γ}}$$

Here $\eta_{0}$ is the intercept, $\beta( \cdot )$ is the functional
coefficient, and $\mathbf{γ}$ is the vector of scalar regression
coefficients. The integral is approximated by a Riemann sum using the
time point matrix `tmat` and the integration weight matrix `lmat` (see
below).

The log-likelihood for this model has the form:

$$\ell\left( \mathbf{Y},{\mathbf{δ}};h_{0},{\mathbf{η}} \right) = \sum\limits_{i = 1}^{n}\left\lbrack \left( 1 - \delta_{i} \right)\left\{ \log h_{0}\left( y_{i} \right) + \eta_{i} - H_{0}\left( y_{i} \right)\exp\left( \eta_{i} \right) \right\} + \delta_{i}\left\{ - H_{0}\left( y_{i} \right)\exp\left( \eta_{i} \right) \right\} \right\rbrack$$

where $H_{0}(t) = \int_{0}^{t}h_{0}(u)\, du$ is the cumulative baseline
hazard function.

### Modeling the Baseline Hazard with M-Splines

Unlike the classical Cox model, which leaves the baseline hazard
unspecified and uses partial likelihood, the Bayesian approach requires
a full likelihood specification. The
[`fcox_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fcox_bayes.md)
function models the baseline hazard function $h_{0}(t)$ using M-splines:

$$h_{0}(t) = \sum\limits_{l = 1}^{L}c_{l}M_{l}(t;\mathbf{k},\tau)$$

where $M_{l}(t;\mathbf{k},\tau)$ denotes the $l$-th M-spline basis
function with knots $\mathbf{k}$ and degree $\tau$, and
$\mathbf{c} = \left( c_{1},\ldots,c_{L} \right)^{t}$ are the spline
coefficients. The constraints $c_{l} \geq 0$ and
$\sum_{l = 1}^{L}c_{l} > 0$ ensure that the baseline hazard is
non-negative and non-trivial. M-splines are a family of non-negative,
piecewise polynomial basis functions that integrate to one over their
support.

The cumulative baseline hazard function is modeled using I-splines,
which are the integrated forms of M-splines:

$$H_{0}(t) = \sum\limits_{l = 1}^{L}c_{l}I_{l}(t;\mathbf{k},\tau)$$

where
$I_{l}(t;\mathbf{k},\tau) = \int_{0}^{t}M_{l}(u;\mathbf{k},\tau)\, du$.
Because M-splines are non-negative, the resulting I-spline cumulative
hazard is non-decreasing by construction.

### Identifiability Constraint

The I-spline coefficients $\mathbf{c}$ and the intercept $\eta_{0}$ are
not simultaneously identifiable: for any $a > 0$, the parameter pairs
$\left( \mathbf{c},\eta_{0} \right)$ and
$\left( a\mathbf{c},\eta_{0} - \log a \right)$ yield the same hazard
function. To resolve this, the model imposes the constraint
$\sum_{l = 1}^{L}c_{l} = 1$ by assigning a non-informative Dirichlet
prior $D(\mathbf{c};{\mathbf{α}})$ with ${\mathbf{α}} = (1,\ldots,1)$ to
the coefficients. Consequently, the intercept defaults to `FALSE` in
[`fcox_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fcox_bayes.md).

### Functional Coefficient via Penalized Splines

As in the SoFR model, the functional coefficient $\beta(s)$ is
represented using $K$ spline basis functions
$\psi_{1}(s),\ldots,\psi_{K}(s)$:

$$\beta(s) = \sum\limits_{k = 1}^{K}b_{k}\psi_{k}(s)$$

The spline coefficients $\mathbf{b}$ are reparametrised using the
spectral decomposition of the penalty matrix $\mathbf{S}$ (the
Wahba-O’Sullivan smoothing penalty) to obtain transformed coefficients
$\widetilde{\mathbf{b}} = \left( {\widetilde{\mathbf{b}}}_{r}^{t},{\widetilde{\mathbf{b}}}_{f}^{t} \right)^{t}$,
where:

- ${\widetilde{\mathbf{b}}}_{r}$ (random effects) receive a
  $N\left( \mathbf{0},\sigma_{b}^{2}\mathbf{I} \right)$ prior,
- ${\widetilde{\mathbf{b}}}_{f}$ (fixed effects) receive non-informative
  priors.

This reparametrisation ensures numerical stability and efficient
sampling in Stan.

#### The Role of `tmat`, `lmat`, and `wmat`

The integral $\int_{\mathcal{T}}W_{i}(s)\beta(s)\, ds$ is approximated
by the Riemann sum
$\sum_{m = 1}^{M}L_{m}W_{i}\left( t_{m} \right)\psi_{k}\left( t_{m} \right)$,
constructed from three user-supplied matrices:

- **`tmat`** (an $n \times M$ matrix): contains the time points $t_{m}$
  at which the functional predictor is observed. The $(i,m)$-th entry
  equals $t_{m}$. The range of values in `tmat` determines the domain of
  integration $\mathcal{T}$. For example, if the functional predictor is
  observed at minutes $1,2,\ldots,1440$ within a day, then `tmat` has
  entries ranging from $1$ to $1440$ and
  $\mathcal{T} = \lbrack 1,1440\rbrack$. There is no requirement that
  the domain be rescaled to $\lbrack 0,1\rbrack$.

- **`lmat`** (an $n \times M$ matrix): contains the integration weights
  $L_{m} = t_{m + 1} - t_{m}$ for the Riemann sum approximation. For
  equally spaced time points with spacing $\Delta t$, every entry of
  `lmat` equals $\Delta t$. These weights, together with `tmat`, fully
  specify how the numerical integration is performed and over what
  domain.

- **`wmat`** (an $n \times M$ matrix): contains the functional predictor
  values. The $i$-th row contains the $M$ observed values
  $W_{i}\left( t_{1} \right),\ldots,W_{i}\left( t_{M} \right)$ for
  subject $i$.

In the formula `s(tmat, by = lmat * wmat, bs = "cc", k = 10)`, the
`mgcv` infrastructure uses `tmat` to construct the spline basis
$\psi_{k}\left( t_{m} \right)$ at the observed time points, and the
`by = lmat * wmat` argument provides the element-wise product
$L_{m} \cdot W_{i}\left( t_{m} \right)$ that enters the Riemann sum. The
basis functions are evaluated on the scale of `tmat`, and the
integration weights in `lmat` ensure that the discrete sum correctly
approximates the integral over the actual domain $\mathcal{T}$ —
regardless of whether it is $\lbrack 0,1\rbrack$,
$\lbrack 1,1440\rbrack$, or any other interval.

Note that for all subjects, the time points are assumed to be identical
so that $t_{im} = t_{m}$ for all $i = 1,\ldots,n$. Thus every row of
`tmat` is the same, and every row of `lmat` is the same. The matrices
are replicated across rows to match the `mgcv` syntax, which expects all
terms in the formula to have the same dimensions.

### Full Bayesian Model

The complete Bayesian Functional Cox Regression model is:

$$\left\{ \begin{array}{l}
{\mathbf{Y} \sim \ell\left( \mathbf{Y},{\mathbf{δ}};h_{0},{\mathbf{η}} \right)} \\
{{\mathbf{η}} = \eta_{0}\mathbf{J}_{n} + {\widetilde{\mathbf{X}}}_{r}^{t}{\widetilde{\mathbf{b}}}_{r} + {\widetilde{\mathbf{X}}}_{f}^{t}{\widetilde{\mathbf{b}}}_{f} + \mathbf{Z}^{t}{\mathbf{γ}}} \\
{{\widetilde{\mathbf{b}}}_{r} \sim N\left( \mathbf{0},\sigma_{b}^{2}\mathbf{I} \right)} \\
{\eta_{0} \sim p\left( \eta_{0} \right);\;{\widetilde{\mathbf{b}}}_{f} \sim p\left( {\widetilde{\mathbf{b}}}_{f} \right);\;{\mathbf{γ}} \sim p({\mathbf{γ}});\;\sigma_{b}^{2} \sim p\left( \sigma_{b}^{2} \right)} \\
{h_{0}(t) = \sum\limits_{l = 1}^{L}c_{l}M_{l}(t;\mathbf{k},\tau)} \\
{\mathbf{c} \sim D(\mathbf{c};{\mathbf{α}})} \\
\end{array} \right.$$

where most components are shared with the SoFR model, except for the Cox
likelihood (first line) and the baseline hazard specification via
M-splines with a Dirichlet prior on $\mathbf{c}$ (last two lines).

### Optional: Joint FPCA Modeling

When functional predictors are observed with measurement error, the
[`fcox_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fcox_bayes.md)
function supports a joint modeling approach that simultaneously
estimates FPCA scores and fits the Cox regression model. In this case,
the model regresses on the latent trajectory $D_{i}(s)$ rather than the
observed (noisy) function $W_{i}(s) = D_{i}(s) + \epsilon_{i}(s)$,
properly accounting for the uncertainty in the score estimates and
avoiding the errors-in-variables attenuation that follows from plugging
$W_{i}(s)$ directly into the integral $\int W_{i}(s)\beta(s)\, ds$.

The joint FPCA option is enabled via the `joint_FPCA` argument and is
shared with
[`sofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/sofr_bayes.md)
and
[`fofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fofr_bayes.md).
The full model specification (FPCA likelihood for the predictor,
FPC-score prior centered on the initial
[`refund::fpca.sc()`](https://rdrr.io/pkg/refund/man/fpca.sc.html)
scores, joint Stan code, and the resulting Cox-with-joint-FPCA program),
together with a worked survival example, is presented in the dedicated
**[Joint FPCA
vignette](https://zirenjiang.github.io/refundBayes/articles/joint_FPCA_vignette.md)**.

## The `fcox_bayes()` Function

### Usage

``` r
fcox_bayes(
  formula,
  data,
  cens,
  joint_FPCA = NULL,
  intercept = FALSE,
  runStan = TRUE,
  niter = 3000,
  nwarmup = 1000,
  nchain = 3,
  ncores = 1
)
```

### Arguments

| Argument     | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|:-------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `formula`    | Functional regression formula, using the same syntax as [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html). The left-hand side should be the observed survival time $Y_{i} = \min\left( T_{i},C_{i} \right)$. Functional predictors are specified using the [`s()`](https://rdrr.io/pkg/mgcv/man/s.html) term with `by = lmat * wmat` to encode the Riemann sum integration.                                                                                                                                                                                                                |
| `data`       | A data frame containing all scalar and functional variables used in the model, as well as the survival time as the response variable.                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| `cens`       | A binary vector indicating censoring status for each subject, following the convention: $\delta_{i} = 0$ if the event is observed and $\delta_{i} = 1$ if censored. Must have the same length as the number of observations in `data`.                                                                                                                                                                                                                                                                                                                                                        |
| `joint_FPCA` | A logical (`TRUE`/`FALSE`) vector of the same length as the number of functional predictors, indicating whether to jointly model FPCA for each functional predictor. When `TRUE`, the observed functional predictor is replaced by an FPCA representation and its FPC scores are sampled jointly with the survival model (errors-in-variables-aware fit). See the [Joint FPCA vignette](https://zirenjiang.github.io/refundBayes/articles/joint_FPCA_vignette.md) for the model specification and a worked Cox example. Default is `NULL`, which sets all entries to `FALSE` (no joint FPCA). |
| `intercept`  | Logical. Whether to include an intercept term in the linear predictor. Default is `FALSE`, because the intercept is not identifiable simultaneously with the M-spline coefficients for the baseline hazard (see the Identifiability Constraint section above).                                                                                                                                                                                                                                                                                                                                |
| `runStan`    | Logical. Whether to run the Stan program. If `FALSE`, the function only generates the Stan code and data without sampling. Default is `TRUE`.                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `niter`      | Total number of Bayesian posterior sampling iterations (including warmup). Default is `3000`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `nwarmup`    | Number of warmup (burn-in) iterations. Default is `1000`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| `nchain`     | Number of Markov chains for posterior sampling. Default is `3`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `ncores`     | Number of CPU cores to use when executing the chains in parallel. Default is `1`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |

### Return Value

The function returns a list of class `"refundBayes"` containing the
following elements:

| Element           | Description                                                                                                                                                                                                  |
|:------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `stanfit`         | The Stan fit object (class `stanfit`). Can be used for convergence diagnostics, traceplots, and additional summaries via the `rstan` package.                                                                |
| `spline_basis`    | Basis functions used to reconstruct the functional coefficients from the posterior samples.                                                                                                                  |
| `stancode`        | A character string containing the generated Stan model code.                                                                                                                                                 |
| `standata`        | A list containing the data passed to the Stan model.                                                                                                                                                         |
| `int`             | A vector of posterior samples for the intercept term $\eta_{0}$. `NULL` by default since `intercept = FALSE`.                                                                                                |
| `scalar_coef`     | A matrix of posterior samples for scalar coefficients $\mathbf{γ}$, where each row is one posterior sample and each column corresponds to one scalar predictor. `NULL` if no scalar predictors are included. |
| `func_coef`       | A list of posterior samples for functional coefficients. Each element is a matrix where each row is one posterior sample and each column corresponds to one location on the functional domain.               |
| `baseline_hazard` | A list containing posterior samples of baseline hazard parameters: `bhaz` (baseline hazard values) and `cbhaz` (cumulative baseline hazard values).                                                          |
| `family`          | The family type: `"Cox"`.                                                                                                                                                                                    |

### Formula Syntax

The formula follows the
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) syntax. The
left-hand side is the observed survival time (not a `Surv` object). The
key component for specifying functional predictors is:

``` r
s(tmat, by = lmat * wmat, bs = "cc", k = 10)
```

where:

- `tmat`: an $n \times M$ matrix of time points defining the domain
  $\mathcal{T}$ over which the functional predictor is observed and the
  integral is computed. The domain adapts to the actual values in `tmat`
  (e.g., $\lbrack 0,1\rbrack$, $\lbrack 1,1440\rbrack$, etc.).
- `lmat`: an $n \times M$ matrix of Riemann sum weights
  ($L_{m} = t_{m + 1} - t_{m}$), controlling the numerical integration
  over $\mathcal{T}$.
- `wmat`: an $n \times M$ matrix of functional predictor values.
- `bs`: the type of spline basis (e.g., `"cr"` for cubic regression
  splines, `"cc"` for cyclic cubic regression splines).
- `k`: the number of basis functions.

Scalar predictors are included as standard formula terms. The censoring
information is provided separately via the `cens` argument, not through
the formula.

## Example: Bayesian Functional Cox Regression

We demonstrate the
[`fcox_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fcox_bayes.md)
function using a simulated example dataset with a time-to-event outcome
and one functional predictor.

### Load and Prepare Data

``` r
## Load the example data
Func_Cox_Data <- readRDS("data/example_data_Cox.rds")

## Prepare the functional predictor and censoring indicator
Func_Cox_Data$wmat <- Func_Cox_Data$MIMS
Func_Cox_Data$cens <- 1 - Func_Cox_Data$event
```

The example dataset contains:

- `survtime`: the observed survival time
  $Y_{i} = \min\left( T_{i},C_{i} \right)$,
- `event`: a binary event indicator ($1$ = event occurred, $0$ =
  censored),
- `X1`: a scalar predictor,
- `tmat`: the $n \times M$ time point matrix (defines the domain
  $\mathcal{T}$),
- `lmat`: the $n \times M$ Riemann sum weight matrix (defines the
  integration weights over $\mathcal{T}$),
- `MIMS`: the $n \times M$ functional predictor matrix (physical
  activity data).

Note that the censoring indicator `cens` is constructed as `1 - event`,
following the convention where $\delta_{i} = 0$ indicates an observed
event and $\delta_{i} = 1$ indicates censoring.

### Fit the Bayesian Functional Cox Model

``` r
library(refundBayes)

refundBayes_FCox <- refundBayes::fcox_bayes(
  survtime ~ X1 + s(tmat, by = lmat * wmat, bs = "cc", k = 10),
  data = Func_Cox_Data,
  cens = Func_Cox_Data$cens,
  runStan = TRUE,
  niter = 1500,
  nwarmup = 500,
  nchain = 3,
  ncores = 3
)
```

In this call:

- The formula specifies the observed survival time `survtime` as the
  response, with one scalar predictor `X1` and one functional predictor
  encoded via `s(tmat, by = lmat * wmat, bs = "cc", k = 10)`.
- The spline basis is evaluated at the time points stored in `tmat`, and
  the integration over the domain $\mathcal{T}$ (determined by `tmat`)
  is approximated using the weights in `lmat`.
- The censoring vector `cens` is supplied separately, with $0$ = event
  observed and $1$ = censored.
- `bs = "cc"` uses cyclic cubic regression splines, appropriate for
  periodic functional predictors (e.g., 24-hour activity patterns).
- `k = 10` specifies 10 basis functions. In practice, 30–40 basis
  functions are recommended for moderately smooth functional data on
  dense grids.
- The `intercept` argument is not specified and defaults to `FALSE`,
  which is appropriate for Cox models due to the identifiability
  constraint with the baseline hazard.
- The sampler runs 3 chains in parallel, each with 1500 total iterations
  (500 warmup + 1000 posterior samples).

### Visualization

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
`refundBayes` objects displays the estimated functional coefficient
$\widehat{\beta}(t)$ along with pointwise 95% credible intervals:

``` r
library(ggplot2)
plot(refundBayes_FCox)
```

### Extracting Posterior Summaries

Posterior summaries of the functional coefficient can be computed
directly from the `func_coef` element:

``` r
## Posterior mean of the functional coefficient
mean_curve <- apply(refundBayes_FCox$func_coef[[1]], 2, mean)

## Pointwise 95% credible interval
upper_curve <- apply(refundBayes_FCox$func_coef[[1]], 2,
                     function(x) quantile(x, prob = 0.975))
lower_curve <- apply(refundBayes_FCox$func_coef[[1]], 2,
                     function(x) quantile(x, prob = 0.025))
```

The posterior samples in `func_coef[[1]]` are stored as a $Q \times M$
matrix, where $Q$ is the number of posterior samples and $M$ is the
number of time points on the functional domain.

### Comparison with Frequentist Results

The Bayesian results can be compared with frequentist estimates obtained
via [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html), which handles
Cox proportional hazards models through its
[`cox.ph()`](https://rdrr.io/pkg/mgcv/man/coxph.html) family using
partial likelihood:

``` r
library(mgcv)

## Fit frequentist functional Cox model using mgcv
fit_freq <- gam(
  survtime ~ s(tmat, by = lmat * wmat, bs = "cc", k = 10) + X1,
  data = Func_Cox_Data,
  family = cox.ph(),
  weights = Func_Cox_Data$event
)

## Extract frequentist estimates
freq_result <- plot(fit_freq)
```

Note that the frequentist approach uses
[`cox.ph()`](https://rdrr.io/pkg/mgcv/man/coxph.html) family with
partial likelihood and a Breslow-type nonparametric estimator for the
baseline hazard, while the Bayesian approach models the baseline hazard
parametrically using M-splines. Despite this difference, simulation
studies in Jiang et al. (2025) show that both approaches achieve
comparable performance in terms of estimation accuracy and coverage.

### Inspecting the Generated Stan Code

Setting `runStan = FALSE` allows you to inspect or modify the Stan code
before running the model:

``` r
## Generate Stan code without running the sampler
fcox_code <- refundBayes::fcox_bayes(
  survtime ~ X1 + s(tmat, by = lmat * wmat, bs = "cc", k = 10),
  data = Func_Cox_Data,
  cens = Func_Cox_Data$cens,
  runStan = FALSE
)

## Print the generated Stan code
cat(fcox_code$stancode)
```

The generated Stan code includes a `functions` block that defines custom
log-likelihood functions for the Cox model (`cox_log_lpdf` for observed
events and `cox_log_lccdf` for censored observations), in addition to
the standard `data`, `parameters`, and `model` blocks.

## Practical Recommendations

- **Number of basis functions (`k`)**: For illustrative purposes,
  `k = 10` is often used. In practice, 30–40 basis functions are
  recommended for moderately smooth functional data observed on dense
  grids.
- **Spline type (`bs`)**: Use `"cr"` (cubic regression splines) for
  general functional predictors. Use `"cc"` (cyclic cubic regression
  splines) when the functional predictor is periodic (e.g., 24-hour
  activity patterns).
- **Intercept**: The intercept defaults to `FALSE` for Cox models
  because it is not simultaneously identifiable with the baseline hazard
  M-spline coefficients. In most cases, this default should not be
  changed.
- **Censoring vector**: Ensure the `cens` vector correctly follows the
  convention $\delta_{i} = 0$ for observed events and $\delta_{i} = 1$
  for censored observations. If your data has an event indicator (1 =
  event, 0 = censored), use `cens = 1 - event`.
- **Number of iterations**: Cox models may require more iterations than
  standard SoFR models. A common setup is `niter = 3000` with
  `nwarmup = 1000`. Watch for divergent transitions in the Stan output;
  if they occur, consider increasing the number of warmup iterations or
  adjusting the `adapt_delta` parameter.
- **Convergence diagnostics**: After fitting, examine traceplots and
  $\widehat{R}$ statistics using the `rstan` package (e.g.,
  `rstan::traceplot(refundBayes_FCox$stanfit)`) to ensure that the
  Markov chains have converged.
- **Joint FPCA**: When functional predictors are measured with
  substantial noise, consider setting `joint_FPCA = TRUE` for the
  relevant predictor to jointly estimate FPCA scores and regression
  coefficients within the Cox model. See the [Joint FPCA
  vignette](https://zirenjiang.github.io/refundBayes/articles/joint_FPCA_vignette.md)
  for details.

## References

- Jiang, Z., Crainiceanu, C., and Cui, E. (2025). Tutorial on Bayesian
  Functional Regression Using Stan. *Statistics in Medicine*, 44(20–22),
  e70265.
- Crainiceanu, C. M., Goldsmith, J., Leroux, A., and Cui, E. (2024).
  *Functional Data Analysis with R*. CRC Press.
- Cui, E., Crainiceanu, C. M., and Leroux, A. (2021). Additive
  Functional Cox Model. *Journal of Computational and Graphical
  Statistics*, 30(3), 780–793.
- Brilleman, S. L., Elci, E. M., Novik, J. B., and Wolfe, R. (2020).
  Bayesian Survival Analysis Using the rstanarm R Package. arXiv
  preprint, arXiv:2002.09633.
- Cox, D. (1972). Regression Models and Life-Tables. *Journal of the
  Royal Statistical Society, Series B*, 34(2), 187–220.
- Carpenter, B., Gelman, A., Hoffman, M. D., et al. (2017). Stan: A
  Probabilistic Programming Language. *Journal of Statistical Software*,
  76(1), 1–32.

# Bayesian Scalar-on-Function Regression with \`refundBayes::sofr_bayes\`

## Introduction

This vignette provides a detailed guide to the
[`sofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/sofr_bayes.md)
function in the `refundBayes` package, which fits Bayesian
Scalar-on-Function Regression (SoFR) models using Stan. The function is
designed with a syntax similar to
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html), making it
accessible to users familiar with the frequentist approach while
providing full Bayesian posterior inference.

The methodology follows the tutorial by Jiang, Crainiceanu, and Cui
(2025), *Tutorial on Bayesian Functional Regression Using Stan*,
published in *Statistics in Medicine*.

## Install the `refundBayes` Package

The `refundBayes` package can be installed from GitHub:

``` r
library(remotes)
remotes::install_github("https://github.com/ZirenJiang/refundBayes")
```

## Statistical Model

### The SoFR Model

Scalar-on-Function Regression (SoFR) models the relationship between a
scalar outcome and one or more functional predictors (curves or
trajectories observed over a continuum), along with optional scalar
covariates.

For subject $i = 1,\ldots,n$, let $Y_{i}$ be the scalar outcome,
$\mathbf{Z}_{i}$ be a $p \times 1$ vector of scalar predictors, and
$\{ W_{i}\left( t_{m} \right),t_{m} \in \mathcal{T}\}$ for
$m = 1,\ldots,M$ be a functional predictor observed at $M$ time points
over a domain $\mathcal{T}$. The SoFR model assumes that the
distribution of $Y_{i}$ belongs to an exponential family with mean
$\mu_{i}$, and the linear predictor $\eta_{i} = g\left( \mu_{i} \right)$
has the following structure:

$$\eta_{i} = \eta_{0} + \int_{\mathcal{T}}W_{i}(t)\beta(t)\, dt + \mathbf{Z}_{i}^{t}{\mathbf{γ}}$$

where:

- $\eta_{0}$ is the overall intercept,
- $\beta( \cdot ) \in L^{2}(\mathcal{T})$ is the unknown functional
  coefficient that characterizes the effect of the functional predictor
  on the outcome,
- $\mathbf{γ}$ is a $p \times 1$ vector of scalar regression
  coefficients,
- $\mathcal{T}$ is the domain of the functional predictor. This domain
  is **not** restricted to $\lbrack 0,1\rbrack$; it is determined by the
  actual time points in the data (e.g.,
  $\mathcal{T} = \lbrack 1,1440\rbrack$ for minute-level 24-hour
  activity data, or $\mathcal{T} = \lbrack 0,100\rbrack$ for a generic
  index set).

The integral $\int_{\mathcal{T}}W_{i}(t)\beta(t)\, dt$ is approximated
using a Riemann sum over the observed time points.

### Basis Expansion and Penalized Splines

The functional coefficient $\beta(t)$ is represented nonparametrically
using a set of $K$ pre-specified spline basis functions
$\psi_{1}(t),\ldots,\psi_{K}(t)$:

$$\beta(t) = \sum\limits_{k = 1}^{K}b_{k}\psi_{k}(t)$$

With this expansion, the linear predictor becomes:

$$\eta_{i} = \eta_{0} + \mathbf{X}_{i}^{t}\mathbf{b} + \mathbf{Z}_{i}^{t}{\mathbf{γ}}$$

where $\mathbf{X}_{i} = \left( X_{i1},\ldots,X_{iK} \right)^{t}$ is a
$K \times 1$ vector with entries:

$$X_{ik} = \sum\limits_{m = 1}^{M}L_{m}W_{i}\left( t_{m} \right)\psi_{k}\left( t_{m} \right)$$

Here $L_{m} = t_{m + 1} - t_{m}$ are the Riemann sum integration weights
and $t_{m}$ are the time points at which the functional predictor is
observed.

#### The Role of `tmat`, `lmat`, and `wmat`

The Riemann sum approximation
$\sum_{m = 1}^{M}L_{m}W_{i}\left( t_{m} \right)\psi_{k}\left( t_{m} \right)$
to the integral $\int_{\mathcal{T}}W_{i}(t)\psi_{k}(t)\, dt$ is
constructed directly from three user-supplied matrices in the `data`
argument:

- **`tmat`** (an $n \times M$ matrix): contains the time points $t_{m}$
  at which the functional predictor is observed. The $(i,m)$-th entry
  equals $t_{m}$. The range of values in `tmat` determines the domain of
  integration $\mathcal{T}$. For example, if the functional predictor is
  observed at minutes $1,2,\ldots,1440$ within a day, then `tmat` has
  entries ranging from $1$ to $1440$ and
  $\mathcal{T} = \lbrack 1,1440\rbrack$. There is no requirement that
  the domain be rescaled to $\lbrack 0,1\rbrack$.

- **`lmat`** (an $n \times M$ matrix): contains the integration weights
  $L_{m} = t_{m + 1} - t_{m}$ for the Riemann sum approximation. The
  $(i,m)$-th entry equals $L_{m}$. For equally spaced time points with
  spacing $\Delta t$, every entry of `lmat` equals $\Delta t$. For
  unevenly spaced time points, `lmat` reflects the varying widths of the
  integration intervals. These weights, together with `tmat`, fully
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
$L_{m} \cdot W_{i}\left( t_{m} \right)$ that enters the Riemann sum.
This means the basis functions are evaluated on the scale of `tmat`, and
the integration weights in `lmat` ensure that the discrete sum correctly
approximates the integral over the actual domain $\mathcal{T}$ —
regardless of whether it is $\lbrack 0,1\rbrack$,
$\lbrack 1,1440\rbrack$, or any other interval.

Note that for all subjects, the time points are assumed to be identical
so that $t_{im} = t_{m}$ for all $i = 1,\ldots,n$. Thus every row of
`tmat` is the same, and every row of `lmat` is the same. The matrices
are replicated across rows to match the `mgcv` syntax, which expects all
terms in the formula to have the same dimensions.

### Smoothness Penalty

To induce smoothness on $\beta(t)$, a quadratic penalty on the spline
coefficients is applied. The penalty is based on the integrated squared
second derivative of $\beta(t)$:

$$\int\{\beta''(t)\}^{2}\, dt = \mathbf{b}^{t}\mathbf{S}\mathbf{b}$$

where $\mathbf{S} = \int{\mathbf{ψ}}''(t)\{{\mathbf{ψ}}''(t)\}^{t}\, dt$
is the penalty matrix. In the Bayesian framework, this penalty is
equivalent to placing a multivariate normal prior on the spline
coefficients:

$$p(\mathbf{b}) \propto \exp\left( - \frac{\mathbf{b}^{t}\mathbf{S}\mathbf{b}}{\sigma_{b}^{2}} \right)$$

where $\sigma_{b}^{2}$ is the smoothing parameter that controls the
smoothness of $\beta(t)$, and is estimated from the data.

### Full Bayesian Model

The complete Bayesian SoFR model is:

$$\left\{ \begin{array}{l}
{\mathbf{Y} \sim \text{Exponential\ Family}({\mathbf{η}},a)} \\
{{\mathbf{η}} = \eta_{0}\mathbf{J}_{n} + {\widetilde{\mathbf{X}}}_{r}^{t}{\widetilde{\mathbf{b}}}_{r} + {\widetilde{\mathbf{X}}}_{f}^{t}{\widetilde{\mathbf{b}}}_{f} + \mathbf{Z}^{t}{\mathbf{γ}}} \\
{{\widetilde{\mathbf{b}}}_{r} \sim N\left( \mathbf{0},\sigma_{b}^{2}\mathbf{I} \right)} \\
{\eta_{0} \sim p\left( \eta_{0} \right);\;{\widetilde{\mathbf{b}}}_{f} \sim p\left( {\widetilde{\mathbf{b}}}_{f} \right);\;{\mathbf{γ}} \sim p({\mathbf{γ}})} \\
{\sigma_{b}^{2} \sim p\left( \sigma_{b}^{2} \right);\; a \sim p(a)} \\
\end{array} \right.$$

where ${\widetilde{\mathbf{X}}}_{r}$ and ${\widetilde{\mathbf{X}}}_{f}$
are the correspondingly transformed design matrices, and $p( \cdot )$
denotes non-informative priors (uniform or weakly informative). The
smoothing variance $\sigma_{b}^{2}$ is assigned an inverse-Gamma prior
$IG(0.001,0.001)$.

## The `sofr_bayes()` Function

### Usage

``` r
sofr_bayes(
  formula,
  data,
  family = gaussian(),
  joint_FPCA = NULL,
  intercept = TRUE,
  runStan = TRUE,
  niter = 3000,
  nwarmup = 1000,
  nchain = 3,
  ncores = 1
)
```

### Arguments

| Argument     | Description                                                                                                                                                                                                                                                                                      |
|:-------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `formula`    | Functional regression formula, using the same syntax as [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html). Functional predictors are specified using the [`s()`](https://rdrr.io/pkg/mgcv/man/s.html) term with `by = lmat * wmat` to encode the Riemann sum integration (see Example below). |
| `data`       | A data frame containing all scalar and functional variables used in the model.                                                                                                                                                                                                                   |
| `family`     | Distribution of the outcome variable. Currently supports [`gaussian()`](https://rdrr.io/r/stats/family.html) and [`binomial()`](https://rdrr.io/r/stats/family.html). Default is [`gaussian()`](https://rdrr.io/r/stats/family.html).                                                            |
| `joint_FPCA` | A logical (`TRUE`/`FALSE`) vector of the same length as the number of functional predictors, indicating whether to jointly model FPCA for each functional predictor. Default is `NULL`, which sets all entries to `FALSE` (no joint FPCA).                                                       |
| `intercept`  | Logical. Whether to include an intercept term in the linear predictor. Default is `TRUE`.                                                                                                                                                                                                        |
| `runStan`    | Logical. Whether to run the Stan program. If `FALSE`, the function only generates the Stan code and data without sampling. This is useful for inspecting or modifying the generated Stan code. Default is `TRUE`.                                                                                |
| `niter`      | Total number of Bayesian posterior sampling iterations (including warmup). Default is `3000`.                                                                                                                                                                                                    |
| `nwarmup`    | Number of warmup (burn-in) iterations. These samples are discarded and not used for inference. Default is `1000`.                                                                                                                                                                                |
| `nchain`     | Number of Markov chains for posterior sampling. Multiple chains help assess convergence. Default is `3`.                                                                                                                                                                                         |
| `ncores`     | Number of CPU cores to use when executing the chains in parallel. Default is `1`.                                                                                                                                                                                                                |

### Return Value

The function returns a list of class `"refundBayes"` containing the
following elements:

| Element        | Description                                                                                                                                                                                                  |
|:---------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `stanfit`      | The Stan fit object (class `stanfit`). Can be used for convergence diagnostics, traceplots, and additional summaries via the `rstan` package.                                                                |
| `spline_basis` | Basis functions used to reconstruct the functional coefficients from the posterior samples.                                                                                                                  |
| `stancode`     | A character string containing the generated Stan model code.                                                                                                                                                 |
| `standata`     | A list containing the data passed to the Stan model.                                                                                                                                                         |
| `int`          | A vector of posterior samples for the intercept term $\eta_{0}$. `NULL` if `intercept = FALSE`.                                                                                                              |
| `scalar_coef`  | A matrix of posterior samples for scalar coefficients $\mathbf{γ}$, where each row is one posterior sample and each column corresponds to one scalar predictor. `NULL` if no scalar predictors are included. |
| `func_coef`    | A list of posterior samples for functional coefficients. Each element is a matrix where each row is one posterior sample and each column corresponds to one location on the functional domain.               |
| `family`       | The distribution family used for the outcome.                                                                                                                                                                |

### Formula Syntax

The formula follows the
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) syntax. The key
component for specifying functional predictors is:

``` r
s(tmat, by = lmat * wmat, bs = "cc", k = 10)
```

where:

- `tmat`: an $n \times M$ matrix of time points. Each row contains the
  same $M$ observation times (replicated across subjects). The values in
  `tmat` define the domain $\mathcal{T}$ over which the functional
  predictor is observed and the integral is computed. For example,
  `tmat` may contain values from $1$ to $100$, from $0$ to $1$, or from
  $1$ to $1440$ — the integration domain adapts accordingly.
- `lmat`: an $n \times M$ matrix of Riemann sum weights. The $(i,m)$-th
  entry equals $L_{m} = t_{m + 1} - t_{m}$. These weights control the
  numerical integration and should be consistent with the spacing of the
  time points in `tmat`. For equally spaced time points, `lmat` is a
  constant matrix; for irregular spacing, it reflects the actual gaps
  between consecutive time points.
- `wmat`: an $n \times M$ matrix of functional predictor values. The
  $i$-th row contains the $M$ observed values for subject $i$.
- `bs`: the type of spline basis. Common choices include `"cr"` (cubic
  regression splines) and `"cc"` (cyclic cubic regression splines,
  suitable for periodic data).
- `k`: the number of basis functions (maximum degrees of freedom for the
  spline).

Scalar predictors are included as standard formula terms (e.g.,
`y ~ X1 + s(tmat, by = lmat * wmat, bs = "cr", k = 10)`).

## Example: Bayesian SoFR with Binary Outcome

We demonstrate the
[`sofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/sofr_bayes.md)
function using a simulated example dataset with a binary outcome and one
functional predictor.

### Load Data

``` r
## Alternative sample data 

# data.SoFR <- readRDS("data/example_data_sofr.rds")

## Load the example data
set.seed(123)
n <- 100
M <- 50
tgrid <- seq(0, 1, length.out = M)
dt    <- tgrid[2] - tgrid[1]
tmat  <- matrix(rep(tgrid, each = n), nrow = n)
lmat  <- matrix(dt, nrow = n, ncol = M)
wmat  <- t(apply(matrix(rnorm(n * M), n, M), 1, cumsum)) / sqrt(M)
beta_true <- sin(2 * pi * tgrid)
X1 <- rnorm(n)
eta <- 0.5 * X1 + wmat %*% (beta_true * dt)
prob <- plogis(eta)
y <- rbinom(n, 1, prob)
data.SoFR <- data.frame(y = y, X1 = X1)
data.SoFR$tmat <- tmat
data.SoFR$lmat <- lmat
data.SoFR$wmat <- wmat
```

The example dataset `data.SoFR` contains:

- `y`: a binary outcome variable,
- `X1`: a scalar predictor,
- `tmat`: the $n \times M$ time point matrix (defines the domain
  $\mathcal{T}$),
- `lmat`: the $n \times M$ Riemann sum weight matrix (defines the
  integration weights over $\mathcal{T}$),
- `wmat`: the $n \times M$ functional predictor matrix.

### Fit the Bayesian SoFR Model

``` r
library(refundBayes)

refundBayes_SoFR <- refundBayes::sofr_bayes(
  y ~ X1 + s(tmat, by = lmat * wmat, bs = "cc", k = 10),
  data = data.SoFR,
  family = binomial(),
  runStan = TRUE,
  niter = 1500,
  nwarmup = 500,
  nchain = 3,
  ncores = 3
)
```

In this call:

- The formula specifies a binary outcome `y` with one scalar predictor
  `X1` and one functional predictor encoded via
  `s(tmat, by = lmat * wmat, bs = "cc", k = 10)`.
- The spline basis is evaluated at the time points stored in `tmat`, and
  the integration over the domain $\mathcal{T}$ (determined by `tmat`)
  is approximated using the weights in `lmat`.
- `bs = "cc"` uses cyclic cubic regression splines, which are
  appropriate when the functional predictor is periodic (e.g., physical
  activity measured over a 24-hour cycle).
- `k = 10` specifies 10 basis functions. In practice, 30–40 basis
  functions are often sufficient for moderately smooth functional data
  on dense grids.
- `family = binomial()` specifies a logistic regression for the binary
  outcome.
- The sampler runs 3 chains in parallel, each with 1500 total iterations
  (500 warmup + 1000 posterior samples).

### Visualization

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
`refundBayes` objects displays the estimated functional coefficient
$\widehat{\beta}(t)$ along with pointwise 95% credible intervals:

``` r
library(ggplot2)
plot(refundBayes_SoFR)
```

### Extracting Posterior Summaries

Posterior summaries of the functional coefficient can be computed
directly from the `func_coef` element:

``` r
## Posterior mean of the functional coefficient
mean_curve <- apply(refundBayes_SoFR$func_coef[[1]], 2, mean)

## Pointwise 95% credible interval
upper_curve <- apply(refundBayes_SoFR$func_coef[[1]], 2,
                     function(x) quantile(x, prob = 0.975))
lower_curve <- apply(refundBayes_SoFR$func_coef[[1]], 2,
                     function(x) quantile(x, prob = 0.025))
```

The posterior samples in `func_coef[[1]]` are stored as a $Q \times M$
matrix, where $Q$ is the number of posterior samples and $M$ is the
number of time points on the functional domain.

### Comparison with Frequentist Results

The Bayesian results can be compared with frequentist estimates obtained
via [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html):

``` r
library(mgcv)

## Fit frequentist SoFR using mgcv
fit_freq <- gam(
  y ~ s(tmat, by = lmat * wmat, bs = "cc", k = 10) + X1,
  data = data.SoFR,
  family = "binomial"
)

## Extract frequentist estimates
freq_result <- plot(fit_freq)
```

The functional coefficient estimates from the Bayesian and frequentist
approaches are generally comparable in shape and magnitude, though the
Bayesian credible intervals tend to be slightly wider due to accounting
for uncertainty in the smoothing parameter.

### Inspecting the Generated Stan Code

Setting `runStan = FALSE` allows you to inspect or modify the Stan code
before running the model:

``` r
## Generate Stan code without running the sampler
sofr_code <- refundBayes::sofr_bayes(
  y ~ X1 + s(tmat, by = lmat * wmat, bs = "cc", k = 10),
  data = data.SoFR,
  family = binomial(),
  runStan = FALSE
)

## Print the generated Stan code
cat(sofr_code$stancode)
```

## Practical Recommendations

- **Number of basis functions (`k`)**: For illustrative purposes,
  `k = 10` is often used. In practice, 30–40 basis functions are
  recommended for moderately smooth functional data observed on dense
  grids.
- **Spline type (`bs`)**: Use `"cr"` (cubic regression splines) for
  general functional predictors. Use `"cc"` (cyclic cubic regression
  splines) when the functional predictor is periodic (e.g., 24-hour
  activity patterns).
- **Number of iterations**: Ensure sufficient posterior samples for
  reliable inference. A common setup is `niter = 3000` with
  `nwarmup = 1000`, but more complex models may require additional
  iterations.
- **Convergence diagnostics**: After fitting, examine traceplots and
  $\widehat{R}$ statistics using the `rstan` package (e.g.,
  `rstan::traceplot(refundBayes_SoFR$stanfit)`) to ensure that the
  Markov chains have converged.
- **Joint FPCA**: When functional predictors are measured with
  substantial noise, consider setting `joint_FPCA = TRUE` for the
  relevant predictor to jointly estimate FPCA scores and regression
  coefficients.

## References

- Jiang, Z., Crainiceanu, C., and Cui, E. (2025). Tutorial on Bayesian
  Functional Regression Using Stan. *Statistics in Medicine*, 44(20–22),
  e70265.
- Crainiceanu, C. M., Goldsmith, J., Leroux, A., and Cui, E. (2024).
  *Functional Data Analysis with R*. CRC Press.
- Crainiceanu, C. M. and Goldsmith, A. J. (2010). Bayesian Functional
  Data Analysis Using WinBUGS. *Journal of Statistical Software*,
  32(11), 1–33.
- Wood, S. (2001). mgcv: GAMs and Generalized Ridge Regression for R. *R
  News*, 1(2), 20–25.
- Carpenter, B., Gelman, A., Hoffman, M. D., et al. (2017). Stan: A
  Probabilistic Programming Language. *Journal of Statistical Software*,
  76(1), 1–32.

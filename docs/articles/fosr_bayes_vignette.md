# Bayesian Function-on-Scalar Regression

## Introduction

This vignette provides a detailed guide to the
[`fosr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fosr_bayes.md)
function in the `refundBayes` package, which fits Bayesian
Function-on-Scalar Regression (FoSR) models using Stan. In contrast to
scalar-on-function regression (SoFR), where the outcome is scalar and
the predictors include functional variables, FoSR reverses this
relationship: the outcome is a functional variable (a curve observed
over a continuum) and the predictors are scalar.

The methodology follows the tutorial by Jiang, Crainiceanu, and Cui
(2025), *Tutorial on Bayesian Functional Regression Using Stan*,
published in *Statistics in Medicine*. The procedure for fitting a
generalized multilevel Bayesian FoSR was introduced by Goldsmith,
Zipunnikov, and Schrack (2015).

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

### The FoSR Model

Function-on-Scalar Regression (FoSR) models the relationship between a
functional outcome and one or more scalar predictors. Two key
differences from traditional regression are: (1) the outcome $Y_{i}(t)$
is multivariate and often high-dimensional, and (2) the residuals
$e_{i}(t)$ are correlated across $t$ (points in the domain).

For subject $i = 1,\ldots,n$, let $Y_{i}(t)$ be the functional response
observed at time points $t = t_{1},\ldots,t_{M}$ over a domain
$\mathcal{T}$. The domain $\mathcal{T}$ is **not** restricted to
$\lbrack 0,1\rbrack$; it is determined by the actual time points in the
data (e.g., $\mathcal{T} = \lbrack 1,1440\rbrack$ for minute-level
24-hour data). The FoSR model assumes:

$$Y_{i}(t) = \sum\limits_{p = 1}^{P}X_{ip}\beta_{p}(t) + e_{i}(t)$$

where:

- $X_{ip}$ with $p = 1,\ldots,P$ are the scalar predictors (the first
  covariate is often the intercept, $X_{i1} = 1$),
- $\beta_{p}(t)$ are the corresponding domain-varying (or functional)
  coefficients,
- $e_{i}(t)$ is the residual process, which is correlated across $t$.

Each functional coefficient $\beta_{p}(t)$ describes how the $p$-th
scalar predictor affects the outcome at each point $t$ in the domain.
This allows the effect of a scalar predictor to vary smoothly over the
domain.

### Modeling the Residual Structure

To account for the within-subject correlation in the residuals, the
model decomposes $e_{i}(t)$ using functional principal components:

$$e_{i}(t) = \sum\limits_{j = 1}^{J}\xi_{ij}\phi_{j}(t) + \epsilon_{i}(t)$$

where $\phi_{1}(t),\ldots,\phi_{J}(t)$ are the eigenfunctions estimated
via FPCA (using
[`refund::fpca.face`](https://rdrr.io/pkg/refund/man/fpca.face.html)),
$\xi_{ij}$ are the subject-specific scores with
$\xi_{ij} \sim N\left( 0,\lambda_{j} \right)$, and
$\epsilon_{i}(t) \sim N\left( 0,\sigma_{\epsilon}^{2} \right)$ is
independent measurement error. The eigenvalues $\lambda_{j}$ and the
error variance $\sigma_{\epsilon}^{2}$ are estimated from the data.

### Functional Coefficients via Penalized Splines

Each functional coefficient $\beta_{p}(t)$ is represented using $K$
spline basis functions $\psi_{1}(t),\ldots,\psi_{K}(t)$:

$$\beta_{p}(t) = \sum\limits_{k = 1}^{K}b_{pk}\psi_{k}(t)$$

Substituting into the model:

$$Y_{i}(t) = \sum\limits_{k = 1}^{K}\left( \sum\limits_{p = 1}^{P}X_{ip}b_{pk} \right)\psi_{k}(t) + \sum\limits_{j = 1}^{J}\xi_{ij}\phi_{j}(t) + \epsilon_{i}(t)$$

Let $B_{ik} = \sum_{p = 1}^{P}X_{ip}b_{pk}$ and denote by
$\mathbf{\Psi}$ the $M \times K$ matrix of spline basis evaluations, by
$\mathbf{\Phi}$ the $M \times J$ matrix of eigenfunctions, and by
$\mathbf{B}_{i} = \left( B_{i1},\ldots,B_{iK} \right)^{t}$. The model in
matrix form becomes:

$$\mathbf{Y}_{i} = \mathbf{\Psi}\mathbf{B}_{i} + \mathbf{\Phi}{\mathbf{ξ}}_{i} + {\mathbf{ϵ}}_{i}$$

where
$\mathbf{Y}_{i} = \{ Y_{i}\left( t_{1} \right),\ldots,Y_{i}\left( t_{M} \right)\}^{t}$
is the $M \times 1$ vector of observed functional data for subject $i$.

### Smoothness Penalty

Smoothness of each $\beta_{p}(t)$ is induced through a quadratic penalty
on the spline coefficients
$\mathbf{b}_{p} = \left( b_{p1},\ldots,b_{pK} \right)^{t}$:

$$p\left( \mathbf{b}_{p} \right) \propto \exp\left( - \frac{\mathbf{b}_{p}^{t}\mathbf{S}\mathbf{b}_{p}}{2\sigma_{p}^{2}} \right)$$

where $\mathbf{S}$ is the penalty matrix and $\sigma_{p}^{2}$ is the
smoothing parameter for the $p$-th functional coefficient. Importantly,
each functional coefficient $\beta_{p}(t)$ has its own smoothing
parameter $\sigma_{p}^{2}$, allowing different levels of smoothness
across predictors.

Unlike the SoFR and FCox models, the FoSR implementation uses the
original spline basis without the spectral reparametrisation. This
choice follows Goldsmith et al. (2015) and simplifies interpretation by
avoiding the need to back-transform the estimated coefficients.

### Full Bayesian Model

The complete Bayesian FoSR model is:

$$\left\{ \begin{array}{l}
{\mathbf{Y}_{i} = \mathbf{\Psi}\mathbf{B}_{i} + \mathbf{\Phi}{\mathbf{ξ}}_{i} + {\mathbf{ϵ}}_{i}} \\
{B_{ik} = \sum\limits_{p = 1}^{P}X_{ip}b_{pk},\quad k = 1,\ldots,K} \\
{p\left( \mathbf{b}_{p} \right) \propto \exp\left( - \mathbf{b}_{p}^{t}\mathbf{S}\mathbf{b}_{p}/2\sigma_{p}^{2} \right),\quad p = 1,\ldots,P} \\
{\xi_{ij} \sim N\left( 0,\lambda_{j} \right),\quad\lambda_{j} \sim p\left( \lambda_{j} \right),\quad j = 1,\ldots,J} \\
{\epsilon_{i}\left( t_{m} \right) \sim N\left( 0,\sigma_{\epsilon}^{2} \right),\quad i = 1,\ldots,n,\; t = t_{1},\ldots,t_{M}} \\
{\sigma_{\epsilon}^{2} \sim p\left( \sigma_{\epsilon}^{2} \right),\quad\sigma_{p}^{2} \sim p\left( \sigma_{p}^{2} \right),\quad p = 1,\ldots,P} \\
\end{array} \right.$$

The priors $p\left( \lambda_{j} \right)$,
$p\left( \sigma_{\epsilon}^{2} \right)$, and
$p\left( \sigma_{p}^{2} \right)$ are non-informative priors on variance
components, using inverse-Gamma priors $IG(0.001,0.001)$.

The model assumes that the functional responses are observed on a common
set of grid points across all subjects. This ensures that the
$\mathbf{\Psi}$ and $\mathbf{\Phi}$ matrices are identical across
subjects, simplifying computation in Stan.

## The `fosr_bayes()` Function

### Usage

``` r
fosr_bayes(
  formula,
  data,
  joint_FPCA = NULL,
  runStan = TRUE,
  niter = 3000,
  nwarmup = 1000,
  nchain = 3,
  ncores = 1,
  spline_type = "bs",
  spline_df = 10
)
```

### Arguments

| Argument      | Description                                                                                                                                                                                                                 |
|:--------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `formula`     | Functional regression formula. The left-hand side should be the name of the functional response variable (an $n \times M$ matrix in `data`). The right-hand side specifies scalar predictors using standard formula syntax. |
| `data`        | A data frame containing all variables used in the model. The functional response should be stored as an $n \times M$ matrix, where each row is one subject and each column is one time point.                               |
| `joint_FPCA`  | A logical (`TRUE`/`FALSE`) vector of the same length as the number of functional predictors, indicating whether to jointly model FPCA. Default is `NULL`, which sets all entries to `FALSE`.                                |
| `runStan`     | Logical. Whether to run the Stan program. If `FALSE`, the function only generates the Stan code and data without sampling. Default is `TRUE`.                                                                               |
| `niter`       | Total number of Bayesian posterior sampling iterations (including warmup). Default is `3000`.                                                                                                                               |
| `nwarmup`     | Number of warmup (burn-in) iterations. Default is `1000`.                                                                                                                                                                   |
| `nchain`      | Number of Markov chains for posterior sampling. Default is `3`.                                                                                                                                                             |
| `ncores`      | Number of CPU cores to use when executing the chains in parallel. Default is `1`.                                                                                                                                           |
| `spline_type` | Type of spline basis used for the functional coefficients. Default is `"bs"` (B-splines).                                                                                                                                   |
| `spline_df`   | Number of degrees of freedom (basis functions) for the spline basis. Default is `10`.                                                                                                                                       |

### Return Value

The function returns a list of class `"refundBayes"` containing the
following elements:

| Element        | Description                                                                                                                                                                                                                                            |
|:---------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `stanfit`      | The Stan fit object (class `stanfit`). Can be used for convergence diagnostics, traceplots, and additional summaries via the `rstan` package.                                                                                                          |
| `spline_basis` | Basis functions used to reconstruct the functional coefficients from the posterior samples.                                                                                                                                                            |
| `stancode`     | A character string containing the generated Stan model code.                                                                                                                                                                                           |
| `standata`     | A list containing the data passed to the Stan model.                                                                                                                                                                                                   |
| `func_coef`    | An array of posterior samples for functional coefficients, with dimensions $Q \times P \times M$, where $Q$ is the number of posterior samples, $P$ is the number of scalar predictors, and $M$ is the number of time points on the functional domain. |
| `family`       | The family type: `"functional"`.                                                                                                                                                                                                                       |

### Formula Syntax

The formula syntax for FoSR differs from SoFR and FCox because the
outcome is functional and the predictors are scalar. The left-hand side
specifies the functional response, and the right-hand side lists scalar
predictors:

``` r
y ~ X
```

where:

- `y`: the name of the functional response variable in `data`. This
  should be an $n \times M$ matrix, where each row contains the
  functional observations for one subject across $M$ time points.
- `X`: scalar predictor(s). Multiple scalar predictors can be included
  using standard formula syntax (e.g., `y ~ X1 + X2 + X3`).

Note that the FoSR formula does **not** use the
[`s()`](https://rdrr.io/pkg/mgcv/man/s.html) term with `tmat`, `lmat`,
and `wmat` that appears in SoFR and FCox models. This is because in
FoSR, the functional variable is the *outcome* (not a predictor), and
the spline basis for the functional coefficients is constructed
internally using the `spline_type` and `spline_df` arguments.

## Example: Bayesian FoSR

We demonstrate the
[`fosr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fosr_bayes.md)
function using a simulated example dataset with a functional response
and one scalar predictor.

### Load and Prepare Data

``` r
## Load the example data
FoSR_exp_data <- readRDS("data/example_data_FoSR.rds")

## Set the functional response
FoSR_exp_data$y <- FoSR_exp_data$MIMS
```

The example dataset contains:

- `MIMS`: an $n \times M$ matrix of functional response values (physical
  activity data), assigned to `y`,
- `X`: a scalar predictor variable.

### Fit the Bayesian FoSR Model

``` r
library(refundBayes)

refundBayes_FoSR <- refundBayes::fosr_bayes(
  y ~ X,
  data = FoSR_exp_data,
  runStan = TRUE,
  niter = 1500,
  nwarmup = 500,
  nchain = 1,
  ncores = 1
)
```

In this call:

- The formula specifies the functional response `y` (an $n \times M$
  matrix) with one scalar predictor `X`.
- The spline basis for the functional coefficients is constructed
  internally using B-splines (`spline_type = "bs"`) with 10 degrees of
  freedom (`spline_df = 10`) by default.
- The eigenfunctions for the residual structure are estimated
  automatically via FPCA using
  [`refund::fpca.face`](https://rdrr.io/pkg/refund/man/fpca.face.html).
- The sampler runs 1 chain with 1500 total iterations (500 warmup + 1000
  posterior samples). For production analyses, using 3 or more chains is
  recommended for convergence assessment.

#### A Note on Computation

FoSR models are computationally more demanding than SoFR models because
the likelihood involves all $M$ time points for each of the $n$
subjects. The example above uses a single chain for demonstration
purposes. In practice, consider the trade-off between the number of
basis functions, the number of FPCA components, and the computational
budget.

### Visualization

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
`refundBayes` objects displays the estimated functional coefficients
with pointwise 95% credible intervals:

``` r
library(ggplot2)
plot(refundBayes_FoSR)
```

### Extracting Posterior Summaries

Posterior summaries of the functional coefficients can be computed from
the `func_coef` element. The `func_coef` object is an array with
dimensions $Q \times P \times M$, where $Q$ is the number of posterior
samples, $P$ is the number of scalar predictors (including the
intercept), and $M$ is the number of time points:

``` r
## Posterior mean of the functional coefficient for the first scalar predictor
mean_curve <- apply(refundBayes_FoSR$func_coef[, 1, ], 2, mean)

## Pointwise 95% credible interval
upper_curve <- apply(refundBayes_FoSR$func_coef[, 1, ], 2,
                     function(x) quantile(x, prob = 0.975))
lower_curve <- apply(refundBayes_FoSR$func_coef[, 1, ], 2,
                     function(x) quantile(x, prob = 0.025))
```

### Comparison with Frequentist Results

The Bayesian results can be compared with frequentist estimates obtained
via [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html). As described
in Crainiceanu et al. (2024), one key distinction is that the
frequentist approach fits the model at each time point independently or
using fast algorithms, while the Bayesian approach jointly models all
time points with a shared smoothness prior and FPCA-based residual
structure:

``` r
library(refund)

## The frequentist approach can be implemented using mgcv or refund::pffr
## See Crainiceanu et al. (2024) for details

fit.freq = pffr(y~-1+X,data=FoSR_exp_data,bs.yindex=list(bs="cc", k=10))
plotfot = plot(fit.freq)
```

Simulation studies in Jiang et al. (2025) show that the Bayesian FoSR
achieves similar estimation accuracy (RISE) to the frequentist approach,
while providing superior coverage of pointwise credible intervals.

### Inspecting the Generated Stan Code

Setting `runStan = FALSE` allows you to inspect or modify the Stan code
before running the model:

``` r
## Generate Stan code without running the sampler
fosr_code <- refundBayes::fosr_bayes(
  y ~ X,
  data = FoSR_exp_data,
  runStan = FALSE
)

## Print the generated Stan code
cat(fosr_code$stancode)
```

## Practical Recommendations

- **Number of basis functions (`spline_df`)**: The default is
  `spline_df = 10`. In practice, the appropriate number depends on the
  complexity of the underlying functional coefficients and the
  resolution of the data. More basis functions allow greater flexibility
  but increase computational cost.
- **Spline type (`spline_type`)**: The default `"bs"` uses B-splines.
  Other basis types supported by `mgcv` may also be used depending on
  the application.
- **Number of FPCA components**: The number of eigenfunctions $J$ used
  to model the residual correlation is determined automatically by
  [`refund::fpca.face`](https://rdrr.io/pkg/refund/man/fpca.face.html).
  If the residual structure is complex, more components may be needed.
- **Number of iterations and chains**: FoSR models are computationally
  intensive. Start with fewer iterations and a single chain for
  exploratory analysis, then increase for final inference. A recommended
  setup for production is `niter = 3000`, `nwarmup = 1000`,
  `nchain = 3`.
- **Convergence diagnostics**: After fitting, examine traceplots and
  $\widehat{R}$ statistics using the `rstan` package (e.g.,
  `rstan::traceplot(refundBayes_FoSR$stanfit)`). Warnings about bulk
  ESS, tail ESS, or $\widehat{R}$ indicate that more iterations or
  chains may be needed.
- **Common grid assumption**: The current implementation assumes that
  the functional responses are observed on a common set of grid points
  across all subjects. Subject-specific observation grids are not yet
  supported.

## References

- Jiang, Z., Crainiceanu, C., and Cui, E. (2025). Tutorial on Bayesian
  Functional Regression Using Stan. *Statistics in Medicine*, 44(20–22),
  e70265.
- Crainiceanu, C. M., Goldsmith, J., Leroux, A., and Cui, E. (2024).
  *Functional Data Analysis with R*. CRC Press.
- Goldsmith, J., Zipunnikov, V., and Schrack, J. (2015). Generalized
  Multilevel Function-on-Scalar Regression and Principal Component
  Analysis. *Biometrics*, 71(2), 344–353.
- Carpenter, B., Gelman, A., Hoffman, M. D., et al. (2017). Stan: A
  Probabilistic Programming Language. *Journal of Statistical Software*,
  76(1), 1–32.

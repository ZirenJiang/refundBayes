# Bayesian Function-on-Function Regression with \`refundBayes::fofr_bayes\`

## Introduction

This vignette provides a detailed guide to the
[`fofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fofr_bayes.md)
function in the `refundBayes` package, which fits Bayesian
Function-on-Function Regression (FoFR) models using Stan. FoFR extends
both the scalar-on-function regression (SoFR) and function-on-scalar
regression (FoSR) frameworks by modeling a **functional response** as a
function of one or more **functional predictors**, along with optional
scalar covariates.

In contrast to SoFR, where the outcome is scalar and the predictors are
functional, and to FoSR, where the outcome is functional and the
predictors are scalar, FoFR allows the effect of a functional predictor
on a functional response to be captured through a **bivariate
coefficient function** $\beta(s,t)$. The model is specified using the
same `mgcv`-style syntax as the other regression functions in the
`refundBayes` package.

The methodology extends the framework described in Jiang, Crainiceanu,
and Cui (2025), *Tutorial on Bayesian Functional Regression Using Stan*,
published in *Statistics in Medicine*.

## Install the `refundBayes` Package

The `refundBayes` package can be installed from GitHub:

``` r
library(remotes)
remotes::install_github("https://github.com/ZirenJiang/refundBayes")
```

## Statistical Model

### The FoFR Model

Function-on-Function Regression (FoFR) models the relationship between a
functional response and one or more functional predictors, optionally
adjusted for scalar covariates. Both the response and the predictors are
curves observed over (potentially different) continua.

For subject $i = 1,\ldots,n$, let $Y_{i}(t)$ be the functional response
observed at time points $t = t_{1},\ldots,t_{M}$ over a response domain
$\mathcal{T}$, and let
$\{ W_{i}\left( s_{l} \right),s_{l} \in \mathcal{S}\}$ for
$l = 1,\ldots,L$ be a functional predictor observed at $L$ points over a
predictor domain $\mathcal{S}$. Let
$\mathbf{X}_{i} = \left( X_{i1},\ldots,X_{iP} \right)^{t}$ be a
$P \times 1$ vector of scalar predictors (the first covariate may be an
intercept, $X_{i1} = 1$). The FoFR model assumes:

$$Y_{i}(t) = \sum\limits_{p = 1}^{P}X_{ip}\,\alpha_{p}(t) + \int_{\mathcal{S}}W_{i}(s)\,\beta(s,t)\, ds + e_{i}(t)$$

where:

- $\alpha_{p}(t)$ are functional coefficients for the scalar predictors,
  each describing how the $p$-th scalar predictor affects the response
  at each point $t \in \mathcal{T}$,
- $\beta(s,t)$ is the **bivariate coefficient function** that
  characterizes how the functional predictor at predictor-domain
  location $s$ affects the response at response-domain location $t$,
- $e_{i}(t)$ is the residual process, which is correlated across $t$.

The integral $\int_{\mathcal{S}}W_{i}(s)\,\beta(s,t)\, ds$ is
approximated using a Riemann sum over the observed predictor-domain grid
points. The domains $\mathcal{S}$ and $\mathcal{T}$ are **not**
restricted to $\lbrack 0,1\rbrack$; they are determined by the actual
observation grids in the data.

When multiple functional predictors are present, the model extends
naturally:

$$Y_{i}(t) = \sum\limits_{p = 1}^{P}X_{ip}\,\alpha_{p}(t) + \sum\limits_{j = 1}^{J}\int_{\mathcal{S}_{j}}W_{ij}(s)\,\beta_{j}(s,t)\, ds + e_{i}(t)$$

where $\beta_{j}(s,t)$ is the bivariate coefficient function for the
$j$-th functional predictor.

### Modeling the Residual Structure

To account for the within-subject correlation in the residuals, the
model decomposes $e_{i}(t)$ using functional principal components,
following the same approach as in FoSR (Goldsmith, Zipunnikov, and
Schrack, 2015):

$$e_{i}(t) = \sum\limits_{r = 1}^{R}\xi_{ir}\,\phi_{r}(t) + \epsilon_{i}(t)$$

where $\phi_{1}(t),\ldots,\phi_{R}(t)$ are the eigenfunctions estimated
via FPCA (using
[`refund::fpca.face`](https://rdrr.io/pkg/refund/man/fpca.face.html)),
$\xi_{ir}$ are the subject-specific FPCA scores with
$\xi_{ir} \sim N\left( 0,\lambda_{r} \right)$, and
$\epsilon_{i}(t) \sim N\left( 0,\sigma_{\epsilon}^{2} \right)$ is
independent measurement error. The eigenvalues $\lambda_{r}$ and the
error variance $\sigma_{\epsilon}^{2}$ are estimated from the data.

### Scalar Predictor Coefficients via Penalized Splines

Each scalar predictor coefficient function $\alpha_{p}(t)$ is
represented using $K$ spline basis functions
$\psi_{1}(t),\ldots,\psi_{K}(t)$ in the response domain:

$$\alpha_{p}(t) = \sum\limits_{k = 1}^{K}a_{pk}\,\psi_{k}(t)$$

Smoothness is induced through a quadratic penalty:

$$p\left( \mathbf{a}_{p} \right) \propto \exp\left( - \frac{\mathbf{a}_{p}^{t}\mathbf{S}\mathbf{a}_{p}}{2\sigma_{p}^{2}} \right)$$

where $\mathbf{S}$ is the penalty matrix derived from the spline basis
and $\sigma_{p}^{2}$ is the smoothing parameter for the $p$-th scalar
predictor, estimated from the data.

### Bivariate Coefficient via Tensor Product Basis

The key feature of FoFR is the bivariate coefficient function
$\beta(s,t)$, which lives on the product domain
$\mathcal{S} \times \mathcal{T}$. This function is represented using a
**tensor product** of two sets of basis functions:

- **Predictor-domain basis**: $B_{1}(s),\ldots,B_{Q}(s)$, constructed
  from the [`s()`](https://rdrr.io/pkg/mgcv/man/s.html) term in the
  formula using the `mgcv` spline basis (e.g., cubic regression
  splines). These are further decomposed into random and fixed effect
  components via the spectral reparametrisation of
  [`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html),
  yielding
  ${\widetilde{B}}_{1}^{r}(s),\ldots,{\widetilde{B}}_{Q_{r}}^{r}(s)$
  (random effects) and
  ${\widetilde{B}}_{1}^{f}(s),\ldots,{\widetilde{B}}_{Q_{f}}^{f}(s)$
  (fixed effects).

- **Response-domain basis**: $\psi_{1}(t),\ldots,\psi_{K}(t)$, the same
  spline basis used for the scalar predictor coefficients.

The bivariate coefficient is then:

$$\beta(s,t) = \sum\limits_{k = 1}^{K}\left\lbrack \sum\limits_{q = 1}^{Q_{r}}\theta_{qk}^{r}\,{\widetilde{B}}_{q}^{r}(s) + \sum\limits_{q = 1}^{Q_{f}}\theta_{qk}^{f}\,{\widetilde{B}}_{q}^{f}(s) \right\rbrack\psi_{k}(t)$$

where $\mathbf{\Theta}^{r}$ is the $Q_{r} \times K$ matrix of random
effect coefficients and $\mathbf{\Theta}^{f}$ is the $Q_{f} \times K$
matrix of fixed effect coefficients.

With this representation, the integral contribution of the functional
predictor becomes:

$$\int_{\mathcal{S}}W_{i}(s)\,\beta(s,t)\, ds \approx \left\lbrack {\widetilde{\mathbf{X}}}_{i}^{r}\mathbf{\Theta}^{r} + {\widetilde{\mathbf{X}}}_{i}^{f}\mathbf{\Theta}^{f} \right\rbrack{\mathbf{ψ}}(t)$$

where ${\widetilde{\mathbf{X}}}_{i}^{r}$ and
${\widetilde{\mathbf{X}}}_{i}^{f}$ are the $1 \times Q_{r}$ and
$1 \times Q_{f}$ row vectors of the transformed predictor-domain design
matrices for subject $i$ (the same matrices used in SoFR).

In matrix notation for all subjects, the functional predictor
contribution to the mean is:

$$\left( {\widetilde{\mathbf{X}}}^{r}\mathbf{\Theta}^{r} + {\widetilde{\mathbf{X}}}^{f}\mathbf{\Theta}^{f} \right)\mathbf{\Psi}$$

which is an $n \times M$ matrix, where $\mathbf{\Psi}$ is the
$K \times M$ matrix of response-domain spline basis evaluations.

### Dual-Direction Smoothness Penalties

Because the bivariate coefficient $\beta(s,t)$ lives on a
two-dimensional domain, smoothness must be enforced in **both**
directions:

#### Predictor-domain smoothness ($s$-direction)

Following the same approach as in SoFR, the random effect coefficients
use a variance-component reparametrisation:

$$\mathbf{\Theta}^{r} = \sigma_{s} \cdot \mathbf{Z}^{r},\quad\text{vec}\left( \mathbf{Z}^{r} \right) \sim N(\mathbf{0},\mathbf{I})$$

where $\sigma_{s}^{2}$ is the predictor-domain smoothing parameter. This
is the standard non-centered parameterisation that separates the scale
($\sigma_{s}$) from the direction ($\mathbf{Z}^{r}$), improving sampling
efficiency in Stan.

#### Response-domain smoothness ($t$-direction)

The response-domain penalty matrix $\mathbf{S}$ is applied row-wise to
the coefficient matrices. For each row $q$ of $\mathbf{\Theta}^{r}$ and
$\mathbf{\Theta}^{f}$:

$$p\left( \mathbf{\Theta}_{q}^{r} \right) \propto \exp\left( - \frac{\mathbf{\Theta}_{q}^{r}\mathbf{S}\,\left( \mathbf{\Theta}_{q}^{r} \right)^{t}}{2\sigma_{t}^{2}} \right),\qquad p\left( \mathbf{\Theta}_{q}^{f} \right) \propto \exp\left( - \frac{\mathbf{\Theta}_{q}^{f}\mathbf{S}\,\left( \mathbf{\Theta}_{q}^{f} \right)^{t}}{2\sigma_{t}^{2}} \right)$$

where $\sigma_{t}^{2}$ is the response-domain smoothing parameter. This
ensures that $\beta(s,t)$ is smooth in the $t$-direction for every fixed
$s$.

The two smoothing parameters $\sigma_{s}^{2}$ and $\sigma_{t}^{2}$ are
estimated from the data with weakly informative inverse-Gamma priors.

### Full Bayesian Model

The complete Bayesian FoFR model combines the mean structure, residual
FPCA decomposition, and all priors:

$$\left\{ \begin{array}{l}
{Y_{i}(t) = {\mathbf{μ}}_{i}(t) + e_{i}(t),\quad i = 1,\ldots,n} \\
{{\mathbf{μ}}_{i}(t) = \sum\limits_{p = 1}^{P}X_{ip}\,\alpha_{p}(t) + \int_{\mathcal{S}}W_{i}(s)\,\beta(s,t)\, ds} \\
{e_{i}(t) = \sum\limits_{r = 1}^{R}\xi_{ir}\,\phi_{r}(t) + \epsilon_{i}(t)} \\
\end{array} \right.$$

In matrix form for all subjects:

$$\mathbf{Y} = \mathbf{X}\,\mathbf{A}^{t}\mathbf{\Psi} + \left( {\widetilde{\mathbf{X}}}^{r}\mathbf{\Theta}^{r} + {\widetilde{\mathbf{X}}}^{f}\mathbf{\Theta}^{f} \right)\mathbf{\Psi} + \mathbf{\Xi}\,\mathbf{\Phi} + \mathbf{E}$$

where $\mathbf{Y}$ is the $n \times M$ matrix of functional responses,
$\mathbf{X}$ is the $n \times P$ scalar design matrix, $\mathbf{A}$ is
the $K \times P$ matrix of scalar predictor spline coefficients,
$\mathbf{\Xi}$ is the $n \times R$ matrix of FPCA scores, and
$\mathbf{\Phi}$ is the $R \times M$ matrix of eigenfunctions.

#### Prior Specification

The full prior specification is:

| Parameter                                                               | Prior                                                                                                                                             | Description                                              |
|:------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------|
| $\mathbf{a}_{p}$ (scalar predictor spline coefs)                        | $p\left( \mathbf{a}_{p} \right) \propto \exp\left( - \frac{\mathbf{a}_{p}^{t}\mathbf{S}\,\mathbf{a}_{p}}{2\sigma_{p}^{2}} \right)$                | Penalized spline prior for smoothness of $\alpha_{p}(t)$ |
| $\sigma_{p}^{2}$ (scalar smoothing parameter)                           | $\sigma_{p}^{2} \sim \text{Inv-Gamma}(0.001,0.001)$                                                                                               | Weakly informative prior on smoothing                    |
| $\text{vec}\left( \mathbf{Z}^{r} \right)$ (standardized random effects) | $\text{vec}\left( \mathbf{Z}^{r} \right) \sim N(\mathbf{0},\mathbf{I})$                                                                           | Non-centered parameterisation for $s$-direction          |
| $\sigma_{s}^{2}$ (predictor-domain smoothing)                           | $\sigma_{s}^{2} \sim \text{Inv-Gamma}(0.0005,0.0005)$                                                                                             | Prior on predictor-domain smoothing                      |
| $\mathbf{\Theta}_{q}^{r},\mathbf{\Theta}_{q}^{f}$ (row-wise)            | $p\left( \mathbf{\Theta}_{q} \right) \propto \exp\left( - \frac{\mathbf{\Theta}_{q}\mathbf{S}\,\mathbf{\Theta}_{q}^{t}}{2\sigma_{t}^{2}} \right)$ | Penalized spline prior for response-domain smoothness    |
| $\sigma_{t}^{2}$ (response-domain smoothing)                            | $\sigma_{t}^{2} \sim \text{Inv-Gamma}(0.001,0.001)$                                                                                               | Prior on response-domain smoothing                       |
| $\xi_{ir}$ (FPCA scores)                                                | $\xi_{ir} \sim N\left( 0,\lambda_{r} \right)$                                                                                                     | FPCA score distribution                                  |
| $\lambda_{r}$ (FPCA eigenvalues)                                        | $\lambda_{r}^{2} \sim \text{Inv-Gamma}(0.001,0.001)$                                                                                              | Weakly informative prior on eigenvalues                  |
| $\sigma_{\epsilon}^{2}$ (residual variance)                             | $\sigma_{\epsilon}^{2} \sim \text{Inv-Gamma}(0.001,0.001)$                                                                                        | Weakly informative prior on noise variance               |

#### Likelihood

The likelihood for the functional response is Gaussian:

$$\log p\left( \mathbf{Y} \mid {\mathbf{μ}},\sigma_{\epsilon}^{2} \right) = - \frac{nM}{2}\log\sigma_{\epsilon} - \frac{1}{2\sigma_{\epsilon}^{2}}\sum\limits_{i = 1}^{n}\sum\limits_{m = 1}^{M}\{ Y_{i}\left( t_{m} \right) - \mu_{i}\left( t_{m} \right)\}^{2}$$

where the mean function $\mu_{i}\left( t_{m} \right)$ includes
contributions from scalar predictors, functional predictors, and FPCA
scores.

### Relationship to SoFR and FoSR

The FoFR model nests both the SoFR and FoSR models as special cases:

| Model    | Response                  | Predictors                                  | Coefficient                                             | Implemented in                                                                         |
|:---------|:--------------------------|:--------------------------------------------|:--------------------------------------------------------|:---------------------------------------------------------------------------------------|
| SoFR     | Scalar $Y_{i}$            | Functional $W_{i}(s)$                       | Univariate $\beta(s)$                                   | [`sofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/sofr_bayes.md)     |
| FoSR     | Functional $Y_{i}(t)$     | Scalar $X_{ip}$                             | Univariate $\alpha_{p}(t)$                              | [`fosr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fosr_bayes.md)     |
| **FoFR** | **Functional $Y_{i}(t)$** | **Functional $W_{i}(s)$ + Scalar $X_{ip}$** | **Bivariate $\beta(s,t)$ + Univariate $\alpha_{p}(t)$** | **[`fofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fofr_bayes.md)** |

The
[`fofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fofr_bayes.md)
function inherits:

- From FoSR: the response-domain spline basis $\mathbf{\Psi}$, the FPCA
  residual structure, and the scalar predictor handling.
- From SoFR: the predictor-domain spectral reparametrisation
  (random/fixed effect decomposition via
  [`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html)),
  the basis extraction, and the functional coefficient reconstruction.

## The `fofr_bayes()` Function

### Usage

``` r
fofr_bayes(
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

| Argument      | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|:--------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `formula`     | Functional regression formula, using the same syntax as [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html). The left-hand side is the functional response (an $n \times M$ matrix in `data`). The right-hand side includes scalar predictors as standard terms and functional predictors via `s(..., by = ...)` terms. At least one functional predictor must be present; otherwise use [`fosr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fosr_bayes.md). |
| `data`        | A data frame containing all variables used in the model. The functional response and functional predictor should be stored as $n \times M$ and $n \times L$ matrices, respectively.                                                                                                                                                                                                                                                                                           |
| `joint_FPCA`  | A logical (`TRUE`/`FALSE`) vector of the same length as the number of functional predictors, indicating whether to jointly model FPCA for each functional predictor. Default is `NULL`, which sets all entries to `FALSE`.                                                                                                                                                                                                                                                    |
| `runStan`     | Logical. Whether to run the Stan program. If `FALSE`, the function only generates the Stan code and data without sampling. This is useful for inspecting or modifying the generated Stan code. Default is `TRUE`.                                                                                                                                                                                                                                                             |
| `niter`       | Total number of Bayesian posterior sampling iterations (including warmup). Default is `3000`.                                                                                                                                                                                                                                                                                                                                                                                 |
| `nwarmup`     | Number of warmup (burn-in) iterations. These samples are discarded and not used for inference. Default is `1000`.                                                                                                                                                                                                                                                                                                                                                             |
| `nchain`      | Number of Markov chains for posterior sampling. Multiple chains help assess convergence. Default is `3`.                                                                                                                                                                                                                                                                                                                                                                      |
| `ncores`      | Number of CPU cores to use when executing the chains in parallel. Default is `1`.                                                                                                                                                                                                                                                                                                                                                                                             |
| `spline_type` | Type of spline basis used for the response-domain component. Default is `"bs"` (B-splines). Other types supported by `mgcv` may also be used.                                                                                                                                                                                                                                                                                                                                 |
| `spline_df`   | Number of degrees of freedom (basis functions) for the response-domain spline basis. Default is `10`.                                                                                                                                                                                                                                                                                                                                                                         |

### Return Value

The function returns a list of class `"refundBayes"` containing the
following elements:

| Element            | Description                                                                                                                                                                                                                                                                                                 |
|:-------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `stanfit`          | The Stan fit object (class `stanfit`). Can be used for convergence diagnostics, traceplots, and additional summaries via the `rstan` package.                                                                                                                                                               |
| `spline_basis`     | Basis functions used to reconstruct the functional coefficients from the posterior samples.                                                                                                                                                                                                                 |
| `stancode`         | A character string containing the generated Stan model code.                                                                                                                                                                                                                                                |
| `standata`         | A list containing the data passed to the Stan model.                                                                                                                                                                                                                                                        |
| `scalar_func_coef` | A 3-d array ($Q \times P \times M$) of posterior samples for scalar predictor coefficient functions $\alpha_{p}(t)$, where $Q$ is the number of posterior samples, $P$ is the number of scalar predictors, and $M$ is the number of response-domain time points. `NULL` if no scalar predictors.            |
| `bivar_func_coef`  | A list of 3-d arrays. Each element corresponds to one functional predictor and is an array of dimension $Q \times L \times M$, representing posterior samples of the bivariate coefficient function $\beta(s,t)$ evaluated on the predictor-domain grid ($L$ points) and response-domain grid ($M$ points). |
| `func_coef`        | Same as `scalar_func_coef`; included for compatibility with the [`plot.refundBayes()`](https://zirenjiang.github.io/refundBayes/reference/plot.refundBayes.md) method.                                                                                                                                      |
| `family`           | The model family: `"fofr"`.                                                                                                                                                                                                                                                                                 |

### Formula Syntax

The formula combines the FoSR syntax (functional response on the
left-hand side) with the SoFR syntax (functional predictors via
[`s()`](https://rdrr.io/pkg/mgcv/man/s.html) terms on the right-hand
side):

``` r
Y_mat ~ X1 + X2 + s(sindex, by = X_func, bs = "cr", k = 10)
```

where:

- `Y_mat`: the name of the functional response variable in `data`. This
  should be an $n \times M$ matrix, where each row contains the
  functional observations for one subject across $M$ response-domain
  time points.
- `X1`, `X2`: scalar predictor(s), included using standard formula
  syntax.
- `s(sindex, by = X_func, bs = "cr", k = 10)`: the functional predictor
  term:
  - `sindex`: an $n \times L$ matrix of predictor-domain grid points.
    Each row contains the same $L$ observation points (replicated across
    subjects).
  - `X_func`: an $n \times L$ matrix of functional predictor values. The
    $i$-th row contains the $L$ observed values for subject $i$.
  - `bs`: the type of spline basis for the predictor domain (e.g.,
    `"cr"` for cubic regression splines).
  - `k`: the number of basis functions in the predictor domain.

The response-domain spline basis is controlled separately via the
`spline_type` and `spline_df` arguments to
[`fofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fofr_bayes.md).
This design separates the two basis specifications: the predictor-domain
basis is specified in the formula (as in SoFR), while the
response-domain basis is specified via function arguments (as in FoSR).

Multiple functional predictors can be included by adding additional
[`s()`](https://rdrr.io/pkg/mgcv/man/s.html) terms.

## Example: Bayesian FoFR with Simulated Data

We demonstrate the
[`fofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fofr_bayes.md)
function using a simulation study with a known bivariate coefficient
function $\beta(s,t)$ and a scalar predictor coefficient function
$\alpha(t)$.

### Simulate Data

``` r
library(refundBayes)

set.seed(42)

# --- Dimensions ---
n  <- 200   # number of subjects
L  <- 30    # number of predictor-domain grid points
M  <- 30    # number of response-domain grid points

sindex <- seq(0, 1, length.out = L)   # predictor domain grid
tindex <- seq(0, 1, length.out = M)   # response domain grid

# --- Functional predictor X(s): smooth random curves ---
X_func <- matrix(0, nrow = n, ncol = L)
for (i in 1:n) {
  X_func[i, ] <- rnorm(1) * sin(2 * pi * sindex) +
                 rnorm(1) * cos(2 * pi * sindex) +
                 rnorm(1) * sin(4 * pi * sindex) +
                 rnorm(1, sd = 0.3)
}

# --- Scalar predictor ---
age <- rnorm(n)

# --- True coefficient functions ---
# Bivariate coefficient: beta(s, t) = sin(2*pi*s) * cos(2*pi*t)
beta_true <- outer(sin(2 * pi * sindex), cos(2 * pi * tindex))

# Scalar coefficient function: alpha(t) = 0.5 * sin(pi*t)
alpha_true <- 0.5 * sin(pi * tindex)

# --- Generate functional response ---
# Y_i(t) = age_i * alpha(t) + integral X_i(s) beta(s,t) ds + epsilon_i(t)
signal_scalar <- outer(age, alpha_true)                    # n x M
signal_func   <- (X_func %*% beta_true) / L               # n x M  (Riemann sum)
epsilon        <- matrix(rnorm(n * M, sd = 0.3), nrow = n) # n x M

Y_mat <- signal_scalar + signal_func + epsilon

# --- Organize data ---
dat <- data.frame(age = age)
dat$Y_mat  <- Y_mat
dat$X_func <- X_func
dat$sindex <- matrix(rep(sindex, n), nrow = n, byrow = TRUE)
```

The simulated dataset `dat` contains:

- `Y_mat`: an $n \times M$ matrix of functional response values,
- `age`: a scalar predictor,
- `X_func`: an $n \times L$ matrix of functional predictor values
  (smooth random curves),
- `sindex`: an $n \times L$ matrix of predictor-domain grid points
  (identical rows).

The true data-generating model is:
$$Y_{i}(t) = \text{age}_{i} \cdot 0.5\sin(\pi t) + \frac{1}{L}\sum\limits_{l = 1}^{L}X_{i}\left( s_{l} \right)\,\sin\left( 2\pi s_{l} \right)\cos(2\pi t) + \epsilon_{i}(t),\quad\epsilon_{i}(t) \sim N\left( 0,0.3^{2} \right)$$

### Fit the Bayesian FoFR Model

``` r
fit_fofr <- fofr_bayes(
  formula     = Y_mat ~ age + s(sindex, by = X_func, bs = "cr", k = 10),
  data        = dat,
  spline_type = "bs",
  spline_df   = 10,
  niter       = 2000,
  nwarmup     = 1000,
  nchain      = 3,
  ncores      = 3
)
```

In this call:

- The formula specifies the functional response `Y_mat` (an $n \times M$
  matrix) with one scalar predictor `age` and one functional predictor
  `X_func`.
- The predictor-domain spline basis uses cubic regression splines
  (`bs = "cr"`) with `k = 10` basis functions.
- The response-domain spline basis uses B-splines (`spline_type = "bs"`)
  with `spline_df = 10` degrees of freedom.
- The eigenfunctions for the residual structure are estimated
  automatically via FPCA using
  [`refund::fpca.face`](https://rdrr.io/pkg/refund/man/fpca.face.html).
- The sampler runs 3 chains in parallel, each with 2000 total iterations
  (1000 warmup + 1000 posterior samples).

#### A Note on Computation

FoFR models are the most computationally demanding among the models in
`refundBayes` because the Stan program estimates bivariate coefficient
matrices (with $Q_{r} \times K + Q_{f} \times K$ parameters per
functional predictor) in addition to the scalar predictor coefficients
and FPCA scores. For exploratory analyses, consider using fewer basis
functions (e.g., `k = 5`, `spline_df = 5`) and a single chain. For final
inference, use the full setup with multiple chains and convergence
diagnostics.

### Visualisation

#### Bivariate Coefficient $\widehat{\beta}(s,t)$

The estimated bivariate coefficient $\widehat{\beta}(s,t)$ is stored as
a 3-d array in `bivar_func_coef`. The posterior mean surface and
comparison with the truth can be visualised using heatmaps:

``` r
# Posterior mean of the bivariate coefficient
beta_est  <- apply(fit_fofr$bivar_func_coef[[1]], c(2, 3), mean)

# Pointwise 95% credible interval bounds
beta_lower <- apply(fit_fofr$bivar_func_coef[[1]], c(2, 3),
                    function(x) quantile(x, 0.025))
beta_upper <- apply(fit_fofr$bivar_func_coef[[1]], c(2, 3),
                    function(x) quantile(x, 0.975))

# Side-by-side heatmaps: true vs estimated vs difference
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
image(sindex, tindex, beta_true,
      xlab = "s (predictor domain)", ylab = "t (response domain)",
      main = expression("True " * beta(s, t)),
      col = hcl.colors(64, "Blue-Red 3"))
image(sindex, tindex, beta_est,
      xlab = "s (predictor domain)", ylab = "t (response domain)",
      main = expression("Estimated " * hat(beta)(s, t)),
      col = hcl.colors(64, "Blue-Red 3"))
image(sindex, tindex, beta_est - beta_true,
      xlab = "s (predictor domain)", ylab = "t (response domain)",
      main = "Difference (Est - True)",
      col = hcl.colors(64, "Blue-Red 3"))
```

For richer 3-d surface visualisations, use
[`fields::image.plot()`](https://rdrr.io/pkg/fields/man/image.plot.html)
or [`plotly::plot_ly()`](https://rdrr.io/pkg/plotly/man/plot_ly.html)
with type `"surface"`.

#### Scalar Coefficient Function $\widehat{\alpha}(t)$

The estimated scalar predictor coefficient function can be plotted with
pointwise credible intervals:

``` r
alpha_est   <- apply(fit_fofr$scalar_func_coef[, 1, ], 2, mean)
alpha_lower <- apply(fit_fofr$scalar_func_coef[, 1, ], 2,
                     function(x) quantile(x, 0.025))
alpha_upper <- apply(fit_fofr$scalar_func_coef[, 1, ], 2,
                     function(x) quantile(x, 0.975))

par(mfrow = c(1, 1))
plot(tindex, alpha_true, type = "l", lwd = 2, col = "black",
     ylim = range(c(alpha_lower, alpha_upper)),
     xlab = "t (response domain)", ylab = expression(alpha(t)),
     main = "Scalar coefficient function: age")
lines(tindex, alpha_est, col = "blue", lwd = 2)
polygon(c(tindex, rev(tindex)),
        c(alpha_lower, rev(alpha_upper)),
        col = rgb(0, 0, 1, 0.2), border = NA)
legend("topright",
       legend = c("Truth", "Posterior mean", "95% CI"),
       col = c("black", "blue", rgb(0, 0, 1, 0.2)),
       lwd = c(2, 2, 10), bty = "n")
```

#### Slices of the Bivariate Coefficient

To examine $\beta(s,t)$ at fixed values of $s$ or $t$, extract slices
from the posterior:

``` r
# Fix s at the midpoint of the predictor domain and plot beta(s_mid, t)
s_mid_idx <- which.min(abs(sindex - 0.5))

beta_slice_est   <- apply(fit_fofr$bivar_func_coef[[1]][, s_mid_idx, ], 2, mean)
beta_slice_lower <- apply(fit_fofr$bivar_func_coef[[1]][, s_mid_idx, ], 2,
                          function(x) quantile(x, 0.025))
beta_slice_upper <- apply(fit_fofr$bivar_func_coef[[1]][, s_mid_idx, ], 2,
                          function(x) quantile(x, 0.975))
beta_slice_true  <- beta_true[s_mid_idx, ]

plot(tindex, beta_slice_true, type = "l", lwd = 2, col = "black",
     ylim = range(c(beta_slice_lower, beta_slice_upper)),
     xlab = "t (response domain)",
     ylab = expression(beta(s[mid], t)),
     main = paste0("Slice at s = ", round(sindex[s_mid_idx], 2)))
lines(tindex, beta_slice_est, col = "red", lwd = 2)
polygon(c(tindex, rev(tindex)),
        c(beta_slice_lower, rev(beta_slice_upper)),
        col = rgb(1, 0, 0, 0.2), border = NA)
legend("topright",
       legend = c("Truth", "Posterior mean", "95% CI"),
       col = c("black", "red", rgb(1, 0, 0, 0.2)),
       lwd = c(2, 2, 10), bty = "n")
```

### Numerical Summary

``` r
# RMSE of the bivariate coefficient surface
cat("RMSE of beta(s,t):", sqrt(mean((beta_est - beta_true)^2)), "\n")

# RMSE of the scalar coefficient function
cat("RMSE of alpha(t): ", sqrt(mean((alpha_est - alpha_true)^2)), "\n")
```

### Inspecting the Generated Stan Code

Setting `runStan = FALSE` allows you to inspect or modify the Stan code
before running the model:

``` r
# Generate Stan code without running the sampler
fofr_code <- fofr_bayes(
  formula     = Y_mat ~ age + s(sindex, by = X_func, bs = "cr", k = 10),
  data        = dat,
  spline_type = "bs",
  spline_df   = 10,
  runStan     = FALSE
)

# Print the generated Stan code
cat(fofr_code$stancode)
```

The generated Stan code includes all five standard blocks (`data`,
`transformed data`, `parameters`, `transformed parameters`, `model`).
The `parameters` block declares matrix-valued parameters for the
bivariate coefficients, and the `model` block includes both
$s$-direction and $t$-direction smoothness priors.

## Practical Recommendations

- **Number of predictor-domain basis functions (`k`)**: Controls the
  flexibility of $\beta(s,t)$ in the $s$-direction. Start with `k = 10`
  for exploration. In practice, 10–20 basis functions are typically
  sufficient, but this depends on the complexity of the true
  $\beta(s,t)$.

- **Number of response-domain basis functions (`spline_df`)**: Controls
  the flexibility in the $t$-direction. The default `spline_df = 10` is
  often adequate for moderately smooth coefficient functions.

- **Spline types (`bs` and `spline_type`)**: Use `"cr"` (cubic
  regression splines) for general functional data. Use `"cc"` (cyclic
  cubic regression splines) when the functional data are periodic. The
  predictor-domain basis (`bs`) and response-domain basis
  (`spline_type`) may use different types.

- **Sample size and grid resolution**: FoFR requires estimating a
  surface $\beta(s,t)$, which demands more data than SoFR or FoSR. As a
  rough guide, ensure $n > k \times \texttt{𝚜𝚙𝚕𝚒𝚗𝚎\_𝚍𝚏}$.

- **Number of iterations and chains**: FoFR models have more parameters
  than SoFR or FoSR. A recommended starting point is `niter = 3000`,
  `nwarmup = 1000`, `nchain = 3`. Increase iterations if convergence
  diagnostics indicate issues.

- **Convergence diagnostics**: After fitting, examine traceplots and
  $\widehat{R}$ statistics using the `rstan` package:

  ``` r
  rstan::traceplot(fit_fofr$stanfit, pars = c("sigma_eps", "sigmabr_1", "sigma_t_1"))
  print(fit_fofr$stanfit, pars = c("sigma_eps", "sigmabr_1", "sigma_t_1"))
  ```

  Warnings about bulk ESS, tail ESS, or $\widehat{R} > 1.01$ indicate
  that more iterations or chains may be needed.

- **Common grid assumption**: The current implementation assumes that
  both the functional response and functional predictors are observed on
  common grids across all subjects. Subject-specific observation grids
  are not yet supported.

- **Multiple functional predictors**: Multiple functional predictors can
  be included by adding additional
  [`s()`](https://rdrr.io/pkg/mgcv/man/s.html) terms in the formula.
  Each functional predictor receives its own pair of smoothing
  parameters ($\sigma_{s,j}^{2}$, $\sigma_{t,j}^{2}$).

## References

- Jiang, Z., Crainiceanu, C., and Cui, E. (2025). Tutorial on Bayesian
  Functional Regression Using Stan. *Statistics in Medicine*, 44(20–22),
  e70265.
- Crainiceanu, C. M., Goldsmith, J., Leroux, A., and Cui, E. (2024).
  *Functional Data Analysis with R*. CRC Press.
- Goldsmith, J., Zipunnikov, V., and Schrack, J. (2015). Generalized
  Multilevel Function-on-Scalar Regression and Principal Component
  Analysis. *Biometrics*, 71(2), 344–353.
- Ramsay, J. O. and Silverman, B. W. (2005). *Functional Data Analysis*,
  2nd Edition. Springer.
- Ivanescu, A. E., Staicu, A.-M., Scheipl, F., and Greven, S. (2015).
  Penalized Function-on-Function Regression. *Computational Statistics*,
  30(2), 539–568.
- Scheipl, F., Staicu, A.-M., and Greven, S. (2015). Functional Additive
  Mixed Models. *Journal of Computational and Graphical Statistics*,
  24(2), 477–501.
- Carpenter, B., Gelman, A., Hoffman, M. D., et al. (2017). Stan: A
  Probabilistic Programming Language. *Journal of Statistical Software*,
  76(1), 1–32.

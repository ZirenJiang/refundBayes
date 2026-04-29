# Joint FPCA Modeling in refundBayes

## Introduction

This vignette describes the **joint FPCA** option that is available
across the `refundBayes` regression functions:

- [`sofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/sofr_bayes.md)
  – Bayesian Scalar-on-Function Regression,
- [`fcox_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fcox_bayes.md)
  – Bayesian Functional Cox Regression,
- [`fofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fofr_bayes.md)
  – Bayesian Function-on-Function Regression.

In every case, the option is exposed via the same `joint_FPCA` argument:
a logical (`TRUE` / `FALSE`) vector of length equal to the number of
functional predictors. When the $i$-th entry is `TRUE`, the observed
functional predictor $W_{i}(s)$ is *not* used directly as a covariate;
instead, the regression model is built on top of a **functional
principal component analysis (FPCA) sub-model** for $W_{i}(s)$, and the
FPC scores are sampled jointly with the regression coefficients.

The methodology mirrors **Section 4 of Jiang, Crainiceanu, and Cui
(2025), *Tutorial on Bayesian Functional Regression Using Stan*,
*Statistics in Medicine* 44(20–22), e70265**. The default option is
(`joint_FPCA = NULL`) which do not adopt joint modelling.

## Why Joint FPCA?

In standard scalar-on-function or functional Cox regression, the
integral $\int W_{i}(s)\beta(s)\, ds$ is plugged in *as if* $W_{i}(s)$
were observed without error. In practice, the curves are sampled at
finitely many points and are typically corrupted by measurement noise.
Treating $W_{i}(s)$ as exact ignores this uncertainty and biases the
regression coefficient $\beta( \cdot )$ toward zero (the classical
errors-in-variables attenuation), and underestimates the posterior
credible-interval width.

The joint FPCA approach replaces $W_{i}(s)$ by a low-rank FPCA
representation

$$W_{i}(s)\; = \;\mu(s)\; + \;\sum\limits_{j = 1}^{J}\xi_{ij}\,\phi_{j}(s)\; + \;\epsilon_{i}(s),$$

and treats the subject-specific scores $\xi_{ij}$ as **parameters** that
are sampled together with the regression coefficients. The FPC
eigenfunctions $\phi_{j}(s)$ are estimated once with
[`refund::fpca.sc()`](https://rdrr.io/pkg/refund/man/fpca.sc.html) and
held fixed, and the initial frequentist FPCA scores
${\widehat{\xi}}_{ij}$ are used as the **prior mean** for $\xi_{ij}$.
This propagates the FPCA uncertainty into the posterior of
$\beta( \cdot )$, so the regression credible bands are correctly
inflated.

## The Joint Bayesian Model We Adopted

The model has two layers that share the FPC scores $\xi_{ij}$:

1.  an **FPCA likelihood** on the observed functional predictor; and
2.  the **outcome regression** (SoFR / FCox / FoFR), in which the scores
    $\xi_{ij}$ enter the linear predictor.

Throughout this section we assume one functional predictor for clarity;
the multi-predictor case is identical with one set of FPCA quantities
per functional predictor.

### FPCA likelihood

For each subject $i = 1,\ldots,n$ and each grid point $s_{m}$,
$m = 1,\ldots,M$,

$$W_{i}\left( s_{m} \right)\; = \;\sum\limits_{j = 1}^{J}\xi_{ij}\,\phi_{j}\left( s_{m} \right)\; + \;\epsilon_{i}\left( s_{m} \right),\qquad\epsilon_{i}\left( s_{m} \right)\overset{\text{iid}}{\sim}N\!\left( 0,\sigma_{e}^{2} \right),$$

so that

$$\log p\left( \mathbf{W} \mid {\mathbf{ξ}},\mathbf{\Phi},\sigma_{e} \right)\; = \; - \, n\, M\,\log\sigma_{e}\; - \;\frac{1}{2\sigma_{e}^{2}}\sum\limits_{i = 1}^{n}\sum\limits_{m = 1}^{M}\!\left\{ W_{i}\left( s_{m} \right) - \sum\limits_{j = 1}^{J}\xi_{ij}\,\phi_{j}\left( s_{m} \right) \right\}^{2}.$$

In Stan code, this corresponds to

``` stan
target += - N_num * M_num_i * log(sigma_e_i)
          - sum((xi_i * Phi_mat_i - M_mat_i)^2) / (2 * sigma_e_i^2);
```

Here `M_mat_i` is the observed $n \times M$ matrix of
$W_{i}\left( s_{m} \right)$ values, `Phi_mat_i` is the $J \times M$
matrix of fixed eigenfunctions, and `xi_i` is the $n \times J$ matrix of
FPC score parameters.

### FPC score prior

Each FPC score is given an independent Gaussian prior centered on the
initial frequentist FPCA score with eigenvalue-determined scale,

$$\xi_{ij}\; \sim \; N\!\left( {\widehat{\xi}}_{ij},\,\lambda_{j}^{2} \right),\qquad i = 1,\ldots,n,\;\; j = 1,\ldots,J.$$

In Stan code, the kernel is implemented as

``` stan
for (nj in 1:J_num_i) {
   target += - N_num * log(lambda_i[nj])
             - sum((xi_i[, nj] - xi_hat_i[, nj])^2) / (2 * lambda_i[nj]^2);
   target += inv_gamma_lpdf(lambda_i[nj]^2 | 0.001, 0.001);
}
target += inv_gamma_lpdf(sigma_e_i^2 | 0.001, 0.001);
```

The eigenvalue scales $\lambda_{j}^{2}$ and the residual scale
$\sigma_{e}^{2}$ both receive weakly informative inverse-Gamma priors.
The prior mean ${\widehat{\xi}}_{ij}$ is computed once via
[`refund::fpca.sc()`](https://rdrr.io/pkg/refund/man/fpca.sc.html) on
the observed functional data.

### Plugging the FPCA into the linear predictor

Because $W_{i}(s) = \sum_{j}\xi_{ij}\,\phi_{j}(s)$ in the FPCA
representation, the functional integral that drives every regression
model in `refundBayes` becomes

$$\int_{\mathcal{S}}W_{i}(s)\,\beta(s)\, ds\; = \;\sum\limits_{j = 1}^{J}\xi_{ij}\,\underset{\, \equiv \,\alpha_{j}}{\underbrace{\int_{\mathcal{S}}\phi_{j}(s)\,\beta(s)\, ds}}.$$

Expanding $\beta(s)$ in the spline basis
$\beta(s) = \sum_{k = 1}^{K}b_{k}\,\psi_{k}(s)$, we have

$$\alpha_{j}\; = \;\sum\limits_{k = 1}^{K}b_{k}\int_{\mathcal{S}}\phi_{j}(s)\,\psi_{k}(s)\, ds\; \equiv \;\sum\limits_{k = 1}^{K}\mathbf{X}_{jk}\, b_{k},$$

where the $J \times K$ FPCA-spline cross-product matrix $\mathbf{X}$ is
approximated by the Riemann sum

$$\mathbf{X}_{jk}\; \approx \;\sum\limits_{m = 1}^{M}L_{m}\,\phi_{j}\left( s_{m} \right)\,\psi_{k}\left( s_{m} \right)$$

using the same integration weights $L_{m}$ as in the standard regression
formulation. After the spectral reparameterization of $\mathbf{X}$ via
$\mathbf{S} = \int{\mathbf{ψ}}''(s){\mathbf{ψ}}''(s)^{t}\, ds$ (the
penalty matrix), the design matrix $\mathbf{X}$ is split into a
penalized random-effect block ${\widetilde{\mathbf{X}}}_{r}$ (size
$K_{r} \times J$) and an unpenalized fixed-effect block
${\widetilde{\mathbf{X}}}_{f}$ (size $K_{f} \times J$). The linear
predictor for SoFR / FCox is then

$$\eta_{i}\; = \;\eta_{0}\; + \;\mathbf{Z}_{i}^{t}{\mathbf{γ}}\; + \;{\mathbf{ξ}}_{i}^{t}\!\left( {\widetilde{\mathbf{X}}}_{r}^{t}\,{\widetilde{\mathbf{b}}}_{r}\; + \;{\widetilde{\mathbf{X}}}_{f}^{t}\,{\widetilde{\mathbf{b}}}_{f} \right),$$

with
${\widetilde{\mathbf{b}}}_{r} \sim N\left( \mathbf{0},\sigma_{b}^{2}\mathbf{I} \right)$
and ${\widetilde{\mathbf{b}}}_{f}$ assigned non-informative priors —
exactly the same prior structure used for the regression coefficients in
the no-joint-FPCA case. This is the chunk in the Stan code:

``` stan
mu += xi_i * (X_mat_r_i' * br_i + X_mat_f_i' * bf_i);
```

For the function-on-function regression model the integration over $s$
is combined with a response-domain expansion in the basis
${\mathbf{ψ}}(t)$ of size $K_{\text{resp}}$, so that the bivariate
coefficient
$\beta(s,t) = \sum_{k = 1}^{K_{\text{resp}}}\beta_{k}(s)\,\psi_{k}(t)$
leads to

$$\int W_{i}(s)\,\beta(s,t)\, ds\; = \;{\mathbf{ξ}}_{i}^{t}\!\left( {\widetilde{\mathbf{X}}}_{r}^{t}\,\mathbf{B}_{r}\; + \;{\widetilde{\mathbf{X}}}_{f}^{t}\,\mathbf{B}_{f} \right){\mathbf{ψ}}(t),$$

where now $\mathbf{B}_{r}$ is a $K_{r} \times K_{\text{resp}}$ matrix
and $\mathbf{B}_{f}$ is a $K_{f} \times K_{\text{resp}}$ matrix of
bivariate coefficients. The resulting Stan model contribution is

``` stan
mu += (xi_i * (X_mat_r_i' * br_i + X_mat_f_i' * bf_i)) * Psi_mat;
```

This is the only structural difference relative to the joint-FPCA SoFR /
FCox model: the spline coefficients become *matrices* indexed by both
the predictor-domain and the response-domain bases, and the same
row-wise response-domain penalty as in the standard FoFR is applied.

### Full Bayesian model

The complete joint model (one functional predictor) reads

$$\left\{ \begin{array}{l}
{\textbf{𝐎𝐮𝐭𝐜𝐨𝐦𝐞\ 𝐫𝐞𝐠𝐫𝐞𝐬𝐬𝐢𝐨𝐧}:\quad\eta_{i} = \eta_{0} + \mathbf{Z}_{i}^{t}{\mathbf{γ}} + {\mathbf{ξ}}_{i}^{t}\!\left( {\widetilde{\mathbf{X}}}_{r}^{t}\,{\widetilde{\mathbf{b}}}_{r} + {\widetilde{\mathbf{X}}}_{f}^{t}\,{\widetilde{\mathbf{b}}}_{f} \right)} \\
{\textbf{𝐅𝐏𝐂𝐀\ 𝐥𝐢𝐤𝐞𝐥𝐢𝐡𝐨𝐨𝐝}:\quad W_{i}\left( s_{m} \right) \mid {\mathbf{ξ}}_{i},\mathbf{\Phi},\sigma_{e} \sim N\!\left( \sum\limits_{j}\xi_{ij}\phi_{j}\left( s_{m} \right),\,\sigma_{e}^{2} \right)} \\
{\textbf{𝐅𝐏𝐂\ 𝐬𝐜𝐨𝐫𝐞\ 𝐩𝐫𝐢𝐨𝐫}:\quad\xi_{ij} \sim N\left( {\widehat{\xi}}_{ij},\,\lambda_{j}^{2} \right)} \\
{\textbf{𝐇𝐲𝐩𝐞𝐫𝐩𝐫𝐢𝐨𝐫𝐬}:\quad{\widetilde{\mathbf{b}}}_{r} \sim N\left( \mathbf{0},\sigma_{b}^{2}\mathbf{I} \right),\;\;\sigma_{b}^{2},\,\sigma_{e}^{2},\,\lambda_{j}^{2} \sim \text{Inv-Gamma}(0.001,0.001)} \\
{\textbf{𝐎𝐮𝐭𝐜𝐨𝐦𝐞\ 𝐝𝐢𝐬𝐭𝐫𝐢𝐛𝐮𝐭𝐢𝐨𝐧}:\quad Y_{i} \sim p\left( Y_{i} \mid \eta_{i} \right)\quad\text{(Gaussian\ /\ binomial\ /\ Cox\ /\ functional)}} \\
\end{array} \right.$$

The four building blocks (outcome regression, FPCA likelihood, score
prior, hyperpriors) are *independent of the family*; only the outcome
distribution in the last line changes. This is precisely why the joint
FPCA option is available across SoFR, FCox, and FoFR with the same
`joint_FPCA` argument.

### Note on the prior mean of $\xi_{ij}$

The prior is centered on the initial frequentist FPCA score
${\widehat{\xi}}_{ij}$ returned by
[`refund::fpca.sc()`](https://rdrr.io/pkg/refund/man/fpca.sc.html),
**not** on $0$. This corresponds to using the standard FPCA fit as a
working informative prior. As $\lambda_{j}^{2}$ becomes large the prior
reverts to a diffuse Gaussian and the joint model reduces (in spirit) to
plugging in the sample FPC scores; for moderate $\lambda_{j}^{2}$ the
posterior trades off the FPCA fit and the regression likelihood. The
observed-data matrix $W_{i}\left( s_{m} \right)$ is passed to Stan in
its original (uncentered) form, exactly as in Section 4 of the Tutorial
supplement.

### Identifying the functional data and the integration weights

The joint FPCA design matrix $\mathbf{X}_{jk}$ is built with the *same*
Riemann-sum integration weights as the no-joint-FPCA design matrix.
Specifically, the convention is:

- in `s(tindex, by = lmat * wmat, ...)`, the **rightmost** variable in
  the product (here `wmat`) is treated as the observed functional data
  $W_{i}( \cdot )$;
- the product of the remaining variables (here `lmat`) is treated as the
  Riemann-sum integration weight matrix.

This matches the Tutorial supplementary code and means that you do not
need to change your formulas when turning joint FPCA on or off; only the
`joint_FPCA` argument changes.

## How to enable joint FPCA

Every regression function in `refundBayes` accepts the same argument:

``` r
joint_FPCA = c(TRUE, FALSE, ...)   # one entry per s(...) functional term
```

The default is `joint_FPCA = NULL`, which is treated as
`rep(FALSE, n_func)` and gives exactly the no-joint-FPCA behavior of the
respective regression vignettes. Below we walk through one example for
each of
[`sofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/sofr_bayes.md),
[`fcox_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fcox_bayes.md),
and
[`fofr_bayes()`](https://zirenjiang.github.io/refundBayes/reference/fofr_bayes.md).

## Example 1: Joint FPCA in `sofr_bayes()`

We reuse the simulation setup from the SoFR vignette and switch on joint
FPCA for the (single) functional predictor.

### Simulate Data

``` r
set.seed(123)
n <- 100
M <- 50
tgrid <- seq(0, 1, length.out = M)
dt    <- tgrid[2] - tgrid[1]
tmat  <- matrix(rep(tgrid, each = n), nrow = n)
lmat  <- matrix(dt, nrow = n, ncol = M)

# Smooth latent trajectory + measurement noise on the functional predictor
D_true <- t(apply(matrix(rnorm(n * M), n, M), 1, cumsum)) / sqrt(M)
wmat   <- D_true + matrix(rnorm(n * M, sd = 0.3), nrow = n) ## noisy

beta_true <- sin(2 * pi * tgrid)
X1        <- rnorm(n)
eta       <- 0.5 * X1 + D_true %*% (beta_true * dt)   ## regression on D, not W
prob      <- plogis(eta)
y         <- rbinom(n, 1, prob)

data.SoFR <- data.frame(y = y, X1 = X1)
data.SoFR$tmat <- tmat
data.SoFR$lmat <- lmat
data.SoFR$wmat <- wmat
```

Notice that the outcome is generated from the **latent** trajectory
$D_{i}(s)$ but only the **noisy**
$W_{i}(s) = D_{i}(s) + \epsilon_{i}(s)$ is observed. This is exactly the
setting that motivates joint FPCA.

### Fit the joint FPCA SoFR model

``` r
library(refundBayes)

fit_sofr_joint <- sofr_bayes(
  formula    = y ~ X1 + s(tmat, by = lmat * wmat, bs = "cc", k = 10),
  data       = data.SoFR,
  family     = binomial(),
  joint_FPCA = c(TRUE),       ## << turn joint FPCA on for the (only) functional term
  niter      = 1500,
  nwarmup    = 500,
  nchain     = 3,
  ncores     = 3
)
```

For comparison, the **no**-joint-FPCA fit uses the same call without the
`joint_FPCA` argument:

``` r
fit_sofr_plain <- sofr_bayes(
  formula = y ~ X1 + s(tmat, by = lmat * wmat, bs = "cc", k = 10),
  data    = data.SoFR,
  family  = binomial(),
  niter   = 1500, nwarmup = 500, nchain = 3, ncores = 3
)
```

The two posterior estimates of $\beta(s)$ can be plotted together:

``` r
plot(fit_sofr_joint)
plot(fit_sofr_plain)
```

Joint-FPCA credible bands are typically wider in the regions where the
measurement noise on the functional predictor is informative.

### Joint-FPCA quantities exposed in the Stan fit

When `joint_FPCA = TRUE` for the $i$-th functional predictor, the Stan
fit additionally exposes:

| Stan parameter | Meaning                                 |
|:---------------|:----------------------------------------|
| `xi_i`         | $n \times J$ matrix of joint FPC scores |
| `lambda_i`     | $J$-vector of FPC eigenvalue SDs        |
| `sigma_e_i`    | scalar, FPC residual SD                 |

These can be inspected with
[`rstan::extract()`](https://mc-stan.org/rstan/reference/stanfit-method-extract.html)
for further analysis (for example, plotting the posterior of the FPC
scores or the eigenvalue SDs).

## Example 2: Joint FPCA in `fcox_bayes()`

The functional Cox setup is identical to the SoFR setup, with one extra
piece (the censoring vector). The joint FPCA model is wired in by
setting `joint_FPCA = c(TRUE)` on a model that has one functional
predictor:

``` r
library(refundBayes)

## Use the example dataset shipped with the FCox vignette
Func_Cox_Data <- readRDS("data/example_data_Cox.rds")
Func_Cox_Data$wmat <- Func_Cox_Data$MIMS
Func_Cox_Data$cens <- 1 - Func_Cox_Data$event

fit_cox_joint <- fcox_bayes(
  formula    = survtime ~ X1 + s(tmat, by = lmat * wmat, bs = "cc", k = 10),
  data       = Func_Cox_Data,
  cens       = Func_Cox_Data$cens,
  joint_FPCA = c(TRUE),
  niter      = 5000,
  nwarmup    = 2000,
  nchain     = 1,
  ncores     = 1
)
```

The Stan code generated under the hood matches Section 4 of the Tutorial
supplement: the linear predictor

$$\eta_{i}\; = \;\mathbf{Z}_{i}^{t}{\mathbf{γ}} + {\mathbf{ξ}}_{i}^{t}\!\left( {\widetilde{\mathbf{X}}}_{r}^{t}{\widetilde{\mathbf{b}}}_{r} + {\widetilde{\mathbf{X}}}_{f}^{t}{\widetilde{\mathbf{b}}}_{f} \right)$$

is plugged into the Cox log-likelihood

$$\ell(\mathbf{Y},{\mathbf{δ}})\; = \;\sum\limits_{i = 1}^{n}\left\lbrack \left( 1 - \delta_{i} \right)\!\left\{ \log h_{0}\left( Y_{i} \right) + \eta_{i} - H_{0}\left( Y_{i} \right)e^{\eta_{i}} \right\} + \delta_{i}\!\left\{ - H_{0}\left( Y_{i} \right)e^{\eta_{i}} \right\} \right\rbrack,$$

with the same M-spline / I-spline baseline-hazard construction as in the
no-joint case. The FPCA likelihood and the FPC-score prior are appended
to the joint Stan target function exactly as described in *The Joint
Bayesian Model We Adopted* above.

### Inspecting the generated Stan code

``` r
fcox_code <- fcox_bayes(
  formula    = survtime ~ X1 + s(tmat, by = lmat * wmat, bs = "cc", k = 10),
  data       = Func_Cox_Data,
  cens       = Func_Cox_Data$cens,
  joint_FPCA = c(TRUE),
  runStan    = FALSE
)
cat(fcox_code$stancode)
```

You will see the FPCA-related declarations in the `data` block
(`Phi_mat_1`, `xi_hat_1`, `M_mat_1`, `J_num_1`, `M_num_1`, `X_mat_r_1`,
`X_mat_f_1`), the FPCA-related parameters in the `parameters` block
(`xi_1`, `lambda_1`, `sigma_e_1`), and the FPCA log-likelihood + score
prior contributions in the `model` block.

## Example 3: Joint FPCA in `fofr_bayes()`

For function-on-function regression, joint FPCA on the functional
predictor is constructed in exactly the same way — only the contribution
to ${\mathbf{μ}}_{i}(t)$ now carries the response-domain basis
$\mathbf{\Psi}$:

$$\int W_{i}(s)\,\beta(s,t)\, ds\; = \;{\mathbf{ξ}}_{i}^{t}\!\left( {\widetilde{\mathbf{X}}}_{r}^{t}\,\mathbf{B}_{r} + {\widetilde{\mathbf{X}}}_{f}^{t}\,\mathbf{B}_{f} \right){\mathbf{ψ}}(t),$$

with $\mathbf{B}_{r},\mathbf{B}_{f}$ as **matrices** of size
$K_{r} \times K_{\text{resp}}$ and $K_{f} \times K_{\text{resp}}$
respectively. Both directions of smoothness ($s$ via the random-effect
reparameterization and $t$ via the response-domain penalty) are imposed
exactly as in the no-joint FoFR model.

### Simulated FoFR with measurement error on the predictor

``` r
library(refundBayes)
set.seed(42)

n <- 200
L <- 30
M <- 30
sindex <- seq(0, 1, length.out = L)
tindex <- seq(0, 1, length.out = M)

# Smooth latent functional predictor + measurement noise
D_true <- matrix(0, nrow = n, ncol = L)
for (i in 1:n) {
  D_true[i, ] <- rnorm(1) * sin(2 * pi * sindex) +
                 rnorm(1) * cos(2 * pi * sindex) +
                 rnorm(1) * sin(4 * pi * sindex)
}
X_func <- D_true + matrix(rnorm(n * L, sd = 0.3), nrow = n)   ## noisy

age <- rnorm(n)

# True coefficient functions
beta_true  <- outer(sin(2 * pi * sindex), cos(2 * pi * tindex))
alpha_true <- 0.5 * sin(pi * tindex)

# Generate response from the latent D_true (not from X_func!)
signal_scalar <- outer(age, alpha_true)
signal_func   <- (D_true %*% beta_true) / L
epsilon       <- matrix(rnorm(n * M, sd = 0.3), nrow = n)
Y_mat         <- signal_scalar + signal_func + epsilon

dat <- data.frame(age = age)
dat$Y_mat  <- Y_mat
dat$X_func <- X_func
dat$sindex <- matrix(rep(sindex, n), nrow = n, byrow = TRUE)
```

### Fit the joint FPCA FoFR model

``` r
fit_fofr_joint <- fofr_bayes(
  formula     = Y_mat ~ age + s(sindex, by = X_func, bs = "cr", k = 10),
  data        = dat,
  joint_FPCA  = c(TRUE),
  spline_type = "bs",
  spline_df   = 10,
  niter       = 2000,
  nwarmup     = 1000,
  nchain      = 3,
  ncores      = 3
)
```

The joint FPCA contributes the additional Stan parameters `xi_1`
($n \times J$), `lambda_1` ($J$-vector), and `sigma_e_1` (scalar), and
adds the FPCA log-likelihood and prior to the model block; the bivariate
coefficient $\widehat{\beta}(s,t)$ is reconstructed from
$\mathbf{B}_{r},\mathbf{B}_{f}$ exactly as in the no-joint FoFR case.

### Visualize the bivariate coefficient

``` r
beta_est <- apply(fit_fofr_joint$bivar_func_coef[[1]], c(2, 3), mean)

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
image(sindex, tindex, beta_true,
      xlab = "s (predictor domain)", ylab = "t (response domain)",
      main = expression("True " * beta(s,t)),
      col = hcl.colors(64, "Blue-Red 3"))
image(sindex, tindex, beta_est,
      xlab = "s (predictor domain)", ylab = "t (response domain)",
      main = expression("Joint-FPCA " * hat(beta)(s,t)),
      col = hcl.colors(64, "Blue-Red 3"))
```

## Practical Recommendations

- **When to use joint FPCA**. Joint FPCA is most useful when the
  observed functional predictor $W_{i}(s)$ contains substantial
  measurement noise relative to the signal in the regression, or when
  the curves are sparse / irregularly observed. For dense, low-noise
  functional data the difference is small.

- **Number of FPC components $J$**. Controlled by the default
  variance-explained criterion in
  [`refund::fpca.sc()`](https://rdrr.io/pkg/refund/man/fpca.sc.html). A
  value such that cumulative variance explained exceeds 95 % is a
  reasonable default. A larger $J$ increases the number of joint FPC
  score parameters ($n \times J$) and may slow down sampling.

- **Convergence**. Joint FPCA introduces $n \times J$ extra latent
  parameters. Multiple chains and the standard `rstan` diagnostics
  (traceplots, $\widehat{R}$, bulk / tail ESS) are recommended:

  ``` r
  rstan::traceplot(fit_sofr_joint$stanfit, pars = c("sigma_e_1", "lambda_1"))
  print(fit_sofr_joint$stanfit, pars = c("sigma_e_1", "lambda_1"))
  ```

- **Multiple functional predictors**. The argument is a vector with one
  entry per `s(...)` functional term. Setting some entries to `TRUE` and
  others to `FALSE` is supported: each functional predictor is treated
  independently.

- **Prior strength**. The default `Inv-Gamma(0.001, 0.001)` priors on
  $\lambda_{j}^{2}$ and $\sigma_{e}^{2}$ are weakly informative. If the
  posterior of $\xi_{i}$ is shrunk too hard toward
  ${\widehat{\xi}}_{i}$, consider loosening the prior by editing the
  Stan code (set `runStan = FALSE`, modify the `inv_gamma_lpdf(...)`
  lines, and pass the resulting program to
  [`rstan::stan()`](https://mc-stan.org/rstan/reference/stan.html)
  directly).

## References

- Jiang, Z., Crainiceanu, C., and Cui, E. (2025). Tutorial on Bayesian
  Functional Regression Using Stan. *Statistics in Medicine*, 44(20–22),
  e70265. (See Section 4 for the joint FPCA formulation that this
  vignette implements.)
- Crainiceanu, C. M., Goldsmith, J., Leroux, A., and Cui, E. (2024).
  *Functional Data Analysis with R*. CRC Press.
- Goldsmith, J., Bobb, J., Crainiceanu, C. M., Caffo, B., and Reich, D.
  (2011). Penalized Functional Regression. *Journal of Computational and
  Graphical Statistics*, 20(4), 830–851.
- Yao, F., Mueller, H.-G., and Wang, J.-L. (2005). Functional Data
  Analysis for Sparse Longitudinal Data. *Journal of the American
  Statistical Association*, 100(470), 577–590.
- Goldsmith, J., Greven, S., and Crainiceanu, C. M. (2013). Corrected
  Confidence Bands for Functional Data Using Principal Components.
  *Biometrics*, 69(1), 41–51.
- Carpenter, B., Gelman, A., Hoffman, M. D., et al. (2017). Stan: A
  Probabilistic Programming Language. *Journal of Statistical Software*,
  76(1), 1–32.

# Bayesian Scalar-on-Function Regression

Fit the Bayesian Scalar-on-Function Regression (SoFR) model using Stan.

## Usage

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

## Arguments

- formula:

  Functional regression formula, with the same syntax as that in the R
  mgcv package.

- data:

  A data frame containing data of all scalar and functional variables
  used in the model.

- family:

  Distribution of the outcome variable. Currently support "gaussian" and
  "binomial".

- joint_FPCA:

  A True/False vector of the same length of the number of functional
  predictors, indicating whether to jointly model the functional
  predictor via FPCA together with the regression model. When the entry
  is `TRUE`, the corresponding observed functional predictor is replaced
  by an FPCA representation, and the FPC scores are sampled jointly with
  the regression coefficients, following Section 4 of Jiang et al.
  (2025). Default to `NULL` (no joint FPCA, equivalent to
  `rep(FALSE, n_func)`).

- intercept:

  True/False variable for whether include an intercept term in the
  linear predictor. Default to TRUE.

- runStan:

  True/False variable for whether to run the Stan program. If False, the
  function only generates the Stan code and data.

- niter:

  Total number of Bayesian iterations.

- nwarmup:

  Number of warmup (burnin) iterations for posterior sampling.

- nchain:

  Number of chains for posterior sampling. Default to 3.

- ncores:

  Number of cores to use when executing the chains in parallel. Default
  to 1.

## Value

A list containing:

- stanfit:

  The Stan fit object.

- spline_basis:

  Basis functions used to reconstruct the functional coefficients from
  posterior samples.

- stancode:

  A character string containing the code to fit the Stan model.

- standate:

  A list containing the data to fit the Stan model.

- int:

  A vector containing posterior samples of the intercept term.

- scalar_coef:

  A matrix containing posterior samples of scalar coefficients, where
  each row is one sample and each column is one variable.

- func_coef:

  A list containing posterior samples of functional coefficients. Each
  element is a matrix, where each row is one sample and each column is
  one location of the functional domain.

- family:

  Distribution of the outcome variable.

## Details

The Bayesian SoFR model is implemented following the tutorial by Jiang
et al., 2025. The model is specified using the same syntax as in the R
mgcv package.

## References

Jiang, Z., Crainiceanu, C., and Cui, E. (2025). Tutorial on Bayesian
Functional Regression Using Stan. *Statistics in Medicine*, 44(20-22),
e70265.

## Author

Erjia Cui <ecui@umn.edu>, Ziren Jiang <jian0746@umn.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data for a Gaussian SoFR model
set.seed(1)
n  <- 100  # number of subjects
L  <- 50   # number of functional domain points
Lindex <- seq(0, 1, length.out = L)       # functional domain grid
X_func <- matrix(rnorm(n * L), nrow = n)  # functional predictor (n x L)
age    <- rnorm(n)                         # scalar predictor
beta_true <- sin(pi * Lindex)         # true functional coefficient
eta <- X_func %*% beta_true / L
Y <- eta + 0.5 * age + rnorm(n, sd = 0.5)

dat <- data.frame(Y = Y, age = age)
dat$X_func  <- X_func
dat$Lindex  <- matrix(rep(Lindex, n), nrow = n, byrow = TRUE)

# Fit Gaussian SoFR
fit_sofr <- sofr_bayes(
  formula = Y ~ age + s(Lindex, by = X_func, bs = "cr", k = 10),
  data    = dat,
  family  = "gaussian",
  niter   = 2000,
  nwarmup = 1000,
  nchain  = 3
)

# Summarise and plot estimated functional coefficient
summary(fit_sofr)
plot(fit_sofr)

# Fit binomial SoFR
prob <- plogis(X_func %*% beta_true / L)
Y_bin <- rbinom(n, 1, prob)
dat$Y_bin <- Y_bin
fit_bin <- sofr_bayes(
  formula = Y_bin ~ s(Lindex, by = X_func, bs = "cr", k = 10),
  data    = dat,
  family  = "binomial",
  niter   = 2000,
  nwarmup = 1000,
  nchain  = 3
)
} # }
```

# Bayesian Functional Cox Regression

Fit the Bayesian Functional Cox Regression model using Stan.

## Usage

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

## Arguments

- formula:

  Functional regression formula, with the same syntax as that in the R
  mgcv package.

- data:

  A data frame containing data of all scalar and functional variables
  used in the model.

- cens:

  A vector indicating censoring status (1 = event observed, 0 =
  censored). Must be the same length as the number of observations.

- joint_FPCA:

  A True/False vector of the same length of the number of functional
  predictors, indicating whether to jointly model the functional
  predictor via FPCA together with the survival model. When the entry is
  `TRUE`, the corresponding observed functional predictor is replaced by
  an FPCA representation, and the FPC scores are sampled jointly with
  the regression coefficients, following Section 4 of Jiang et al.
  (2025). Default to `NULL` (no joint FPCA, equivalent to
  `rep(FALSE, n_func)`).

- intercept:

  True/False variable for whether include an intercept term in the
  linear predictor. Default to FALSE.

- runStan:

  True/False variable for whether to run the Stan program. If False, the
  function only generates the Stan code and data.

- niter:

  Total number of Bayesian iterations. Default to 3000.

- nwarmup:

  Number of warmup (burnin) iterations for posterior sampling. Default
  to 1000.

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

- standata:

  A list containing the data to fit the Stan model.

- int:

  A vector containing posterior samples of the intercept term (NULL for
  Cox models by default).

- scalar_coef:

  A matrix containing posterior samples of scalar coefficients, where
  each row is one sample and each column is one variable.

- func_coef:

  A list containing posterior samples of functional coefficients. Each
  element is a matrix, where each row is one sample and each column is
  one location of the functional domain.

- baseline_hazard:

  Posterior samples of baseline hazard parameters.

- family:

  Family type: "Cox".

## Details

The Bayesian Functional Cox model extends the scalar-on-function
regression framework to survival outcomes with right censoring. The
model is specified using similar syntax as in the R mgcv package.

## References

Jiang, Z., Crainiceanu, C., and Cui, E. (2025). Tutorial on Bayesian
Functional Regression Using Stan. *Statistics in Medicine*, 44(20-22),
e70265.

## Author

Erjia Cui <ecui@umn.edu>, Ziren Jiang <jian0746@umn.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate survival data with a functional predictor
set.seed(1)
n  <- 150  # number of subjects
L  <- 50   # number of functional domain points
Lindex <- seq(0, 1, length.out = L)       # functional domain grid
X_func <- matrix(rnorm(n * L), nrow = n)  # functional predictor (n x L)
age    <- rnorm(n)                         # scalar predictor

# True functional effect and linear predictor
beta_true <- cos(pi * Lindex)
lp <- X_func %*% beta_true / L + age

# Generate survival times from an exponential baseline hazard
time  <- rexp(n, rate = exp(lp))
cens_time <- runif(n, min = 0.5, max = 3)
obs_time  <- pmin(time, cens_time)
cens_ind  <- as.integer(time <= cens_time)  # 1 = event, 0 = censored

dat <- data.frame(obs_time = obs_time, age = age)
dat$X_func <- X_func
dat$Lindex <- matrix(rep(Lindex, n), nrow = n, byrow = TRUE)

# Fit the Bayesian Functional Cox model
fit_cox <- fcox_bayes(
  formula = obs_time ~ age + s(Lindex, by = X_func, bs = "cr", k = 10),
  data    = dat,
  cens    = cens_ind,
  niter   = 2000,
  nwarmup = 1000,
  nchain  = 3
)

# Summarise scalar coefficients and plot functional coefficient
summary(fit_cox)
plot(fit_cox)
} # }
```

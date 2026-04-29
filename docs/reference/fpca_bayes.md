# Bayesian Functional Principal Component Analysis

Fit the Bayesian Functional Principal Component Analysis (FPCA) model
using Stan.

## Usage

``` r
fpca_bayes(
  formula,
  data,
  npc = NULL,
  runStan = TRUE,
  niter = 3000,
  nwarmup = 1000,
  nchain = 3,
  ncores = 1,
  spline_type = "bs",
  spline_df = 10
)
```

## Arguments

- formula:

  Functional formula of the form `Y_mat ~ 1`, where `Y_mat` is the
  functional response variable stored as a matrix column of `data`.

- data:

  A data frame containing the functional response variable used in the
  model.

- npc:

  Number of functional principal components. If NULL, it is selected
  automatically by the initial
  [`refund::fpca.sc()`](https://rdrr.io/pkg/refund/man/fpca.sc.html) fit
  based on the percentage of variance explained. Default to NULL.

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

- spline_type:

  Type of spline basis for modelling the population mean function.
  Default to "bs".

- spline_df:

  Degrees of freedom for the spline basis for modelling the population
  mean function. Default to 10.

## Value

A list containing:

- stanfit:

  The Stan fit object.

- stancode:

  A character string containing the code to fit the Stan model.

- standata:

  A list containing the data to fit the Stan model.

- mu:

  A matrix of posterior samples of the population mean function, where
  each row is one sample and each column is one location of the
  functional domain.

- efunctions:

  A matrix of the (fixed) eigenfunctions from the initial FPCA used as
  the FPC basis. Each column is one eigenfunction.

- scores:

  A 3-d array of posterior samples of FPC scores with dimensions
  (n_samples x n_subjects x npc).

- evalues:

  A matrix of posterior samples of FPC eigenvalue standard deviations,
  where each row is one sample and each column is one principal
  component.

- sigma:

  A vector of posterior samples of the residual standard deviation.

- family:

  Family type: "fpca".

## Details

The Bayesian FPCA model is implemented following the tutorial by Jiang
et al., 2025. The model decomposes a dense functional response
\\Y_i(t)\\ into a smooth population mean function \\\mu(t)\\ and a sum
of functional principal components, \$\$Y_i(t) = \mu(t) +
\sum\_{k=1}^{K} \xi\_{ik} \phi_k(t) + \epsilon_i(t),\$\$ where
\\\phi_k(t)\\ are orthonormal eigenfunctions obtained from an initial
frequentist FPCA fit (used as a fixed basis), and the mean function
\\\mu(t)\\, the FPC scores \\\xi\_{ik}\\, the eigenvalues \\\lambda_k\\,
and the residual variance are estimated via posterior sampling. The
population mean function is modelled with a penalized spline using the
same syntax as in the R mgcv package.

## References

Jiang, Z., Crainiceanu, C., and Cui, E. (2025). Tutorial on Bayesian
Functional Regression Using Stan. *Statistics in Medicine*, 44(20-22),
e70265.

## Author

Erjia Cui <ecui@umn.edu>, Ziren Jiang <jian0746@umn.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate functional data with two underlying principal components
set.seed(1)
n  <- 100   # number of subjects
M  <- 50    # number of functional observation points
tindex <- seq(0, 1, length.out = M)

# True mean function and eigenfunctions
mu_true <- sin(pi * tindex)
phi1    <- sqrt(2) * sin(2 * pi * tindex)
phi2    <- sqrt(2) * cos(2 * pi * tindex)

# Simulate scores and noisy observations
xi1 <- rnorm(n, 0, sqrt(2))
xi2 <- rnorm(n, 0, sqrt(0.5))
Y_mat <- matrix(rep(mu_true, n), nrow = n, byrow = TRUE) +
         outer(xi1, phi1) + outer(xi2, phi2) +
         matrix(rnorm(n * M, sd = 0.3), nrow = n)

dat <- data.frame(inx = 1:n)
dat$Y_mat <- Y_mat

# Fit the Bayesian FPCA model
fit_fpca <- fpca_bayes(
  formula     = Y_mat ~ 1,
  data        = dat,
  spline_type = "bs",
  spline_df   = 10,
  niter       = 2000,
  nwarmup     = 1000,
  nchain      = 3
)

# Posterior mean of the mean function
mu_est <- apply(fit_fpca$mu, 2, mean)
plot(tindex, mu_est, type = "l", ylab = expression(hat(mu)(t)))

# Posterior means of the FPC scores and eigenvalues
scores_est  <- apply(fit_fpca$scores, c(2, 3), mean)
evalues_est <- apply(fit_fpca$evalues, 2, mean)
} # }
```

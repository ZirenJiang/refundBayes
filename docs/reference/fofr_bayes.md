# Bayesian Function-on-Function Regression

Fit the Bayesian Function-on-Function Regression (FoFR) model using
Stan.

## Usage

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

## Arguments

- formula:

  Functional regression formula, with the same syntax as that in the R
  mgcv package. The response must be a matrix (functional response) and
  at least one `s(..., by = ...)` term must be present for the
  functional predictor(s).

- data:

  A data frame containing data of all scalar and functional variables
  used in the model.

- joint_FPCA:

  A True/False vector of the same length of the number of functional
  predictors, indicating whether jointly modeling FPCA for the
  functional predictors. Default to NULL.

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

  Type of spline basis for modelling the residual process and the
  response-domain component of the bivariate coefficient. Default to
  "bs".

- spline_df:

  Degrees of freedom for the spline basis for modelling the
  response-domain component. Default to 10.

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

- scalar_func_coef:

  A 3-d array (n_samples x P x M) containing posterior samples of scalar
  predictor coefficient functions. Each slice `[,p,]` is the coefficient
  function for the p-th scalar predictor. NULL if no scalar predictors.

- bivar_func_coef:

  A list of 3-d arrays. Each element corresponds to one functional
  predictor and is an array of dimension (n_samples x S_grid x M),
  representing posterior samples of the bivariate coefficient function
  \\\beta(s,t)\\.

- func_coef:

  A 3-d array for the scalar predictor coefficient functions, stored for
  compatibility with the `plot.refundBayes` method. Same as
  `scalar_func_coef`.

- family:

  Model family: "fofr".

## Details

The Bayesian FoFR model extends the function-on-scalar regression (FoSR)
framework by allowing functional predictors in addition to (optional)
scalar predictors. The bivariate coefficient function \\\beta(s,t)\\ is
represented using a tensor product of the predictor-domain spline basis
and the response-domain spline basis. The model is specified using the
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
# Simulate data for a Function-on-Function Regression model
set.seed(1)
n  <- 100  # number of subjects
L  <- 30   # number of functional predictor domain points
M  <- 30   # number of functional response domain points
sindex <- seq(0, 1, length.out = L)  # predictor domain grid
tindex <- seq(0, 1, length.out = M)  # response domain grid

# Functional predictor
X_func <- matrix(rnorm(n * L), nrow = n)

# Scalar predictor
age <- rnorm(n)

# True bivariate coefficient beta(s, t)
beta_true <- outer(sin(2 * pi * sindex), cos(2 * pi * tindex))

# True scalar coefficient function
alpha_true <- sin(pi * tindex)

# Generate functional response: Y(t) = age * alpha(t) + integral X(s) beta(s,t) ds + error
Y_mat <- outer(age, alpha_true) + X_func %*% beta_true / L +
          matrix(rnorm(n * M, sd = 0.3), nrow = n)

dat <- data.frame(age = age)
dat$Y_mat  <- Y_mat
dat$X_func <- X_func
dat$sindex <- matrix(rep(sindex, n), nrow = n, byrow = TRUE)

# Fit the Bayesian FoFR model
fit_fofr <- fofr_bayes(
  formula    = Y_mat ~ age + s(sindex, by = X_func, bs = "cr", k = 10),
  data       = dat,
  spline_type = "bs",
  spline_df  = 10,
  niter      = 2000,
  nwarmup    = 1000,
  nchain     = 3
)

# Examine the bivariate coefficient beta(s, t) (posterior mean)
beta_est <- apply(fit_fofr$bivar_func_coef[[1]], c(2, 3), mean)
image(sindex, tindex, beta_est,
      xlab = "s (predictor domain)", ylab = "t (response domain)",
      main = expression(hat(beta)(s, t)))

# Fit a FoFR model with functional predictors only (no scalar predictors)
fit_fofr2 <- fofr_bayes(
  formula    = Y_mat ~ s(sindex, by = X_func, bs = "cr", k = 10),
  data       = dat,
  niter      = 2000,
  nwarmup    = 1000,
  nchain     = 3
)
} # }
```

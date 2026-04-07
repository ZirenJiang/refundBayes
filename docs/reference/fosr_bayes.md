# Bayesian Function-on-Scalar Regression

Fit the Bayesian Function-on-Scalar Regression (FOSR) model using Stan.

## Usage

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

## Arguments

- formula:

  Functional regression formula, with the same syntax as that in the R
  mgcv package.

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

  Type of spline basis for modelling the residual process.

- spline_df:

  Degrees of freedom for the spline basis for modelling the residual
  process.

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

The Bayesian FOSR model is implemented following the tutorial by Jiang
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
# Simulate data for a Function-on-Scalar Regression model
set.seed(1)
n  <- 100   # number of subjects
M  <- 50    # number of functional response observation points
tindex <- seq(0, 1, length.out = M)  # response functional domain grid

# Scalar predictors
age <- rnorm(n)
sex <- rbinom(n, 1, 0.5)

# True coefficient functions
beta_age <- sin(2 * pi * tindex)
beta_sex <- cos(2 * pi * tindex)

# Generate functional response (n x M matrix)
epsilon  <- matrix(rnorm(n * M, sd = 0.3), nrow = n)
Y_mat    <- outer(age, beta_age) + outer(sex, beta_sex) + epsilon

dat <- data.frame(age = age, sex = sex)
dat$Y_mat <- Y_mat

# Fit the Bayesian FoSR model
fit_fosr <- fosr_bayes(
  formula    = Y_mat ~ age + sex,
  data       = dat,
  spline_type = "bs",
  spline_df  = 10,
  niter      = 2000,
  nwarmup    = 1000,
  nchain     = 3
)

# Plot estimated coefficient functions
plot(fit_fosr)
} # }
```

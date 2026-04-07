#----------------------------------------------------------------------------
#' Bayesian Function-on-Scalar Regression
#'
#' Fit the Bayesian Function-on-Scalar Regression (FOSR) model using Stan.
#' 
#' The Bayesian FOSR model is implemented following the tutorial by Jiang et al., 2025.
#' The model is specified using the same syntax as in the R mgcv package.
#'
#' @param formula Functional regression formula, with the same syntax as that in the R mgcv package.
#' @param data A data frame containing data of all scalar and functional variables used in the model.
#' @param joint_FPCA A True/False vector of the same length of the number of functional predictors, indicating whether jointly modeling FPCA for the functional predictors. Default to NULL.
#' @param spline_type Type of spline basis for modelling the residual process.
#' @param spline_df Degrees of freedom for the spline basis for modelling the residual process.
#' @param runStan True/False variable for whether to run the Stan program. If False, the function only generates the Stan code and data.
#' @param niter Total number of Bayesian iterations.
#' @param nwarmup Number of warmup (burnin) iterations for posterior sampling.
#' @param nchain Number of chains for posterior sampling. Default to 3.
#' @param ncores Number of cores to use when executing the chains in parallel. Default to 1.
#'
#' @return A list containing:
#' \item{stanfit}{The Stan fit object.}
#' \item{spline_basis}{Basis functions used to reconstruct the functional coefficients from posterior samples.}
#' \item{stancode}{A character string containing the code to fit the Stan model.}
#' \item{standate}{A list containing the data to fit the Stan model.}
#' \item{int}{A vector containing posterior samples of the intercept term.}
#' \item{scalar_coef}{A matrix containing posterior samples of scalar coefficients, where each row is one sample and each column is one variable.}
#' \item{func_coef}{A list containing posterior samples of functional coefficients. Each element is a matrix, where each row is one sample and each column is one location of the functional domain.}
#' \item{family}{Distribution of the outcome variable.}
#' 
#' @author Erjia Cui \email{ecui@@umn.edu}, Ziren Jiang \email{jian0746@@umn.edu}
#'
#' @references Jiang, Z., Crainiceanu, C., and Cui, E. (2025). Tutorial on Bayesian Functional Regression Using Stan. \emph{Statistics in Medicine}, 44(20-22), e70265.
#'
#' @examples
#' \dontrun{
#' # Simulate data for a Function-on-Scalar Regression model
#' set.seed(1)
#' n  <- 100   # number of subjects
#' M  <- 50    # number of functional response observation points
#' tindex <- seq(0, 1, length.out = M)  # response functional domain grid
#'
#' # Scalar predictors
#' age <- rnorm(n)
#' sex <- rbinom(n, 1, 0.5)
#'
#' # True coefficient functions
#' beta_age <- sin(2 * pi * tindex)
#' beta_sex <- cos(2 * pi * tindex)
#'
#' # Generate functional response (n x M matrix)
#' epsilon  <- matrix(rnorm(n * M, sd = 0.3), nrow = n)
#' Y_mat    <- outer(age, beta_age) + outer(sex, beta_sex) + epsilon
#'
#' dat <- data.frame(age = age, sex = sex)
#' dat$Y_mat <- Y_mat
#'
#' # Fit the Bayesian FoSR model
#' fit_fosr <- fosr_bayes(
#'   formula    = Y_mat ~ age + sex,
#'   data       = dat,
#'   spline_type = "bs",
#'   spline_df  = 10,
#'   niter      = 2000,
#'   nwarmup    = 1000,
#'   nchain     = 3
#' )
#'
#' # Plot estimated coefficient functions
#' plot(fit_fosr)
#' }
#' 
#' @import mgcv
#' @import splines2
#' @import rstan
#' @importFrom brms brmsformula
#' @importFrom refund fpca.face
#' @importFrom stats gaussian
#' @export fosr_bayes

fosr_bayes <- function(formula, data, 
                       joint_FPCA = NULL, runStan = TRUE, 
                       niter = 3000, nwarmup = 1000, nchain = 3, ncores = 1, 
                       spline_type = "bs", spline_df = 10){
  
  
  ## Set family to be "functional"
  family <- "functional"
  
  func_parameter = list(type = spline_type,
                        df = spline_df)
  intercept = TRUE
  
  # Organize the formula as the input for refundBayesdata() and refundBayesmodel()
  formula_use <- refundBayesformula(formula)
  if(is.null(joint_FPCA)){
    if(length(formula_use$func_var) == 0){
      func_comp <- NULL
    }else{
      func_comp <- rep(FALSE,length(formula_use$func_var))
    }
  }else{
    func_comp <- joint_FPCA
  }
  
  if (length(formula_use$func_var) != length(func_comp)) {
    stop("Length of func_comp must be the total number of functional predictors!")
  }
  
  # Construct the data and code (data, transformed data, parameters, and transformed parameters blocks) as input for Stan
  datacode_use <- refundBayesdata(formula_use, data, family, func_comp, intercept, func_parameter = func_parameter)
  Standata_use <- datacode_use$standata
  Stancode_data <- datacode_use$stancode_data
  Stancode_transdata <- datacode_use$stancode_transdata
  Stancode_para <- datacode_use$stancode_para
  Stancode_transpara <- datacode_use$stancode_transpara
  
  # Construct the model block as input for Stan
  if(family!="functional"){
    Stancode_model <- refundBayesmodel(formula_use, data, family, func_comp, intercept)
  }else{
    Stancode_model <- refundBayesmodel_functional(formula_use, data, family, func_comp, intercept)
  }
  
  
  # Combine all blocks of the Stan code
  Stancode_function <- NULL ## no manual function needed for SoFR
  Stancode_use <- paste0(Stancode_function, "\n",
                         Stancode_data, "\n",
                         Stancode_transdata, "\n",
                         Stancode_para, "\n",
                         Stancode_transpara, "\n",
                         Stancode_model)
  
  
  if(runStan){
    # Run the Stan program
    fit <- stan(model_code = Stancode_use,
                data = Standata_use,
                iter = niter,
                warmup = nwarmup,
                chains = nchain,
                cores = ncores)
    
    # Extract the basis used to reconstruct estimated functional coefficients
    
    if(!is.null(func_comp)){ ## Have functional predictor
      func_effect_basis <- extract_basis(formula_use, data, func_comp)
    }else{
      func_effect_basis <- NA
    }
    # Extract posterior samples from the Stan fit object
    fit.samp <- rstan::extract(fit, permuted = TRUE)
    
    func_coef_pos = array(dim = c(dim(fit.samp$b)[1], Standata_use$P_num, dim(Standata_use$Psi_mat)[2]))
    for(iter in 1:dim(fit.samp$b)[1]){
      func_coef_pos[iter,,] = t(fit.samp$beta[iter,,]) %*% Standata_use$Psi_mat
    }
    func_coef_pos_mean = apply(func_coef_pos, c(2,3), mean)
    
    
    
    # Reconstruct the estimated functional coefficients based on the posterior sampling
    # Each functional predictor sample is stored as a matrix, where each each row is one sample, each column is one location of the functional domain
    if(!is.null(func_comp)){ ## Have functional predictor
      func_coef <- recon_fun_coef(basis = func_effect_basis, ## basis from extract_basis()
                                  fit.samp = fit.samp, ## samples from stan() and extract()
                                  func_comp = func_comp, 
                                  trans.mat = datacode_use$trans.mat)
    }else{
      func_coef = func_coef_pos
    }
    
    
    
    # Extract the estimated intercept as a vector
    # if(intercept == TRUE){
    #   int <- fit.samp$eta_0
    # }else{
    #   int <- NULL
    # }
  }else{
    # Do not run the Stan program
    fit <- NA ## no sampling
    if(!is.null(func_comp)){ ## Have functional predictor
      func_effect_basis <- extract_basis(formula_use, data, func_comp)
    }else{
      func_effect_basis = NA
    }
    fit.samp <- NA ## no sampling
    func_coef <- NA ## no sampling
    if (length(formula_use$scalar_var) > 0){
      scalar_coef <- NA ## no sampling
    }else{
      scalar_coef <- NULL
    }
    # if(intercept == TRUE){
    #   int <- NA ## no sampling
    # }else{
    #   int <- NULL
    # }
  }
  
  # Organize the results
  res.fit = list(stanfit = fit,
                 spline_basis = func_effect_basis,
                 stancode = Stancode_use,
                 standata = Standata_use,
                 func_coef = func_coef,
                 family = family)
  class(res.fit) = "refundBayes"
  return(res.fit)
}

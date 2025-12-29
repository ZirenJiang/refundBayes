#----------------------------------------------------------------------------
#' Bayesian Functional Cox Regression
#'
#' Fit the Bayesian Functional Cox Regression model using Stan.
#' 
#' The Bayesian Functional Cox model extends the scalar-on-function regression
#' framework to survival outcomes with right censoring. The model is specified 
#' using similar syntax as in the R mgcv package.
#'
#' @param formula Functional regression formula, with the same syntax as that in the R mgcv package.
#' @param data A data frame containing data of all scalar and functional variables used in the model.
#' @param cens A vector indicating censoring status (1 = event observed, 0 = censored). Must be the same length as the number of observations.
#' @param joint_FPCA A True/False vector of the same length of the number of functional predictors, indicating whether jointly modeling FPCA for the functional predictors. Default to NULL.
#' @param intercept True/False variable for whether include an intercept term in the linear predictor. Default to FALSE.
#' @param runStan True/False variable for whether to run the Stan program. If False, the function only generates the Stan code and data.
#' @param niter Total number of Bayesian iterations. Default to 3000.
#' @param nwarmup Number of warmup (burnin) iterations for posterior sampling. Default to 1000.
#' @param nchain Number of chains for posterior sampling. Default to 3.
#' @param ncores Number of cores to use when executing the chains in parallel. Default to 1.
#'
#' @return A list containing:
#' \item{stanfit}{The Stan fit object.}
#' \item{spline_basis}{Basis functions used to reconstruct the functional coefficients from posterior samples.}
#' \item{stancode}{A character string containing the code to fit the Stan model.}
#' \item{standata}{A list containing the data to fit the Stan model.}
#' \item{int}{A vector containing posterior samples of the intercept term (NULL for Cox models by default).}
#' \item{scalar_coef}{A matrix containing posterior samples of scalar coefficients, where each row is one sample and each column is one variable.}
#' \item{func_coef}{A list containing posterior samples of functional coefficients. Each element is a matrix, where each row is one sample and each column is one location of the functional domain.}
#' \item{baseline_hazard}{Posterior samples of baseline hazard parameters.}
#' \item{family}{Family type: "Cox".}
#' 
#' @author Erjia Cui \email{ecui@@umn.edu}, Ziren Jiang \email{jian0746@@umn.edu}
#'
#' @references Jiang, Z., Crainiceanu, C., and Cui, E. (2025). Tutorial on Bayesian Functional Regression Using Stan. \emph{Statistics in Medicine}, 44(20-22), e70265.
#'
#' @import mgcv
#' @import splines2
#' @import rstan
#' @importFrom brms brmsformula
#' @importFrom refund fpca.face
#' @export fcox_bayes

fcox_bayes <- function(formula, data, cens, 
                       joint_FPCA = NULL, intercept = FALSE, runStan = TRUE, 
                       niter = 3000, nwarmup = 1000, nchain = 3, ncores = 1){
  
  # Family is fixed to Cox for this function
  family <- "Cox"
  
  # Check the input arguments
  if(is.null(cens)){
    stop("Censoring indicator 'cens' must be provided for Cox regression!")
  }
  
  if(length(cens) != nrow(data)){
    stop("Length of 'cens' must equal the number of observations in 'data'!")
  }
  
  if(!all(cens %in% c(0, 1))){
    stop("'cens' must be a binary vector with values 0 (censored) and 1 (event)!")
  }
  
  # Organize the formula as the input for refundBayesdata() and refundBayesmodel()
  formula_use <- refundBayesformula(formula)
  
  if(is.null(joint_FPCA)){
    if(length(formula_use$func_var) == 0){
      func_comp <- NULL
    }else{
      func_comp <- rep(FALSE, length(formula_use$func_var))
    }
  }else{
    func_comp <- joint_FPCA
  }
  
  if (length(formula_use$func_var) != length(func_comp)) {
    stop("Length of joint_FPCA must be the total number of functional predictors!")
  }
  
  # Construct the data and code (data, transformed data, parameters, and transformed parameters blocks) as input for Stan
  datacode_use <- refundBayesdata(formula_use, data, family, func_comp, intercept, cens)
  Standata_use <- datacode_use$standata
  Stancode_data <- datacode_use$stancode_data
  Stancode_transdata <- datacode_use$stancode_transdata
  Stancode_para <- datacode_use$stancode_para
  Stancode_transpara <- datacode_use$stancode_transpara
  
  # Construct the model block as input for Stan
  Stancode_model <- refundBayesmodel(formula_use, data, family, func_comp, intercept)
  
  # Construct the functions block for Cox likelihood
  Stancode_function <- paste0("functions {\n",
                              "  real cox_log_lhaz(real y, real log_mu, real bhaz, real cbhaz) {\n",
                              "    return log(bhaz) + log_mu;\n",
                              "  }\n",
                              "  real cox_log_lccdf(real y, real log_mu, real bhaz, real cbhaz) {\n",
                              "    return - cbhaz * exp(log_mu);\n",
                              "  }\n",
                              "  real cox_log_lcdf(real y, real log_mu, real bhaz, real cbhaz) {\n",
                              "    return log1m_exp(cox_log_lccdf(y | log_mu, bhaz, cbhaz));\n",
                              "  }\n",
                              "  real cox_log_lpdf(real y, real log_mu, real bhaz, real cbhaz) {\n",
                              "    return cox_log_lhaz(y, log_mu, bhaz, cbhaz) + cox_log_lccdf(y | log_mu, bhaz, cbhaz);\n",
                              "  }\n",
                              "}\n")
  
  # Combine all blocks of the Stan code
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
    
    # Reconstruct the estimated functional coefficients based on the posterior sampling
    if(!is.null(func_comp)){ ## Have functional predictor
      func_coef <- recon_fun_coef(basis = func_effect_basis,
                                  fit.samp = fit.samp,
                                  func_comp = func_comp, 
                                  trans.mat = datacode_use$trans.mat)
    }else{
      func_coef <- NA
    }
    
    # Extract the estimated scalar coefficients as a matrix
    if (length(formula_use$scalar_var) > 0){
      scalar_coef <- fit.samp$gamma
      colnames(scalar_coef) <- colnames(datacode_use$X_scalar)
    }else{
      scalar_coef <- NULL
    }
    
    # Extract the estimated intercept as a vector (typically NULL for Cox)
    if(intercept == TRUE){
      int <- fit.samp$eta_0
    }else{
      int <- NULL
    }
    
    # Extract baseline hazard parameters
    baseline_hazard <- list(
      bhaz = fit.samp$bhaz,
      cbhaz = fit.samp$cbhaz
    )
    
  }else{
    # Do not run the Stan program
    fit <- NA
    if(!is.null(func_comp)){ ## Have functional predictor
      func_effect_basis <- extract_basis(formula_use, data, func_comp)
    }else{
      func_effect_basis <- NA
    }
    fit.samp <- NA
    func_coef <- NA
    if (length(formula_use$scalar_var) > 0){
      scalar_coef <- NA
    }else{
      scalar_coef <- NULL
    }
    if(intercept == TRUE){
      int <- NA
    }else{
      int <- NULL
    }
    baseline_hazard <- NA
  }
  
  # Organize the results
  res.fit <- list(stanfit = fit,
                  spline_basis = func_effect_basis,
                  stancode = Stancode_use,
                  standata = Standata_use,
                  int = int,
                  scalar_coef = scalar_coef,
                  func_coef = func_coef,
                  baseline_hazard = baseline_hazard,
                  family = family)
  class(res.fit) <- "refundBayes"
  return(res.fit)
}


#----------------------------------------------------------------------------
#' Bayesian Scalar-on-Function Regression
#'
#' Fit the Bayesian Scalar-on-Function Regression (SoFR) model using Stan.
#' 
#' The Bayesian SoFR model is implemented following the tutorial by Jiang et al., 2025.
#' The model is specified using the same syntax as in the R mgcv package.
#'
#' @param formula Functional regression formula, with the same syntax as that in the R mgcv package.
#' @param data A data frame containing data of all scalar and functional variables used in the model.
#' @param family Distribution of the outcome variable. Currently support "gaussian" and "binomial".
#' @param joint_FPCA A True/False vector of the same length of the number of functional predictors, indicating whether jointly modeling FPCA for the functional predictors. Default to NULL.
#' @param intercept True/False variable for whether include an intercept term in the linear predictor. Default to TRUE.
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
#' @import mgcv
#' @import splines2
#' @import rstan
#' @importFrom brms brmsformula
#' @importFrom refund fpca.face
#' @export sofr_bayes

sofr_bayes <- function(formula, data, family = gaussian(), 
                       joint_FPCA = NULL, intercept = TRUE, runStan = TRUE, 
                       niter = 3000, nwarmup = 1000, nchain = 3, ncores = 1){
  
  # Currently support family
  family_list <- c("gaussian", "binomial")
  
  # Check the input arguments
  ## func_comp is a binary vector of length of the number of functional predictors 
  ## indicating whether to do a joint FPCA for each included functional predictor 
  if(is.null(joint_FPCA)){
    func_comp <- F ## by default FALSE
  }else{
    func_comp <- joint_FPCA
  }
  if(class(family) == "family"){
    family <- family$family
  }
  if(!family %in% family_list){
    stop(paste0("Family should be one of the following arguments: gaussian, binomial."))
  }
  
  # Organize the formula as the input for refundBayesdata() and refundBayesmodel()
  formula_use <- refundBayesformula(formula)
  if (length(formula_use$func_var) != length(func_comp)) {
    stop("Length of func_comp must be the total number of functional predictors!")
  }
  
  # Construct the data and code (data, transformed data, parameters, and transformed parameters blocks) as input for Stan
  datacode_use <- refundBayesdata(formula_use, data, family, func_comp, intercept)
  Standata_use <- datacode_use$standata
  Stancode_data <- datacode_use$stancode_data
  Stancode_transdata <- datacode_use$stancode_transdata
  Stancode_para <- datacode_use$stancode_para
  Stancode_transpara <- datacode_use$stancode_transpara
  
  # Construct the model block as input for Stan
  Stancode_model <- refundBayesmodel(formula_use, data, family, func_comp, intercept)
  
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
                chain = nchain,
                cores = ncores)
    
    # Extract the basis used to reconstruct estimated functional coefficients
    func_effect_basis <- extract_basis(formula_use, data, func_comp)
    
    # Extract posterior samples from the Stan fit object
    fit.samp <- rstan::extract(fit, permuted = TRUE)
    
    # Reconstruct the estimated functional coefficients based on the posterior sampling
    # Each functional predictor sample is stored as a matrix, where each each row is one sample, each column is one location of the functional domain
    func_coef <- recon_fun_coef(basis = func_effect_basis, ## basis from extract_basis()
                                fit.samp = fit.samp, ## samples from stan() and extract()
                                func_comp = func_comp, 
                                trans.mat = datacode_use$trans.mat)
    
    # Extract the estimated scalar coefficients as a matrix, where each each row is one sample, each column is one variable
    if (length(formula_use$scalar_var) > 0){
      scalar_coef <- fit.samp$gamma
      colnames(scalar_coef) <- colnames(datacode_use$X_scalar)
    }else{
      scalar_coef <- NULL
    }
    
    # Extract the estimated intercept as a vector
    if(intercept == TRUE){
      int <- fit.samp$eta_0
    }else{
      int <- NULL
    }
  }else{
    # Do not run the Stan program
    fit <- NA ## no sampling
    func_effect_basis <- extract_basis(formula_use, data, func_comp)
    fit.samp <- NA ## no sampling
    func_coef <- NA ## no sampling
    if (length(formula_use$scalar_var) > 0){
      scalar_coef <- NA ## no sampling
    }else{
      scalar_coef <- NULL
    }
    if(intercept == TRUE){
      int <- NA ## no sampling
    }else{
      int <- NULL
    }
  }
  
  # Organize the results
  return(list(stanfit = fit,
              spline_basis = func_effect_basis,
              stancode = Stancode_use,
              standata = Standata_use,
              int = int,
              scalar_coef = scalar_coef,
              func_coef = func_coef,
              family = family))
}

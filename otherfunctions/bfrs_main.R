#----------------------------------------------------------------------------
#' Bayesian functional regression model with Stan
#'
#' Fit the Bayesian functional regression model using R Stan. At this moment, bfrs supports the Scalar-on-Function Regression (SoFR),
#' functional Cox regression, and Function-on-Scalar Regression (FoSR).
#'
#' @param formula functional regression formula
#' @param data data.frame containing data of all variables used in the model
#' @param family type of outcome variable, bfrs currently support "guassian", "binomial", "Cox", and "functional"
#' @param joint_FPCA A vector of the same length of the functional predictors, indicating whether jointly modeling FPCA for the functional predictors.
#' @param intercept True/False variable for whether include an intercept term in the linear predictor
#' @param n.iter Total number of Bayesian iterations
#' @param n.warmup Number of warmup (burnin) iterations for posterior sampling
#' @param n.knots Number of cores to use when executing the chains in parallel
#' @param cens Indicator variable for censoring if family = "Cox"
#' @param func_parameter Specify the spline parameter for the FoSR.
#' @param runStan True/False variable for whether to run the Stan program. If false, only generate the Stan code and data.
#'
#' @import parallel
#' @import methods
#' @import stats
#' @import Rcpp
#' @import rstan
#' @import mgcv
#' @import refund
#' @import scam
#' @import splines2
#' @export bfrs
#'


bfrs <- function(formula, data, family = gaussian(), joint_FPCA=NULL,
                 intercept = TRUE,n.iter=5000,n.warmup=2000,n.knots=1, cens=NULL,func_parameter=NULL,runStan=TRUE,
                ...){

  # Currently support family
  family_list=c("gaussian", "binomial", "Cox", "functional")

  # Check the input arguments
  if(is.null(joint_FPCA)){
    func_comp=F
  }else{
    func_comp=joint_FPCA
  }
  if(class(family)=="family"){
    family=family$family
  }
  if(!family %in% family_list){
    stop(paste0("Family should be one of the following arguments: gaussian, binomial, Cox, functional!"))
  }
  formula_use = brfs_formula(formula)
  if ((length(formula_use$func_var)!=length(func_comp))&(family!="functional")) {
    stop("Length of func_comp must be the total number of functional predictors!")
  }

  # Construct the Stan data and code
  brfs_datafit=brfs_data(formula_use,data,family,func_comp,intercept,cens=cens,func_parameter=func_parameter)
  Standata_use=brfs_datafit$standata
  Stancode_data=brfs_datafit$stancode_data
  Stancode_para=brfs_datafit$stancode_para
  Stancode_transpara=brfs_datafit$stancode_transpara

  if(family!="functional"){
    if(length(formula_use$scalar_var)>0){
      Stancode_transdata="transformed data {\n   matrix[N_num, K_num] X_sc;\n   vector[K_num] mean_Xs;\n   "
      Stancode_transdata=paste0(Stancode_transdata, "for (i in 1:K_num) {\n      mean_Xs[i] = mean(Z_mat[, i]);
     X_sc[, i] = Z_mat[, i] - mean_Xs[i];\n   }\n}")
    }else{
      Stancode_transdata="transformed data {\n}"
    }
    Stancode_model=brfs_code_model(formula_use,data,family,func_comp,intercept)
  }else{
    Stancode_model=bfrs_code_functional(formula_use,Standata_use,family,func_comp,intercept)
    Stancode_transdata="transformed data {\n}"
  }
  if(family=="Cox"){
    Stancode_function=paste0("functions {\n   real cox_log_lhaz(real y, real log_mu, real bhaz, real cbhaz) {\n   ")
    Stancode_function=paste0(Stancode_function, "   return log(bhaz) + log_mu;\n   }\n")
    Stancode_function=paste0(Stancode_function, "   real cox_log_lccdf(real y, real log_mu, real bhaz, real cbhaz) {\n   ")
    Stancode_function=paste0(Stancode_function, "   return - cbhaz * exp(log_mu);\n   }\n")
    Stancode_function=paste0(Stancode_function, "   real cox_log_lcdf(real y, real log_mu, real bhaz, real cbhaz) {\n   ")
    Stancode_function=paste0(Stancode_function, "   return log1m_exp(cox_log_lccdf(y | log_mu, bhaz, cbhaz));\n   }\n")
    Stancode_function=paste0(Stancode_function, "   real cox_log_lpdf(real y, real log_mu, real bhaz, real cbhaz) {\n   ")
    Stancode_function=paste0(Stancode_function, "   return cox_log_lhaz(y, log_mu, bhaz, cbhaz) +
             cox_log_lccdf(y | log_mu, bhaz, cbhaz);\n   }\n}\n")
  }else{
    Stancode_function=NULL
  }
  Stancode_use=paste0(Stancode_function,
                      Stancode_data,"\n",
                      Stancode_transdata,"\n",
                      Stancode_para, "\n",
                      Stancode_transpara, "\n",
                      Stancode_model)

  # Run the Stan program
  if(runStan){
    fit <- stan(model_code=Stancode_use,data=Standata_use,iter = n.iter,warmup = n.warmup,chain=n.knots,cores = n.knots)
    if(family!="functional"){
      func_effect_basis=brfs_extract_basis(formula_use,data,func_comp)
      fit.samp = rstan::extract(fit, permuted = TRUE)
      func_effect=brfs_effect(basis = func_effect_basis,
                              fit.samp=fit.samp,
                              func_comp=func_comp,
                              trans.mat=brfs_datafit$trans.mat)
    }else{
      func_effect_basis=Standata_use$phi
      fit.samp = rstan::extract(fit, permuted = TRUE)
      func_effect_est=apply(fit.samp$b,c(1,2),function(x){
        x%*%func_effect_basis
      })
      func_effect=func_effect_est
    }
  }else{
    fit = NA
    if(family!="functional"){
      func_effect_basis=brfs_extract_basis(formula_use,data,func_comp)
      fit.samp = NA
      func_effect=NA
    }else{
      func_effect_basis=Standata_use$phi
      fit.samp = NA
      func_effect=NA

    }
  }

  return(list(stanfit = fit,
              spline_basis = func_effect_basis,
              Stancode = Stancode_use,
              Standata = Standata_use,
              func_effect = func_effect,
              family = family,
              scalar_pred = colnames(brfs_datafit$X_scalar)))
}






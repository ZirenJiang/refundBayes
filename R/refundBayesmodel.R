#' Functions to build the model block of the Stan code.
#'

#----------------------------------------------------------------------------
#' The main function for building the model block for SoFR and functional Cox model.
#----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
#' 
refundBayesmodel = function(formula, data, family, func_comp, intercept){
  result_code_model <- paste0("model{ \n")
  result_code_model <- paste0(result_code_model, "   vector[N_num] mu = rep_vector(0.0, N_num);\n")
  result_code_model <- paste0(result_code_model, "   //Linear predictor\n")
  result_code_model <- paste0(result_code_model, "   mu += ")

  if(intercept){ ## add intercept
    result_code_model <- paste0(result_code_model, "eta_0 ")
  }
  if(length(formula$scalar_var) > 0){ ## add scalar predictors
    result_code_model <- paste0(result_code_model, "+ X_sc * gamma")
  }
  if(is.null(func_comp)){ ## Do not have any functional predictor

  }else{
    for(i in 1:length(func_comp)){
      if(func_comp[i]){
        ## Joint FPCA: with X_mat_r_i (Kr x J_num) and X_mat_f_i (Kf x J_num),
        ## the contribution to the linear predictor for subject n is
        ##   sum_j xi_{n,j} * (X_mat_r_i' * br_i + X_mat_f_i' * bf_i)_j
        ## which in matrix form is xi_i * (X_mat_r_i' * br_i + X_mat_f_i' * bf_i).
        result_code_model <- paste0(result_code_model,
                                    "+ ", paste0("xi_", i), " * (",
                                    paste0("X_mat_r_", i), "' * ", paste0("br_", i),
                                    " + ", paste0("X_mat_f_", i), "' * ", paste0("bf_", i),
                                    ")")
      }
      if(!func_comp[i]){ ## add functional predictors
        result_code_model <- paste0(result_code_model, "+ ",paste0("X_mat_r_",i)," * ",paste0("br_",i),"+ ",paste0("X_mat_f_",i)," * ",paste0("bf_",i))
      }
    }
  }
  
  
  
  
  result_code_model <- paste0(result_code_model, ";\n")
  
  ## specify the likelihood for different outcome distributions
  if(family == "gaussian"){
    result_code_model <- paste0(result_code_model, "   for (n in 1:N_num) {\n     //Gaussian log-likelihood\n     target += normal_lpdf(Y[n]|mu[n],sigma);\n   }\n")
  }
  if(family == "binomial"){
    result_code_model <- paste0(result_code_model, "   for (n in 1:N_num) {\n     //Binomial log-likelihood\n     target += bernoulli_logit_lpmf(Y[n]|mu[n]);\n   }\n")
  }
  if(family == "Cox"){
    result_code_model <- paste0(result_code_model, "   //Fit the baseline hazard function\n")
    result_code_model <- paste0(result_code_model, "   vector[N_num] bhaz = Mbasis * c;\n")
    result_code_model <- paste0(result_code_model, "   //Fit the cumulative baseline hazard function\n")
    result_code_model <- paste0(result_code_model, "   vector[N_num] cbhaz = Ibasis * c;\n")
    result_code_model <- paste0(result_code_model, "   for (n in 1:N_num) {\n    if (cens[n] == 0) {\n   ")
    result_code_model <- paste0(result_code_model, "   target += cox_log_lpdf(Y[n] | mu[n], bhaz[n], cbhaz[n]);\n")
    result_code_model <- paste0(result_code_model, "    } else if (cens[n] == 1) {;   \n   ")
    result_code_model <- paste0(result_code_model, "   target += cox_log_lccdf(Y[n] | mu[n], bhaz[n], cbhaz[n]);\n")
    result_code_model <- paste0(result_code_model, "    } else if (cens[n] == -1) {   \n   ")
    result_code_model <- paste0(result_code_model, "   target += cox_log_lcdf(Y[n] | mu[n], bhaz[n], cbhaz[n]);\n    }\n   }\n")
  }
  
  result_code_model <- paste0(result_code_model, set_prior(formula, data, family, func_comp, intercept))
  result_code_model <- paste0(result_code_model, "}")
  return(result_code_model)
}

#----------------------------------------------------------------------------
#' The main function for building the model block for FoSR.
#----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
#' 
refundBayesmodel_functional = function(formula, standata, family, func_comp, intercept){
  result_code_model <- paste0("model{ \n")
  result_code_model <- paste0(result_code_model, "   matrix[N_num, M_num] mu;\n")
  result_code_model <- paste0(result_code_model, "   // Fitted mean matrix \n")
  result_code_model <- paste0(result_code_model, "   mu = X_mat * beta' * Psi_mat + xi * Phi_mat;\n")
  result_code_model <- paste0(result_code_model, "   // Log-likelihood for functional response\n")
  result_code_model <- paste0(result_code_model, "   target += - N_num * M_num * log(sigma_eps) / 2 - sum((mu - Y_mat)^2) / (2 * sigma_eps^2);\n")
  
  result_code_model <- paste0(result_code_model, "   // Prior for the penalized spline coefficients \n")
  result_code_model <- paste0(result_code_model, "   for(np in 1:P_num){\n")
  result_code_model <- paste0(result_code_model, "        target += (- beta[,np]' * S_mat * beta[,np]) / (2 * sigma[np]^2);\n")
  result_code_model <- paste0(result_code_model, "        target += inv_gamma_lpdf(sigma[np]^2|0.001,0.001);\n")
  result_code_model <- paste0(result_code_model, "   }\n")
  
  result_code_model <- paste0(result_code_model, "    // Prior for the FPCA scores \n")  
  result_code_model <- paste0(result_code_model, "   for(nj in 1:J_num){\n")
  result_code_model <- paste0(result_code_model, "        target += - N_num * log(lambda[nj]) / 2 - sum((xi[,nj])^2) / (2 * lambda[nj]^2);\n")
  result_code_model <- paste0(result_code_model, "        target += inv_gamma_lpdf(lambda[nj]^2|0.001,0.001);\n")
  result_code_model <- paste0(result_code_model, "   }\n")
  result_code_model <- paste0(result_code_model, "   target += inv_gamma_lpdf(sigma_eps^2|0.001,0.001);\n")
  result_code_model <- paste0(result_code_model, "}\n")
  return(result_code_model)
}

#----------------------------------------------------------------------------
#' A convenient function for setting the priors.
#----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
#' 
set_prior = function(formula, data, family, func_comp, intercept){
  result_code_model = ""
  if(intercept){ ## priors for the intercept term eta_0
    if(family == "gaussian"){
      result_code_model <- paste0("   target += normal_lpdf(eta_0 | plocation, pscale);\n")
    }
    if(family == "binomial"){
      result_code_model <- paste0("   target += student_t_lpdf(eta_0 | 3, plocation, pscale);\n")
    }
    if(family == "Cox"){
      result_code_model <- paste0("   target += normal_lpdf(eta_0 | real_inter, 10);\n")
      result_code_model <- paste0(result_code_model, "   target += dirichlet_lpdf(c | con_sbhaz);\n")
    }
  }
  if(is.null(func_comp)){ ## Do not have any functional predictor

  }else{
    for(i in 1:length(func_comp)){
      if(func_comp[i]){
        ## Joint FPCA likelihood and priors: this introduces a working FPCA
        ## sub-model for the functional predictor and shrinks the FPC scores
        ## toward the initial fpca.sc estimates xi_hat_i (Tutorial Section 4).
        result_code_model <- paste0(result_code_model,
                                    "   //Joint FPCA likelihood for functional predictor ", i, "\n")
        result_code_model <- paste0(result_code_model,
                                    "   target += -N_num * ", paste0("M_num_", i), " * log(", paste0("sigma_e_", i), ")",
                                    " - sum((", paste0("xi_", i), " * ", paste0("Phi_mat_", i), " - ", paste0("M_mat_", i), ")^2) / (2 * ", paste0("sigma_e_", i), "^2);\n")
        result_code_model <- paste0(result_code_model,
                                    "   //Prior for the joint FPCA scores for functional predictor ", i, "\n")
        result_code_model <- paste0(result_code_model,
                                    "   for (nj in 1:", paste0("J_num_", i), ") {\n")
        result_code_model <- paste0(result_code_model,
                                    "     target += -N_num * log(", paste0("lambda_", i), "[nj]) - sum((",
                                    paste0("xi_", i), "[, nj] - ", paste0("xi_hat_", i), "[, nj])^2) / (2 * ",
                                    paste0("lambda_", i), "[nj]^2);\n")
        result_code_model <- paste0(result_code_model,
                                    "     target += inv_gamma_lpdf(", paste0("lambda_", i), "[nj]^2 | 0.001, 0.001);\n")
        result_code_model <- paste0(result_code_model, "   }\n")
        result_code_model <- paste0(result_code_model,
                                    "   target += inv_gamma_lpdf(", paste0("sigma_e_", i), "^2 | 0.001, 0.001);\n")
        ## Random/fixed effect priors for the functional regression coefficient
        result_code_model <- paste0(result_code_model,
                                    "   target += std_normal_lpdf(", paste0("zbr_", i), ");\n")
        result_code_model <- paste0(result_code_model,
                                    "   target += inv_gamma_lpdf(", paste0("sigmabr_", i), "^2 | 0.0005, 0.0005);\n")
      }
      if(!func_comp[i]){ ## priors for the random effect coefficients zbr and sigmabr
        result_code_model <- paste0(result_code_model,"   target += std_normal_lpdf(", paste0("zbr_",i),");\n")
        result_code_model <- paste0(result_code_model,"   target += inv_gamma_lpdf(", paste0("sigmabr_",i),"^2 | 0.0005,0.0005);\n")
      }
    }
  }
  return(result_code_model)
}
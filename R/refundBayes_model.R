#' This document describes functions for building the Stan code.
#'
#'


#----------------------------------------------------------------------------
#' Building the Stan code for SoFR.
#----------------------------------------------------------------------------



brfs_code_model=function(formula,data,family,func_comp,intercept){

  result_code_model=paste0("model{ \n")
  result_code_model=paste0(result_code_model, "   vector[N_num] mu = rep_vector(0.0, N_num);\n")


  result_code_model=paste0(result_code_model, "   //Linear predictor\n")
  result_code_model=paste0(result_code_model, "   mu += ")

  if(intercept){
    result_code_model=paste0(result_code_model, "eta_0 ")
  }

  if(length(formula$scalar_var)>0){
    result_code_model=paste0(result_code_model, "+ X_sc * gamma")
  }



  for(i in 1:length(func_comp)){
    if(func_comp[i]){

    }
    if(!func_comp[i]){
      result_code_model=paste0(result_code_model, "+ ",paste0("X_mat_r_",i)," * ",paste0("br_",i),"+ ",paste0("X_mat_f_",i)," * ",paste0("bf_",i))
    }
  }
  result_code_model=paste0(result_code_model, ";\n")


  if(family=="gaussian"){
    result_code_model=paste0(result_code_model, "   for (n in 1:N_num) {\n     //Gaussian log-likelihood\n     target += normal_lpdf(Y[n]|mu[n],sigma);\n   }\n")
  }

  if(family=="binomial"){
    result_code_model=paste0(result_code_model, "   for (n in 1:N_num) {\n     //Binomial log-likelihood\n     target += bernoulli_logit_lpmf(Y[n]|mu[n]);\n   }\n")
  }

  if(family=="Cox"){
    result_code_model=paste0(result_code_model, "   //Fit the baseline hazard function\n")
    result_code_model=paste0(result_code_model, "   vector[N_num] bhaz = Mbasis * c;\n")

    result_code_model=paste0(result_code_model, "   //Fit the cumulative baseline hazard function\n")
    result_code_model=paste0(result_code_model, "   vector[N_num] cbhaz = Ibasis * c;\n")

    result_code_model=paste0(result_code_model, "   for (n in 1:N_num) {\n    if (cens[n] == 0) {\n   ")
    result_code_model=paste0(result_code_model, "   target += cox_log_lpdf(Y[n] | mu[n], bhaz[n], cbhaz[n]);\n")
    result_code_model=paste0(result_code_model, "    } else if (cens[n] == 1) {;   \n   ")
    result_code_model=paste0(result_code_model, "   target += cox_log_lccdf(Y[n] | mu[n], bhaz[n], cbhaz[n]);\n")
    result_code_model=paste0(result_code_model, "    } else if (cens[n] == -1) {   \n   ")
    result_code_model=paste0(result_code_model, "   target += cox_log_lcdf(Y[n] | mu[n], bhaz[n], cbhaz[n]);\n    }\n   }\n")
  }

  result_code_model=paste0(result_code_model, brfs_prior(formula,data,family,func_comp,intercept))


  result_code_model=paste0(result_code_model, "}")
  return(result_code_model)
}


#----------------------------------------------------------------------------
#' Building the Stan code for FoSR.
#----------------------------------------------------------------------------


bfrs_code_functional=function(formula,standata,family,func_comp,intercept){

  result_code_model=paste0("model{ \n")
  result_code_model=paste0(result_code_model, "   //Linear predictor\n")
  result_code_model=paste0(result_code_model, "   matrix[N_num, T] mu;\n")


  if(TRUE){
  #if(standata[["K_num"]]>1){
    result_code_model=paste0(result_code_model, "   matrix[K_num, T] betaT;\n")
  }else{
    result_code_model=paste0(result_code_model, "   row_vector[T] betaT;\n")
  }

  result_code_model=paste0(result_code_model, "   betaT = b *phi;\n")
  result_code_model=paste0(result_code_model, "   mu = Z_mat * betaT + zxi * psi_mat;\n")
  result_code_model=paste0(result_code_model, "   target += -N_num*T*log(sigma_epis)/2- sum((mu-Y)^2)/(2*sigma_epis^2);\n")
  result_code_model=paste0(result_code_model, "   for(nj in 1:ncol_psi_mat){\n")
  result_code_model=paste0(result_code_model, "     target += -N_num*log(sigma_k[nj])/2- sum((zxi[,nj])^2)/(2*sigma_k[nj]^2);\n")
  result_code_model=paste0(result_code_model, "     target += inv_gamma_lpdf(sigma_k[nj]^2|0.005,0.005);\n")
  result_code_model=paste0(result_code_model, "   }\n")

  result_code_model=paste0(result_code_model, "   for(nk in 1:K_num){\n")
  result_code_model=paste0(result_code_model, "     target += (- b[nk,]*S_beta* b[nk,]') / (2*sigma_b[nk]^2);\n")
  result_code_model=paste0(result_code_model, "     target += inv_gamma_lpdf(sigma_b[nk]^2|0.005,0.005);;\n")
  result_code_model=paste0(result_code_model, "   };\n")


  result_code_model=paste0(result_code_model, "   target +=  inv_gamma_lpdf(sigma_epis^2|0.005,0.005);\n")
  result_code_model=paste0(result_code_model, "}\n")
  return(result_code_model)
}



#----------------------------------------------------------------------------
#' Stan code for setting the priors.
#----------------------------------------------------------------------------




brfs_prior=function(formula,data,family,func_comp,intercept){


  if(intercept){
    if(family=="gaussian"){
      result_code_model=paste0("   target += normal_lpdf(eta_0 |plocation, pscale);\n")
    }
    if(family=="binomial"){
      result_code_model=paste0("   target += student_t_lpdf(eta_0 | 3, plocation, pscale);\n")
    }
    if(family=="Cox"){
      result_code_model=paste0("   target += normal_lpdf(eta_0 | real_inter, 10);\n")
      result_code_model=paste0(result_code_model,"   target += dirichlet_lpdf(c | con_sbhaz);\n")

    }
  }



  for(i in 1:length(func_comp)){
    if(func_comp[i]){

    }
    if(!func_comp[i]){
      result_code_model=paste0(result_code_model,"   target += std_normal_lpdf(",paste0("zbr_",i),");\n")
      result_code_model=paste0(result_code_model,"   target += inv_gamma_lpdf(",paste0("sigmabr_",i),"|0.0005,0.0005);\n")

      #result_code_model=paste0(result_code_model, "+ ",paste0("X_mat_r_",i)," * ",paste0("br_",i),"+ ",paste0("X_mat_f_",i)," * ",paste0("bf_",i))
    }
  }
  return(result_code_model)
}

#' Functions for reconstructing the estimated functional coefficient based on the posterior sampling
#'

#----------------------------------------------------------------------------
#' The function for reconstructing the SoFR fitted functional coefficient according to Stan posterior sampling
#----------------------------------------------------------------------------

recon_fun_coef = function(basis, fit.samp, func_comp, trans.mat){
  sp_samp <- list()
  for(inx in 1:length(func_comp)){
    eigendecomp <- basis[[inx]][["eigendecomp"]]
    base_mat <- basis[[inx]][["base_mat"]]
    E <- basis[[inx]][["E"]]
    para <- cbind(fit.samp[[paste0("br_", inx)]], fit.samp[[paste0("bf_", inx)]])
    para.trans <- para %*% trans.mat[[inx]]
    para.trans.untilde <- para.trans
    for(i in 1:nrow(para.trans)){
      para.trans.untilde[i,] <- (eigendecomp$vectors %*% diag(1/E)) %*% matrix(para.trans[i,], ncol=1)
    }
    curve_est <- matrix(ncol = nrow(base_mat), nrow = nrow(para.trans.untilde))
    for(i in 1:nrow(curve_est)){
      curve_est[i,] <- base_mat %*% para.trans.untilde[i,]
    }
    sp_samp[[inx]] <- curve_est
  }
  return(sp_samp)
}

#----------------------------------------------------------------------------
#' The function for reconstructing the FoSR fitted functional coefficient according to Stan posterior sampling
#----------------------------------------------------------------------------

recon_fun_coef_fosr = function(basis, b.samp){
  sp_samp <- list()
  for(inx in 1:length(func_comp)){
    eigendecomp <- basis[[inx]][["eigendecomp"]]
    base_mat <- basis[[inx]][["base_mat"]]
    E <- basis[[inx]][["E"]]
    para <- cbind(fit.samp[[paste0("br_",inx)]], fit.samp[[paste0("bf_",inx)]])
    para.trans <- para %*% trans.mat[[inx]]
    para.trans.untilde <- para.trans
    for(i in 1:nrow(para.trans)){
      para.trans.untilde[i,] <- (eigendecomp$vectors %*% diag(1/E)) %*% matrix(para.trans[i,], ncol=1)
    }
    curve_est <- matrix(ncol = nrow(base_mat), nrow = nrow(para.trans.untilde))
    for(i in 1:nrow(curve_est)){
      curve_est[i,] <- base_mat %*% para.trans.untilde[i,]
    }
    sp_samp[[inx]] <- curve_est
  }
  return(sp_samp)
}


#----------------------------------------------------------------------------
#' Generate the summary table for the Bayesian model
#' 
#' @param object A fitted object returned by \code{sofr_bayes()}.
#' @param prob Coverage probability for the reported confidence intervals. Defaults to 0.95. 
#' @param ... Other parameters
#' 
#' @return A list of two objects, the first is the summary table for the estimated scalar coefficients, the second is the plots for the estimated functional coefficients.
#' 
#' @export 
#' @method summary refundBayes


summary.refundBayes = function(object = NULL,...,prob = 0.95){
  res = list()
  res[[1]] = summary_scalar.bfrs(object ,prob)
  if(is.na(object$func_coef)){
    res[[2]] = NA
  }else{
    res[[2]] = plot.refundBayes(object ,prob , include = "both")
  }
  res
}

#' @keywords internal
#' @noRd
#' 
summary_scalar.bfrs = function(brfs.fit=NULL,prob = 0.95){
  q_upper = (1+prob)/2
  q_lower = (1-prob)/2
  if(is.null(brfs.fit$scalar_coef)){
    return(NULL)
  }else{
    scalar_samp = rstan::extract(brfs.fit$stanfit, permuted = TRUE)
    res_table = matrix(nrow = ncol(brfs.fit$scalar_coef),ncol=3)
    for(i in 1:ncol(brfs.fit$scalar_coef)){
      res_table[i,1]=mean(scalar_samp$gamma[,i])
      res_table[i,2]=stats::quantile(scalar_samp$gamma[,i],prob = q_lower)
      res_table[i,3]=stats::quantile(scalar_samp$gamma[,i],prob = q_upper)
    }
    res = as.data.frame(res_table)
    res = cbind(colnames(brfs.fit$scalar_coef),res)
    colnames(res) = c("Scalar Predictor", "mean", paste0(q_lower," quantile"), paste0(q_upper," quantile"))
    return(res)
  }
}

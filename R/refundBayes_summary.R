#----------------------------------------------------------------------------
#' Generate the summary table for the Bayesian model
#'
#' @export summary.bfrs
#' @export summary_scalar.bfrs


summary.bfrs = function(brfs.fit=NULL,prob = 0.95){

  res = list()

  res[[1]] = summary_scalar.bfrs(brfs.fit ,prob)

  res[[2]] = plot.bfrs(brfs.fit ,prob , include = "both")

  res
}



summary_scalar.bfrs = function(brfs.fit=NULL,prob = 0.95){

  q_upper = (1+prob)/2
  q_lower = (1-prob)/2

  if(length(brfs.fit$scalar_pred)==0){
    return(NULL)
  }else{
    scalar_samp = rstan::extract(brfs.fit$stanfit, permuted = TRUE)

    res_table = matrix(nrow = length(brfs.fit$scalar_pred),ncol=3)

    for(i in 1:length(brfs.fit$scalar_pred)){
      res_table[i,1]=mean(scalar_samp$gamma[,i])
      res_table[i,2]=quantile(scalar_samp$gamma[,i],prob = q_lower)
      res_table[i,3]=quantile(scalar_samp$gamma[,i],prob = q_upper)
    }

    res = as.data.frame(res_table)

    res = cbind(brfs.fit$scalar_pred,res)

    colnames(res) = c("Scalar Predictor", "mean", paste0(q_lower," quantile"), paste0(q_upper," quantile"))

    return(res)

  }

}

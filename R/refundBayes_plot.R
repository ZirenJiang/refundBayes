#----------------------------------------------------------------------------
#' Plot the estimated functional coefficients with the corresponding credible interval(s). 
#' @param x A fitted object returned by \code{sofr_bayes()}.
#' @param prob Coverage probability for the credible interval(s). Defaults to 0.95.
#' @param include Type of interval to include. \code{"pointwise"} produces pointwise credible intervals;
#'   \code{"CMA"} produces the CMA credible band; \code{"both"} produces both. Defaults to \code{"both"}.
#' @param ... Other parameters
#' 
#' @return A list of \code{ggplot} objects, one for each functional coefficient.
#' 
#' @import ggplot2
#' @import dplyr
#' @importFrom stats as.formula qnorm quantile sd
#' @export 
#' @method plot refundBayes

plot.refundBayes=function(x = NULL,...,prob = 0.95,include = "both"){

  if(x$family=="functional"){
    ## for functional outcomes (FoSR)
    plot.res=list()
    for(inx.effect in 1:dim(x$func_effect)[3]){
      curve_est=x$func_effect[,,inx.effect]
      mean.curve.est=apply(curve_est,1,mean)
      upper.curve.est.quantile=apply(curve_est,1,function(xx){stats::quantile(xx,probs = (1+prob)/2)})
      lower.curve.est.quantile=apply(curve_est,1,function(xx){stats::quantile(xx,probs = (1-prob)/2)})
      upper.curve.wald=upper.curve.est.quantile
      lower.curve.wald=lower.curve.est.quantile

      for(i in 1:length(upper.curve.wald)){
        upper.curve.wald[i] <- mean.curve.est[i]+qnorm((1+prob)/2)*sd(curve_est[,i])
        lower.curve.wald[i] <- mean.curve.est[i]-qnorm((1+prob)/2)*sd(curve_est[,i])
      }
      
      # We now only use the quantile approach to calculate the credible interval 
      upper.curve.est = upper.curve.est.quantile
      lower.curve.est = lower.curve.est.quantile

      plotdata <- data.frame(value = c(mean.curve.est,
                                  upper.curve.est,
                                  lower.curve.est),
                          xmat = c(1: length(upper.curve.est),
                                 1: length(upper.curve.est),
                                 1: length(upper.curve.est)),
                          type = c(rep("mean",length(upper.curve.est)),
                                 rep("upper",length(upper.curve.est)),
                                 rep("lower",length(upper.curve.est))))
      plot.res[[inx.effect]] <- ggplot2::ggplot(plotdata,aes(y = .data$value, x = .data$xmat))+
        ggplot2::geom_line(aes(type=.data$type,color=.data$type))+
        ggplot2::ylab(colnames(x$Standata[["X_s"]])[inx.effect])
    }
  }else{
    
    ## for scalar/survival outcomes (SoFR / functional Cox regression)
    n.func.coeff=length(x$func_coef)
    plot.res=list()
    for(inx.effect in 1:n.func.coeff){
      curve_est=x$func_coef[[inx.effect]]
      mean.curve.est=apply(curve_est,2,mean)
      upper.curve.est.quantile=apply(curve_est,2,function(xx){stats::quantile(xx,probs = (1+prob)/2)})
      lower.curve.est.quantile=apply(curve_est,2,function(xx){stats::quantile(xx,probs = (1-prob)/2)})

      if(include=="pointwise"){
        upper.curve.est = upper.curve.est.quantile
        lower.curve.est = lower.curve.est.quantile
        CMA.upper.curve.est = NA
        CMA.lower.curve.est = NA
      }
      if(include=="CMA"){
        func_coeff_sd=apply(curve_est,2,sd)
        func_coeff_extr=1:dim(curve_est)[1]
        func_coeff_est=apply(curve_est,2,mean)
        for(i in 1:dim(curve_est)[1]){
          func_coeff_extr[i]=max(abs(curve_est[i,]-func_coeff_est)/func_coeff_sd)
        }
        cutpoint=stats::quantile(func_coeff_extr,probs = prob)
        CMA.upper.curve.est = mean.curve.est + cutpoint*func_coeff_sd
        CMA.lower.curve.est = mean.curve.est - cutpoint*func_coeff_sd
        upper.curve.est = NA
        lower.curve.est = NA
      }
      if(include=="both"){
        func_coeff_sd=apply(curve_est,2,sd)
        func_coeff_extr=1:dim(curve_est)[1]
        func_coeff_est=apply(curve_est,2,mean)
        for(i in 1:dim(curve_est)[1]){
          func_coeff_extr[i]=max(abs(curve_est[i,]-func_coeff_est)/func_coeff_sd)
        }
        cutpoint=stats::quantile(func_coeff_extr,probs = prob)
        CMA.upper.curve.est = mean.curve.est + cutpoint*func_coeff_sd
        CMA.lower.curve.est = mean.curve.est - cutpoint*func_coeff_sd
        upper.curve.est = upper.curve.est.quantile
        lower.curve.est = lower.curve.est.quantile
      }

      plotdata2 <- data.frame( value = c(mean.curve.est),
                               CI.upper = c(upper.curve.est),
                               CI.lower = c(lower.curve.est),
                               CMA.upper = c(CMA.upper.curve.est),
                               CMA.lower = c(CMA.lower.curve.est),
                               xmat = c(1:length(upper.curve.est)),
                               Method = c(rep("Bayesian",length(mean.curve.est))))

      plot.res[[inx.effect]] <- ggplot2::ggplot(plotdata2, aes(y = .data$value, x = .data$xmat, color = .data$Method)) +
        ggplot2::geom_line(size = 1) +
        ggplot2::geom_ribbon(aes(ymin = .data$CI.lower, ymax = .data$CI.upper, fill = "Pointwise CI"),
                    alpha = 0.2, color = "black", linetype = "dashed") +
        ggplot2::geom_ribbon(aes(ymin = .data$CMA.lower, ymax = .data$CMA.upper, fill = "CMA CI"),
                    alpha = 0.2, linetype = "dotted", color = "black") +
        ggplot2::ylab("Functional Effect") + ggplot2::facet_wrap(~ .data$Method)+
        ggplot2::scale_fill_manual(values = c("Pointwise CI" = "darkgrey", "CMA CI" = "lightgrey"), name = "Interval Type") +
        ggplot2::theme_minimal()
    }
  }
  return(plot.res)
}


#----------------------------------------------------------------------------
#' Plot the estimated functional coefficients with the corresponding credible interval
#'
#' @import ggplot2
#' @import dplyr
#' @export plot.bfrs

plot.bfrs=function(brfs.fit=NULL,prob = 0.95,include = "both"){

  if(brfs.fit$family=="functional"){
    plot.res=list()
    for(inx.effect in 1:dim(brfs.fit$func_effect)[3]){
      curve_est=brfs.fit$func_effect[,,inx.effect]
      mean.curve.est=apply(curve_est,1,mean)
      upper.curve.est.quantile=apply(curve_est,1,function(x){quantile(x,probs = (1+prob)/2)})
      lower.curve.est.quantile=apply(curve_est,1,function(x){quantile(x,probs = (1-prob)/2)})
      upper.curve.wald=upper.curve.est.quantile
      lower.curve.wald=lower.curve.est.quantile

      for(i in 1:length(upper.curve.wald)){
        upper.curve.wald[i]=mean.curve.est[i]+qnorm((1+prob)/2)*sd(curve_est[,i])
        lower.curve.wald[i]=mean.curve.est[i]-qnorm((1+prob)/2)*sd(curve_est[,i])
      }
      if(type=="quantile"){
        upper.curve.est = upper.curve.est.quantile
        lower.curve.est = lower.curve.est.quantile
      }else if(type=="Wald"){
        upper.curve.est = upper.curve.wald
        lower.curve.est = lower.curve.wald
      }else{
        stop("Type should be either quantile or Wald!")
      }

      plotdata=data.frame(value=c(mean.curve.est,
                                  upper.curve.est,
                                  lower.curve.est),
                          xmat=c(1: length(upper.curve.est),
                                 1: length(upper.curve.est),
                                 1: length(upper.curve.est)),
                          type=c(rep("mean",length(upper.curve.est)),
                                 rep("upper",length(upper.curve.est)),
                                 rep("lower",length(upper.curve.est))))
      plot.res[[inx.effect]]=ggplot(plotdata,aes(y=value,x=xmat))+geom_line(aes(type=type,color=type))+
        ylab(colnames(brfs.fit$Standata[["X_s"]])[inx.effect])
    }
  }else{
    n.func.coeff=length(brfs.fit$func_effect)
    plot.res=list()
    for(inx.effect in 1:n.func.coeff){
      curve_est=brfs.fit$func_effect[[inx.effect]]
      mean.curve.est=apply(curve_est,2,mean)
      upper.curve.est.quantile=apply(curve_est,2,function(x){quantile(x,probs = (1+prob)/2)})
      lower.curve.est.quantile=apply(curve_est,2,function(x){quantile(x,probs = (1-prob)/2)})

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
        cutpoint=quantile(func_coeff_extr,probs = prob)
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
        cutpoint=quantile(func_coeff_extr,probs = prob)
        CMA.upper.curve.est = mean.curve.est + cutpoint*func_coeff_sd
        CMA.lower.curve.est = mean.curve.est - cutpoint*func_coeff_sd
        upper.curve.est = upper.curve.est.quantile
        lower.curve.est = lower.curve.est.quantile
      }

      plotdata2=data.frame(value=c(mean.curve.est),
                           CI.upper=c(upper.curve.est),
                           CI.lower=c(lower.curve.est),
                           CMA.upper = c(CMA.upper.curve.est),
                           CMA.lower = c(CMA.lower.curve.est),
                           xmat=c(1:length(upper.curve.est)),
                           Method=c(rep("Bayesian",dim(nhanes_lite_use$MIMS)[2])))

      plot.res[[inx.effect]]=ggplot(plotdata2, aes(y = value, x = xmat, color = Method)) +
        geom_line(size = 1) +
        geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper, fill = "Pointwise CI"),
                    alpha = 0.2, color = "black", linetype = "dashed") +
        geom_ribbon(aes(ymin = CMA.lower, ymax = CMA.upper, fill = "CMA CI"),
                    alpha = 0.2, linetype = "dotted", color = "black") +
        ylab("Functional Effect") +facet_wrap(~ Method)+
        scale_fill_manual(values = c("Pointwise CI" = "darkgrey", "CMA CI" = "lightgrey"), name = "Interval Type") +
        theme_minimal()
    }
  }
  return(plot.res)
}


## Fit a Bayesian SoFR using Stan and compare it with mgcv::gam

# library(rstan)
library(mgcv)
library(dplyr)
# library(MASS)
# library(splines)
# library(scam)
# library(refund)
library(refundBayes)

# ## mgcv internal function
# ExtractData <- function(object,data,knots) {
#   ## `data' and `knots' contain the data needed to evaluate the `terms', `by'
#   ## and `knots' elements of `object'. This routine does so, and returns
#   ## a list with element `data' containing just the evaluated `terms', 
#   ## with the by variable as the last column. If the `terms' evaluate matrices, 
#   ## then a check is made of whether repeat evaluations are being made, 
#   ## and if so only the unique evaluation points are returned in data, along 
#   ## with the `index' attribute required to re-assemble the full dataset.
#   knt <- dat <- list()
#   ## should data be processed as for summation convention with matrix arguments?
#   vecMat <- if (!is.list(object$xt)||is.null(object$xt$sumConv)) TRUE else object$xt$sumConv
#   for (i in 1:length(object$term)) { 
#     dat[[object$term[i]]] <- get.var(object$term[i],data,vecMat=vecMat)
#     knt[[object$term[i]]] <- get.var(object$term[i],knots,vecMat=vecMat)
#     
#   }
#   names(dat) <- object$term; m <- length(object$term)
#   if (!is.null(attr(dat[[1]],"matrix")) && vecMat) { ## strip down to unique covariate combinations
#     n <- length(dat[[1]])
#     #X <- matrix(unlist(dat),n,m) ## no use for factors!
#     X <- data.frame(dat)
#     #if (is.numeric(X)) {
#     X <- uniquecombs(X)
#     if (nrow(X)<n*.9) { ## worth the hassle
#       for (i in 1:m) dat[[i]] <- X[,i]     ## return only unique rows
#       attr(dat,"index") <- attr(X,"index") ## index[i] is row of dat[[i]] containing original row i
#     }
#     #} ## end if(is.numeric(X))
#   }    
#   if (object$by!="NA") {
#     by <- get.var(object$by,data) 
#     if (!is.null(by))
#     { dat[[m+1]] <- by 
#     names(dat)[m+1] <- object$by
#     }
#   }
#   return(list(data=dat,knots=knt))
# }

## generate data
n_num = 2000 ## number of subjects
T_num = 100 ## number of time points
K_num = 6 ## number of scalar predictors
basis_use = readRDS("./../../../../../Downloads/Supplementary materials/example_data_basis.rds")
set.seed(1234)
eigenfunc=basis_use[2:(K_num+1),]
scores = matrix(rnorm(n_num*K_num,0,1),nrow = n_num, ncol = K_num)
mu_func = basis_use[1,]*1
W_mat = matrix(mu_func,nrow = n_num,ncol = T_num,byrow = T) + scores %*% eigenfunc
true.beta.func=function(x){
  (-(-0.5+(x-0.5)^2)-0.416)*10
}
true.beta=true.beta.func((1:T_num)/T_num)
p_num = 1
data.X = rnorm(n_num,sd = 0.1)
true.X.effect = 0.1
linearPred.func = W_mat%*%matrix(true.beta,ncol=1)/length(true.beta)
linearPred.X = data.X%*%matrix(true.X.effect,ncol=1)
linearPred = linearPred.func+linearPred.X
logit.linearPred=exp(linearPred)/(1+exp(linearPred))
y.use = ifelse(runif(length(logit.linearPred))<logit.linearPred,1,0)
nt = ncol(W_mat)
tind = seq(0, 1, length.out = nt)
data.SoFR = data.frame(wmat = I(W_mat),
                       lmat = I(matrix(1/nt, ncol = nt, nrow = n_num)),
                       tmat = I(matrix(tind, ncol = nt, nrow = n_num, byrow = TRUE)))
data.SoFR$y=c(y.use)
data.SoFR$X1=data.X

# fit_freq2 = gam(y~ s(tmat, by=lmat*wmat, bs="cc", k=10)+X1, data=data.SoFR, family="binomial", method = "REML")
# plot(fit_freq2, unconditional = TRUE)

## center the functional predictor
data.SoFR$wmat_c <- I(scale(data.SoFR$wmat, center = TRUE, scale = FALSE))

fit_freq3 = gam(y~ s(tmat, by=lmat*wmat_c, bs="cc", k=10), data=data.SoFR, family="binomial", method = "REML") ## the intercept is more interpretable
summary(fit_freq3)
plot(fit_freq3, unconditional = TRUE)


## Bayesian SoFR
fit_bayes <- sofr_bayes(formula = y ~ s(tmat, by=lmat*wmat_c, bs="cc", k=10),
                        data = data.SoFR,
                        family = "binomial")

gamma <- apply(fit_bayes$scalar_coef, 2, mean)
Intercept <- mean(fit_bayes$int)

beta.post <- fit_bayes$func_coef[[1]]
mean.curve.est=apply(beta.post,2,mean)
upper.curve.est=apply(beta.post,2,function(x){quantile(x,probs = 0.975)})
lower.curve.est=apply(beta.post,2,function(x){quantile(x,probs = 0.025)})  

plotdata=data.frame(value=c(mean.curve.est,
                            upper.curve.est,
                            lower.curve.est),
                    xmat=c(rep(1: dim(data.SoFR$wmat)[2],3)),
                    Method=c(rep("Bayesian",dim(data.SoFR$wmat)[2]),
                             rep("Bayesian",dim(data.SoFR$wmat)[2]),
                             rep("Bayesian",dim(data.SoFR$wmat)[2])),
                    type=c(rep("Estimate",dim(data.SoFR$wmat)[2]),
                           rep("CI_upper",dim(data.SoFR$wmat)[2]),
                           rep("CI_lower",dim(data.SoFR$wmat)[2])))

library(ggplot2)
ggplot(plotdata,aes(y=value,x=xmat))+geom_line(aes(linetype=type,color=Method))+
  ylab("Functional effect")+xlab("Time (hour)")+ 
  #scale_x_continuous(breaks=seq(0,1,by=3))+
  scale_linetype_manual(values=c("twodash", "longdash","solid"),name="Line Type")+
  theme_minimal()



#----------------------------------------------------------------------------
#'
#' This document describes ancillary functions for refundBayes package that is used for data process.
#'
#'


#----------------------------------------------------------------------------
#' Construct the spline basis according to the spline type and functional data
#' For now, only used for FoSR.
#----------------------------------------------------------------------------

brfs_extract_basis=function(formula,data,func_comp){
  sp_basis=list()
  for(i in 1:length(func_comp)){
    sp_type=formula$func_var[[i]][[4]]

    if(sp_type=="cr"){
      sp_basis[[i]]=brfs_extract_basis_cr(formula$func_var[[i]],data)
    }

    if(sp_type=="cc"){
      sp_basis[[i]]=brfs_extract_basis_cc(formula$func_var[[i]],data)
    }
  }

  return(sp_basis)
}


#----------------------------------------------------------------------------
#' Extract basis for cubic regression spline
#----------------------------------------------------------------------------


brfs_extract_basis_cr=function(term,data){

  bymat=data[[term[["by"]][[3]]]]
  curr_term=term[["by"]]
  while(class(curr_term[[2]])!="name"){

    curr_term=curr_term[[2]]
    bymat=bymat*data[[curr_term[[3]]]]
    print(curr_term[[3]])
  }
  bymat=bymat*data[[curr_term[[2]]]]


  eval(parse(text = paste0("object=",format(term))))

  #object=s(tmat, by=lmat*MIMS, bs="cc", k=10)
  #data = data_sim
  knots = NULL

  dk <- ExtractData(object,data,knots)


  crspline.cons <- smooth.construct.cr.smooth.spec(object,dk$data,dk$knots)



  beta.base=crspline.cons$X


  ### for cr:
  beta.base.new=matrix(1,nrow = prod(dim(bymat)), ncol=NCOL(beta.base))

  for(i in 1:NROW(beta.base)){
    for(inx in ((i-1)*dim(bymat)[1]+1):(i*dim(bymat)[1])){
      beta.base.new[inx,]=beta.base[i,]
    }
  }

  beta.base=beta.base.new
  beta.S=crspline.cons$S



  modCon=0
  sm <-crspline.cons
  knots=NULL
  absorb.cons=TRUE
  scale.penalty=TRUE
  n=nrow(data)
  q=dim(bymat)[2]
  dataX = NULL
  null.space.penalty = FALSE
  sparse.cons=0
  diagonal.penalty=TRUE
  apply.by=TRUE



  ############# Scale the S matrix


  beta.S.scale <- rep(1,length(beta.S))
  maXX <- norm(beta.base,type="I")^2 ##mean(abs(t(sm$X)%*%sm$X)) # `size' of X'X

  maS <- norm(beta.S[[1]])/maXX  ## mean(abs(sm$S[[i]])) / maXX
  beta.S[[1]] <- beta.S[[1]] / maS
  beta.S.scale[1] <- maS ## multiply S[[i]] by this to get original S[[i]]



  eigendecomp=eigen(beta.S[[1]], symmetric = TRUE)




  ############# take product of the by matrix


  smlX <- as.numeric(bymat)*beta.base


  ind <- 1:n
  Xmat <- smlX[ind,,drop=FALSE]
  for (i in 2:q) {
    ind <- ind + n
    Xmat <- Xmat +smlX[ind,,drop=FALSE]
  }
  smlX <- Xmat




  ############# some transportation




  rank=crspline.cons$rank

  E <- rep(1, ncol(Xmat))
  E[1:rank] <- sqrt(eigendecomp$value[1:rank])
  Xmat <- Xmat %*% eigendecomp$vectors
  col.norm <- colSums(Xmat^2)
  col.norm <- col.norm/E^2
  av.norm <- mean(col.norm[1:rank])
  for (i in (rank + 1):ncol(Xmat)) {
    E[i] <- sqrt(col.norm[i]/av.norm)
  }

  Xmat <- t(t(Xmat)/E)



  base_mat=matrix(nrow = dim(bymat)[2],ncol = 10)

  for(i in 1:nrow(base_mat)){
    base_mat[i,]=beta.base[dim(bymat)[1]*(i-1)+1,]
  }

  return(list(E=E,
              eigendecomp=eigendecomp,
              base_mat=base_mat))
}



#----------------------------------------------------------------------------
#' Extract basis for cyclic cubic regression spline
#----------------------------------------------------------------------------


brfs_extract_basis_cc=function(term,data){

  bymat=data[[term[["by"]][[3]]]]
  curr_term=term[["by"]]
  while(class(curr_term[[2]])!="name"){

    curr_term=curr_term[[2]]
    bymat=bymat*data[[curr_term[[3]]]]
    print(curr_term[[3]])
  }
  bymat=bymat*data[[curr_term[[2]]]]


  eval(parse(text = paste0("object=",format(term))))

  #object=s(tmat, by=lmat*MIMS, bs="cc", k=10)
  #data = data_sim
  knots = NULL

  dk <- ExtractData(object,data,knots)


  crspline.cons <- smooth.construct.cc.smooth.spec(object,dk$data,dk$knots)



  beta.base=crspline.cons$X


  ### for cc:
  beta.base.new=matrix(1,nrow = prod(dim(bymat)), ncol=NCOL(beta.base))

  for(i in 1:NROW(beta.base)){
    for(inx in ((i-1)*dim(bymat)[1]+1):(i*dim(bymat)[1])){
      beta.base.new[inx,]=beta.base[i,]
    }
  }

  beta.base=beta.base.new
  beta.S=crspline.cons$S



  modCon=0
  sm <-crspline.cons
  knots=NULL
  absorb.cons=TRUE
  scale.penalty=TRUE
  n=nrow(data)
  q=dim(bymat)[2]
  dataX = NULL
  null.space.penalty = FALSE
  sparse.cons=0
  diagonal.penalty=TRUE
  apply.by=TRUE



  ############# Scale the S matrix


  beta.S.scale <- rep(1,length(beta.S))
  maXX <- norm(beta.base,type="I")^2 ##mean(abs(t(sm$X)%*%sm$X)) # `size' of X'X

  maS <- norm(beta.S[[1]])/maXX  ## mean(abs(sm$S[[i]])) / maXX
  beta.S[[1]] <- beta.S[[1]] / maS
  beta.S.scale[1] <- maS ## multiply S[[i]] by this to get original S[[i]]



  eigendecomp=eigen(beta.S[[1]], symmetric = TRUE)




  ############# take product of the by matrix


  smlX <- as.numeric(bymat)*beta.base


  ind <- 1:n
  Xmat <- smlX[ind,,drop=FALSE]
  for (i in 2:q) {
    ind <- ind + n
    Xmat <- Xmat +smlX[ind,,drop=FALSE]
  }
  smlX <- Xmat




  ############# some transportation




  rank=crspline.cons$rank

  E <- rep(1, ncol(Xmat))
  E[1:rank] <- sqrt(eigendecomp$value[1:rank])
  Xmat <- Xmat %*% eigendecomp$vectors
  col.norm <- colSums(Xmat^2)
  col.norm <- col.norm/E^2
  av.norm <- mean(col.norm[1:rank])
  for (i in (rank + 1):ncol(Xmat)) {
    E[i] <- sqrt(col.norm[i]/av.norm)
  }

  Xmat <- t(t(Xmat)/E)



  base_mat=matrix(nrow = dim(bymat)[2],ncol = NCOL(beta.base))

  for(i in 1:nrow(base_mat)){
    base_mat[i,]=beta.base[dim(bymat)[1]*(i-1)+1,]
  }

  return(list(E=E,
              eigendecomp=eigendecomp,
              base_mat=base_mat))
}



#----------------------------------------------------------------------------
#' Construct the SoFR fitted functional coefficient according to Stan posterior sampling
#----------------------------------------------------------------------------



brfs_effect=function(basis,fit.samp,func_comp,trans.mat){
  sp_samp=list()
  for(inx in 1:length(func_comp)){
    eigendecomp=basis[[inx]][["eigendecomp"]]
    base_mat=basis[[inx]][["base_mat"]]
    E=basis[[inx]][["E"]]

    para=cbind(fit.samp[[paste0("br_",inx)]],fit.samp[[paste0("bf_",inx)]])
    para.trans=para %*% trans.mat[[inx]]
    para.trans.untilde=para.trans

    for(i in 1:nrow(para.trans)){
      para.trans.untilde[i,]=(eigendecomp$vectors %*% diag(1/E)) %*% matrix(para.trans[i,],ncol=1)
    }

    curve_est=matrix(ncol = nrow(base_mat),nrow = nrow(para.trans.untilde))


    for(i in 1:nrow(curve_est)){
      curve_est[i,]=base_mat%*%para.trans.untilde[i,]
    }

    sp_samp[[inx]]=curve_est
  }

  return(sp_samp)
}

#----------------------------------------------------------------------------
#' Construct the FoSR fitted functional coefficient according to Stan posterior sampling
#----------------------------------------------------------------------------


bfrs_effect_functional=function(basis,b.samp){
  sp_samp=list()
  for(inx in 1:length(func_comp)){
    eigendecomp=basis[[inx]][["eigendecomp"]]
    base_mat=basis[[inx]][["base_mat"]]
    E=basis[[inx]][["E"]]

    para=cbind(fit.samp[[paste0("br_",inx)]],fit.samp[[paste0("bf_",inx)]])
    para.trans=para %*% trans.mat[[inx]]
    para.trans.untilde=para.trans

    for(i in 1:nrow(para.trans)){
      para.trans.untilde[i,]=(eigendecomp$vectors %*% diag(1/E)) %*% matrix(para.trans[i,],ncol=1)
    }

    curve_est=matrix(ncol = nrow(base_mat),nrow = nrow(para.trans.untilde))


    for(i in 1:nrow(curve_est)){
      curve_est[i,]=base_mat%*%para.trans.untilde[i,]
    }

    sp_samp[[inx]]=curve_est
  }

  return(sp_samp)
}


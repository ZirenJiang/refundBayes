#' Functions to extract the spline basis from data based on formula
#'

#----------------------------------------------------------------------------
#' The main function for extracting the spline basis from data based on formula
#----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
#' 
extract_basis = function(formula, data, func_comp){ ## currently only support "cc" and "cr"
  sp_basis <- list()
  for(i in 1:length(func_comp)){
    sp_type <- formula$func_var[[i]][[4]]
    # if(sp_type == "cr"){
    #   sp_basis[[i]] <- extract_basis_cr(formula$func_var[[i]], data)
    # }
    # if(sp_type == "cc"){
    #   sp_basis[[i]] <- extract_basis_cc(formula$func_var[[i]], data)
    # }
    if(isTRUE(func_comp[i])){
      sp_basis[[i]] <- extract_basis_fpca(formula$func_var[[i]], data)
    }else{
      sp_basis[[i]] <- extract_basis_all(formula$func_var[[i]], data)
    }
  }
  return(sp_basis)
}

#----------------------------------------------------------------------------
#' Extract the basis info needed to reconstruct the functional coefficient
#' \eqn{\beta(s)} when the functional predictor is jointly modelled with
#' FPCA (i.e. when \code{func_comp[i] == TRUE}).
#'
#' Mirrors the logic of \code{ext_func_pred_fpca()} (in refundBayesdata.R) so
#' that the same \code{Psi_mat}, eigendecomposition, and \code{E} scaling are
#' used for posterior reconstruction. The Stan parameters \code{br_i} and
#' \code{bf_i} are already in the (eigendecomp, E)-transformed space, so the
#' reconstruction is:
#'   1. spline coef = (eigendecomp$vectors %*% diag(1 / E)) %*% c(br_i, bf_i)
#'   2. beta(s)     = base_mat %*% spline coef
#'
#' @keywords internal
#' @noRd
#'
extract_basis_fpca = function(term, data){
  knots <- NULL
  obj <- 1
  eval(parse(text = paste0("obj = ", format(term))))
  dk <- ExtractData(obj, data, knots)
  splinecons <- mgcv::smooth.construct(obj, dk$data, dk$knots)
  Psi_mat <- splinecons$X
  S_mat <- splinecons$S[[1]]
  rank <- splinecons$rank
  K_num <- ncol(Psi_mat)

  ## Scale the penalty matrix to match ext_func_pred_fpca()
  maXX <- norm(Psi_mat, type = "I") ^ 2
  maS <- norm(S_mat) / maXX
  S_mat_scaled <- S_mat / maS
  eigendecomp <- eigen(S_mat_scaled, symmetric = TRUE)

  ## Recompute X_mat_t to derive the same E that was used at data-prep time.
  parsed <- parse_by_func(term, data)
  funcmat <- parsed$funcmat
  lmat <- parsed$lmat
  fpca_fit <- refund::fpca.sc(Y = unclass(funcmat))
  Phi <- fpca_fit$efunctions
  if(is.null(dim(Phi))){
    Phi <- matrix(Phi, ncol = 1)
  }
  l_vec <- as.numeric(lmat[1, ])
  X_mat_t <- t(Phi) %*% (l_vec * Psi_mat)

  E <- rep(1, K_num)
  E[1:rank] <- sqrt(eigendecomp$value[1:rank])
  X_mat_t <- X_mat_t %*% eigendecomp$vectors
  if(rank < K_num){
    col.norm <- colSums(X_mat_t ^ 2)
    col.norm <- col.norm / E ^ 2
    av.norm <- mean(col.norm[1:rank])
    for (kk in (rank + 1):K_num) {
      E[kk] <- sqrt(col.norm[kk] / av.norm)
    }
  }

  return(list(E = E,
              eigendecomp = eigendecomp,
              base_mat = Psi_mat,
              joint_FPCA = TRUE))
}

#----------------------------------------------------------------------------
#' Extract basis for all types of regression spline
#----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
#' 
extract_basis_all = function(term, data){
  if (is.name(term[["by"]])) {
    bymat <- data[[as.character(term[["by"]])]]
  } else {
    bymat <- data[[term[["by"]][[3]]]]
    curr_term <- term[["by"]]
    while(!is.name(curr_term[[2]])) {
      curr_term <- curr_term[[2]]
      bymat <- bymat * data[[curr_term[[3]]]]
    }
    bymat <- bymat * data[[curr_term[[2]]]]
  }
  object = 1
  eval(parse(text = paste0("object = ", format(term))))
  knots <- NULL
  dk <- ExtractData(object, data, knots) ## organize data as input for smooth.construct-type functions
  crspline.cons <- mgcv::smooth.construct(object, dk$data, dk$knots)
  beta.base <- crspline.cons$X
  
  ## for cr:
  beta.base.new <- matrix(1, nrow = prod(dim(bymat)), ncol = NCOL(beta.base))
  for(i in 1:NROW(beta.base)){
    for(inx in ((i-1)*dim(bymat)[1]+1):(i*dim(bymat)[1])){
      beta.base.new[inx,] <- beta.base[i,]
    }
  }
  beta.base <- beta.base.new
  beta.S <- crspline.cons$S
  modCon <- 0
  sm <- crspline.cons
  knots <- NULL
  absorb.cons <- TRUE
  scale.penalty <- TRUE
  n <- nrow(data)
  q <- dim(bymat)[2]
  dataX <- NULL
  null.space.penalty <- FALSE
  sparse.cons <- 0
  diagonal.penalty <- TRUE
  apply.by <- TRUE
  
  ### Scale the S matrix
  beta.S.scale <- rep(1,length(beta.S))
  maXX <- norm(beta.base, type = "I") ^ 2 ##mean(abs(t(sm$X)%*%sm$X)) # `size' of X'X
  maS <- norm(beta.S[[1]]) / maXX  ## mean(abs(sm$S[[i]])) / maXX
  beta.S[[1]] <- beta.S[[1]] / maS
  beta.S.scale[1] <- maS ## multiply S[[i]] by this to get original S[[i]]
  eigendecomp <- eigen(beta.S[[1]], symmetric = TRUE)
  
  ### take product of the by matrix
  smlX <- as.numeric(bymat) * beta.base
  ind <- 1:n
  Xmat <- smlX[ind,,drop=FALSE]
  for (i in 2:q) {
    ind <- ind + n
    Xmat <- Xmat + smlX[ind,,drop=FALSE]
  }
  smlX <- Xmat
  
  ### some transformation
  rank <- crspline.cons$rank
  E <- rep(1, ncol(Xmat))
  E[1:rank] <- sqrt(eigendecomp$value[1:rank])
  Xmat <- Xmat %*% eigendecomp$vectors
  col.norm <- colSums(Xmat ^ 2)
  col.norm <- col.norm / E ^ 2
  av.norm <- mean(col.norm[1:rank])
  for (i in (rank + 1):ncol(Xmat)) {
    E[i] <- sqrt(col.norm[i] / av.norm)
  }
  Xmat <- t(t(Xmat) / E)
  base_mat <- matrix(nrow = dim(bymat)[2], ncol = NCOL(beta.base))
  for(i in 1:nrow(base_mat)){
    base_mat[i,] <- beta.base[dim(bymat)[1]*(i-1)+1,]
  }
  return(list(E = E,
              eigendecomp = eigendecomp,
              base_mat = base_mat))
}



#----------------------------------------------------------------------------
#' Extract basis for cubic regression spline ("cr")
#----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
#' 
extract_basis_cr = function(term, data){
  if (is.name(term[["by"]])) {
    bymat <- data[[as.character(term[["by"]])]]
  } else {
    bymat <- data[[term[["by"]][[3]]]]
    curr_term <- term[["by"]]
    while(!is.name(curr_term[[2]])) {
      curr_term <- curr_term[[2]]
      bymat <- bymat * data[[curr_term[[3]]]]
    }
    bymat <- bymat * data[[curr_term[[2]]]]
  }
  object = 1
  eval(parse(text = paste0("object = ", format(term))))
  knots <- NULL
  dk <- ExtractData(object, data, knots) ## organize data as input for smooth.construct-type functions
  crspline.cons <- mgcv::smooth.construct.cr.smooth.spec(object, dk$data, dk$knots)
  beta.base <- crspline.cons$X
  
  ## for cr:
  beta.base.new <- matrix(1, nrow = prod(dim(bymat)), ncol = NCOL(beta.base))
  for(i in 1:NROW(beta.base)){
    for(inx in ((i-1)*dim(bymat)[1]+1):(i*dim(bymat)[1])){
      beta.base.new[inx,] <- beta.base[i,]
    }
  }
  beta.base <- beta.base.new
  beta.S <- crspline.cons$S
  modCon <- 0
  sm <- crspline.cons
  knots <- NULL
  absorb.cons <- TRUE
  scale.penalty <- TRUE
  n <- nrow(data)
  q <- dim(bymat)[2]
  dataX <- NULL
  null.space.penalty <- FALSE
  sparse.cons <- 0
  diagonal.penalty <- TRUE
  apply.by <- TRUE
  
  ### Scale the S matrix
  beta.S.scale <- rep(1,length(beta.S))
  maXX <- norm(beta.base, type = "I") ^ 2 ##mean(abs(t(sm$X)%*%sm$X)) # `size' of X'X
  maS <- norm(beta.S[[1]]) / maXX  ## mean(abs(sm$S[[i]])) / maXX
  beta.S[[1]] <- beta.S[[1]] / maS
  beta.S.scale[1] <- maS ## multiply S[[i]] by this to get original S[[i]]
  eigendecomp <- eigen(beta.S[[1]], symmetric = TRUE)
  
  ### take product of the by matrix
  smlX <- as.numeric(bymat) * beta.base
  ind <- 1:n
  Xmat <- smlX[ind,,drop=FALSE]
  for (i in 2:q) {
    ind <- ind + n
    Xmat <- Xmat + smlX[ind,,drop=FALSE]
  }
  smlX <- Xmat
  
  ### some transformation
  rank <- crspline.cons$rank
  E <- rep(1, ncol(Xmat))
  E[1:rank] <- sqrt(eigendecomp$value[1:rank])
  Xmat <- Xmat %*% eigendecomp$vectors
  col.norm <- colSums(Xmat ^ 2)
  col.norm <- col.norm / E ^ 2
  av.norm <- mean(col.norm[1:rank])
  for (i in (rank + 1):ncol(Xmat)) {
    E[i] <- sqrt(col.norm[i] / av.norm)
  }
  Xmat <- t(t(Xmat) / E)
  base_mat <- matrix(nrow = dim(bymat)[2], ncol = 10)
  for(i in 1:nrow(base_mat)){
    base_mat[i,] <- beta.base[dim(bymat)[1]*(i-1)+1,]
  }
  return(list(E = E,
              eigendecomp = eigendecomp,
              base_mat = base_mat))
}

#----------------------------------------------------------------------------
#' Extract basis for cyclic cubic regression spline ("cc")
#----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
#' 
extract_basis_cc = function(term, data){
  if (is.name(term[["by"]])) {
    bymat <- data[[as.character(term[["by"]])]]
  } else {
    bymat <- data[[term[["by"]][[3]]]]
    curr_term <- term[["by"]]
    while(!is.name(curr_term[[2]])) {
      curr_term <- curr_term[[2]]
      bymat <- bymat * data[[curr_term[[3]]]]
    }
    bymat <- bymat * data[[curr_term[[2]]]]
  }
  object = 1
  eval(parse(text = paste0("object = ", format(term))))
  knots <- NULL
  dk <- ExtractData(object,data,knots)
  crspline.cons <- smooth.construct.cc.smooth.spec(object, dk$data, dk$knots)
  beta.base <- crspline.cons$X
  
  ### for cc:
  beta.base.new <- matrix(1, nrow = prod(dim(bymat)), ncol=NCOL(beta.base))
  for(i in 1:NROW(beta.base)){
    for(inx in ((i-1)*dim(bymat)[1]+1):(i*dim(bymat)[1])){
      beta.base.new[inx,] <- beta.base[i,]
    }
  }
  beta.base <- beta.base.new
  beta.S <- crspline.cons$S
  modCon <- 0
  sm <- crspline.cons
  knots <- NULL
  absorb.cons <- TRUE
  scale.penalty <- TRUE
  n <- nrow(data)
  q <- dim(bymat)[2]
  dataX <- NULL
  null.space.penalty <- FALSE
  sparse.cons <- 0
  diagonal.penalty <- TRUE
  apply.by <- TRUE
  
  ### Scale the S matrix
  beta.S.scale <- rep(1,length(beta.S))
  maXX <- norm(beta.base, type = "I") ^ 2 ##mean(abs(t(sm$X)%*%sm$X)) # `size' of X'X
  maS <- norm(beta.S[[1]]) / maXX  ## mean(abs(sm$S[[i]])) / maXX
  beta.S[[1]] <- beta.S[[1]] / maS
  beta.S.scale[1] <- maS ## multiply S[[i]] by this to get original S[[i]]
  eigendecomp <- eigen(beta.S[[1]], symmetric = TRUE)
  
  ### take product of the by matrix
  smlX <- as.numeric(bymat) * beta.base
  ind <- 1:n
  Xmat <- smlX[ind,,drop=FALSE]
  for (i in 2:q) {
    ind <- ind + n
    Xmat <- Xmat +smlX[ind,,drop=FALSE]
  }
  smlX <- Xmat
  
  ### some transformation
  rank <- crspline.cons$rank
  E <- rep(1, ncol(Xmat))
  E[1:rank] <- sqrt(eigendecomp$value[1:rank])
  Xmat <- Xmat %*% eigendecomp$vectors
  col.norm <- colSums(Xmat^2)
  col.norm <- col.norm/E^2
  av.norm <- mean(col.norm[1:rank])
  for (i in (rank + 1):ncol(Xmat)) {
    E[i] <- sqrt(col.norm[i] / av.norm)
  }
  Xmat <- t(t(Xmat) / E)
  base_mat <- matrix(nrow = dim(bymat)[2], ncol = NCOL(beta.base))
  for(i in 1:nrow(base_mat)){
    base_mat[i,] <- beta.base[dim(bymat)[1]*(i-1)+1,]
  }
  return(list(E = E,
              eigendecomp = eigendecomp,
              base_mat = base_mat))
}

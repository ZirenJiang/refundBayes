
ExtractData <- function(object,data,knots) {
  ## `data' and `knots' contain the data needed to evaluate the `terms', `by'
  ## and `knots' elements of `object'. This routine does so, and returns
  ## a list with element `data' containing just the evaluated `terms',
  ## with the by variable as the last column. If the `terms' evaluate matrices,
  ## then a check is made of whether repeat evaluations are being made,
  ## and if so only the unique evaluation points are returned in data, along
  ## with the `index' attribute required to re-assemble the full dataset.
  knt <- dat <- list()
  ## should data be processed as for summation convention with matrix arguments?
  vecMat <- if (!is.list(object$xt)||is.null(object$xt$sumConv)) TRUE else object$xt$sumConv
  for (i in 1:length(object$term)) {
    dat[[object$term[i]]] <- get.var(object$term[i],data,vecMat=vecMat)
    knt[[object$term[i]]] <- get.var(object$term[i],knots,vecMat=vecMat)

  }
  names(dat) <- object$term; m <- length(object$term)
  if (!is.null(attr(dat[[1]],"matrix")) && vecMat) { ## strip down to unique covariate combinations
    n <- length(dat[[1]])
    #X <- matrix(unlist(dat),n,m) ## no use for factors!
    X <- data.frame(dat)
    #if (is.numeric(X)) {
    X <- uniquecombs(X)
    if (nrow(X)<n*.9) { ## worth the hassle
      for (i in 1:m) dat[[i]] <- X[,i]     ## return only unique rows
      attr(dat,"index") <- attr(X,"index") ## index[i] is row of dat[[i]] containing original row i
    }
    #} ## end if(is.numeric(X))
  }
  if (object$by!="NA") {
    by <- get.var(object$by,data)
    if (!is.null(by))
    { dat[[m+1]] <- by
    names(dat)[m+1] <- object$by
    }
  }
  return(list(data=dat,knots=knt))
}




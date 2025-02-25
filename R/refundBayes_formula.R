#' Functions for building the data that is used for Stan.
#'

#----------------------------------------------------------------------------
#' General functional for recongnizing the formula input of refundBayes
#----------------------------------------------------------------------------

brfs_formula=function(formula){
  formula_val=brms::brmsformula(formula)
  y_var=formula_val$formula[[2]]
  if(class(formula_val$formula[[3]])=="name"){
    return(list(scalar_var=list(formula_val$formula[[3]]),
                func_var=list(),
                y_var=y_var))
  }else{
    var_extract=extract_var(formula_val$formula[[3]])
    return(list(scalar_var=var_extract$scalar_var,
                func_var=var_extract$func_var,
                y_var=y_var))
  }
}

#----------------------------------------------------------------------------
#' Whether the object is a functional object with smoothness
#----------------------------------------------------------------------------

is.func.cov=function(term){
  func=FALSE
  if(term[[1]]=="s"){
    func=TRUE
  }else{
    if(class(term[[2]])!="name"){
      func1=is.func.cov(term[[2]])
    }else{
      func1=FALSE
    }
    if(class(term[[3]])!="name"){
      func2=is.func.cov(term[[3]])
    }else{
      func2=FALSE
    }
    func=func1|func2
  }
  return(func)
}

#----------------------------------------------------------------------------
#' Identify the scalar predictors and functional predictors
#----------------------------------------------------------------------------

extract_var=function(term){
  scalar_var=list()
  func_var=list()

  if(term[[1]]=="+"){
    if(class(term[[2]])!="name"){
      scalar_var=extract_var(term[[2]])$scalar_var
      func_var=extract_var(term[[2]])$func_var
    }else{
      scalar_var=term[[2]]
    }
    if(class(term[[3]])!="name"){
      scalar_var=c(scalar_var,extract_var(term[[3]])$scalar_var)
      func_var=c(func_var,extract_var(term[[3]])$func_var)
    }else{
      scalar_var=c(scalar_var,term[[3]])
    }
  }

  if(term[[1]]=="*"){
    if(is.func.cov(term)){
      if(class(term[[2]])!="name"){
        func_var=c(func_var,unlist(extract_var(term[[2]])$func_var))
      }else{
        func_var=c(func_var,unlist(term[[2]]))
      }
      if(class(term[[3]])!="name"){
        func_var=c(func_var,unlist(extract_var(term[[3]])$func_var))
      }else{
        func_var=c(func_var,unlist(term[[3]]))
      }
      func_var=unlist(func_var)
      func_var=list(func_var)
    }else{
      if(class(term[[2]])!="name"){
        scalar_var=c(scalar_var,unlist(extract_var(term[[2]])$scalar_var))
      }else{
        scalar_var=c(scalar_var,unlist(term[[2]]))
      }
      if(class(term[[3]])!="name"){
        scalar_var=c(scalar_var,unlist(extract_var(term[[3]])$scalar_var))
      }else{
        scalar_var=c(scalar_var,unlist(term[[3]]))
      }
      scalar_var=unlist(scalar_var)
      scalar_var=list(scalar_var)
    }
  }

  if(term[[1]]=="s"){
    func_var[[1]]=term
  }
  if(term[[1]]=="$"){
    scalar_var=term[[3]]
  }
  return(list(scalar_var=scalar_var,
              func_var=func_var))
}





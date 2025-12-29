#' Functions to generate the formula for subsequent Stan modeling.
#'

#----------------------------------------------------------------------------
#' The main function for generating the formula as input for other functions
#----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
#' 
refundBayesformula = function(formula){
  formula_val <- brms::brmsformula(formula)
  y_var <- formula_val$formula[[2]]
  #if(class(formula_val$formula[[3]]) == "name"){ ## when there is no functional variable
    if(is.name(formula_val$formula[[3]])){
      
    return(list(scalar_var = list(formula_val$formula[[3]]),
                func_var = list(),
                y_var = y_var))
  }else{ ## when there are functional variables
    var_extract <- extract_var(formula_val$formula[[3]])
    return(list(scalar_var = var_extract$scalar_var,
                func_var = var_extract$func_var,
                y_var = y_var))
  }
}

#----------------------------------------------------------------------------
#' A function for deciding whether the object is a functional object with smoothness
#----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
#' 
is.func.cov = function(term){
  func <- FALSE
  if(term[[1]] == "s"){ #!!! in mgcv, the functional term is specified using s(), but s() is also used to specify additive terms
    func <- TRUE
  }else{
    #if(class(term[[2]]) != "name"){
      if(!is.name(term[[2]])){
      func1 <- is.func.cov(term[[2]])
    }else{
      func1 <- FALSE
    }
    #if(class(term[[3]]) != "name"){
      if(!is.name(term[[3]])){
      func2 <- is.func.cov(term[[3]])
    }else{
      func2 <- FALSE
    }
    func <- func1 | func2
  }
  return(func)
}

#----------------------------------------------------------------------------
#' A function for identifying the scalar and functional predictors
#----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
#' 
extract_var = function(term){
  scalar_var <- list()
  func_var <- list()
  
  if(term[[1]] == "+"){
    #if(class(term[[2]]) != "name"){
      if(!is.name(term[[2]])){
      scalar_var <- extract_var(term[[2]])$scalar_var
      func_var <- extract_var(term[[2]])$func_var
    }else{
      scalar_var <- term[[2]]
    }
    #if(class(term[[3]]) != "name"){
      if(!is.name(term[[3]])){
      scalar_var <- c(scalar_var, extract_var(term[[3]])$scalar_var)
      func_var <- c(func_var, extract_var(term[[3]])$func_var)
    }else{
      scalar_var <- c(scalar_var, term[[3]])
    }
  }
  
  if(term[[1]] == "*"){
    if(is.func.cov(term)){
      #if(class(term[[2]]) != "name"){
        if(!is.name(term[[2]])){
        func_var <- c(func_var, unlist(extract_var(term[[2]])$func_var))
      }else{
        func_var <- c(func_var, unlist(term[[2]]))
      }
      #if(class(term[[3]]) != "name"){
        if(!is.name(term[[3]])){
        func_var <- c(func_var, unlist(extract_var(term[[3]])$func_var))
      }else{
        func_var <- c(func_var, unlist(term[[3]]))
      }
      func_var <- unlist(func_var)
      func_var <- list(func_var)
    }else{
      #if(class(term[[2]]) != "name"){
        if(!is.name(term[[2]])){
        scalar_var <- c(scalar_var, unlist(extract_var(term[[2]])$scalar_var))
      }else{
        scalar_var <- c(scalar_var, unlist(term[[2]]))
      }
      #if(class(term[[3]]) != "name"){
        if(!is.name(term[[3]])){
        scalar_var <- c(scalar_var, unlist(extract_var(term[[3]])$scalar_var))
      }else{
        scalar_var <- c(scalar_var, unlist(term[[3]]))
      }
      scalar_var <- unlist(scalar_var)
      scalar_var <- list(scalar_var)
    }
  }
  
  if(term[[1]] == "s"){
    func_var[[1]] <- term
  }
  if(term[[1]] == "$"){
    scalar_var <- term[[3]]
  }
  return(list(scalar_var = scalar_var,
              func_var = func_var))
}





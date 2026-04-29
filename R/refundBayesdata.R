#' Functions to build the data and code for subsequent Stan modeling
#'

#----------------------------------------------------------------------------
#' The main function for constructing the data and code used as input for Stan
#----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
refundBayesdata = function(formula, data, family, 
                           func_comp, intercept, 
                           cens = NULL, func_parameter = NULL){
  
  trans.mat <- list() ## to store the transformation matrix U
  result_code_data <- paste0("data{ \n")
  result_code_para <- paste0("parameters{ \n")
  result_code_transpara <- paste0("transformed parameters { \n real lprior = 0;\n")
  result <- list() ## a list storing all data input for Stan
  
  X_scalar <- ext_scalar_pred(formula$scalar_var, data) ## organize scalar predictors as input for Stan
  
  if(intercept){ ## if include an intercept
    result_code_para <- paste0(result_code_para, "   //Linear predictor intercept  \n")
    result_code_para <- paste0(result_code_para, "   real eta_0;\n")
  }
  result[["N_num"]] <- NROW(data)
  result_code_data <- paste0(result_code_data, "   //Total number of observations  \n")
  result_code_data <- paste0(result_code_data, "   int<lower=1> N_num;\n")
  
  ## for functional Cox model, need to model hazard functions
  if(family == "Cox"){
    args <- list()
    y <- data[[formula$y_var]]
    min_y <- min(y, na.rm = TRUE)
    max_y <- max(y, na.rm = TRUE)
    diff_y <- max_y - min_y
    lower_knot <- max(min_y - diff_y / 50, 0)
    upper_knot <- max_y + diff_y / 50
    args$Boundary.knots <- c(lower_knot, upper_knot)
    args$x <- y
    
    Mbasis <- splines2::mSpline(x = args$x, Boundary.knots = args$Boundary.knots, df=5L, intercept=TRUE)
    Ibasis <- splines2::iSpline(x = args$x, Boundary.knots = args$Boundary.knots, df=5L)
    
    result[["real_inter"]] <- 0
    result_code_data <- paste0(result_code_data, "   real real_inter;\n")
    
    result[["cens"]] <- cens
    result_code_data <- paste0(result_code_data, "   //Indicator variable for censoring  \n")
    result_code_data <- paste0(result_code_data, "   array[N_num] int<lower=-1,upper=2> cens;\n")
    
    result[["L_num"]] <- NCOL(Mbasis)
    result_code_data <- paste0(result_code_data, "   //Number of baseline hazard function basis  \n")
    result_code_data <- paste0(result_code_data, "   int L_num;\n")
    
    result[["con_sbhaz"]] <- rep(1, NCOL(Mbasis))
    result_code_data <- paste0(result_code_data, "   vector<lower=0>[L_num] con_sbhaz;\n")
    
    result[["Mbasis"]] <- Mbasis
    result_code_data <- paste0(result_code_data, "   //Baseline hazard function basis  \n")
    result_code_data <- paste0(result_code_data, "   matrix[N_num, L_num] Mbasis;\n")
    
    result[["Ibasis"]] <- Ibasis
    result_code_data <- paste0(result_code_data, "   //Cumulative baseline hazard function basis  \n")
    result_code_data <- paste0(result_code_data, "   matrix[N_num, L_num] Ibasis;\n")
    
    result[["Y"]] <- data[[formula$y_var]]
    result_code_data <- paste0(result_code_data, "   //Time-to-event outcome variable  \n")
    result_code_data <- paste0(result_code_data, "   real Y[N_num];\n")
    
    result_code_para <- paste0(result_code_para, "   //Baseline hazard spline coefficient   \n")
    result_code_para <- paste0(result_code_para, "   simplex[L_num] c;\n")
  }
  
  ## for SoFR with Gaussian outcomes
  if(family == "gaussian"){
    result_code_para <- paste0(result_code_para, "   real<lower=0> sigma;\n")
    result[["Y"]] <- data[[formula$y_var]]
    result_code_data <- paste0(result_code_data, "   //Outcome variable   \n")
    result_code_data <- paste0(result_code_data, "   real Y[N_num];\n")
    result[["plocation"]] <- mean(data[[formula$y_var]]) #!!! check do we really need priors for eta_0? (Line 78-81)
    result[["pscale"]] <- 10
    result_code_data <- paste0(result_code_data, "   real plocation;\n")
    result_code_data <- paste0(result_code_data, "   real pscale;\n")
  }
  
  ## for SoFR with binomial outcomes
  if(family == "binomial"){
    result[["Y"]] <- data[[formula$y_var]]
    result_code_data <- paste0(result_code_data, "   //Outcome variable   \n")
    result_code_data <- paste0(result_code_data, "   int Y[N_num];\n")
    result[["plocation"]] <- log((mean(data[[formula$y_var]]))/(1-mean(data[[formula$y_var]])))
    result[["pscale"]] <- 10
    result_code_data <- paste0(result_code_data, "   real plocation;\n")
    result_code_data <- paste0(result_code_data, "   real pscale;\n")
  }
  
  if(family == "functional"){
    ## for FoSR
    result[["M_num"]] <- dim(data[[formula$y_var]])[2]
    result_code_data <- paste0(result_code_data, "   int<lower=1> M_num; // Total number of observed functional time point\n")
    
    result[["Y_mat"]] <- data[[formula$y_var]]
    result_code_data <- paste0(result_code_data, "   matrix[N_num, M_num] Y_mat; // Functional response\n")
    
    #freq.fpca <- refund::fpca.face(unclass(data[[formula$y_var]]))
    freq.fpca <- refund::fpca.sc(unclass(data[[formula$y_var]]))
    
    data_temp <- data.frame(inx = 1:NROW(data))
    data_temp[["yindex.vec"]] <- matrix(rep(1:dim(data[[formula$y_var]])[2], NROW(data)), nrow = NROW(data))
    object = 1
    obj_exp <- paste0("object=s(yindex.vec, bs = \"",func_parameter[["type"]],"\", k =",func_parameter[["df"]],")")
    
    eval(parse(text = obj_exp))
    dk <- ExtractData(object, data_temp, NULL)
    crspline.cons <- mgcv::smooth.construct(object, dk$data, dk$knots)
    beta.base <- crspline.cons$X
    
    beta.S <- crspline.cons$S[[1]]
    maXX <- norm(beta.base, type = "I") ^ 2
    maS <- norm(beta.S) / maXX
    beta.S = beta.S / maS
    
    result[["J_num"]] <- ifelse(is.null(dim(freq.fpca$efunctions)[2]), 1, dim(freq.fpca$efunctions)[2])
    result_code_data <- paste0(result_code_data, "    int<lower=1> J_num; // Number of FPCA eigenfunctions \n")
    
    result[["Phi_mat"]] <- t(freq.fpca$efunctions)
    result_code_data <- paste0(result_code_data, "   matrix[J_num, M_num] Phi_mat; // Matrix of FPCA eigenfunctions \n")
    
    #result[["K_num"]] <- ifelse(is.null(dim(freq.fpca$efunctions)[2]), 1, dim(freq.fpca$efunctions)[2])
    result[["K_num"]] <- dim(crspline.cons$X)[2]
    result_code_data <- paste0(result_code_data, "   int<lower=1> K_num; // Number of spline basis \n")
    
    result[["Psi_mat"]] <- t(crspline.cons$X)
    result_code_data <- paste0(result_code_data, "   matrix[K_num, M_num] Psi_mat; // Matrix of spline basis  \n")
    
    result[["S_mat"]] <- beta.S
    result_code_data <- paste0(result_code_data, "   matrix[K_num, K_num] S_mat; // Penalty matrix  \n")
    
    result[["P_num"]] <- NCOL(X_scalar)
    result_code_data <- paste0(result_code_data, "   //Number of scalar predictors\n")
    result_code_data <- paste0(result_code_data, "   int<lower=0> P_num;\n")
    
    if(!is.null(X_scalar)){
      result[["X_mat"]] <- X_scalar
      result_code_data <- paste0(result_code_data, "   matrix[N_num, P_num] X_mat; // Design matrix for the scalar predictor\n")
      
      result_code_para <- paste0(result_code_para, "   matrix[K_num, P_num] beta; // Spline coefficients\n")
      result_code_para <- paste0(result_code_para, "   matrix[N_num, J_num] xi; // FPCA Scores\n")
      result_code_para <- paste0(result_code_para, "   real<lower=0> sigma_eps; // Standard deviation of independent error \n")
      result_code_para <- paste0(result_code_para, "   vector<lower=0>[P_num] sigma; // Smoothing parameter\n")
      result_code_para <- paste0(result_code_para, "   vector<lower=0>[J_num] lambda; // FPCA eigenvalues\n")
      
    }
  }else{
    # For SoFR and Functional Cox
    result[["K_num"]] <- NCOL(X_scalar)
    result_code_data <- paste0(result_code_data, "   //Number of scalar predictors   \n")
    result_code_data <- paste0(result_code_data, "   int<lower=0> K_num;\n")
    if(!is.null(X_scalar)){
      result[["Z_mat"]] <- X_scalar
      result_code_data <- paste0(result_code_data, "   //Matrix of scalar predictors   \n")
      result_code_data <- paste0(result_code_data, "   matrix[N_num,K_num] Z_mat;\n")
      result_code_para <- paste0(result_code_para, "   vector[K_num] gamma;\n")
    }
  }
  if(length(formula$func_var) > 0){ ## if we have functional predictors
    for(i in 1:length(func_comp)){
      if(func_comp[i]){
        ## Joint FPCA: replace the observed functional predictor with its FPCA
        ## representation and treat the FPC scores as parameters that are
        ## sampled jointly with the regression coefficients (Tutorial Section 4).
        func_data <- ext_func_pred_fpca(term = formula$func_var[[i]], data = data)
        ## In the joint FPCA case the random/fixed effect coefficients are
        ## already in the (eigendecomp, E)-reparameterized space, so no
        ## transformation is applied at reconstruction time. We use the identity
        ## as a placeholder so downstream bookkeeping stays consistent.
        trans.mat[[i]] <- diag(NROW(func_data$Xr) + NROW(func_data$Xf))

        ## Number of FPCA components and functional time points for predictor i
        result[[paste0("J_num_", i)]] <- func_data$J_num
        result_code_data <- paste0(result_code_data, "   //Number of FPCA components for functional predictor ", i, "\n")
        result_code_data <- paste0(result_code_data, "   int<lower=1> ", paste0("J_num_", i), ";\n")

        result[[paste0("M_num_", i)]] <- func_data$M_num
        result_code_data <- paste0(result_code_data, "   //Number of functional time points for functional predictor ", i, "\n")
        result_code_data <- paste0(result_code_data, "   int<lower=1> ", paste0("M_num_", i), ";\n")

        ## FPCA eigenfunctions, initial scores, and centered functional data
        result[[paste0("Phi_mat_", i)]] <- t(func_data$Phi)
        result_code_data <- paste0(result_code_data, "   //FPCA eigenfunctions for functional predictor ", i, "\n")
        result_code_data <- paste0(result_code_data, "   matrix[", paste0("J_num_", i), ", ", paste0("M_num_", i), "] ", paste0("Phi_mat_", i), ";\n")

        result[[paste0("xi_hat_", i)]] <- func_data$xi_hat
        result_code_data <- paste0(result_code_data, "   //Initial FPCA scores (prior mean) for functional predictor ", i, "\n")
        result_code_data <- paste0(result_code_data, "   matrix[N_num, ", paste0("J_num_", i), "] ", paste0("xi_hat_", i), ";\n")

        result[[paste0("M_mat_", i)]] <- func_data$M_mat
        result_code_data <- paste0(result_code_data, "   //Observed functional data for functional predictor ", i, "\n")
        result_code_data <- paste0(result_code_data, "   matrix[N_num, ", paste0("M_num_", i), "] ", paste0("M_mat_", i), ";\n")

        ## Xr (penalized) and Xf (unpenalized) parts of the spline basis
        ## reparameterized into the FPC eigenfunction space. In the joint FPCA
        ## case these have shape (Kr x J_num) and (Kf x J_num).
        result[[paste0("Kr_", i)]] <- NROW(func_data$Xr)
        result_code_data <- paste0(result_code_data, "   int<lower=0> ", paste0("Kr_", i), ";\n")

        result[[paste0("X_mat_r_", i)]] <- func_data$Xr
        result_code_data <- paste0(result_code_data, "   matrix[", paste0("Kr_", i), ", ", paste0("J_num_", i), "] ", paste0("X_mat_r_", i), ";\n")

        result[[paste0("Kf_", i)]] <- NROW(func_data$Xf)
        result_code_data <- paste0(result_code_data, "   int<lower=0> ", paste0("Kf_", i), ";\n")

        result[[paste0("X_mat_f_", i)]] <- func_data$Xf
        result_code_data <- paste0(result_code_data, "   matrix[", paste0("Kf_", i), ", ", paste0("J_num_", i), "] ", paste0("X_mat_f_", i), ";\n")

        ## Parameters block: FPC scores, FPC eigenvalue SDs, FPC residual SD,
        ## and the spline coefficients (in the random/fixed reparameterization).
        result_code_para <- paste0(result_code_para, "   //Joint FPCA scores for functional predictor ", i, "\n")
        result_code_para <- paste0(result_code_para, "   matrix[N_num, ", paste0("J_num_", i), "] ", paste0("xi_", i), ";\n")
        result_code_para <- paste0(result_code_para, "   //Joint FPCA eigenvalue SDs for functional predictor ", i, "\n")
        result_code_para <- paste0(result_code_para, "   vector<lower=0>[", paste0("J_num_", i), "] ", paste0("lambda_", i), ";\n")
        result_code_para <- paste0(result_code_para, "   //Joint FPCA residual SD for functional predictor ", i, "\n")
        result_code_para <- paste0(result_code_para, "   real<lower=0> ", paste0("sigma_e_", i), ";\n")
        result_code_para <- paste0(result_code_para, "   vector[", paste0("Kr_", i), "] ", paste0("zbr_", i), ";\n")
        result_code_para <- paste0(result_code_para, "   real<lower=0> ", paste0("sigmabr_", i), ";\n")
        result_code_para <- paste0(result_code_para, "   vector[", paste0("Kf_", i), "] ", paste0("bf_", i), ";\n")

        ## Transformed parameters
        result_code_transpara <- paste0(result_code_transpara, "   vector[", paste0("Kr_", i), "] ", paste0("br_", i), ";\n")
        result_code_transpara <- paste0(result_code_transpara, "   ", paste0("br_", i), " = ", paste0("sigmabr_", i), " * ", paste0("zbr_", i), ";\n")
      }
      if(!func_comp[i]){
        func_data <- ext_func_pred_obs(term = formula$func_var[[i]], data = data) ## organize functional predictors as input for Stan
        trans.mat[[i]] <- func_data$trans.mat

        ## Xr
        result[[paste0("Kr_",i)]] <- NCOL(func_data$Xr)
        result_code_data <- paste0(result_code_data, "   int<lower=0> ", paste0("Kr_",i), ";\n")

        result[[paste0("X_mat_r_",i)]] <- func_data$Xr
        result_code_data <- paste0(result_code_data, "   matrix[N_num, ", paste0("Kr_",i),"] ", paste0("X_mat_r_",i), ";\n")
        result_code_para <- paste0(result_code_para, "   vector[", paste0("Kr_",i),"] ", paste0("zbr_",i), ";\n")
        result_code_para <- paste0(result_code_para, "   real<lower=0>", paste0("sigmabr_",i), ";\n")
        result_code_transpara <- paste0(result_code_transpara, "   vector[", paste0("Kr_",i),"] ", paste0("br_",i),";\n")
        result_code_transpara <- paste0(result_code_transpara, "   ", paste0("br_",i), " = ", paste0("sigmabr_",i), " * ", paste0("zbr_",i), ";\n")

        ## Xf
        result[[paste0("Kf_",i)]] <- NCOL(func_data$Xf)
        result_code_data <- paste0(result_code_data, "   int<lower=0> ", paste0("Kf_",i), ";\n")

        result[[paste0("X_mat_f_",i)]] <- func_data$Xf
        result_code_data <- paste0(result_code_data, "   matrix[N_num, ", paste0("Kf_",i),"] ", paste0("X_mat_f_",i), ";\n")
        result_code_para <- paste0(result_code_para, "   vector[", paste0("Kf_",i),"] ", paste0("bf_",i), ";\n")
      }
    }
  }
  result_code_data <- paste0(result_code_data, "}")
  result_code_para <- paste0(result_code_para, "}")
  result_code_transpara <- paste0(result_code_transpara, "}")
  
  ## Center all scalar variables before putting them to Stan #!!! how would you handle categorical variables?
  if(family != "functional"){
    if(length(formula$scalar_var) > 0){
      result_code_transdata <- "transformed data {\n   matrix[N_num, K_num] X_sc;\n   vector[K_num] mean_Xs;\n   "
      result_code_transdata <- paste0(result_code_transdata, "for (i in 1:K_num) {\n      mean_Xs[i] = mean(Z_mat[, i]);\n
                                 X_sc[, i] = Z_mat[, i] - mean_Xs[i];\n   }\n}")
    }else{
      result_code_transdata <- "transformed data {\n}"
    }
  }else{
    result_code_transdata <- "transformed data {\n}"
  }
  
  return(list(standata = result, ## data as input for Stan
              stancode_data = result_code_data, ## code data block
              stancode_transdata = result_code_transdata, ## code transformed data block
              stancode_para = result_code_para, ## code parameters block
              stancode_transpara = result_code_transpara, ## code transformed parameters block
              trans.mat = trans.mat, ## transformation matrix used to recover functional coefficients
              X_scalar = X_scalar) ## scalar variables data as input for Stan
  )
}

#----------------------------------------------------------------------------
#' A function for extracting the data for the scalar predictors and organize them as input for Stan
#----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
#' 
ext_scalar_pred = function(term, data){
  if(length(term) == 0){
    return(NULL)
  }else{
    X_use <- matrix(1, nrow = NROW(data), ncol = 1)
    for(i in 1:length(term)){
      if(length(term[[i]]) > 1){ ## if there are interaction terms
        term_now <- term[[i]][[1]]
        for(k in 2:length(term[[i]])){
          term_now <- paste0(term_now, "*", term[[i]][[k]])
        }
        X_itm <- stats::model.matrix(stats::as.formula(paste0("~", term_now)), data = data)
      }else{ ## if there are just additive terms
        X_itm <- stats::model.matrix(stats::as.formula(paste0("~", term[[i]])), data = data)
        term_now <- term[[i]]
      }
      
      colsave <- colnames(X_itm)
      X_itm <- as.matrix(X_itm[, -1]) ## remove intercept
      colnames(X_itm) <- colsave[-1]
      X_use <- cbind(X_use,X_itm) ## add the variable as a new column of the design matrix
    }
    X_use <- X_use[, -1, drop = FALSE] ## remove intercept
    good_inx <- c()
    for(j in 1:ncol(X_use)){
      if(sum(abs(X_use[, j])) == 0){ ## avoid the situation when all values are 0
        
      }else{
        good_inx <- c(good_inx,j)
      }
    }
    X_use <- X_use[, good_inx, drop = FALSE]
    return(X_use = X_use)
  }
}

#----------------------------------------------------------------------------
#' A function for extracting the data for the functional predictors and organize them as input for Stan
#----------------------------------------------------------------------------
#' @keywords internal
#' @noRd
#'
ext_func_pred_obs = function(term, data){
  if(term[[1]] == "s"){ ## currently only support s() as a functional term input
    obj <- format(term)
    sm = 1
    eval(parse(text = paste0("sm = mgcv:::smoothCon(",obj,", data = data, absorb.cons = TRUE, diagonal.penalty = TRUE)")))
    re <- mgcv::smooth2random(sm[[1]], names(data), type = 2)
    Xr <- re$rand$Xr
    Xf <- re$Xf
    trans.mat <- re$trans.U
  }else{

  }
  return(list(Xr = Xr,
              Xf = Xf,
              trans.mat = trans.mat))
}

#----------------------------------------------------------------------------
#' A helper to identify the functional data matrix and the integration weights
#' from a smooth-term \code{by=} expression.
#'
#' The convention follows the tutorial supplementary code: when \code{by} is a
#' product expression such as \code{lmat * wmat} (or more generally
#' \code{a * b * ... * wmat}), the rightmost variable (\code{wmat}) is treated
#' as the actual functional data matrix, and the product of all other
#' variables is treated as the integration weight matrix.
#'
#' @keywords internal
#' @noRd
#'
parse_by_func = function(term, data){
  if (is.name(term[["by"]])) {
    funcmat <- data[[as.character(term[["by"]])]]
    M_num <- ncol(funcmat)
    lmat <- matrix(1 / M_num, nrow = nrow(funcmat), ncol = M_num)
  } else {
    funcmat <- data[[as.character(term[["by"]][[3]])]]
    M_num <- ncol(funcmat)
    rest_term <- term[["by"]][[2]]
    if (is.name(rest_term)) {
      lmat <- data[[as.character(rest_term)]]
    } else {
      lmat <- data[[as.character(rest_term[[3]])]]
      curr_term <- rest_term
      while (!is.name(curr_term[[2]])) {
        curr_term <- curr_term[[2]]
        lmat <- lmat * data[[as.character(curr_term[[3]])]]
      }
      lmat <- lmat * data[[as.character(curr_term[[2]])]]
    }
  }
  return(list(funcmat = funcmat, lmat = lmat, M_num = M_num))
}

#----------------------------------------------------------------------------
#' A function for extracting the FPCA-based design matrices for a functional
#' predictor that is to be modeled jointly with the regression model.
#'
#' Given a functional predictor \eqn{W_i(s)}, fpca.sc is run to obtain
#' the eigenfunctions \eqn{\phi_j(s)}, the initial estimated scores
#' \eqn{\hat\xi_{ij}}, and the mean function \eqn{\mu(s)}. The functional
#' coefficient \eqn{\beta(s)} is represented in the spline basis
#' \eqn{\beta(s) = \sum_k \tilde\beta_k \psi_k(s)}, after the same eigen-based
#' reparameterization used in the Tutorial supplementary code (Section 4).
#'
#' Returns the design-matrix pieces \code{Xr} (\code{rank} x \code{J_num}) and
#' \code{Xf} (\code{(K_num - rank)} x \code{J_num}), the FPCA pieces
#' \code{Phi}, \code{xi_hat}, the observed functional data \code{M_mat},
#' and the basis info needed for reconstructing \eqn{\beta(s)}.
#'
#' @keywords internal
#' @noRd
#'
ext_func_pred_fpca = function(term, data, npc = NULL){
  if(term[[1]] != "s"){
    stop("Joint FPCA currently only supports s() functional terms.")
  }

  ## Identify the functional data matrix and integration weights
  parsed <- parse_by_func(term, data)
  funcmat <- parsed$funcmat
  lmat <- parsed$lmat
  M_num <- parsed$M_num
  N_num <- nrow(funcmat)

  ## Initial frequentist FPCA on the functional data (uncentered)
  fpca_fit <- refund::fpca.sc(Y = unclass(funcmat), npc = npc)
  Phi <- fpca_fit$efunctions  ## M_num x J_num
  xi_hat <- fpca_fit$scores   ## N_num x J_num
  if(is.null(dim(Phi))){
    Phi <- matrix(Phi, ncol = 1)
  }
  if(is.null(dim(xi_hat))){
    xi_hat <- matrix(xi_hat, ncol = 1)
  }
  J_num <- ncol(Phi)

  ## We follow the Tutorial supplementary code (Section 4) and pass the
  ## observed functional data directly to the joint Stan model. The prior
  ## centers xi on the initial FPCA scores xi_hat.
  M_mat_for_stan <- unclass(funcmat)

  ## Build the spline basis for beta(s) following the Tutorial supplementary
  ## code (Section 4): use smooth.construct directly (no absorb.cons), so the
  ## same Psi_mat / S_mat / eigendecomp will be used in both data preparation
  ## and the reconstruction step.
  knots <- NULL
  obj <- 1
  eval(parse(text = paste0("obj = ", format(term))))
  dk <- ExtractData(obj, data, knots)
  splinecons <- mgcv::smooth.construct(obj, dk$data, dk$knots)
  Psi_mat <- splinecons$X     ## M_num x K_num
  S_mat <- splinecons$S[[1]]  ## K_num x K_num
  rank <- splinecons$rank
  K_num <- ncol(Psi_mat)

  ## Scale the penalty matrix for numerical stability (matches the rest of
  ## the package and the Tutorial supplementary code)
  maXX <- norm(Psi_mat, type = "I") ^ 2
  maS <- norm(S_mat) / maXX
  S_mat_scaled <- S_mat / maS
  eigendecomp <- eigen(S_mat_scaled, symmetric = TRUE)

  ## X_mat_t[j,k] = integral phi_j(s) psi_k(s) ds
  ##              ~= sum_t l_t * phi_j(t) * psi_k(t)
  ## Use the first row of lmat as the integration weights (assumes all
  ## subjects share the same functional grid / integration weights).
  l_vec <- as.numeric(lmat[1, ])
  X_mat_t <- t(Phi) %*% (l_vec * Psi_mat)  ## J_num x K_num

  ## Apply the eigenbasis transformation
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
  X_mat_t <- t(t(X_mat_t) / E)

  ## Split into penalized (Xr) and unpenalized (Xf) parts. For the joint FPCA
  ## design these are stored as (Kr x J_num) and (Kf x J_num), matching
  ## the Tutorial supplementary code (Section 4).
  Xr <- t(X_mat_t[, 1:rank, drop = FALSE])                      ## rank x J_num
  if (rank < K_num) {
    Xf <- t(X_mat_t[, (rank + 1):K_num, drop = FALSE])          ## (K_num - rank) x J_num
  } else {
    Xf <- matrix(0, nrow = 0, ncol = J_num)
  }

  return(list(Xr = Xr,
              Xf = Xf,
              Phi = Phi,
              xi_hat = xi_hat,
              M_mat = M_mat_for_stan,
              J_num = J_num,
              M_num = M_num,
              eigendecomp = eigendecomp,
              E = E,
              base_mat = Psi_mat,
              rank = rank,
              K_num = K_num))
}
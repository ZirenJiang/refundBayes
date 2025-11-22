#' Functions to build the data and code for subsequent Stan modeling
#'

#----------------------------------------------------------------------------
#' The main function for constructing the data and code used as input for Stan
#----------------------------------------------------------------------------

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
    
    result[["cens"]] <- data[[cens]]
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
    result[["T"]] <- dim(data[[formula$y_var]])[2]
    result_code_data <- paste0(result_code_data, "   //Number of observed time points\n")
    result_code_data <- paste0(result_code_data, "   int<lower=1> T;\n")
    
    result[["Y"]] <- data[[formula$y_var]]
    result_code_data <- paste0(result_code_data, "   //Functional outcome for FoSR\n")
    result_code_data <- paste0(result_code_data, "   matrix[N_num, T] Y;\n")
    freq.fpca <- refund::fpca.face(unclass(data[[formula$y_var]]))
    
    data_temp <- data.frame(inx = 1:NROW(data))
    data_temp[["yindex.vec"]] <- matrix(rep(1:dim(data[[formula$y_var]])[2], NROW(data)), nrow = NROW(data))
    obj_exp <- paste0("object=s(yindex.vec, bs = \"",func_parameter[["type"]],"\", k =",func_parameter[["df"]],")")
    
    eval(parse(text = obj_exp))
    dk <- ExtractData(object, data_temp, NULL)
    crspline.cons <- mgcv::smooth.construct(object, dk$data, dk$knots)
    beta.base <- crspline.cons$X
    
    beta.S <- crspline.cons$S[[1]]
    maXX <- norm(beta.base, type = "I") ^ 2
    maS <- norm(beta.S) / maXX
    beta.S = beta.S / maS
    
    result[["ncol_phi"]] <- dim(beta.base)[2]
    result_code_data <- paste0(result_code_data, "   //Number of FPCA eigenfunctions\n")
    result_code_data <- paste0(result_code_data, "   int ncol_phi; \n")
    
    result[["phi"]] <- t(beta.base)
    result_code_data <- paste0(result_code_data, "   //Matrix of FPCA eigenfunctions\n")
    result_code_data <- paste0(result_code_data, "   matrix[ncol_phi, T] phi; \n")
    
    result[["ncol_psi_mat"]] <- ifelse(is.null(dim(freq.fpca$efunctions)[2]), 1, dim(freq.fpca$efunctions)[2])
    result_code_data <- paste0(result_code_data, "   //Number of spline basis\n")
    result_code_data <- paste0(result_code_data, "   int ncol_psi_mat; \n")
    
    result[["psi_mat"]] <- t(freq.fpca$efunctions)
    result_code_data <- paste0(result_code_data, "   //Matrix of spline basis\n")
    result_code_data <- paste0(result_code_data, "   matrix[ncol_psi_mat, T] psi_mat;  \n")
    
    result[["S_beta"]] <- beta.S
    result_code_data <- paste0(result_code_data, "   //Penalty matrix\n")
    result_code_data <- paste0(result_code_data, "   matrix[ncol_phi, ncol_phi] S_beta;  \n")
    
    result[["K_num"]] <- NCOL(X_scalar)
    result_code_data <- paste0(result_code_data, "   //Number of scalar predictors\n")
    result_code_data <- paste0(result_code_data, "   int<lower=0> K_num;\n")
    
    if(!is.null(X_scalar)){
      result[["Z_mat"]] <- X_scalar
      result_code_data <- paste0(result_code_data, "   //Matrix of scalar predictors\n")
      result_code_data <- paste0(result_code_data, "   matrix[N_num,K_num] Z_mat;\n")
      
      result_code_para <- paste0(result_code_para, "   //FPCA scores \n")
      result_code_para <- paste0(result_code_para, "   matrix[N_num,ncol_psi_mat] zxi;\n")
      result_code_para <- paste0(result_code_para, "   real<lower=0> sigma_epis;\n")
      result_code_para <- paste0(result_code_para, "   //Smoothing parameters\n")
      result_code_para <- paste0(result_code_para, "   vector<lower=0>[K_num] sigma_b;\n")
      result_code_para <- paste0(result_code_para, "   vector<lower=0>[ncol_psi_mat] sigma_k;\n")
      if(FALSE){
        result_code_para <- paste0(result_code_para, "   row_vector[ncol_phi] b;\n")
      }else{
        result_code_para <- paste0(result_code_para, "   //Functional effect spline coefficients\n")
        result_code_para <- paste0(result_code_para, "   matrix[K_num,ncol_phi] b;\n")
      }
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
    if(sum(func_comp) > 0){
      # For now, we have not implemented the joint model of FPCA
    }
    for(i in 1:length(func_comp)){
      if(func_comp[i]){
        
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
  if(length(formula$scalar_var) > 0){
    result_code_transdata <- "transformed data {\n   matrix[N_num, K_num] X_sc;\n   vector[K_num] mean_Xs;\n   "
    result_code_transdata <- paste0(result_code_transdata, "for (i in 1:K_num) {\n      mean_Xs[i] = mean(Z_mat[, i]);\n
                                 X_sc[, i] = Z_mat[, i] - mean_Xs[i];\n   }\n}")
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
        X_itm <- stats::model.matrix(as.formula(paste0("~", term_now)), data = data)
      }else{ ## if there are just additive terms
        X_itm <- stats::model.matrix(as.formula(paste0("~", term[[i]])), data = data)
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

ext_func_pred_obs = function(term, data){
  if(term[[1]] == "s"){ ## currently only support s() as a functional term input
    obj <- format(term)
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

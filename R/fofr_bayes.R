#----------------------------------------------------------------------------
#' Bayesian Function-on-Function Regression
#'
#' Fit the Bayesian Function-on-Function Regression (FoFR) model using Stan.
#'
#' The Bayesian FoFR model extends the function-on-scalar regression (FoSR)
#' framework by allowing functional predictors in addition to (optional) scalar
#' predictors. The bivariate coefficient function \eqn{\beta(s,t)} is represented
#' using a tensor product of the predictor-domain spline basis and the
#' response-domain spline basis. The model is specified using the same syntax
#' as in the R mgcv package.
#'
#' @param formula Functional regression formula, with the same syntax as that in the R mgcv package.
#'   The response must be a matrix (functional response) and at least one \code{s(..., by = ...)}
#'   term must be present for the functional predictor(s).
#' @param data A data frame containing data of all scalar and functional variables used in the model.
#' @param joint_FPCA A True/False vector of the same length of the number of functional predictors,
#'   indicating whether jointly modeling FPCA for the functional predictors. Default to NULL.
#' @param spline_type Type of spline basis for modelling the residual process and the response-domain
#'   component of the bivariate coefficient. Default to "bs".
#' @param spline_df Degrees of freedom for the spline basis for modelling the response-domain
#'   component. Default to 10.
#' @param runStan True/False variable for whether to run the Stan program. If False, the function
#'   only generates the Stan code and data.
#' @param niter Total number of Bayesian iterations.
#' @param nwarmup Number of warmup (burnin) iterations for posterior sampling.
#' @param nchain Number of chains for posterior sampling. Default to 3.
#' @param ncores Number of cores to use when executing the chains in parallel. Default to 1.
#'
#' @return A list containing:
#' \item{stanfit}{The Stan fit object.}
#' \item{spline_basis}{Basis functions used to reconstruct the functional coefficients from posterior samples.}
#' \item{stancode}{A character string containing the code to fit the Stan model.}
#' \item{standata}{A list containing the data to fit the Stan model.}
#' \item{scalar_func_coef}{A 3-d array (n_samples x P x M) containing posterior samples of scalar predictor
#'   coefficient functions. Each slice \code{[,p,]} is the coefficient function for the p-th scalar predictor.
#'   NULL if no scalar predictors.}
#' \item{bivar_func_coef}{A list of 3-d arrays. Each element corresponds to one functional predictor
#'   and is an array of dimension (n_samples x S_grid x M), representing posterior samples of the
#'   bivariate coefficient function \eqn{\beta(s,t)}.}
#' \item{func_coef}{A 3-d array for the scalar predictor coefficient functions, stored for compatibility
#'   with the \code{plot.refundBayes} method. Same as \code{scalar_func_coef}.}
#' \item{family}{Model family: "fofr".}
#'
#' @author Erjia Cui \email{ecui@@umn.edu}, Ziren Jiang \email{jian0746@@umn.edu}
#'
#' @references Jiang, Z., Crainiceanu, C., and Cui, E. (2025). Tutorial on Bayesian Functional Regression Using Stan. \emph{Statistics in Medicine}, 44(20-22), e70265.
#'
#' @examples
#' \dontrun{
#' # Simulate data for a Function-on-Function Regression model
#' set.seed(1)
#' n  <- 100  # number of subjects
#' L  <- 30   # number of functional predictor domain points
#' M  <- 30   # number of functional response domain points
#' sindex <- seq(0, 1, length.out = L)  # predictor domain grid
#' tindex <- seq(0, 1, length.out = M)  # response domain grid
#'
#' # Functional predictor
#' X_func <- matrix(rnorm(n * L), nrow = n)
#'
#' # Scalar predictor
#' age <- rnorm(n)
#'
#' # True bivariate coefficient beta(s, t)
#' beta_true <- outer(sin(2 * pi * sindex), cos(2 * pi * tindex))
#'
#' # True scalar coefficient function
#' alpha_true <- sin(pi * tindex)
#'
#' # Generate functional response: Y(t) = age * alpha(t) + integral X(s) beta(s,t) ds + error
#' Y_mat <- outer(age, alpha_true) + X_func %*% beta_true / L +
#'           matrix(rnorm(n * M, sd = 0.3), nrow = n)
#'
#' dat <- data.frame(age = age)
#' dat$Y_mat  <- Y_mat
#' dat$X_func <- X_func
#' dat$sindex <- matrix(rep(sindex, n), nrow = n, byrow = TRUE)
#'
#' # Fit the Bayesian FoFR model
#' fit_fofr <- fofr_bayes(
#'   formula    = Y_mat ~ age + s(sindex, by = X_func, bs = "cr", k = 10),
#'   data       = dat,
#'   spline_type = "bs",
#'   spline_df  = 10,
#'   niter      = 2000,
#'   nwarmup    = 1000,
#'   nchain     = 3
#' )
#'
#' # Examine the bivariate coefficient beta(s, t) (posterior mean)
#' beta_est <- apply(fit_fofr$bivar_func_coef[[1]], c(2, 3), mean)
#' image(sindex, tindex, beta_est,
#'       xlab = "s (predictor domain)", ylab = "t (response domain)",
#'       main = expression(hat(beta)(s, t)))
#'
#' # Fit a FoFR model with functional predictors only (no scalar predictors)
#' fit_fofr2 <- fofr_bayes(
#'   formula    = Y_mat ~ s(sindex, by = X_func, bs = "cr", k = 10),
#'   data       = dat,
#'   niter      = 2000,
#'   nwarmup    = 1000,
#'   nchain     = 3
#' )
#' }
#'
#' @import mgcv
#' @import splines2
#' @import rstan
#' @importFrom refund fpca.face
#' @importFrom stats gaussian
#' @export fofr_bayes

fofr_bayes <- function(formula, data,
                       joint_FPCA = NULL, runStan = TRUE,
                       niter = 3000, nwarmup = 1000, nchain = 3, ncores = 1,
                       spline_type = "bs", spline_df = 10){

  ## Set family to be "functional" for data preparation (functional response)
  family <- "functional"

  func_parameter <- list(type = spline_type,
                         df = spline_df)
  intercept <- TRUE

  # Organize the formula as the input for refundBayesdata() and refundBayesmodel()
  formula_use <- refundBayesformula(formula)

  # Validate that functional predictors are present
  if(length(formula_use$func_var) == 0){
    stop("No functional predictors found in the formula. Use fosr_bayes() for function-on-scalar regression.")
  }

  if(is.null(joint_FPCA)){
    func_comp <- rep(FALSE, length(formula_use$func_var))
  }else{
    func_comp <- joint_FPCA
  }

  if(length(formula_use$func_var) != length(func_comp)){
    stop("Length of joint_FPCA must equal the total number of functional predictors!")
  }

  # Construct the data using refundBayesdata
  # The standata (actual data values) are correct for all parts:
  #   - Functional response: Y_mat, M_num, J_num, Phi_mat, K_num, Psi_mat, S_mat
  #   - Scalar predictors: P_num, X_mat
  #   - Functional predictors: Kr_i, X_mat_r_i, Kf_i, X_mat_f_i
  datacode_use <- refundBayesdata(formula_use, data, family, func_comp, intercept,
                                  func_parameter = func_parameter)
  Standata_use <- datacode_use$standata

  # Use the data block from refundBayesdata (all declarations are correct)
  Stancode_data <- datacode_use$stancode_data

  # Use the transformed data block from refundBayesdata
  Stancode_transdata <- datacode_use$stancode_transdata

  has_scalar <- !is.null(datacode_use$X_scalar)
  n_func <- length(func_comp)

  ##--------------------------------------------------------------------------
  ## Build the parameters block
  ## For FoFR, functional predictor coefficients need to be MATRICES
  ## (K_pred_basis x K_resp_basis) to represent the bivariate coefficient beta(s,t),
  ## unlike SoFR where they are vectors.
  ##--------------------------------------------------------------------------
  Stancode_para <- paste0("parameters{ \n")

  # Scalar predictor parameters (same as FoSR)
  if(has_scalar){
    Stancode_para <- paste0(Stancode_para,
                            "   matrix[K_num, P_num] beta; // Spline coefficients for scalar predictor effects\n",
                            "   vector<lower=0>[P_num] sigma; // Smoothing parameters for scalar predictors\n")
  }

  # FPCA and residual parameters (always present)
  Stancode_para <- paste0(Stancode_para,
                          "   matrix[N_num, J_num] xi; // FPCA Scores\n",
                          "   real<lower=0> sigma_eps; // Standard deviation of independent error\n",
                          "   vector<lower=0>[J_num] lambda; // FPCA eigenvalues\n")

  # Functional predictor parameters (bivariate coefficient representation)
  for(i in 1:n_func){
    if(!func_comp[i]){
      Stancode_para <- paste0(Stancode_para,
                              "   // Functional predictor ", i, " (bivariate coefficient)\n",
                              "   matrix[", paste0("Kr_", i), ", K_num] ", paste0("zr_", i),
                              "; // Standardized random effects\n",
                              "   real<lower=0> ", paste0("sigmabr_", i),
                              "; // s-direction smoothing parameter\n",
                              "   matrix[", paste0("Kf_", i), ", K_num] ", paste0("bf_", i),
                              "; // Fixed effects\n",
                              "   real<lower=0> ", paste0("sigma_t_", i),
                              "; // t-direction smoothing parameter\n")
    }
  }
  Stancode_para <- paste0(Stancode_para, "}")

  ##--------------------------------------------------------------------------
  ## Build the transformed parameters block
  ##--------------------------------------------------------------------------
  Stancode_transpara <- paste0("transformed parameters { \n real lprior = 0;\n")
  for(i in 1:n_func){
    if(!func_comp[i]){
      Stancode_transpara <- paste0(Stancode_transpara,
                                   "   matrix[", paste0("Kr_", i), ", K_num] ", paste0("br_", i), ";\n",
                                   "   ", paste0("br_", i), " = ", paste0("sigmabr_", i), " * ",
                                   paste0("zr_", i), ";\n")
    }
  }
  Stancode_transpara <- paste0(Stancode_transpara, "}")

  ##--------------------------------------------------------------------------
  ## Build the model block
  ##--------------------------------------------------------------------------
  Stancode_model <- paste0("model{ \n")
  Stancode_model <- paste0(Stancode_model, "   matrix[N_num, M_num] mu;\n")
  Stancode_model <- paste0(Stancode_model, "   // Fitted mean matrix \n")

  # Mean function: scalar predictor contributions + FPCA residual process
  if(has_scalar){
    Stancode_model <- paste0(Stancode_model,
                             "   mu = X_mat * beta' * Psi_mat + xi * Phi_mat;\n")
  }else{
    Stancode_model <- paste0(Stancode_model,
                             "   mu = xi * Phi_mat;\n")
  }

  # Add functional predictor contributions
  # Each functional predictor contributes:
  #   (X_mat_r_i * br_i + X_mat_f_i * bf_i) * Psi_mat
  # where br_i and bf_i are now matrices (K_pred x K_resp),
  # giving an N x K_resp intermediate that is projected onto the response domain via Psi_mat.
  for(i in 1:n_func){
    if(!func_comp[i]){
      Stancode_model <- paste0(Stancode_model,
                               "   // Functional predictor ", i, " contribution\n",
                               "   mu += (", paste0("X_mat_r_", i), " * ", paste0("br_", i),
                               " + ", paste0("X_mat_f_", i), " * ", paste0("bf_", i),
                               ") * Psi_mat;\n")
    }
  }

  # Log-likelihood for functional response
  Stancode_model <- paste0(Stancode_model,
                           "   // Log-likelihood for functional response\n",
                           "   target += - N_num * M_num * log(sigma_eps) / 2 - sum((mu - Y_mat)^2) / (2 * sigma_eps^2);\n")

  # Priors for scalar predictor spline coefficients (same as FoSR)
  if(has_scalar){
    Stancode_model <- paste0(Stancode_model,
                             "   // Prior for the penalized spline coefficients (scalar predictors) \n",
                             "   for(np in 1:P_num){\n",
                             "        target += (- beta[,np]' * S_mat * beta[,np]) / (2 * sigma[np]^2);\n",
                             "        target += inv_gamma_lpdf(sigma[np]^2|0.001,0.001);\n",
                             "   }\n")
  }

  # Priors for FPCA scores
  Stancode_model <- paste0(Stancode_model,
                           "    // Prior for the FPCA scores \n",
                           "   for(nj in 1:J_num){\n",
                           "        target += - N_num * log(lambda[nj]) / 2 - sum((xi[,nj])^2) / (2 * lambda[nj]^2);\n",
                           "        target += inv_gamma_lpdf(lambda[nj]^2|0.001,0.001);\n",
                           "   }\n")

  # Priors for functional predictor bivariate coefficients
  for(i in 1:n_func){
    if(!func_comp[i]){
      Stancode_model <- paste0(Stancode_model,
                               "   // Priors for functional predictor ", i, " bivariate coefficient\n",
                               "   // s-direction smoothness via random effects reparameterization\n",
                               "   to_vector(", paste0("zr_", i), ") ~ std_normal();\n",
                               "   target += inv_gamma_lpdf(", paste0("sigmabr_", i), "^2 | 0.0005, 0.0005);\n",
                               "   // t-direction smoothness via response-domain penalty matrix\n",
                               "   for(kk in 1:", paste0("Kr_", i), "){\n",
                               "        target += (- ", paste0("br_", i), "[kk,] * S_mat * ",
                               paste0("br_", i), "[kk,]') / (2 * ", paste0("sigma_t_", i), "^2);\n",
                               "   }\n",
                               "   for(kk in 1:", paste0("Kf_", i), "){\n",
                               "        target += (- ", paste0("bf_", i), "[kk,] * S_mat * ",
                               paste0("bf_", i), "[kk,]') / (2 * ", paste0("sigma_t_", i), "^2);\n",
                               "   }\n",
                               "   target += inv_gamma_lpdf(", paste0("sigma_t_", i), "^2 | 0.001, 0.001);\n")
    }
  }

  # Prior for residual variance
  Stancode_model <- paste0(Stancode_model,
                           "   target += inv_gamma_lpdf(sigma_eps^2|0.001,0.001);\n",
                           "}\n")

  # Combine all blocks of the Stan code
  Stancode_function <- NULL ## no manual function needed for FoFR
  Stancode_use <- paste0(Stancode_function, "\n",
                         Stancode_data, "\n",
                         Stancode_transdata, "\n",
                         Stancode_para, "\n",
                         Stancode_transpara, "\n",
                         Stancode_model)


  if(runStan){
    # Run the Stan program
    fit <- stan(model_code = Stancode_use,
                data = Standata_use,
                iter = niter,
                warmup = nwarmup,
                chains = nchain,
                cores = ncores)

    # Extract the basis used to reconstruct estimated functional coefficients
    func_effect_basis <- extract_basis(formula_use, data, func_comp)

    # Extract posterior samples from the Stan fit object
    fit.samp <- rstan::extract(fit, permuted = TRUE)

    # Reconstruct scalar predictor coefficient functions (same as FoSR)
    # Each scalar predictor effect is a matrix: n_samples x M
    if(has_scalar){
      scalar_func_coef <- array(dim = c(dim(fit.samp$beta)[1],
                                        Standata_use$P_num,
                                        dim(Standata_use$Psi_mat)[2]))
      for(iter in 1:dim(fit.samp$beta)[1]){
        scalar_func_coef[iter, , ] <- t(fit.samp$beta[iter, , ]) %*% Standata_use$Psi_mat
      }
    }else{
      scalar_func_coef <- NULL
    }

    # Reconstruct bivariate functional coefficients beta_i(s,t)
    # For each functional predictor i, the bivariate coefficient is:
    #   beta(s,t) = sum_l curve_l(s) * Psi_l(t)
    # where curve_l(s) is reconstructed from the posterior samples of br_i[,l] and bf_i[,l]
    # using the same transformation as in SoFR (trans.mat, eigendecomp, base_mat).
    Psi_mat <- Standata_use$Psi_mat  # K_num x M
    K_num_resp <- Standata_use$K_num
    M_resp <- Standata_use$M_num

    bivar_func_coef <- list()
    for(i in 1:n_func){
      if(!func_comp[i]){
        br_samp <- fit.samp[[paste0("br_", i)]]  # n_samp x Kr_i x K_num
        bf_samp <- fit.samp[[paste0("bf_", i)]]  # n_samp x Kf_i x K_num

        eigendecomp <- func_effect_basis[[i]][["eigendecomp"]]
        base_mat <- func_effect_basis[[i]][["base_mat"]]  # S_grid x n_basis
        E <- func_effect_basis[[i]][["E"]]
        trans_mat_i <- datacode_use$trans.mat[[i]]

        n_samp <- dim(br_samp)[1]
        S_grid <- nrow(base_mat)

        bivar_coef_i <- array(dim = c(n_samp, S_grid, M_resp))

        for(iter in 1:n_samp){
          # For each response basis column l, reconstruct the s-domain curve
          # using the same procedure as recon_fun_coef() in SoFR
          curve_matrix <- matrix(nrow = S_grid, ncol = K_num_resp)

          for(l in 1:K_num_resp){
            # Combine random and fixed effects for column l
            para_l <- c(br_samp[iter, , l], bf_samp[iter, , l])
            # Apply the transformation matrix U
            para_l_trans <- para_l %*% trans_mat_i
            # Undo the scaling transformation
            para_l_untilde <- (eigendecomp$vectors %*% diag(1/E)) %*% matrix(para_l_trans, ncol = 1)
            # Evaluate on the predictor-domain grid
            curve_matrix[, l] <- base_mat %*% para_l_untilde
          }

          # Combine with response-domain basis: beta(s,t) = curve_matrix %*% Psi_mat
          bivar_coef_i[iter, , ] <- curve_matrix %*% Psi_mat
        }

        bivar_func_coef[[i]] <- bivar_coef_i
      }
    }

  }else{
    # Do not run the Stan program
    fit <- NA ## no sampling
    func_effect_basis <- extract_basis(formula_use, data, func_comp)
    fit.samp <- NA ## no sampling
    scalar_func_coef <- NA ## no sampling
    bivar_func_coef <- NA ## no sampling
  }

  # Organize the results
  res.fit <- list(stanfit = fit,
                  spline_basis = func_effect_basis,
                  stancode = Stancode_use,
                  standata = Standata_use,
                  scalar_func_coef = scalar_func_coef,
                  bivar_func_coef = bivar_func_coef,
                  func_coef = scalar_func_coef,  ## for plot.refundBayes compatibility
                  family = "fofr")
  class(res.fit) <- "refundBayes"
  return(res.fit)
}

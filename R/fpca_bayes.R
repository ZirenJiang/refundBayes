#----------------------------------------------------------------------------
#' Bayesian Functional Principal Component Analysis
#'
#' Fit the Bayesian Functional Principal Component Analysis (FPCA) model using Stan.
#'
#' The Bayesian FPCA model is implemented following the tutorial by Jiang et al., 2025.
#' The model decomposes a dense functional response \eqn{Y_i(t)} into a smooth
#' population mean function \eqn{\mu(t)} and a sum of functional principal components,
#' \deqn{Y_i(t) = \mu(t) + \sum_{k=1}^{K} \xi_{ik} \phi_k(t) + \epsilon_i(t),}
#' where \eqn{\phi_k(t)} are orthonormal eigenfunctions obtained from an initial
#' frequentist FPCA fit (used as a fixed basis), and the mean function \eqn{\mu(t)},
#' the FPC scores \eqn{\xi_{ik}}, the eigenvalues \eqn{\lambda_k}, and the residual
#' variance are estimated via posterior sampling. The population mean function is
#' modelled with a penalized spline using the same syntax as in the R mgcv package.
#'
#' @param formula Functional formula of the form \code{Y_mat ~ 1}, where \code{Y_mat}
#'   is the functional response variable stored as a matrix column of \code{data}.
#' @param data A data frame containing the functional response variable used in the model.
#' @param npc Number of functional principal components. If NULL, it is selected
#'   automatically by the initial \code{refund::fpca.sc()} fit based on the
#'   percentage of variance explained. Default to NULL.
#' @param spline_type Type of spline basis for modelling the population mean function. Default to "bs".
#' @param spline_df Degrees of freedom for the spline basis for modelling the population mean function. Default to 10.
#' @param runStan True/False variable for whether to run the Stan program. If False, the function only generates the Stan code and data.
#' @param niter Total number of Bayesian iterations.
#' @param nwarmup Number of warmup (burnin) iterations for posterior sampling.
#' @param nchain Number of chains for posterior sampling. Default to 3.
#' @param ncores Number of cores to use when executing the chains in parallel. Default to 1.
#'
#' @return A list containing:
#' \item{stanfit}{The Stan fit object.}
#' \item{stancode}{A character string containing the code to fit the Stan model.}
#' \item{standata}{A list containing the data to fit the Stan model.}
#' \item{mu}{A matrix of posterior samples of the population mean function, where each row is one sample and each column is one location of the functional domain.}
#' \item{efunctions}{A matrix of the (fixed) eigenfunctions from the initial FPCA used as the FPC basis. Each column is one eigenfunction.}
#' \item{scores}{A 3-d array of posterior samples of FPC scores with dimensions (n_samples x n_subjects x npc).}
#' \item{evalues}{A matrix of posterior samples of FPC eigenvalue standard deviations, where each row is one sample and each column is one principal component.}
#' \item{sigma}{A vector of posterior samples of the residual standard deviation.}
#' \item{family}{Family type: "fpca".}
#'
#' @author Erjia Cui \email{ecui@@umn.edu}, Ziren Jiang \email{jian0746@@umn.edu}
#'
#' @references Jiang, Z., Crainiceanu, C., and Cui, E. (2025). Tutorial on Bayesian Functional Regression Using Stan. \emph{Statistics in Medicine}, 44(20-22), e70265.
#'
#' @examples
#' \dontrun{
#' # Simulate functional data with two underlying principal components
#' set.seed(1)
#' n  <- 100   # number of subjects
#' M  <- 50    # number of functional observation points
#' tindex <- seq(0, 1, length.out = M)
#'
#' # True mean function and eigenfunctions
#' mu_true <- sin(pi * tindex)
#' phi1    <- sqrt(2) * sin(2 * pi * tindex)
#' phi2    <- sqrt(2) * cos(2 * pi * tindex)
#'
#' # Simulate scores and noisy observations
#' xi1 <- rnorm(n, 0, sqrt(2))
#' xi2 <- rnorm(n, 0, sqrt(0.5))
#' Y_mat <- matrix(rep(mu_true, n), nrow = n, byrow = TRUE) +
#'          outer(xi1, phi1) + outer(xi2, phi2) +
#'          matrix(rnorm(n * M, sd = 0.3), nrow = n)
#'
#' dat <- data.frame(inx = 1:n)
#' dat$Y_mat <- Y_mat
#'
#' # Fit the Bayesian FPCA model
#' fit_fpca <- fpca_bayes(
#'   formula     = Y_mat ~ 1,
#'   data        = dat,
#'   spline_type = "bs",
#'   spline_df   = 10,
#'   niter       = 2000,
#'   nwarmup     = 1000,
#'   nchain      = 3
#' )
#'
#' # Posterior mean of the mean function
#' mu_est <- apply(fit_fpca$mu, 2, mean)
#' plot(tindex, mu_est, type = "l", ylab = expression(hat(mu)(t)))
#'
#' # Posterior means of the FPC scores and eigenvalues
#' scores_est  <- apply(fit_fpca$scores, c(2, 3), mean)
#' evalues_est <- apply(fit_fpca$evalues, 2, mean)
#' }
#'
#' @import mgcv
#' @import splines2
#' @import rstan
#' @importFrom refund fpca.sc
#' @importFrom stats gaussian
#' @export fpca_bayes

fpca_bayes <- function(formula, data,
                       npc = NULL, runStan = TRUE,
                       niter = 3000, nwarmup = 1000, nchain = 3, ncores = 1,
                       spline_type = "bs", spline_df = 10){

  ## Set family to "fpca"
  family <- "fpca"

  # Parse the formula to get the functional response variable name
  formula_val <- stats::as.formula(formula)
  y_var <- formula_val[[2]]

  # Extract the functional response matrix
  Y_mat <- data[[y_var]]
  if(is.null(Y_mat)){
    stop("Response variable '", as.character(y_var), "' not found in data!")
  }
  Y_mat <- unclass(Y_mat)
  if(!is.matrix(Y_mat)){
    stop("Response must be a matrix (functional response)!")
  }

  N_num <- nrow(Y_mat)
  M_num <- ncol(Y_mat)

  ##--------------------------------------------------------------------------
  ## Construct Stan data
  ##--------------------------------------------------------------------------

  # Initial frequentist FPCA to obtain the fixed eigenfunction basis
  init.fpca <- refund::fpca.sc(Y_mat, npc = npc)
  J_num <- ifelse(is.null(dim(init.fpca$efunctions)[2]), 1, dim(init.fpca$efunctions)[2])
  Phi_mat <- t(init.fpca$efunctions)  # J_num x M_num

  # Construct the spline basis and penalty matrix for the population mean function
  data_temp <- data.frame(inx = 1:N_num)
  data_temp[["yindex.vec"]] <- matrix(rep(1:M_num, N_num), nrow = N_num)
  object <- 1
  obj_exp <- paste0("object=s(yindex.vec, bs = \"", spline_type, "\", k =", spline_df, ")")
  eval(parse(text = obj_exp))
  dk <- ExtractData(object, data_temp, NULL)
  crspline.cons <- mgcv::smooth.construct(object, dk$data, dk$knots)
  mu.base <- crspline.cons$X
  mu.S <- crspline.cons$S[[1]]
  maXX <- norm(mu.base, type = "I") ^ 2
  maS <- norm(mu.S) / maXX
  mu.S <- mu.S / maS
  K_num <- dim(crspline.cons$X)[2]
  Psi_mat <- t(crspline.cons$X)  # K_num x M_num

  Standata_use <- list(
    N_num = N_num,
    M_num = M_num,
    J_num = J_num,
    K_num = K_num,
    Y_mat = Y_mat,
    Phi_mat = Phi_mat,
    Psi_mat = Psi_mat,
    S_mat = mu.S
  )

  ##--------------------------------------------------------------------------
  ## Build the Stan code
  ##--------------------------------------------------------------------------

  # Data block
  Stancode_data <- paste0(
    "data{ \n",
    "   int<lower=1> N_num; // Total number of subjects\n",
    "   int<lower=1> M_num; // Total number of observed functional time points\n",
    "   int<lower=1> J_num; // Number of FPCA eigenfunctions\n",
    "   int<lower=1> K_num; // Number of spline basis for the mean function\n",
    "   matrix[N_num, M_num] Y_mat; // Functional response\n",
    "   matrix[J_num, M_num] Phi_mat; // Matrix of FPCA eigenfunctions\n",
    "   matrix[K_num, M_num] Psi_mat; // Matrix of spline basis for the mean function\n",
    "   matrix[K_num, K_num] S_mat; // Penalty matrix for the mean function\n",
    "}"
  )

  # Transformed data block
  Stancode_transdata <- "transformed data {\n}"

  # Parameters block
  Stancode_para <- paste0(
    "parameters{ \n",
    "   vector[K_num] alpha; // Spline coefficients for the mean function\n",
    "   matrix[N_num, J_num] xi; // FPCA Scores\n",
    "   real<lower=0> sigma_eps; // Standard deviation of independent error\n",
    "   vector<lower=0>[J_num] lambda; // FPCA eigenvalues\n",
    "   real<lower=0> sigma_mu; // Smoothing parameter for the mean function\n",
    "}"
  )

  # Transformed parameters block
  Stancode_transpara <- paste0(
    "transformed parameters { \n",
    " real lprior = 0;\n",
    "}"
  )

  # Model block
  Stancode_model <- paste0(
    "model{ \n",
    "   matrix[N_num, M_num] mu;\n",
    "   row_vector[M_num] mu_t;\n",
    "   // Population mean function\n",
    "   mu_t = alpha' * Psi_mat;\n",
    "   // Fitted mean matrix: population mean + subject-specific FPC deviations\n",
    "   mu = rep_matrix(mu_t, N_num) + xi * Phi_mat;\n",
    "   // Log-likelihood for functional response\n",
    "   target += - N_num * M_num * log(sigma_eps) / 2 - sum((mu - Y_mat)^2) / (2 * sigma_eps^2);\n",
    "   // Prior for the penalized spline coefficients of the mean function\n",
    "   target += (- alpha' * S_mat * alpha) / (2 * sigma_mu^2);\n",
    "   target += inv_gamma_lpdf(sigma_mu^2|0.001,0.001);\n",
    "    // Prior for the FPCA scores \n",
    "   for(nj in 1:J_num){\n",
    "        target += - N_num * log(lambda[nj]) / 2 - sum((xi[,nj])^2) / (2 * lambda[nj]^2);\n",
    "        target += inv_gamma_lpdf(lambda[nj]^2|0.001,0.001);\n",
    "   }\n",
    "   target += inv_gamma_lpdf(sigma_eps^2|0.001,0.001);\n",
    "}\n"
  )

  # Combine all blocks of the Stan code
  Stancode_function <- NULL ## no manual function needed for FPCA
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

    # Extract posterior samples from the Stan fit object
    fit.samp <- rstan::extract(fit, permuted = TRUE)

    # Reconstruct posterior samples of the population mean function
    # Each row is one sample; each column is one location of the functional domain
    mu <- fit.samp$alpha %*% Psi_mat

    # FPC scores (n_samples x N_num x J_num)
    scores <- fit.samp$xi

    # FPC eigenvalue SDs (n_samples x J_num)
    evalues <- fit.samp$lambda

    # Residual standard deviation (n_samples)
    sigma <- fit.samp$sigma_eps

    # Fixed eigenfunctions (M_num x J_num)
    efunctions <- init.fpca$efunctions
  }else{
    # Do not run the Stan program
    fit <- NA ## no sampling
    fit.samp <- NA ## no sampling
    mu <- NA ## no sampling
    scores <- NA ## no sampling
    evalues <- NA ## no sampling
    sigma <- NA ## no sampling
    efunctions <- init.fpca$efunctions
  }

  # Organize the results
  res.fit <- list(stanfit = fit,
                  stancode = Stancode_use,
                  standata = Standata_use,
                  mu = mu,
                  efunctions = efunctions,
                  scores = scores,
                  evalues = evalues,
                  sigma = sigma,
                  family = "fpca")
  class(res.fit) <- "refundBayes"
  return(res.fit)
}

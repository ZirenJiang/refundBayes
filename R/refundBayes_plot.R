#----------------------------------------------------------------------------
#' Plot the estimated functional coefficients with the corresponding credible interval(s).
#'
#' Produces coefficient plots tailored to the model family:
#' \itemize{
#'   \item \code{sofr_bayes()}, \code{fcox_bayes()}: one curve plot per
#'     functional predictor coefficient \eqn{\beta(s)}, with pointwise and/or
#'     CMA credible bands.
#'   \item \code{fosr_bayes()}: one curve plot per scalar predictor
#'     coefficient function \eqn{\alpha_p(t)}.
#'   \item \code{fofr_bayes()}: curve plots for scalar predictor coefficient
#'     functions \eqn{\alpha_p(t)} (if any), followed by heatmap plots for
#'     each bivariate coefficient surface \eqn{\beta_q(s, t)} (posterior
#'     mean) from the functional predictors.
#'   \item \code{fpca_bayes()}: posterior mean function \eqn{\mu(t)} with a
#'     pointwise credible band; a combined plot of the (fixed) FPC
#'     eigenfunctions \eqn{\phi_j(t)}; a point-and-error-bar plot of the
#'     posterior of the eigenvalue SDs \eqn{\lambda_j}; and a histogram of the
#'     residual-SD posterior \eqn{\sigma_\epsilon}.
#' }
#'
#' @param x A fitted object returned by \code{sofr_bayes()},
#'   \code{fosr_bayes()}, \code{fofr_bayes()}, \code{fcox_bayes()}, or
#'   \code{fpca_bayes()}.
#' @param prob Coverage probability for the credible interval(s). Defaults to 0.95.
#' @param include Type of interval to include. \code{"pointwise"} produces pointwise credible intervals;
#'   \code{"CMA"} produces the CMA credible band; \code{"both"} produces both. Defaults to \code{"both"}.
#'   Only used for \code{sofr_bayes()} / \code{fcox_bayes()} curve plots.
#' @param ... Other parameters
#'
#' @return A named list of \code{ggplot} objects. For FoFR, scalar-predictor
#'   curves are named \code{scalar_<p>} and bivariate-predictor heatmaps are
#'   named \code{bivar_<q>}. For FPCA, the plots are named \code{mu},
#'   \code{efunctions}, \code{evalues}, and \code{sigma}.
#'
#' @import ggplot2
#' @importFrom stats as.formula qnorm quantile sd
#' @export
#' @method plot refundBayes

plot.refundBayes=function(x = NULL,...,prob = 0.95,include = "both"){

  if(x$family=="functional"){
    ## for functional outcomes (FoSR)
    plot.res=list()
    for(inx.effect in 1:dim(x$func_coef)[2]){
      curve_est=x$func_coef[,inx.effect,]
      mean.curve.est=apply(curve_est,2,mean)
      upper.curve.est.quantile=apply(curve_est,2,function(xx){stats::quantile(xx,probs = (1+prob)/2)})
      lower.curve.est.quantile=apply(curve_est,2,function(xx){stats::quantile(xx,probs = (1-prob)/2)})
      upper.curve.wald=upper.curve.est.quantile
      lower.curve.wald=lower.curve.est.quantile

      for(i in 1:length(upper.curve.wald)){
        upper.curve.wald[i] <- mean.curve.est[i]+qnorm((1+prob)/2)*sd(curve_est[,i])
        lower.curve.wald[i] <- mean.curve.est[i]-qnorm((1+prob)/2)*sd(curve_est[,i])
      }
      
      # We now only use the quantile approach to calculate the credible interval 
      upper.curve.est = upper.curve.est.quantile
      lower.curve.est = lower.curve.est.quantile

      plotdata <- data.frame(value = c(mean.curve.est,
                                  upper.curve.est,
                                  lower.curve.est),
                          xmat = c(1: length(upper.curve.est),
                                 1: length(upper.curve.est),
                                 1: length(upper.curve.est)),
                          type = c(rep("mean",length(upper.curve.est)),
                                 rep("upper",length(upper.curve.est)),
                                 rep("lower",length(upper.curve.est))))
      plot.res[[inx.effect]] <- ggplot2::ggplot(plotdata,aes(y = .data$value, x = .data$xmat))+
        ggplot2::geom_line(aes(type=.data$type,color=.data$type))+
        ggplot2::ylab(colnames(x$Standata[["X_s"]])[inx.effect])
    }
  }else if(x$family == "fofr"){

    ## for FoFR (function-on-function regression):
    ##   - Scalar predictors (optional): functional coefficients alpha_p(t)
    ##     plotted as curves with pointwise credible ribbons.
    ##   - Functional predictors: bivariate coefficients beta_q(s, t) plotted
    ##     as heatmaps of the posterior mean.
    plot.res <- list()

    # ---- Helper: nice names for scalar predictors (fall back to generic) ----
    scalar_names <- NULL
    if(!is.null(x$scalar_func_coef)){
      P_num_plot <- dim(x$scalar_func_coef)[2]
      scalar_names <- colnames(x$standata$X_mat)
      if(is.null(scalar_names) || length(scalar_names) != P_num_plot){
        scalar_names <- paste0("Scalar predictor ", seq_len(P_num_plot))
      }
    }

    # ==========================================================
    # (1) Curve plots for scalar predictor coefficients alpha_p(t)
    # ==========================================================
    if(!is.null(x$scalar_func_coef)){
      P_num_plot <- dim(x$scalar_func_coef)[2]
      for(p in seq_len(P_num_plot)){
        curve_est <- x$scalar_func_coef[, p, , drop = TRUE]  # n_samp x M

        mean.curve.est  <- apply(curve_est, 2, mean)
        upper.curve.est <- apply(curve_est, 2,
                                 function(xx){stats::quantile(xx, probs = (1+prob)/2)})
        lower.curve.est <- apply(curve_est, 2,
                                 function(xx){stats::quantile(xx, probs = (1-prob)/2)})

        plotdata <- data.frame(
          value  = mean.curve.est,
          upper  = upper.curve.est,
          lower  = lower.curve.est,
          tindex = seq_along(mean.curve.est)
        )

        plot.res[[paste0("scalar_", p)]] <-
          ggplot2::ggplot(plotdata,
                          ggplot2::aes(x = .data$tindex, y = .data$value)) +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower,
                                            ymax = .data$upper,
                                            fill = "Pointwise CI"),
                               alpha = 0.25,
                               colour = "black", linetype = "dashed") +
          ggplot2::geom_line(colour = "#0072B2", linewidth = 1) +
          ggplot2::scale_fill_manual(values = c("Pointwise CI" = "darkgrey"),
                                     name = "Interval") +
          ggplot2::labs(x = "Response-domain index (t)",
                        y = paste0("alpha(t) for ", scalar_names[p]),
                        title = paste0("Scalar predictor coefficient: ",
                                       scalar_names[p])) +
          ggplot2::theme_minimal()
      }
    }

    # ==========================================================
    # (2) Heatmap plots for bivariate coefficients beta_q(s, t)
    # ==========================================================
    if(!is.null(x$bivar_func_coef)){
      bivar_names <- names(x$bivar_func_coef)
      if(is.null(bivar_names) ||
         length(bivar_names) != length(x$bivar_func_coef)){
        bivar_names <- paste0("Functional predictor ",
                              seq_along(x$bivar_func_coef))
      }

      for(q in seq_along(x$bivar_func_coef)){
        bivar_q <- x$bivar_func_coef[[q]]
        if(is.null(bivar_q) || length(dim(bivar_q)) != 3) next

        beta_mean  <- apply(bivar_q, c(2, 3), mean)                # L x M
        L_q        <- nrow(beta_mean)
        M_q        <- ncol(beta_mean)

        heat_df <- expand.grid(s_idx = seq_len(L_q),
                               t_idx = seq_len(M_q))
        heat_df$beta <- as.vector(beta_mean)

        # Diverging colour scale centered at 0 so sign is interpretable
        max_abs <- max(abs(beta_mean), na.rm = TRUE)

        plot.res[[paste0("bivar_", q)]] <-
          ggplot2::ggplot(heat_df,
                          ggplot2::aes(x = .data$s_idx,
                                       y = .data$t_idx,
                                       fill = .data$beta)) +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_gradient2(
            low = "#2166AC", mid = "white", high = "#B2182B",
            midpoint = 0, limits = c(-max_abs, max_abs),
            name = expression(hat(beta)(s, t))
          ) +
          ggplot2::coord_fixed() +
          ggplot2::labs(
            x = "Predictor domain (s)",
            y = "Response domain (t)",
            title = paste0("Bivariate coefficient: ", bivar_names[q]),
            subtitle = "Posterior mean"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(panel.grid = ggplot2::element_blank())
      }
    }

  }else if(x$family == "fpca"){

    ## for FPCA (functional principal component analysis):
    ##   (1) Posterior of the mean function mu(t) with pointwise credible band
    ##   (2) All fixed FPC eigenfunctions phi_j(t) overlaid in one plot
    ##   (3) Posterior of the eigenvalue SDs lambda_j (point + error bar)
    ##   (4) Posterior of the residual SD sigma_eps (histogram)
    plot.res <- list()

    # ---------- (1) Mean function mu(t) ----------
    mu_samp  <- x$mu                                              # n_samp x M
    mu_mean  <- apply(mu_samp, 2, mean)
    mu_upper <- apply(mu_samp, 2,
                      function(xx){stats::quantile(xx, probs = (1+prob)/2)})
    mu_lower <- apply(mu_samp, 2,
                      function(xx){stats::quantile(xx, probs = (1-prob)/2)})

    plotdata_mu <- data.frame(
      value  = mu_mean,
      upper  = mu_upper,
      lower  = mu_lower,
      tindex = seq_along(mu_mean)
    )

    plot.res[["mu"]] <-
      ggplot2::ggplot(plotdata_mu,
                      ggplot2::aes(x = .data$tindex, y = .data$value)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower,
                                        ymax = .data$upper,
                                        fill = "Pointwise CI"),
                           alpha = 0.25, colour = "black", linetype = "dashed") +
      ggplot2::geom_line(colour = "#0072B2", linewidth = 1) +
      ggplot2::scale_fill_manual(values = c("Pointwise CI" = "darkgrey"),
                                 name = "Interval") +
      ggplot2::labs(x = "Functional-domain index (t)",
                    y = expression(mu(t)),
                    title = "Posterior mean function") +
      ggplot2::theme_minimal()

    # ---------- (2) Fixed FPC eigenfunctions phi_j(t) ----------
    Phi        <- x$efunctions                                    # M x J
    J_num_plot <- ncol(Phi)
    M_num_plot <- nrow(Phi)

    plotdata_phi <- data.frame(
      value     = as.vector(Phi),
      tindex    = rep(seq_len(M_num_plot), times = J_num_plot),
      component = factor(rep(paste0("PC", seq_len(J_num_plot)),
                             each = M_num_plot),
                         levels = paste0("PC", seq_len(J_num_plot)))
    )

    plot.res[["efunctions"]] <-
      ggplot2::ggplot(plotdata_phi,
                      ggplot2::aes(x = .data$tindex,
                                   y = .data$value,
                                   colour = .data$component)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(x = "Functional-domain index (t)",
                    y = expression(phi[j](t)),
                    colour = "Component",
                    title = "FPC eigenfunctions (fixed from initial FPCA)") +
      ggplot2::theme_minimal()

    # ---------- (3) Eigenvalue SDs lambda_j ----------
    lambda_samp  <- x$evalues                                     # n_samp x J
    lambda_mean  <- apply(lambda_samp, 2, mean)
    lambda_upper <- apply(lambda_samp, 2,
                          function(xx){stats::quantile(xx, probs = (1+prob)/2)})
    lambda_lower <- apply(lambda_samp, 2,
                          function(xx){stats::quantile(xx, probs = (1-prob)/2)})

    plotdata_lambda <- data.frame(
      component = seq_along(lambda_mean),
      value     = lambda_mean,
      upper     = lambda_upper,
      lower     = lambda_lower
    )

    plot.res[["evalues"]] <-
      ggplot2::ggplot(plotdata_lambda,
                      ggplot2::aes(x = .data$component, y = .data$value)) +
      ggplot2::geom_point(colour = "#0072B2", size = 2.5) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower,
                                          ymax = .data$upper),
                             width = 0.15, colour = "#0072B2") +
      ggplot2::scale_x_continuous(breaks = seq_along(lambda_mean)) +
      ggplot2::labs(x = "Component j",
                    y = expression(lambda[j]),
                    title = "Posterior of FPC eigenvalue SDs") +
      ggplot2::theme_minimal()

    # ---------- (4) Residual SD sigma_eps ----------
    sigma_samp     <- x$sigma
    plotdata_sigma <- data.frame(sigma = sigma_samp)
    sigma_mean     <- mean(sigma_samp)

    plot.res[["sigma"]] <-
      ggplot2::ggplot(plotdata_sigma,
                      ggplot2::aes(x = .data$sigma)) +
      ggplot2::geom_histogram(bins = 30, fill = "steelblue", colour = "white") +
      ggplot2::geom_vline(xintercept = sigma_mean, colour = "red",
                          linewidth = 1) +
      ggplot2::labs(x = expression(sigma[epsilon]),
                    y = "Frequency",
                    title = "Posterior of residual standard deviation",
                    subtitle = "Red line = posterior mean") +
      ggplot2::theme_minimal()

  }else{

    ## for scalar/survival outcomes (SoFR / functional Cox regression)
    n.func.coeff=length(x$func_coef)
    plot.res=list()
    for(inx.effect in 1:n.func.coeff){
      curve_est=x$func_coef[[inx.effect]]
      mean.curve.est=apply(curve_est,2,mean)
      upper.curve.est.quantile=apply(curve_est,2,function(xx){stats::quantile(xx,probs = (1+prob)/2)})
      lower.curve.est.quantile=apply(curve_est,2,function(xx){stats::quantile(xx,probs = (1-prob)/2)})

      if(include=="pointwise"){
        upper.curve.est = upper.curve.est.quantile
        lower.curve.est = lower.curve.est.quantile
        CMA.upper.curve.est = NA
        CMA.lower.curve.est = NA
      }
      if(include=="CMA"){
        func_coeff_sd=apply(curve_est,2,sd)
        func_coeff_extr=1:dim(curve_est)[1]
        func_coeff_est=apply(curve_est,2,mean)
        for(i in 1:dim(curve_est)[1]){
          func_coeff_extr[i]=max(abs(curve_est[i,]-func_coeff_est)/func_coeff_sd)
        }
        cutpoint=stats::quantile(func_coeff_extr,probs = prob)
        CMA.upper.curve.est = mean.curve.est + cutpoint*func_coeff_sd
        CMA.lower.curve.est = mean.curve.est - cutpoint*func_coeff_sd
        upper.curve.est = NA
        lower.curve.est = NA
      }
      if(include=="both"){
        func_coeff_sd=apply(curve_est,2,sd)
        func_coeff_extr=1:dim(curve_est)[1]
        func_coeff_est=apply(curve_est,2,mean)
        for(i in 1:dim(curve_est)[1]){
          func_coeff_extr[i]=max(abs(curve_est[i,]-func_coeff_est)/func_coeff_sd)
        }
        cutpoint=stats::quantile(func_coeff_extr,probs = prob)
        CMA.upper.curve.est = mean.curve.est + cutpoint*func_coeff_sd
        CMA.lower.curve.est = mean.curve.est - cutpoint*func_coeff_sd
        upper.curve.est = upper.curve.est.quantile
        lower.curve.est = lower.curve.est.quantile
      }

      plotdata2 <- data.frame( value = c(mean.curve.est),
                               CI.upper = c(upper.curve.est),
                               CI.lower = c(lower.curve.est),
                               CMA.upper = c(CMA.upper.curve.est),
                               CMA.lower = c(CMA.lower.curve.est),
                               xmat = c(1:length(upper.curve.est)),
                               Method = c(rep("Bayesian",length(mean.curve.est))))

      plot.res[[inx.effect]] <- ggplot2::ggplot(plotdata2, aes(y = .data$value, x = .data$xmat, color = .data$Method)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_ribbon(aes(ymin = .data$CI.lower, ymax = .data$CI.upper, fill = "Pointwise CI"),
                    alpha = 0.2, color = "black", linetype = "dashed") +
        ggplot2::geom_ribbon(aes(ymin = .data$CMA.lower, ymax = .data$CMA.upper, fill = "CMA CI"),
                    alpha = 0.2, linetype = "dotted", color = "black") +
        ggplot2::ylab("Functional Effect") + ggplot2::facet_wrap(~ .data$Method)+
        ggplot2::scale_fill_manual(values = c("Pointwise CI" = "darkgrey", "CMA CI" = "lightgrey"), name = "Interval Type") +
        ggplot2::theme_minimal()
    }
  }
  return(plot.res)
}


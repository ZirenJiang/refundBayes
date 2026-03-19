# ============================================================
#  Targeted fix: diagnose & directly patch fcox_bayes.Rd
#  (and the other two Rd files) so pkgdown can parse them.
#
#  Run from the package root:
#    setwd("D:/Code/R/18. FDA Survival/refundBayes")
#    source("fix_rd_direct.R")
# ============================================================

pkg_root <- getwd()
stopifnot(file.exists(file.path(pkg_root, "DESCRIPTION")))

# ── Helper: rewrite an Rd file with a correct \examples{} block ──────────────
patch_rd <- function(rd_file, examples_text) {
  lines <- readLines(rd_file, warn = FALSE, encoding = "UTF-8")
  
  # Strip any existing \examples{...} block (could be empty or malformed)
  in_ex  <- FALSE
  depth  <- 0
  keep   <- logical(length(lines))
  for (i in seq_along(lines)) {
    ln <- lines[i]
    if (!in_ex && grepl("^\\\\examples\\s*\\{", ln)) {
      in_ex <- TRUE
      depth <- 0
    }
    if (in_ex) {
      depth <- depth +
        nchar(gsub("[^{]", "", ln)) -
        nchar(gsub("[^}]", "", ln))
      keep[i] <- FALSE
      if (depth <= 0) { in_ex <- FALSE }
    } else {
      keep[i] <- TRUE
    }
  }
  
  # Also strip any stray \donttest lines that leaked into \details{}
  cleaned <- lines[keep]
  cleaned <- cleaned[!grepl("^\\\\donttest", cleaned)]
  
  # Remove trailing blank lines, then append the new \examples block
  while (length(cleaned) > 0 && trimws(tail(cleaned, 1)) == "")
    cleaned <- head(cleaned, -1)
  
  out <- c(cleaned, "", examples_text, "")
  writeLines(out, rd_file, useBytes = FALSE)
  cat("  ✔ Patched", basename(rd_file), "\n")
}


# ── Inspect what is currently in fcox_bayes.Rd ───────────────────────────────
cat("── Diagnosing fcox_bayes.Rd ────────────────────────────\n")
rd_fcox <- file.path(pkg_root, "man", "fcox_bayes.Rd")
raw     <- readLines(rd_fcox, warn = FALSE)

# Show lines from \details onward so we can see if \donttest leaked in
detail_start <- grep("^\\\\(details|examples|donttest)", raw)
if (length(detail_start)) {
  cat("  Lines", detail_start[1], "onwards:\n")
  cat(paste(raw[detail_start[1]:min(length(raw), detail_start[1]+30)],
            collapse = "\n"), "\n\n")
} else {
  cat("  No \\details / \\examples / \\donttest found — file may be empty\n\n")
}


# ── Also fix the R source typo more robustly (handles CRLF on Windows) ───────
cat("── Re-fixing R source @examples tags ──────────────────\n")
fix_r_source <- function(rfile) {
  raw_bytes <- readBin(rfile, "raw", file.info(rfile)$size)
  txt       <- rawToChar(raw_bytes)
  # Replace #'#' @examples  (with or without \r before the newline)
  fixed     <- gsub("(#')#'( @examples)", "\\1\\2", txt, perl = TRUE)
  if (!identical(txt, fixed)) {
    writeBin(charToRaw(fixed), rfile)
    cat("  ✔ Fixed @examples tag in", basename(rfile), "\n")
  } else {
    cat("  ℹ No #'#' typo found in", basename(rfile),
        "— checking current tag:\n")
    hits <- regmatches(txt, gregexpr("(?m)^.{0,5}@examples.{0,20}",
                                     txt, perl = TRUE))[[1]]
    cat("   ", paste(hits, collapse = "\n    "), "\n")
  }
}
fix_r_source(file.path(pkg_root, "R", "sofr_bayes.R"))
fix_r_source(file.path(pkg_root, "R", "fosr_bayes.R"))
fix_r_source(file.path(pkg_root, "R", "fcox_bayes.R"))


# ── Patch the three Rd files directly ────────────────────────────────────────
cat("\n── Patching man/*.Rd files directly ───────────────────\n")

patch_rd(
  file.path(pkg_root, "man", "sofr_bayes.Rd"),
  c(
    "\\examples{",
    "\\donttest{",
    "# Simulate data for a Gaussian SoFR model",
    "set.seed(1)",
    "n  <- 100  # number of subjects",
    "L  <- 50   # number of functional domain points",
    "Lindex    <- seq(0, 1, length.out = L)",
    "X_func    <- matrix(rnorm(n * L), nrow = n)",
    "age       <- rnorm(n)",
    "beta_true <- sin(2 * pi * Lindex)",
    "Y <- X_func \\%*\\% beta_true / L + 0.5 * age + rnorm(n, sd = 0.5)",
    "",
    "dat         <- data.frame(Y = Y, age = age)",
    "dat$X_func  <- X_func",
    "dat$Lindex  <- matrix(rep(Lindex, n), nrow = n, byrow = TRUE)",
    "",
    "# Generate Stan code and data without running the sampler",
    "fit_sofr <- sofr_bayes(",
    "  formula = Y ~ age + s(Lindex, by = X_func, bs = \"cr\", k = 10),",
    "  data    = dat,",
    "  family  = \"gaussian\",",
    "  runStan = FALSE",
    ")",
    "}",
    "}"
  )
)

patch_rd(
  file.path(pkg_root, "man", "fosr_bayes.Rd"),
  c(
    "\\examples{",
    "\\donttest{",
    "# Simulate data for a Function-on-Scalar Regression model",
    "set.seed(1)",
    "n      <- 100",
    "M      <- 50",
    "tindex <- seq(0, 1, length.out = M)",
    "age    <- rnorm(n)",
    "sex    <- rbinom(n, 1, 0.5)",
    "beta_age <- sin(2 * pi * tindex)",
    "beta_sex <- cos(2 * pi * tindex)",
    "epsilon  <- matrix(rnorm(n * M, sd = 0.3), nrow = n)",
    "Y_mat    <- outer(age, beta_age) + outer(sex, beta_sex) + epsilon",
    "",
    "dat       <- data.frame(age = age, sex = sex)",
    "dat$Y_mat <- Y_mat",
    "",
    "# Generate Stan code and data without running the sampler",
    "fit_fosr <- fosr_bayes(",
    "  formula     = Y_mat ~ age + sex,",
    "  data        = dat,",
    "  spline_type = \"bs\",",
    "  spline_df   = 10,",
    "  runStan     = FALSE",
    ")",
    "}",
    "}"
  )
)

patch_rd(
  file.path(pkg_root, "man", "fcox_bayes.Rd"),
  c(
    "\\examples{",
    "\\donttest{",
    "# Simulate survival data with a functional predictor",
    "set.seed(1)",
    "n      <- 150",
    "L      <- 50",
    "Lindex <- seq(0, 1, length.out = L)",
    "X_func <- matrix(rnorm(n * L), nrow = n)",
    "age    <- rnorm(n)",
    "",
    "beta_true <- cos(2 * pi * Lindex)",
    "lp        <- X_func \\%*\\% beta_true / L + 0.3 * age",
    "",
    "time      <- rexp(n, rate = exp(lp))",
    "cens_time <- runif(n, min = 0.5, max = 3)",
    "obs_time  <- pmin(time, cens_time)",
    "cens_ind  <- as.integer(time <= cens_time)",
    "",
    "dat        <- data.frame(obs_time = obs_time, age = age)",
    "dat$X_func <- X_func",
    "dat$Lindex <- matrix(rep(Lindex, n), nrow = n, byrow = TRUE)",
    "",
    "# Generate Stan code and data without running the sampler",
    "fit_cox <- fcox_bayes(",
    "  formula = obs_time ~ age + s(Lindex, by = X_func, bs = \"cr\", k = 10),",
    "  data    = dat,",
    "  cens    = cens_ind,",
    "  runStan = FALSE",
    ")",
    "}",
    "}"
  )
)


# ── Verify the patched files parse cleanly ────────────────────────────────────
cat("\n── Verifying Rd files parse cleanly ───────────────────\n")
for (rd in c("fcox_bayes.Rd", "fosr_bayes.Rd", "sofr_bayes.Rd")) {
  path <- file.path(pkg_root, "man", rd)
  tryCatch({
    tools::parse_Rd(path)
    cat("  ✔", rd, "parses OK\n")
  }, error = function(e) {
    cat("  ✖", rd, "STILL HAS PARSE ERROR:\n   ", conditionMessage(e), "\n")
  })
}


# ── Rebuild site ──────────────────────────────────────────────────────────────
cat("\n── Rebuilding pkgdown site ─────────────────────────────\n")
pkgdown::build_site()
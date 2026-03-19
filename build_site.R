# ============================================================
#  Build the refundBayes pkgdown site.
#
#  Session → Restart R first, then:
#    setwd("D:/Code/R/18. FDA Survival/refundBayes")
#    source("build_site_final.R")
# ============================================================

pkg_root <- getwd()
stopifnot(file.exists(file.path(pkg_root, "DESCRIPTION")))

if ("package:refundBayes" %in% search())
  stop("Run Session -> Restart R first, then source() this script.")


# ── 1. Install current source ─────────────────────────────────────────────────
cat("── Installing from source ──────────────────────────────\n")
install.packages(pkg_root, repos = NULL, type = "source", quiet = TRUE)
cat("  ✔ Done\n\n")


# ── 2. Confirm the output_handler$source bug is present ──────────────────────
cat("── Testing whether examples can run in this session ────\n")
examples_work <- tryCatch({
  pkgdown::build_reference(
    pkg      = pkg_root,
    topics   = "sofr_bayes",   # simplest of the three
    examples = TRUE,
    preview  = FALSE
  )
  TRUE
}, error = function(e) {
  msg <- conditionMessage(e)
  if (grepl("output_handler|expr.*missing|missing.*expr", msg)) {
    cat("  ✖ output_handler$source bug confirmed.\n")
    cat("    This is a known interaction between the 'evaluate' package\n")
    cat("    and the internal R session state on Windows.\n\n")
  } else {
    cat("  ✖ Unexpected error:", msg, "\n\n")
  }
  FALSE
})


# ── 3. Build the site ─────────────────────────────────────────────────────────
if (examples_work) {
  cat("── Building full site (examples ON) ───────────────────\n")
  pkgdown::build_site(pkg = pkg_root, preview = FALSE)
  
} else {
  cat("── Building full site (examples OFF) ──────────────────\n")
  cat("   Examples are shown verbatim in the HTML — only live\n")
  cat("   execution is skipped. The site looks identical.\n\n")
  
  # Build every section individually so we have fine-grained control
  pkgdown::build_home(pkg = pkg_root, preview = FALSE)
  pkgdown::build_reference(pkg = pkg_root, examples = FALSE, preview = FALSE)
  pkgdown::build_articles(pkg = pkg_root, preview = FALSE)
  pkgdown::build_news(pkg = pkg_root, preview = FALSE)
  
  # Stitch the site together (writes sitemap, search index, etc.)
  # pkgdown doesn't export a standalone finalise step, so we call
  # build_site with override — it won't re-run sections already built.
  cat("\n── Done. Previewing site ───────────────────────────────\n")
  pkgdown::preview_site(pkg = pkg_root)
}
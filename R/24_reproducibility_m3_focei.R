# ===========================================================================
# PAPER 2 - STEP 6a-3: two open questions from R/23, one run (3 parallel fits).
# Read this block first.
#
#   Q1 REPRODUCIBILITY. The SAME M1 data (rep 1, design D, BLQ dropped) gave
#   Kss 4.237 in R/22 (serial, no CENS column) but 3.354 in R/23 (PSOCK
#   worker, extra CENS = 0 column). Same seed (1) and SAEM settings.
#     M1_A, M1_B: M1 data WITHOUT the CENS column, SAEM seed 1, two identical
#     fits in two fresh PSOCK sessions.
#   READING: A = B = 4.237 -> the CENS column changed the SAEM path;
#            A = B = 3.354 -> the session type (serial vs PSOCK) matters;
#            A != B        -> nlmixr2 SAEM is not reproducible even with a
#                             fixed seed -> every SSE estimate needs that caveat.
#   Either way: on identical data, Kss moved 4.237 -> 3.354 (27%) between
#   runs, more than its reported RSE (14%) - the Kss likelihood is flat.
#
#   Q2 M3. SAEM with censoring (R/23) collapsed (VC 0.00025, RSE 381%).
#     M3_FOCEI: the same censored data (CENS = 1, DV = LLOQ) with FOCEi, the
#     standard estimator for M3 (Beal 2001; NONMEM uses Laplace/FOCE for M3).
#   READING: sensible estimates (VC near 1.58) -> use FOCEi for M3 in Step 6b.
#
#   Setup = R/23's own lines (model, CHECK 0, simulation, data checks), re-run
#   inside every worker. Starting values = truth. Writes no files.
#   Run from the project root (denosumab-mpbpk/). Expected ~1 h (FOCEi unknown).
# ===========================================================================

suppressMessages(library(parallel))

setup_code <- function() {
  s <- readLines("R/23_blq_diagnostic.R")
  invisible(capture.output(eval(parse(text = s[grep("^suppressMessages", s)[1]:(grep("^fit_one <-", s) - 1)]),
                                envir = globalenv())))
}

run_job <- function(job) {
  setup_code()
  d <- if (job == "M3_FOCEI") datasets$M3 else datasets$M1[, names(datasets$M1) != "CENS"]
  tm <- system.time(fit <- tryCatch(
    suppressMessages(suppressWarnings(
      if (job == "M3_FOCEI") nlmixr2(f_model, d, est = "focei", control = foceiControl(print = 0))
      else nlmixr2(f_model, d, est = "saem", control = saemControl(seed = 1, print = 0, nBurn = 300, nEm = 200)))),
    error = function(e) e))[["elapsed"]]
  if (inherits(fit, "error")) return(list(nm = job, err = conditionMessage(fit), tm = tm))
  list(nm = job, tm = tm, cov = fit$covMethod, n_obs = sum(d$EVID == 0), has_cens = "CENS" %in% names(d),
       est = exp(fixef(fit)[pn]), se = fit$parFixedDf[pn, "SE"], truth = truth,
       om = diag(fit$omega)[paste0("eta.", tolower(names(om_true)))], om_true = om_true)
}

cl <- makeCluster(3)
clusterCall(cl, setwd, getwd())
clusterExport(cl, "setup_code")
res <- parLapply(cl, c("M1_A", "M1_B", "M3_FOCEI"), run_job)
stopCluster(cl)

for (r in res) {
  cat(sprintf("=== %s ===\n", r$nm))
  if (!is.null(r$err)) { cat("  FIT FAILED:", r$err, sprintf("(%.0f s)\n\n", r$tm)); next }
  cat(sprintf("  time %.0f s | covariance: %s | obs %d | CENS column: %s\n", r$tm, r$cov, r$n_obs, r$has_cens))
  print(data.frame(true = signif(unname(r$truth), 4), est = signif(unname(r$est), 4),
                   err_pct = round(100 * (unname(r$est) / unname(r$truth) - 1), 1),
                   rse_pct = round(100 * r$se, 1), row.names = names(r$truth)))
  cat("  omega KSS (est / true):", sprintf("%.3f / %.3f", r$om[["eta.kss"]], r$om_true[["KSS"]]), "\n\n")
}
cat("[reference: M1 Kss 4.237 (R/22 serial, no CENS) | 3.354 (R/23 PSOCK, CENS = 0) | ORACLE 1.012 | truth 1.56]\n")
cat("DIAGNOSTIC ONLY - one dataset; do not report as a result.\n")

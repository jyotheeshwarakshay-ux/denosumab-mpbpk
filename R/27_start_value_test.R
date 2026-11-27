# ===========================================================================
# PAPER 2 - STEP 6b-check: do the M3 FOCEi fits stop at their starting values?
# Read this block first.
#
#   EVIDENCE (R/25, 18 fits, all started at the truth): D_M3 reps 14 and 20
#   returned Kss 1.56000 and R0 15.2300 exactly (= starting values), and
#   7 of 11 M3 fits report SE(log Kss) 0.0002-0.005, while every M1 fit
#   (no CENS column) moved and has SE ~0.16-0.24.
#
#   TEST on ONE trial (rep 14; same simulated data as R/25), 4 FOCEi fits:
#     1 D_M3_low   design D, M3, ALL structural starts x0.7, Kss x0.5
#     2 A_M1_truth design A, BLQ dropped (A has ~0.1% BLQ, so nearly the
#                  same data as A_M3), NO CENS column, start = truth
#                  (R/25 A_M3 rep 14: Kss 1.58589, SE 0.00045)
#     3 D_M3_high  design D, M3, ALL structural starts x1.3, Kss x2
#     4 D_M3_bobyqa design D, M3, start = truth, derivative-free outer
#                  optimiser (bobyqa) instead of the default
#   READING:
#     1 and 3 end at the SAME Kss (not their starts), SE ~0.1-0.4
#        -> truth-start fits stopped early; restart M3 jobs from perturbed
#           starts (R/25 change) and re-fit.
#     1 and 3 stay near their own starts -> M3 likelihood/optimiser fails to
#        move; check 4 (bobyqa) as the alternative.
#     2 moves (SE ~0.1-0.4) while A_M3 rep 14 did not -> the censoring
#        (CENS) code path is what stops the fit.
#   Per fit: all 8 estimates to 6 digits, SE(log Kss), objf, covariance,
#   outer iterations (rows of parHistData) and optimiser message.
#
#   2 PSOCK workers (less heat), jobs in order (chunk.size = 1).
#   ~2 rounds x ~1.5-2 h. Writes no files. Run from the project root.
# ===========================================================================

suppressMessages(library(parallel))

setup_code <- function() {
  s <- readLines("R/23_blq_diagnostic.R")
  invisible(capture.output(eval(parse(text = s[grep("^suppressMessages", s)[1]:(grep("^fit_one <-", s) - 1)]),
                                envir = globalenv())))
}

make_data <- function(dn, rep, meth) {
  obs <- sim_obs(designs[[dn]], rep)
  dose <- expand.grid(ID = 1:80, TIME = dose_times)
  dr <- data.frame(ID = dose$ID, TIME = dose$TIME, EVID = 1, AMT = dose_nmol, CMT = "depot", DV = NA)
  if (meth == "M1") {
    o <- obs[!obs$BLQ, ]
    orows <- data.frame(ID = o$ID, TIME = o$TIME, EVID = 0, AMT = 0, CMT = "ctot", DV = o$DV)
  } else {
    dr$CENS <- 0
    orows <- data.frame(ID = obs$ID, TIME = obs$TIME, EVID = 0, AMT = 0, CMT = "ctot",
                        DV = ifelse(obs$BLQ, LLOQ, obs$DV), CENS = as.integer(obs$BLQ))
  }
  d <- rbind(dr, orows); d[order(d$ID, d$TIME, -d$EVID), ]
}

run_job <- function(job) {
  setup_code()
  b <- bpop
  start <- function(f, fk) ui |> ini(lka = log(b[["KA"]] * f), lvc = log(b[["VC"]] * f), lvp = log(b[["VP"]] * f),
                                     lcl = log(b[["CL"]] * f), lkint = log(b[["KINT"]] * f),
                                     lkss = log(b[["KSS"]] * fk), lksyn = log(b[["KSYN"]] * f), lr0 = log(b[["R0"]] * f))
  spec <- switch(job,
    D_M3_low    = list(d = make_data("D", 14, "M3"), m = start(0.7, 0.5), ctl = foceiControl(print = 0)),
    A_M1_truth  = list(d = make_data("A", 14, "M1"), m = ui, ctl = foceiControl(print = 0)),
    D_M3_high   = list(d = make_data("D", 14, "M3"), m = start(1.3, 2), ctl = foceiControl(print = 0)),
    D_M3_bobyqa = list(d = make_data("D", 14, "M3"), m = ui, ctl = foceiControl(print = 0, outerOpt = "bobyqa")))
  st <- exp(spec$m$theta[pn])
  tm <- system.time(fit <- tryCatch(
    suppressMessages(suppressWarnings(nlmixr2(spec$m, spec$d, est = "focei", control = spec$ctl))),
    error = function(e) e))[["elapsed"]]
  if (inherits(fit, "error")) return(list(job = job, err = conditionMessage(fit), tm = tm, start = st))
  list(job = job, tm = tm, start = st, est = exp(fixef(fit)[pn]),
       se_kss = fit$parFixedDf["lkss", "SE"], objf = fit$objf, cov = fit$covMethod,
       n_obs = sum(spec$d$EVID == 0), cens = if (is.null(spec$d$CENS)) "none" else sum(spec$d$CENS),
       iter = tryCatch(nrow(fit$parHistData), error = function(e) NA),
       msg = tryCatch(paste(fit$message, collapse = " "), error = function(e) NA))
}

jobs <- c("D_M3_low", "A_M1_truth", "D_M3_high", "D_M3_bobyqa")
cl <- makeCluster(2)
clusterCall(cl, setwd, getwd())
clusterExport(cl, c("setup_code", "make_data"))
res <- parLapplyLB(cl, jobs, run_job, chunk.size = 1)
stopCluster(cl)

for (r in res) {
  cat(sprintf("=== %s ===\n", r$job))
  if (!is.null(r$err)) { cat("  FIT FAILED:", r$err, sprintf("(%.0f s)\n\n", r$tm)); next }
  cat(sprintf("  time %.0f s | obs %d | censored %s | cov %s | objf %.2f | outer iter %s\n",
              r$tm, r$n_obs, r$cens, r$cov, r$objf, r$iter))
  cat("  optimiser message:", r$msg, "\n")
  print(data.frame(start = signif(unname(r$start), 6), est = signif(unname(r$est), 6),
                   row.names = toupper(sub("^l", "", names(r$est)))))
  cat(sprintf("  SE(log Kss) %.4g\n\n", r$se_kss))
}
cat("[R/25 rep 14 (truth start): D_M3 Kss 1.56000 SE 0.00056 objf 4309.09 | A_M3 Kss 1.58589 SE 0.00045 objf 5107.22 | D_M1 Kss 4.07653]\n")
cat("DIAGNOSTIC ONLY - one trial; do not report as a result.\n")
